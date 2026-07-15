#!/bin/bash
#SBATCH --account=XXXXX
#SBATCH --partition=default_free
#SBATCH --mem=50G
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=40
#SBATCH --job-name=index_prep
#SBATCH --output=slurm_log_%j.out

# Script to prepare RNAseq mutation analysis environment including downloading references and build STAR indexes
# Ruth Cranston 2026

[ $# -ne 0 ] && { echo -en \
"\nRuth Cranston 2026\n\n
*** Script to build RNAseq mutation analysis environment including downloading references and build STAR indexes. 
Runs in current directory *** \n\n" ; exit 1; }

# Set variables
BASE_DIR="$PWD"
REFERENCE_DIR=${BASE_DIR}/"References"
STAR_INDEX_DIR=${BASE_DIR}/"STAR_indexes"

# Load modules
echo -en "Loading modules...\n"
module --force purge
module load STAR
module load BEDTools
module load VEP
echo -en "Environment set up.\n"

# Load Python >=3.10 for gcloud CLI compatibility
module load Python/3.11.5-GCCcore-13.2.0
export CLOUDSDK_PYTHON=$(which python3)

# install google cloud sdk locally for download of files from google repo
GCLOUD_DIR="$(pwd)/google-cloud-sdk"

if [ ! -d "${GCLOUD_DIR}" ]; then
    echo "Installing Google Cloud SDK..."
    curl -s https://sdk.cloud.google.com | bash -s -- --disable-prompts --install-dir="$(pwd)"
fi

# Add to PATH for this script's session
export PATH="${GCLOUD_DIR}/bin:${PATH}"

# Confirm it's available
which gsutil
gsutil --version

# Configure anonymous access for public buckets
gcloud config set auth/disable_credentials True
gcloud config unset project 2>/dev/null || true

# checking if reference files are present - if not then download these
echo "Detecting if Reference directory is present"
if [[ -d ${REFERENCE_DIR} ]]
then
    echo -en " * ${REFERENCE_DIR} exists, no need to re-download\n\n"
else
    echo -en " * ${REFERENCE_DIR} not found, downloading now...\n"

    # downloads: gatk specific genome build, indexes and dict files,  gencodev44 gtf files, BQSR files, Mutect2 resources
    mkdir ${REFERENCE_DIR}

    # gatk genome files
    gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta ${REFERENCE_DIR}/
    gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai ${REFERENCE_DIR}/
    gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict ${REFERENCE_DIR}/

    # BQSR files
    curl -L "https://storage.googleapis.com/storage/v1/b/gcp-public-data--broad-references/o/hg38%2Fv0%2FHomo_sapiens_assembly38.dbsnp138.vcf?alt=media" -o ${REFERENCE_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf
    curl -L "https://storage.googleapis.com/storage/v1/b/gcp-public-data--broad-references/o/hg38%2Fv0%2FHomo_sapiens_assembly38.dbsnp138.vcf.idx?alt=media" -o ${REFERENCE_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf.idx
    gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz ${REFERENCE_DIR}/
    gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi ${REFERENCE_DIR}/
    gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz ${REFERENCE_DIR}/
    gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi ${REFERENCE_DIR}/

    # Mutect2 files
    gsutil cp gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz "${REFERENCE_DIR}/"
    gsutil cp gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi "${REFERENCE_DIR}/"
    gsutil cp gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz "${REFERENCE_DIR}/"
    gsutil cp gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi "${REFERENCE_DIR}/"
    gsutil cp gs://gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz "${REFERENCE_DIR}/"
    gsutil cp gs://gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz.tbi "${REFERENCE_DIR}/"

    # gencode gtf
    wget "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz" -O "${REFERENCE_DIR}/gencode.v44.primary_assembly.annotation.gtf.gz"
    gunzip "${REFERENCE_DIR}/gencode.v44.primary_assembly.annotation.gtf.gz"

    
    # Download and prepare the reference file for filtering of RNA editing events from Mutect2 detected variants
    wget http://rediportal.cloud.ba.infn.it/download/TABLE1_hg38_v3.txt.gz -O ${REFERENCE_DIR}/TABLE1_hg38_v3.txt.gz

    # unzip and convert to BED (skip header row, convert 1-based pos -> 0-based BED start)
    zcat ${REFERENCE_DIR}/TABLE1_hg38_v3.txt.gz | \
	awk 'BEGIN{OFS="\t"} $3 ~ /^[0-9]+$/ {gsub(/\r$/,""); print $2, $3-1, $3}' > ${REFERENCE_DIR}/REDIportal_GRCh38.unsorted.bed

    # sort - needed if you want to use bedtools' fast -sorted mode (recommended, ~16M sites)
    bedtools sort -i ${REFERENCE_DIR}/REDIportal_GRCh38.unsorted.bed \
	     -g ${REFERENCE_DIR}/GRCh38.genome > ${REFERENCE_DIR}/REDIportal_GRCh38.bed
    
    # one-off: generate a chrom-sizes file from your reference dict
    cut -f1,2 ${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta.fai > ${REFERENCE_DIR}/GRCh38.genome
    echo "Reference download complete."

    # Download files for vep cache
    mkdir -p ${REFERENCE_DIR}/vep_cache
    wget https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/homo_sapiens_vep_113_GRCh38.tar.gz \
	 -O ${REFERENCE_DIR}/vep_cache/homo_sapiens_vep_113_GRCh38.tar.gz
    tar xzf ${REFERENCE_DIR}/vep_cache/homo_sapiens_vep_113_GRCh38.tar.gz -C ${REFERENCE_DIR}/vep_cache

fi


# STAR Reference genome Preparation:
# checking if STAR indexes are present - if not then create these
echo "Detecting if STAR indexes are present"
if [[ -d ${STAR_INDEX_DIR} ]]
then
    echo -en " * ${STAR_INDEX_DIR} exists, no need to re-create\n\n"
else
    echo -en " * ${STAR_INDEX_DIR} not found, preparing STAR reference genomes...\n"

    mkdir ${STAR_INDEX_DIR}

    STAR --runMode genomeGenerate \
         --genomeDir ${STAR_INDEX_DIR}/STAR_GRCh38 \
         --genomeFastaFiles ${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta \
         --sjdbGTFfile ${REFERENCE_DIR}/gencode.v44.primary_assembly.annotation.gtf \
         --runThreadN 8

fi

echo -en "*** All done ***"

