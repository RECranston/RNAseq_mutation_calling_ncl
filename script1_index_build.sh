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
REF38=${REFERENCE_DIR}/"GRCh38"
REF37=${REFERENCE_DIR}/"GRCh37"
STAR_INDEX_DIR=${BASE_DIR}/"STAR_indexes"

# Load modules
echo -en "Loading modules...\n"
module --force purge
module load STAR
module load BEDTools
module load GATK/4.6.0.0-GCCcore-13.2.0-Java-17
module load VEP/113.3-GCC-13.3.0
module load Python/3.11.5-GCCcore-13.2.0
set -euo pipefail
echo -en "Environment set up.\n"

# make directories
mkdir -p ${REF38}
mkdir -p ${REF37}
mkdir -p ${STAR_INDEX_DIR}
mkdir -p logs

# For gcloud CLI compatibility
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

######################
### GRCh38 section ###
######################


# checking if reference files are present - if not then download these
echo "Detecting if GRCh38 references are present"
if [[ -d ${REF38}/Homo_sapiens_assembly38.fasta ]];
then
    echo -en " * GRCh38 references already exist in ${REF38}, no need to re-download\n\n"
else
    echo -en " * Downloading GRCh38 references to ${REF38} now...\n"

    # downloads: gatk specific genome build, indexes and dict files,  gencodev44 gtf files, BQSR files, Mutect2 resources

    # gatk genome files
    gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta ${REF38}/
    gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.fasta.fai ${REF38}/
    gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.dict ${REF38}/

    # BQSR files
    curl -L "https://storage.googleapis.com/storage/v1/b/gcp-public-data--broad-references/o/hg38%2Fv0%2FHomo_sapiens_assembly38.dbsnp138.vcf?alt=media" -o ${REF38}/Homo_sapiens_assembly38.dbsnp138.vcf
    curl -L "https://storage.googleapis.com/storage/v1/b/gcp-public-data--broad-references/o/hg38%2Fv0%2FHomo_sapiens_assembly38.dbsnp138.vcf.idx?alt=media" -o ${REF38}/Homo_sapiens_assembly38.dbsnp138.vcf.idx

    gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz ${REF38}/
    gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Homo_sapiens_assembly38.known_indels.vcf.gz.tbi ${REF38}/
    gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz ${REF38}/
    gsutil cp gs://gcp-public-data--broad-references/hg38/v0/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz.tbi ${REF38}/

    # Mutect2 files
    gsutil cp gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz "${REF38}/"
    gsutil cp gs://gatk-best-practices/somatic-hg38/af-only-gnomad.hg38.vcf.gz.tbi "${REF38}/"
    gsutil cp gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz "${REF38}/"
    gsutil cp gs://gatk-best-practices/somatic-hg38/1000g_pon.hg38.vcf.gz.tbi "${REF38}/"
    gsutil cp gs://gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz "${REF38}/"
    gsutil cp gs://gatk-best-practices/somatic-hg38/small_exac_common_3.hg38.vcf.gz.tbi "${REF38}/"

    # gencode gtf
    wget "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.primary_assembly.annotation.gtf.gz" -O "${REF38}/gencode.v44.primary_assembly.annotation.gtf.gz"
    gunzip "${REF38}/gencode.v44.primary_assembly.annotation.gtf.gz"

    # make chromosome sizes file
    cut -f1,2 ${REF38}/Homo_sapiens_assembly38.fasta.fai > ${REF38}/GRCh38.genome

    # Download and prepare the reference file for filtering of RNA editing events from Mutect2 detected variants
    wget http://rediportal.cloud.ba.infn.it/download/TABLE1_hg38_v3.txt.gz -O ${REF38}/TABLE1_hg38_v3.txt.gz

    # unzip and convert to BED (skip header row, convert 1-based pos -> 0-based BED start)
    zcat ${REF38}/TABLE1_hg38_v3.txt.gz | \
	awk 'BEGIN{OFS="\t"} $3 ~ /^[0-9]+$/ {gsub(/\r$/,""); print $2, $3-1, $3}' > ${REF38}/REDIportal_GRCh38.unsorted.bed

    # sort with bedtools
    bedtools sort -i ${REF38}/REDIportal_GRCh38.unsorted.bed \
	     -g ${REF38}/GRCh38.genome > ${REF38}/REDIportal_GRCh38.bed
    rm ${REF38}/REDIportal_GRCh38.unsorted.bed
    
    # one-off: generate a chrom-sizes file from your reference dict
    cut -f1,2 ${REF38}/Homo_sapiens_assembly38.fasta.fai > ${REF38}/GRCh38.genome
    echo "GRCh38 reference download complete."

    # Download files for vep cache
    mkdir -p ${REF38}/vep_cache
    wget https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/homo_sapiens_vep_113_GRCh38.tar.gz \
	 -O ${REF38}/vep_cache/homo_sapiens_vep_113_GRCh38.tar.gz
    tar xzf ${REF38}/vep_cache/homo_sapiens_vep_113_GRCh38.tar.gz -C ${REF38}/vep_cache

    # STAR Reference genome preparation:
    echo "Detecting if STAR indexes are present"

    # STAR GRCh38
    if [[ -d ${STAR_INDEX_DIR}/STAR_GRCh38 ]]
    then
	echo -en " * ${STAR_INDEX_DIR}/STAR_GRCh38 exists, no need to re-create\n\n"
    else
	echo -en " * ${STAR_INDEX_DIR}/STAR_GRCh38 not found, preparing STAR reference genomes...\n"

	mkdir ${STAR_INDEX_DIR}/STAR_GRCh38

	STAR --runMode genomeGenerate \
             --genomeDir ${STAR_INDEX_DIR}/STAR_GRCh38 \
             --genomeFastaFiles ${REF38}/Homo_sapiens_assembly38.fasta \
             --sjdbGTFfile ${REF38}/gencode.v44.primary_assembly.annotation.gtf \
             --runThreadN 8

    fi

fi


######################
### GRCh37 section ###
######################

echo "Detecting if GRCh37 references are present"
if [[ -d ${REF37}/Homo_sapiens_assembly19.fasta ]];
then
    echo -en " * GRCh37 references already exist in ${REF37}, no need to re-download\n\n"
else
    echo -en " * Downloading GRCh38 references to ${REF37} now...\n"

    # downloads: gatk specific genome build, indexes and dict files,  gencodev44 gtf files, BQSR files, Mutect2 resources

    # gatk genome files
    gsutil cp gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta ${REF37}/
    gsutil cp gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.fasta.fai ${REF37}/
    gsutil cp gs://gcp-public-data--broad-references/hg19/v0/Homo_sapiens_assembly19.dict ${REF37}/

    # BQSR files
    gsutil cp gs://gcp-public-data--broad-references/hg19/v0/dbsnp_138.b37.vcf.gz ${REF37}/
    gsutil cp gs://gcp-public-data--broad-references/hg19/v0/dbsnp_138.b37.vcf.gz.tbi ${REF37}/
    gsutil cp gs://gatk-legacy-bundles/b37/Mills_and_1000G_gold_standard.indels.b37.vcf ${REF37}/
    gsutil cp gs://gatk-legacy-bundles/b37/Mills_and_1000G_gold_standard.indels.b37.vcf.idx ${REF37}/
    gsutil cp gs://gatk-legacy-bundles/b37/1000G_phase1.indels.b37.vcf ${REF37}/
    gsutil cp gs://gatk-legacy-bundles/b37/1000G_phase1.indels.b37.vcf.idx ${REF37}/

    # Mutect2 files
    gsutil cp gs://gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf ${REF37}/
    gsutil cp gs://gatk-best-practices/somatic-b37/af-only-gnomad.raw.sites.vcf.idx ${REF37}/
    gsutil cp gs://gatk-best-practices/somatic-b37/Mutect2-WGS-panel-b37.vcf ${REF37}/
    gsutil cp gs://gatk-best-practices/somatic-b37/Mutect2-WGS-panel-b37.vcf.idx ${REF37}/
    gsutil cp gs://gatk-best-practices/somatic-b37/small_exac_common_3.vcf ${REF37}/
    gsutil cp gs://gatk-best-practices/somatic-b37/small_exac_common_3.vcf.idx ${REF37}/
    
    # Download GRCh37 release 87 GTF
    # Release 87 = last native Ensembl GRCh37 release. No chr prefix
    # Using Ensembl GTF rather than GENCODE because GENCODE uses chr prefix
    wget "https://ftp.ensembl.org/pub/grch37/release-87/gtf/homo_sapiens/Homo_sapiens.GRCh37.87.gtf.gz" \
         -O "${REF37}/Homo_sapiens.GRCh37.87.gtf.gz"
    gunzip "${REF37}/Homo_sapiens.GRCh37.87.gtf.gz"
     
    # make chromosome sizes file
    cut -f1,2 ${REF37}/Homo_sapiens_assembly19.fasta.fai > ${REF37}/GRCh37.genome

    # Download and prepare the reference file for filtering of RNA editing events from Mutect2 detected variants
    wget https://rediportal.cloud.ba.infn.it/download/TABLE1_hg19_v3.txt.gz \
	 -O ${REF37}/TABLE1_hg19.txt.gz

    echo "GRCh37 reference download complete."

    # REDIPortal bed conversion
    if [[ -f ${REF37}/REDIportal_GRCh37.bed ]]; then
	echo " * GRCh37 REDIportal BED already present — skipping"
    else
	echo " * Converting REDIportal hg19 to BED..."
	zcat ${REF37}/TABLE1_hg19.txt.gz | \
	    awk 'BEGIN{OFS="\t"} $3 ~ /^[0-9]+$/ {gsub(/\r$/,""); print $2, $3-1, $3}' | \
	    sed 's/^chr//' | \
	    awk 'NR==FNR{valid[$1]=1; next} $1 in valid' ${REF37}/GRCh37.genome - \
		> ${REF37}/REDIportal_GRCh37.unsorted.bed

	# sort bed file
	bedtools sort -i ${REF37}/REDIportal_GRCh37.unsorted.bed \
		 -g ${REF37}/GRCh37.genome \
		 > ${REF37}/REDIportal_GRCh37.bed
	rm ${REF37}/REDIportal_GRCh37.unsorted.bed
	echo "GRCh37 REDIportal BED: $(wc -l < ${REF37}/REDIportal_GRCh37.bed) sites"
	head -3 ${REF37}/REDIportal_GRCh37.bed
    fi

    # Download files for vep cache
    mkdir -p ${REF37}/vep_cache
    wget https://ftp.ensembl.org/pub/release-113/variation/indexed_vep_cache/homo_sapiens_vep_113_GRCh37.tar.gz \
         -O ${REF37}/vep_cache/homo_sapiens_vep_113_GRCh37.tar.gz
    tar xzf ${REF37}/vep_cache/homo_sapiens_vep_113_GRCh37.tar.gz \
         -C ${REF37}/vep_cache

    # STAR index build
    echo "Detecting if STAR indexes are present"
    if [[ -d ${STAR_INDEX_DIR}/STAR_GRCh37 ]]
    then
	echo -en " * ${STAR_INDEX_DIR}/STAR_GRCh37 exists, no need to re-create\n\n"
    else
	echo -en " * ${STAR_INDEX_DIR}/STAR_GRCh37 not found, preparing STAR reference genomes...\n"

	mkdir ${STAR_INDEX_DIR}/STAR_GRCh37

	STAR --runMode genomeGenerate \
             --genomeDir ${STAR_INDEX_DIR}/STAR_GRCh37 \
             --genomeFastaFiles ${REF37}/Homo_sapiens_assembly19.fasta \
             --sjdbGTFfile ${REF37}/Homo_sapiens.GRCh37.87.gtf \
             --runThreadN 8

    fi
 
fi
 
echo -en "*** All done ***"

