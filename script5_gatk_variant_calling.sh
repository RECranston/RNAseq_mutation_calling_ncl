#!/bin/bash
#SBATCH --account=XXXXX
#SBATCH --partition=default_free
#SBATCH --mem=100G
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=8
#SBATCH --job-name=gatk_var_calling
#SBATCH --output=logs/gatk_var_calling_%A_%a.out
#SBATCH --array=1-100%100

# Script to run an array of gatk variant calling on file ids from a paired list of fastq files
# Ruth Cranston 2026

[ $# -ne 3 ] && { echo -en \
"\nRuth Cranston 2026\n\n
*** Script to run gatk mutation calling on a list of sample ids from the original fastq file sample sheet [sample name] [fastq1] [fastq2] (tab delimited sheet). 
Runs in current directory. Input dir is location of gatk preprocessed files. Output directory is created.
<sample sheet> <input dir (relative)> <output dir (relative)>
example run: sbatch ./script5_gatk_variant_calling.sh sample_sheet.txt output_preprocessing/ output_mutation_calling/ *** \n\n" ; exit 1; }

# --array=1-5%10 means run array job IDs 1-5 with a maximum of 10 running at once

# Set variables
BASE_DIR="$PWD"
STAR_INDEX_DIR=${BASE_DIR}/"STAR_indexes/STAR_GRCh38"
SAMPLE_SHEET=$1
INPUT_DIR=${BASE_DIR}/$2
OUTPUT_DIR=${BASE_DIR}/$3
REFERENCE_DIR=${BASE_DIR}/"References"


# Load modules
echo -en " * Loading modules...\n"
module --force purge
module load GATK/4.6.0.0-GCCcore-13.2.0-Java-17
echo -en " * Environment set up.\n"

set -euo pipefail

# Remove output dir file if exists
mkdir -p ${OUTPUT_DIR}
mkdir -p logs

# Get the correct row for this array task
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_SHEET})

SAMPLE_ID=$(echo $LINE | awk '{print $1}')
FILE1=$(echo $LINE | awk '{print $2}')
FILE2=$(echo $LINE | awk '{print $3}')

echo "Processing sample: ${SAMPLE_ID}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"

# Mutect2 (tumor-only mode)
# --dont-use-soft-clipped-bases is important for RNA processing
gatk Mutect2 \
     --java-options "-Xmx90g" \
     -R ${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta \
     -I ${INPUT_DIR}${SAMPLE_ID}_recal.bam \
     --tumor-sample ${SAMPLE_ID} \
     --germline-resource ${REFERENCE_DIR}/af-only-gnomad.hg38.vcf.gz \
     --panel-of-normals ${REFERENCE_DIR}/1000g_pon.hg38.vcf.gz \
     --dont-use-soft-clipped-bases \
     --f1r2-tar-gz ${OUTPUT_DIR}${SAMPLE_ID}_f1r2.tar.gz \
     -O ${OUTPUT_DIR}${SAMPLE_ID}_tumor_raw.vcf.gz

echo -ne "*** Mutect2 finished! ***"

# Learn orientation bias (important for RNA oxidation artifacts)
gatk LearnReadOrientationModel \
     --java-options "-Xmx90g" \
     -I ${OUTPUT_DIR}${SAMPLE_ID}_f1r2.tar.gz \
     -O ${OUTPUT_DIR}${SAMPLE_ID}_artifact_prior.tar.gz

echo -ne "*** LearnReadOrientationModel finished! ***"

# Get contamination estimate
gatk GetPileupSummaries \
    --java-options "-Xmx90g" \
    -I ${INPUT_DIR}${SAMPLE_ID}_recal.bam \
    -V ${REFERENCE_DIR}/small_exac_common_3.hg38.vcf.gz \
    -L ${REFERENCE_DIR}/small_exac_common_3.hg38.vcf.gz \
    -O ${OUTPUT_DIR}${SAMPLE_ID}_pileup_summaries.table

gatk CalculateContamination \
    --java-options "-Xmx90g" \
    -I ${OUTPUT_DIR}${SAMPLE_ID}_pileup_summaries.table \
    -O ${OUTPUT_DIR}${SAMPLE_ID}_contamination.table

echo -ne "*** CalculateContamination finished! ***"

# Apply all filters
gatk FilterMutectCalls \
     --java-options "-Xmx90g" \
     -R ${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta \
     -V ${OUTPUT_DIR}${SAMPLE_ID}_tumor_raw.vcf.gz \
     --ob-priors ${OUTPUT_DIR}${SAMPLE_ID}_artifact_prior.tar.gz \
     --contamination-table ${OUTPUT_DIR}${SAMPLE_ID}_contamination.table \
     --min-allele-fraction 0.05 \
     -O ${OUTPUT_DIR}${SAMPLE_ID}_tumor_filtered.vcf.gz

echo -ne "*** Mutect call filtering finished! ***"

# Extract PASS variants only
gatk SelectVariants \
    --java-options "-Xmx90g" \
     -V ${OUTPUT_DIR}${SAMPLE_ID}_tumor_filtered.vcf.gz \
    --exclude-filtered \
    -O ${OUTPUT_DIR}${SAMPLE_ID}_tumor_filtered_PASS.vcf.gz

echo -ne "*** Extract PASS variants only finished! ***"

echo -ne "*** All done! ***"

