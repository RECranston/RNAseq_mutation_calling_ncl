#!/bin/bash
#SBATCH --account=XXXXX
#SBATCH --partition=default_free
#SBATCH --mem=100G
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=8
#SBATCH --job-name=rna_vep
#SBATCH --output=logs/vep_%A_%a.out
#SBATCH --array=1-100%100

# Script to run an array of vep annotation jobs on file ids from a paired list of fastq files
# Ruth Cranston 2026

[ $# -ne 3 ] && { echo -en \
"\nRuth Cranston 2026\n\n
*** Script to run vep annotation jobs on a list of sample ids from the original fastq file sample sheet [sample name] [fastq1] [fastq2] (tab delimited sheet). 
Runs in current directory. Input dir is location of rna editing filtered files. Output directory is created.
<sample sheet> <input dir (relative)> <output dir (relative)>
example run: sbatch ./script7_vep.sh sample_sheet.txt output_rnaed_filtering/ output_vep/ *** \n\n" ; exit 1; }

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
module load VEP/113.3-GCC-13.3.0
echo -en " * Environment set up.\n"

set -euo pipefail

# Remove output dir file if exists
mkdir -p ${OUTPUT_DIR}
mkdir -p logs
mkdir -p ${REFERENCE_DIR}/vep_cache

# Get the correct row for this array task
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_SHEET})

SAMPLE_ID=$(echo $LINE | awk '{print $1}')
FILE1=$(echo $LINE | awk '{print $2}')
FILE2=$(echo $LINE | awk '{print $3}')

echo "Processing sample: ${SAMPLE_ID}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"

# vep annotation of filtered variants
vep \
    --input_file ${INPUT_DIR}${SAMPLE_ID}_tumor_no_editing.vcf \
    --output_file ${OUTPUT_DIR}${SAMPLE_ID}_tumor_annotated.vcf \
    --format vcf --vcf \
    --cache --offline \
    --dir_cache ${REFERENCE_DIR}/vep_cache \
    --assembly GRCh38 \
    --everything \
    --fork 8 \
    --fasta ${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta

echo -ne "*** vep annotation complete! ***"

echo -ne "*** All done! ***"

