#!/bin/bash
#SBATCH --account=XXXXX
#SBATCH --partition=default_free
#SBATCH --mem=100G
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=8
#SBATCH --job-name=rna_ed_filter
#SBATCH --output=logs/rna_ed_filter_%A_%a.out
#SBATCH --array=1-100%100

# Script to run an array of rna editing filter jobs on file ids from a paired list of fastq files
# Ruth Cranston 2026

[ $# -ne 3 ] && { echo -en \
"\nRuth Cranston 2026\n\n
*** Script to run rna editing filter jobs on a list of sample ids from the fastq file sample sheet [sample name] [fastq1] [fastq2] (tab delimited sheet). 
Runs in current directory. Input dir is location of mutect2 called and filtered files. Output directory is created.
<sample sheet> <input dir (relative)> <output dir (relative)>
example run: sbatch ./script6_rna_editing_filter.sh sample_sheet.txt output_mutation_calling/ output_rnaed_filtering/ *** \n\n" ; exit 1; }

# --array=1-5%10 means run array job IDs 1-5 with a maximum of 10 running at once

# Set variables
BASE_DIR="$PWD"
ASSEMBLY="GRCh37"
SAMPLE_SHEET=$1
INPUT_DIR=${BASE_DIR}/$2
OUTPUT_DIR=${BASE_DIR}/$3
REFERENCE_DIR=${BASE_DIR}/References/${ASSEMBLY}


# Load modules
echo -en " * Loading modules...\n"
module --force purge
module load BEDTools
echo -en " * Environment set up.\n"

set -euo pipefail

# Remove output dir file if exists
mkdir -p ${OUTPUT_DIR}
mkdir -p logs

# Get the correct row for this array task
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_SHEET})

SAMPLE_ID=$(echo $LINE | awk '{print $1}')

echo "Processing sample: ${SAMPLE_ID}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"

# Convert to exclude known RNA editing sites
bedtools intersect \
    -a ${INPUT_DIR}${SAMPLE_ID}_tumor_filtered_PASS.vcf.gz \
    -b ${REFERENCE_DIR}/REDIportal_${ASSEMBLY}.bed \
    -sorted \
    -g ${REFERENCE_DIR}/${ASSEMBLY}.genome \
    -v \
    -header \
    > ${OUTPUT_DIR}${SAMPLE_ID}_tumor_no_editing.vcf

echo -ne "*** Filter for RNA editing events finished! ***"

echo -ne "*** All done! ***"

