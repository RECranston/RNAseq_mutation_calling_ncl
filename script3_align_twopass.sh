#!/bin/bash
#SBATCH --account=XXXXX
#SBATCH --partition=default_free
#SBATCH --mem=100G
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=8
#SBATCH --job-name=align
#SBATCH --output=logs/star_%A_%a.out
#SBATCH --array=1-100%100

# Script to run an array of two pass STAR alignment jobs on a paired list of fastq files
# Ruth Cranston 2026

[ $# -ne 3 ] && { echo -en \
"\nRuth Cranston 2026\n\n
*** Script to run two pass STAR alignment for a paired list of trimmed fastq files [sample name] [fastq1] [fastq2] (tab delimited sheet). 
Runs in current directory. Output directory is created. 
<sample sheet> <input dir (relative)> <output dir (relative)>
example run: sbatch ./script3_align_twopass.sh trimmed_sample_sheet.txt test_trimmed_fastq/ test_aligned_array/ *** \n\n" ; exit 1; }

# example run 
# sbatch --array=1-5 rna_seq_align_arr.sh test_trimmed_sample_sheet.txt test_trimmed_fastq/ test_aligned_array/
# --array=1-5%10 means run array job IDs 1-5 with a maximum of 10 running at once

# Set variables
BASE_DIR="$PWD"
STAR_INDEX_DIR=${BASE_DIR}/"STAR_indexes/STAR_GRCh38"
SAMPLE_SHEET=$1
INPUT_DIR=${BASE_DIR}/$2
OUTPUT_DIR=${BASE_DIR}/$3
PASS1_DIR=${OUTPUT_DIR}"PASS1/"
PASS2_DIR=${OUTPUT_DIR}"PASS2/"

# Load modules
echo -en " * Loading modules...\n"
module --force purge
module load STAR
module load SAMtools
echo -en " * Environment set up.\n"

# Remove output dir file if exists
mkdir -p ${OUTPUT_DIR}
mkdir -p logs
mkdir -p ${PASS1_DIR}
mkdir -p ${PASS2_DIR}

# Get the correct row for this array task
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_SHEET})

SAMPLE_ID=$(echo $LINE | awk '{print $1}')
FILE1=$(echo $LINE | awk '{print $2}')
FILE2=$(echo $LINE | awk '{print $3}')

echo "Processing sample: ${SAMPLE_ID}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"

# Run STAR pass 1
STAR --runMode alignReads \
     --runThreadN ${SLURM_CPUS_PER_TASK} \
     --genomeDir ${STAR_INDEX_DIR} \
     --outFileNamePrefix ${PASS1_DIR}${SAMPLE_ID}_\
     --outTmpDir ${PASS1_DIR}TmpDir_${SAMPLE_ID} \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --outSAMattributes NH HI AS NM MD \
     --outBAMcompression 9 \
     --genomeLoad NoSharedMemory \
     --readFilesIn ${INPUT_DIR}${FILE1} ${INPUT_DIR}${FILE2}

echo -ne "*** STAR pass 1 done! ***"

# Run STAR pass 2
STAR --runMode alignReads \
     --runThreadN ${SLURM_CPUS_PER_TASK} \
     --genomeDir ${STAR_INDEX_DIR} \
     --sjdbFileChrStartEnd ${PASS1_DIR}/${SAMPLE_ID}_SJ.out.tab \
     --readFilesIn ${INPUT_DIR}${FILE1} ${INPUT_DIR}${FILE2} \
     --readFilesCommand zcat \
     --outSAMtype BAM SortedByCoordinate \
     --outFileNamePrefix ${PASS2_DIR}/${SAMPLE_ID}_ \
     --outSAMattributes NH HI AS NM MD \
     --outSAMattrRGline ID:${SAMPLE_ID} SM:${SAMPLE_ID} PL:ILLUMINA LB:lib1 PU:unit1


echo -ne "*** STAR pass 2 done! ***"
