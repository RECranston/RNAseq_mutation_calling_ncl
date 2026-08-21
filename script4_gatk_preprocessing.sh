#!/bin/bash
#SBATCH --account=XXXXX
#SBATCH --partition=default_free
#SBATCH --mem=100G
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=8
#SBATCH --job-name=gatk_preprocess
#SBATCH --output=logs/gatk_preprocess_%A_%a.out
#SBATCH --array=1-100%100

# Script to run an array of gatk preprocessing jobs on file ids from a paired list of fastq files
# Ruth Cranston 2026

[ $# -ne 3 ] && { echo -en \
"\nRuth Cranston 2026\n\n
*** Script to run gatk preprocessing of a list of sample ids from the original fastq file sample sheet [sample name] [fastq1] [fastq2] (tab delimited sheet). 
Runs in current directory. Output directory is created.
<sample sheet> <input dir (relative)> <output dir (relative)>
example run: sbatch ./script4_gatk_preprocessing.sh sample_sheet.txt star_output/ output_preprocessing/ *** \n\n" ; exit 1; }


# Set variables
BASE_DIR="$PWD"
ASSEMBLY="GRCh37"
STAR_INDEX_DIR=${BASE_DIR}/"STAR_indexes/STAR_GRCh38"
STAR_INDEX_DIR=${BASE_DIR}/STAR_indexes/STAR_${ASSEMBLY}
REFERENCE_DIR=${BASE_DIR}/References/${ASSEMBLY}
TMPDIR=${BASE_DIR}/tmp

SAMPLE_SHEET=$1
OUTPUT_DIR=${BASE_DIR}/$3
STAR_DIR=${BASE_DIR}/$2
PASS1_DIR=${STAR_DIR}"PASS1/"
PASS2_DIR=${STAR_DIR}"PASS2/"

# Load modules
echo -en " * Loading modules...\n"
module --force purge
module load GATK/4.6.0.0-GCCcore-13.2.0-Java-17
echo -en " * Environment set up.\n"

# Remove output dir file if exists & make tmp dir
mkdir -p ${OUTPUT_DIR}
mkdir -p logs
mkdir -p ${TMPDIR}

# Get the correct row for this array task
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_SHEET})

SAMPLE_ID=$(echo $LINE | awk '{print $1}')

echo "Processing sample: ${SAMPLE_ID}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"

# Set genome reference files
if [ "${ASSEMBLY}" == "GRCh38" ]; then
    REF_FASTA=${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta
    DBSNP=${REFERENCE_DIR}/Homo_sapiens_assembly38.dbsnp138.vcf
    KNOWN_INDELS=${REFERENCE_DIR}/Homo_sapiens_assembly38.known_indels.vcf.gz
else
    REF_FASTA=${REFERENCE_DIR}/Homo_sapiens_assembly19.fasta
    DBSNP=${REFERENCE_DIR}/dbsnp_138.b37.vcf.gz
    KNOWN_INDELS=${REFERENCE_DIR}/Mills_and_1000G_gold_standard.indels.b37.vcf
fi

# Mark duplicates
gatk MarkDuplicates \
     --java-options "-Xmx90g -Djava.io.tmpdir=${TMPDIR}" \
     -I ${PASS2_DIR}/${SAMPLE_ID}_Aligned.sortedByCoord.out.bam \
     -O ${OUTPUT_DIR}${SAMPLE_ID}_marked_dup.bam \
     -M ${OUTPUT_DIR}${SAMPLE_ID}_marked_dup_metrics.txt \
     --CREATE_INDEX true \
     --TMP_DIR ${TMPDIR}

echo -ne "*** Mark duplicates done! ***\n"

# Split reads at splice junctions
gatk SplitNCigarReads \
     --java-options "-Xmx90g -Djava.io.tmpdir=${TMPDIR}" \
     -R ${REF_FASTA} \
     -I ${OUTPUT_DIR}${SAMPLE_ID}_marked_dup.bam \
     -O ${OUTPUT_DIR}${SAMPLE_ID}_split.bam \
     --tmp-dir ${TMPDIR}

# catch here for if SplitNCigarReads fails
if [ ! -s ${OUTPUT_DIR}${SAMPLE_ID}_split.bam ]; then
    echo "ERROR: SplitNCigarReads produced empty output for ${SAMPLE_ID}" >&2
    exit 1
fi

echo -ne "*** Split reads done! ***\n"


# BQSR (use dbSNP + known indels)
gatk BaseRecalibrator \
     --java-options "-Xmx90g" \
     -R ${REF_FASTA} \
     -I ${OUTPUT_DIR}${SAMPLE_ID}_split.bam \
     --known-sites ${DBSNP} \
     --known-sites ${KNOWN_INDELS} \
     -O ${OUTPUT_DIR}${SAMPLE_ID}_recal.table

echo -ne "*** BaseRecalibrator done! ***\n"


gatk ApplyBQSR \
     --java-options "-Xmx90g -Djava.io.tmpdir=${TMPDIR}" \
     -R ${REF_FASTA} \
     -I ${OUTPUT_DIR}${SAMPLE_ID}_split.bam \
     --bqsr-recal-file ${OUTPUT_DIR}${SAMPLE_ID}_recal.table \
     -O ${OUTPUT_DIR}${SAMPLE_ID}_recal.bam \
     --tmp-dir ${TMPDIR}

echo -ne "*** BQSR done! ***\n"

echo -ne "*** All done! ***\n"
