#!/bin/bash
#SBATCH --account=XXXXX
#SBATCH --partition=default_free
#SBATCH --mem=100G
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=8
#SBATCH --job-name=gatk_var_calling
#SBATCH --output=logs/gatk_var_calling_%A_%a.out
#SBATCH --array=1-100%100

# Script to run an array of GATK somatic variant calling jobs (tumor-only) with rescue of clustered/germline filtered variants
# Ruth Cranston 2026

[ $# -ne 3 ] && { echo -en \
"\nRuth Cranston 2026\n\n
*** Script to run gatk mutation calling on a list of sample ids from the original fastq file
sample sheet [sample name] [fastq1] [fastq2] (tab delimited sheet).
Runs in current directory. Input dir is location of gatk preprocessed files. Output directory is created.
<sample sheet> <input dir (relative)> <output dir (relative)>
example run: sbatch ./script5_gatk_variant_calling.sh sample_sheet.txt output_preprocessing/ output_mutation_calling/ *** \n\n" ; exit 1; }

# Set variables
BASE_DIR="$PWD"
ASSEMBLY="GRCh37"
REFERENCE_DIR=${BASE_DIR}/References/${ASSEMBLY}
TMPDIR=${BASE_DIR}/tmp

SAMPLE_SHEET=$1
INPUT_DIR=${BASE_DIR}/$2
OUTPUT_DIR=${BASE_DIR}/$3

# Load modules
echo -en " * Loading modules...\n"
module --force purge
module load GATK/4.6.0.0-GCCcore-13.2.0-Java-17
module load Python/3.11.5-GCCcore-13.2.0
echo -en " * Environment set up.\n"

set -euo pipefail

# Make output, log and tmp dirs
mkdir -p ${OUTPUT_DIR} logs ${TMPDIR}

# Get the correct row for this array task
LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_SHEET})
SAMPLE_ID=$(echo $LINE | awk '{print $1}')

echo "Processing sample: ${SAMPLE_ID}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Assembly: ${ASSEMBLY}"

# Set genome reference files (b37 resources have no chr prefix)
if [ "${ASSEMBLY}" == "GRCh38" ]; then
    REF_FASTA=${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta
    GNOMAD=${REFERENCE_DIR}/af-only-gnomad.hg38.vcf.gz
    PON=${REFERENCE_DIR}/1000g_pon.hg38.vcf.gz
    EXAC=${REFERENCE_DIR}/small_exac_common_3.hg38.vcf.gz
else
    REF_FASTA=${REFERENCE_DIR}/Homo_sapiens_assembly19.fasta
    GNOMAD=${REFERENCE_DIR}/af-only-gnomad.raw.sites.vcf
    PON=${REFERENCE_DIR}/Mutect2-WGS-panel-b37.vcf
    EXAC=${REFERENCE_DIR}/small_exac_common_3.vcf
fi

# Mutect2 tumor-only (soft clips excluded as splice artefacts)
gatk Mutect2 \
     --java-options "-Xmx90g" \
     -R ${REF_FASTA} \
     -I ${INPUT_DIR}${SAMPLE_ID}_recal.bam \
     --tumor-sample ${SAMPLE_ID} \
     --germline-resource ${GNOMAD} \
     --panel-of-normals ${PON} \
     --dont-use-soft-clipped-bases \
     --f1r2-tar-gz ${OUTPUT_DIR}${SAMPLE_ID}_f1r2.tar.gz \
     -O ${OUTPUT_DIR}${SAMPLE_ID}_tumor_raw.vcf.gz

echo -ne "*** Mutect2 done! ***\n"

# Orientation model - output kept for reference only, not applied (invalid for strand-specific libraries)
gatk LearnReadOrientationModel \
     --java-options "-Xmx90g" \
     -I ${OUTPUT_DIR}${SAMPLE_ID}_f1r2.tar.gz \
     -O ${OUTPUT_DIR}${SAMPLE_ID}_artifact_prior.tar.gz

echo -ne "*** LearnReadOrientationModel done! ***\n"

# Contamination estimation
gatk GetPileupSummaries \
    --java-options "-Xmx90g" \
    -I ${INPUT_DIR}${SAMPLE_ID}_recal.bam \
    -V ${EXAC} \
    -L ${EXAC} \
    -O ${OUTPUT_DIR}${SAMPLE_ID}_pileup_summaries.table

gatk CalculateContamination \
    --java-options "-Xmx90g" \
    -I ${OUTPUT_DIR}${SAMPLE_ID}_pileup_summaries.table \
    -O ${OUTPUT_DIR}${SAMPLE_ID}_contamination.table

echo -ne "*** CalculateContamination done! ***\n"

# Filter calls (--ob-priors intentionally omitted)
gatk FilterMutectCalls \
     --java-options "-Xmx90g" \
     -R ${REF_FASTA} \
     -V ${OUTPUT_DIR}${SAMPLE_ID}_tumor_raw.vcf.gz \
     --contamination-table ${OUTPUT_DIR}${SAMPLE_ID}_contamination.table \
     --min-allele-fraction 0.05 \
     -O ${OUTPUT_DIR}${SAMPLE_ID}_tumor_filtered.vcf.gz

echo -ne "*** FilterMutectCalls done! ***\n"

# Select PASS variants and rescue variants filtered only as clustered_events/haplotype/germline
# Rescue rules (validated against matched DNA):
#   clustered_events;haplotype  -> TLOD > 10
#   any with germline           -> TLOD > 10 and POPAF > 6
# Rescued variants keep their original FILTER label
export RESCUE_INPUT_VCF="${OUTPUT_DIR}${SAMPLE_ID}_tumor_filtered.vcf.gz"
export RESCUE_OUTPUT_VCF="${OUTPUT_DIR}${SAMPLE_ID}_tumor_filtered_PASS.vcf.gz"
export RESCUE_SAMPLE_ID="${SAMPLE_ID}"

python3 << 'PYEOF'
import gzip
import re
import os

def get_info_value(info_field, key):
    match = re.search(r'(?:^|;)' + re.escape(key) + r'=([^;]+)', info_field)
    if match:
        try:
            return float(match.group(1))
        except ValueError:
            return None
    return None

def should_rescue(filter_tag, info_field):
    filter_components = set(filter_tag.split(';'))

    # only rescue if all filters are in the allowed set
    if not filter_components.issubset({'clustered_events', 'germline', 'haplotype'}):
        return False

    tlod = get_info_value(info_field, 'TLOD')
    if tlod is None or tlod <= 10:
        return False

    if 'germline' not in filter_components:
        return True

    popaf = get_info_value(info_field, 'POPAF')
    if popaf is None or popaf <= 6:
        return False

    return True


input_vcf  = os.environ['RESCUE_INPUT_VCF']
output_vcf = os.environ['RESCUE_OUTPUT_VCF']
sample_id  = os.environ['RESCUE_SAMPLE_ID']

pass_count     = 0
rescued_count  = 0
excluded_count = 0

opener     = gzip.open if input_vcf.endswith('.gz') else open
out_opener = gzip.open if output_vcf.endswith('.gz') else open

with opener(input_vcf, 'rt') as infile, out_opener(output_vcf, 'wt') as outfile:
    for line in infile:

        if line.startswith('#'):
            outfile.write(line)
            continue

        fields = line.strip().split('\t')
        if len(fields) < 10:
            continue

        filter_tag = fields[6]
        info_field = fields[7]

        if filter_tag == 'PASS':
            pass_count += 1
            outfile.write(line)
        elif should_rescue(filter_tag, info_field):
            rescued_count += 1
            outfile.write(line)
        else:
            excluded_count += 1

print(f"\nVariant selection complete for {sample_id}:")
print(f"  PASS variants kept:         {pass_count}")
print(f"  Rescued variants:           {rescued_count}")
print(f"  Excluded (failed filters):  {excluded_count}")
print(f"  Total output:               {pass_count + rescued_count}")
print(f"  Output file: {output_vcf}")

PYEOF

echo -ne "*** Variant selection done! ***\n"

echo -ne "*** All done! ***\n"
