#!/bin/bash
#SBATCH --account=XXXXX
#SBATCH --partition=default_free
#SBATCH --mem=100G
#SBATCH --time=23:00:00
#SBATCH --cpus-per-task=8
#SBATCH --job-name=gatk_var_calling
#SBATCH --output=logs/gatk_var_calling_%A_%a.out
#SBATCH --array=1-100%100

# Script to run GATK somatic variant calling on RNA-seq data (tumor-only mode)
# Includes rescue-aware variant selection for KBTBD4-type clustered/germline-filtered
# variants that are genuine somatic mutations suppressed by low RNA-seq depth
# Ruth Cranston 2026

[ $# -ne 3 ] && { echo -en \
"\nRuth Cranston 2026\n\n
*** Script to run gatk mutation calling on a list of sample ids from the original fastq file
sample sheet [sample name] [fastq1] [fastq2] (tab delimited sheet).
Runs in current directory. Input dir is location of gatk preprocessed files. Output directory is created.
<sample sheet> <input dir (relative)> <output dir (relative)>
example run: sbatch ./script5_gatk_variant_calling.sh sample_sheet.txt output_preprocessing/ output_mutation_calling/ *** \n\n" ; exit 1; }

# Define variables
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

# Setup
mkdir -p ${OUTPUT_DIR} logs ${TMPDIR}

LINE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" ${SAMPLE_SHEET})
SAMPLE_ID=$(echo $LINE | awk '{print $1}')

echo "Processing sample: ${SAMPLE_ID}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "Assembly: ${ASSEMBLY}"

# Reference file selection
if [ "${ASSEMBLY}" == "GRCh38" ]; then
    REF_FASTA=${REFERENCE_DIR}/Homo_sapiens_assembly38.fasta
    GNOMAD=${REFERENCE_DIR}/af-only-gnomad.hg38.vcf.gz
    PON=${REFERENCE_DIR}/1000g_pon.hg38.vcf.gz
    EXAC=${REFERENCE_DIR}/small_exac_common_3.hg38.vcf.gz
else
    # GRCh37/b37 — all resources use no chr prefix
    REF_FASTA=${REFERENCE_DIR}/Homo_sapiens_assembly19.fasta
    GNOMAD=${REFERENCE_DIR}/af-only-gnomad.raw.sites.vcf
    PON=${REFERENCE_DIR}/Mutect2-WGS-panel-b37.vcf
    EXAC=${REFERENCE_DIR}/small_exac_common_3.vcf
fi

# Step 1: Mutect2 (tumor-only mode)
# --dont-use-soft-clipped-bases: RNA-seq specific — soft clips are splice artefacts
# --f1r2-tar-gz: collects strand orientation data for LearnReadOrientationModel
echo -en "\n--- Running Mutect2 ---\n"
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

echo -ne "*** Mutect2 finished! ***\n"

# Step 2: LearnReadOrientationModel
# Models strand-specific sequencing artefacts from the f1r2 data
# (--ob-priors is intentionally omitted below). The orientation bias model is
# invalid for strand-specific RNA-seq libraries
# LearnReadOrientationModel still runs to completion; its output is retained for
# reference but not applied to filtering.
echo -en "\n--- Running LearnReadOrientationModel ---\n"
gatk LearnReadOrientationModel \
     --java-options "-Xmx90g" \
     -I ${OUTPUT_DIR}${SAMPLE_ID}_f1r2.tar.gz \
     -O ${OUTPUT_DIR}${SAMPLE_ID}_artifact_prior.tar.gz

echo -ne "*** LearnReadOrientationModel finished! ***\n"

# Step 3: Contamination estimation
echo -en "\n--- Running GetPileupSummaries ---\n"
gatk GetPileupSummaries \
    --java-options "-Xmx90g" \
    -I ${INPUT_DIR}${SAMPLE_ID}_recal.bam \
    -V ${EXAC} \
    -L ${EXAC} \
    -O ${OUTPUT_DIR}${SAMPLE_ID}_pileup_summaries.table

echo -en "\n--- Running CalculateContamination ---\n"
gatk CalculateContamination \
    --java-options "-Xmx90g" \
    -I ${OUTPUT_DIR}${SAMPLE_ID}_pileup_summaries.table \
    -O ${OUTPUT_DIR}${SAMPLE_ID}_contamination.table

echo -ne "*** CalculateContamination finished! ***\n"

# Step 4: FilterMutectCalls
# --ob-priors intentionally omitted (see LearnReadOrientationModel note above)
# --contamination-table: applies sample-specific contamination correction
echo -en "\n--- Running FilterMutectCalls ---\n"
gatk FilterMutectCalls \
     --java-options "-Xmx90g" \
     -R ${REF_FASTA} \
     -V ${OUTPUT_DIR}${SAMPLE_ID}_tumor_raw.vcf.gz \
     --contamination-table ${OUTPUT_DIR}${SAMPLE_ID}_contamination.table \
     --min-allele-fraction 0.05 \
     -O ${OUTPUT_DIR}${SAMPLE_ID}_tumor_filtered.vcf.gz

echo -ne "*** FilterMutectCalls finished! ***\n"

# Step 5: Rescue-aware variant selection
# Replaces SelectVariants --exclude-filtered
#
# Background: RNA-seq has substantially lower depth than DNA panel sequencing.
# At sites with genuine somatic mutations, low depth causes two GATK filters to
# fire incorrectly:
#
#   1. clustered_events;haplotype: fires when multiple variants occur within
#      a short window. In KBTBD4 and similar hotspot genes, these are genuine
#      compound somatic mutations, not misassembly artefacts. GATK cannot
#      distinguish this from a true artefact at low depth.
#
#   2. germline: fires when depth is insufficient to statistically distinguish
#      a somatic mutation (~35-50% VAF) from a heterozygous germline variant.
#      At RNA-seq depths of 15-30x, GERMQ collapses to 1 even when population
#      evidence (POPAF) strongly argues against germline status.
#
# Rescue rules (validated against confirmed positive controls with matched DNA):
#   Filter tag                            Rescue condition
#   clustered_events;haplotype         →  TLOD > 10
#   germline                           →  TLOD > 10 AND POPAF > 6
#   clustered_events;germline;haplotype→  TLOD > 10 AND POPAF > 6
#
#   TLOD > 10:  well above GATK's own PASS threshold (3.0) and FAIL threshold (5.3)
#   POPAF > 6:  population AF < 10^-6 — effectively absent from gnomAD,
#               distinguishing genuine somatic from common germline variants
#
# Any variant carrying additional filter tags (contamination, weak_evidence,
# strand_bias etc.) is excluded normally — those represent independent concerns.


echo -en "\n--- Running rescue-aware variant selection ---\n"

# Pass bash variables to Python via environment variables
export RESCUE_INPUT_VCF="${OUTPUT_DIR}${SAMPLE_ID}_tumor_filtered.vcf.gz"
export RESCUE_OUTPUT_VCF="${OUTPUT_DIR}${SAMPLE_ID}_tumor_filtered_PASS.vcf.gz"
export RESCUE_SAMPLE_ID="${SAMPLE_ID}"


python3 << 'PYEOF'
import gzip
import sys
import re
import os

def get_info_value(info_field, key):
    """Extract numeric value from INFO field by key."""
    match = re.search(r'(?:^|;)' + re.escape(key) + r'=([^;]+)', info_field)
    if match:
        try:
            return float(match.group(1))
        except ValueError:
            return None
    return None

def should_rescue(filter_tag, info_field):
    """
    Return True if a filtered variant should be rescued based on validated rules.

    Rescue conditions — all validated against confirmed positive controls
    with matched high-depth DNA panel data:

      clustered_events;haplotype:
        Fires on genuine compound somatic hotspot mutations (e.g. KBTBD4).
        Rescue if TLOD > 10.

      germline (alone, or combined with clustered_events;haplotype):
        Fires at low RNA-seq depth where model cannot distinguish somatic from
        germline at ~35-50% VAF. POPAF gate ensures variant is absent from
        gnomAD before rescue.
        Rescue if TLOD > 10 AND POPAF > 6 (population AF < 10^-6).

    Any variant with additional filter tags beyond these three is NOT rescued.
    """
    # Normalise filter components — sort so order of tags does not matter
    filter_components = set(filter_tag.split(';'))

    # Permitted rescue components
    allowed_components = {'clustered_events', 'germline', 'haplotype'}

    # Reject immediately if any unexpected filter tags are present
    if not filter_components.issubset(allowed_components):
        return False

    # Must have at least one of our target filters
    if not filter_components.intersection({'clustered_events', 'germline', 'haplotype'}):
        return False

    # TLOD gate applies to all rescue cases
    tlod = get_info_value(info_field, 'TLOD')
    if tlod is None or tlod <= 10:
        return False

    # clustered_events;haplotype without germline — TLOD gate is sufficient
    if 'germline' not in filter_components:
        return True

    # germline present — add POPAF gate to distinguish somatic from germline
    popaf = get_info_value(info_field, 'POPAF')
    if popaf is None or popaf <= 6:
        return False

    return True


# Read configuration from environment
input_vcf  = os.environ['RESCUE_INPUT_VCF']
output_vcf = os.environ['RESCUE_OUTPUT_VCF']
sample_id  = os.environ['RESCUE_SAMPLE_ID']

# Counters for reporting
pass_count     = 0
rescued_count  = 0
excluded_count = 0

opener     = gzip.open if input_vcf.endswith('.gz') else open
out_opener = gzip.open if output_vcf.endswith('.gz') else open

with opener(input_vcf, 'rt') as infile, out_opener(output_vcf, 'wt') as outfile:
    for line in infile:

        # Always write header lines unchanged
        if line.startswith('#'):
            outfile.write(line)
            continue

        fields = line.strip().split('\t')
        if len(fields) < 10:
            continue

        filter_tag   = fields[6]
        info_field   = fields[7]

        # Keep PASS variants
        if filter_tag == 'PASS':
            pass_count += 1
            outfile.write(line)
            continue

        # Attempt rescue of specific filtered variants
        if should_rescue(filter_tag, info_field):
            # Promote to PASS so downstream tools (REDIportal filter, VEP)
            # treat rescued variants identically to genuine PASS calls
            fields[6] = 'PASS'
            outfile.write('\t'.join(fields) + '\n')
            rescued_count += 1
            continue

        # All other filtered variants are excluded
        excluded_count += 1

print(f"\nVariant selection complete for {sample_id}:")
print(f"  PASS variants kept:         {pass_count}")
print(f"  Rescued variants:           {rescued_count}")
print(f"    (clustered/germline rescue — confirmed against matched DNA)")
print(f"  Excluded (failed filters):  {excluded_count}")
print(f"  Total output:               {pass_count + rescued_count}")
print(f"  Output file: {output_vcf}")

PYEOF

echo -ne "*** Variant rescue and selection finished! ***\n"
echo -ne "*** All done! ***\n"
