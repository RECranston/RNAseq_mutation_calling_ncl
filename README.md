# RNAseq_mutation_calling_ncl

Here are seven sequential bash scripts to run an end-to-end tumour-only analysis of RNA-sequencing mutation detection using STAR, [gatk](https://github.com/broadinstitute/gatk) and [vep](https://github.com/ensembl/ensembl-vep) on the Newcastle University server (Comet) with slurm scheduler.

The scripts include:  
* `script1_index_build.sh`:
    * Downloads reference files for GRCh37 and GRCh38. Including files for use with gatk, BQSR, mutect2 and builds indexes/prepares reference files where needed
    * Creates vep caches
    * Builds STAR indexes
* `script2_trimming.sh`:
    * Trims sequencing adapters from paired fastq.gz files using Trim Galore
    * Runs FastQC
* `script3_align_twopass.sh`:
    * Performs STAR two-pass alignment on a paired list of trimmed fastq.gz files
* `script4_gatk_preprocessing.sh`:
    * Performs gatk preprocessing on named files (names derived from paired fastq.gz file sample sheet)
    * Preprocessing includes: Mark duplicates, Split reads at splice junctions (SplitNCigarReads), Base Quality Score Recalibration (BQSR) model building and application
* `script5_gatk_variant_calling.sh`:
    * Performs mutation detection on named files (names derived from paired fastq.gz file sample sheet) using Mutect2 in tumour-only mode, with comparison to the reference genome. Includes `--dont-use-soft-clipped-bases` which is important for excluding alignment noise in RNA-sequencing data.
    * Considers strand/orientation bias, sample contamination and filters variants accordingly. Also filters for PASS mutations (gatk SelectVariants)
* `script6_rna_editing_filter.sh`:
    * Excludes RNA variants detected in RNA editing sites
* `script7_vep.sh`:
    * Annotation of filtered variants by vep.
  
### Setup

* Create a new directory and move into it.
* Git clone this repository.
```
git clone https://github.com/RECranston/RNAseq_mutation_calling_ncl.git
```
* Change into the cloned directory `cd RNAseq_mutation_calling_ncl`. Make all shell scripts executable
```
chmod 777 *.sh
```
* Please edit the script header of all scripts to assign the correct account name to the sbatch run.
* Run the setup script. References, indexes and required files will be downloaded and built as required for all analysis stages.
```
sbatch ./script1_index_build.sh
```
* Move all fastq.gz files for analysis (or symlink using `ln -s`) to a sub-directory within the current directory.
* Create a tab-separated sample sheet of paired fastq files and sample identifiers and save as a `.txt` file. Sample sheet should be three columns including: sample identifier followed by paired files as per the example below:
```
sample_name  sample1_L001_R1.fastq.gz  sample1_L001_R2.fastq.gz
```

### Run the pipeline
Parameters required for each script can be checked by running `./script_name.sh` in the terminal:
* `<tab delimited sample sheet>` is the name of the sample sheet (`.txt` file) of fastq.gz file locations and associated sample names created during the setup stage.
* `<input dir (relative)>` is the location of the directory containing input files relative to the current directory e.g. `trimmed_fastq/` (note the trailing “/”). 
* Similarly `<output dir (relative)>` is the location of the directory where the output data is to be stored, relative to the current directory e.g. `vep_output/` (note the trailing “/”).
* Please edit the script header to assign the correct account name to the sbatch run and the correct number of jobs in the array, and jobs to be simultaneously performed.  
  E.g. this example runs samples 1-10 from the sample sheet, running two samples at a time.
```
#SBATCH --array=1-10%2
```
This can usually be defined by the number of rows in the sample sheet e.g. `cat sample_sheet.txt | wc -l`
* Run scripts in sequential order.
* After `script2_trimming.sh` create a new tab delimited sample sheet including paired trimmed fastq files and sample identifiers and save as a `.txt` file.
* Sample sheet should be three columns including: sample identifier followed by paired files as per the example below:
```
sample_name  sample1_L001_R1_val_2.fastq.gz sample1_L001_R2_val_2.fastq.gz
```
* Resulting filtered, vep-annotated mutation data is saved to the defined output directory as `.vcf` files.
