#!/bin/sh
#SBATCH -t 6:00:00
#SBATCH --mem=30G
#SBATCH -J cellranger-atac-count
#SBATCH -p all
#SBATCH -c 12
#SBATCH -N 1
#SBATCH -o %x-%j.out
# ==============================================================================
# Purpose: SLURM batch submission script to run 10x Genomics 'cellranger-atac count'
#          for single-cell ATAC-seq FASTQ processing (alignment, barcode counting,
#          peak calling, and fragment file generation against GRCh38 reference).
# Referenced in: Guilhamon_preprocessing.Rmd
#
# Positional Arguments:
#   $1 = Sample / output folder ID (e.g., GBM_4218)
#   $2 = Path to directory containing sample FASTQ files (e.g., fastq/GBM_4218)
# ==============================================================================

## uncomment the following line if encountering 'resource unavailable' errors
## despite using --localcores and --localmem
# ulimit -u 4096


# load cellranger-atac
module load cellranger-atac/2.1.0

cellranger-atac count --id $1 \
  --fastqs $2 \
  --reference=/cluster/tools/data/commondata/cellranger/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
  --localcores=$SLURM_CPUS_PER_TASK --localmem=30