#!/bin/sh
#SBATCH -t 6:00:00
#SBATCH --mem=30G
#SBATCH -J cellranger-atac-count
#SBATCH -p all
#SBATCH -c 12
#SBATCH -N 1
#SBATCH -o %x-%j.out

## uncomment the following line if encountering 'resource unavailable' errors
## despite using --localcores and --localmem
# ulimit -u 4096

# Arguments
# $1 = the id name of the output folder
# $2 = folder where the fastq files are for the given sample


# load cellranger-atac
module load cellranger-atac/2.1.0

cellranger-atac count --id $1 \
  --fastqs $2 \
  --reference=/cluster/tools/data/commondata/cellranger/refdata-cellranger-arc-GRCh38-2020-A-2.0.0 \
  --localcores=$SLURM_CPUS_PER_TASK --localmem=30