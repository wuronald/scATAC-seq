#!/bin/bash
#SBATCH --job-name=homer_motifs_v5
#SBATCH --output=homer_motifs_v5_%j.out
#SBATCH --error=homer_motifs_v5_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=58G
#SBATCH --time=18:00:00


#Identify Motifs from ATAC-seq data via Homer

# ** New in Version v4 **
## added -FDR 100 to compute 100 randomizations and outputs the qvalue/FDR value
## added -p 12 to utilize 12 processors
## set to himem to prevent oom

# new in version 5.1

module load homer/5.1
module load weblogo/2.8.2

#1=full path to treatment bed/peak file ie. /cluster/projects/wouterslab/RW492/Genrich/bed/EXCLUSIVE/RW471_H_24_EXCLUSIVE_peaks.bed
#2=full path to output directory ie /cluster/projects/wouterslab/RW492/Genrich/HOMER/RW492-1
#3=full path to preparsed directory ie /cluster/projects/wouterslab/RW221/diffpeaks/motifs_preparsed
#4=full path to background bed/peak file ie. /cluster/projects/wouterslab/RW492/Genrich/bed/catalogue/RW471_RW473_RW478_RW486_normoxia_catalogue_peaks.sort.merged.bed

mkdir -p $2;
mkdir -p $3;

#-S 15: number of motifs to find, default is 25; -p 12 processors
# findMotifsGenome.pl $1 hg19 $2 -preparsedDir $3 -S 20 -p 12

# findMotifsGenome.pl $1 hg19 $2 -mask -p 12 -fdr 100 -size 368 -preparsedDir $3 -bg $4

findMotifsGenome.pl $1 hg19 $2 -mask -p 12 -fdr 100 -size 368 -preparsedDir $3 -bg $4


