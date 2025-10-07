#!/bin/bash
#SBATCH --job-name=run_07c_human_multiome_projectBulkATACseq
#SBATCH --output=run_07c_human_multiome_projectBulkATACseq_%j.out
#SBATCH --error=run_07c_human_multiome_projectBulkATACseq_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=16G
#SBATCH --time=00:30:00

# Parse command line arguments

# Load necessary modules (adjust as needed for your system)
module load R/4.4.1

# Run R script with parameters
Rscript - <<'EOF'

# load libraries
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(here)
set.seed(1)

# Set the number of threads for ArchR
addArchRThreads(threads = 18)

# set genome to hg38
addArchRGenome("hg38")

# Load the project
proj <- loadArchRProject(path = "human_multiome_harmony_merged_malig_peak")

# Load Bulk ATAC-seq SummarizedExperiment
seBulk <- loadRDS("/cluster/projects/wouterslab/RW555/DiffBind/2023-03-27_seBulk_hg38_all_normoxia_vs_hypoxia_DiffBind.rds")

# project bulk ATAC-seq to scATAC-seq
bulkPro <- projectBulkATAC(ArchRProj = proj,
                           seATAC = seBulk,
                           reducedDims = "IterativeLSI",
                           embedding = "UMAP",
                           n = 250
                           )

head(bulkPro$simulatedBulkUMAP)
head(bulkPro$singleCellUMAP)

#concatenate single-cell and pseudocell embedding positions
pro_df <- rbind(bulkPro$singleCellUMAP, bulkPro$simulatedBulkUMAP)

#create a color palette and force the scATAC cells to be grey to enable visualization of the project bulk ATAC data
pal <- paletteDiscrete(values = unique(as.character(pro_df$Type)), set = "stallion")
pal["scATAC"] <- "#BABABA"

#plot using ggPoint
p <- ggPoint(x = pro_df$UMAP1,
        y = pro_df$UMAP2,
        discrete = TRUE,
        color = as.character(pro_df$Type),
        pal = pal,
        xlabel = "UMAP Dimension 1",
        ylabel = "UMAP Dimension 2",
        title = "Bulk ATAC-seq Projection")

EOF