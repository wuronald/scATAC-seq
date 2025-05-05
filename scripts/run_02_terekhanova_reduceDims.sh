#!/bin/bash
#SBATCH --job-name=run_02_terekhanova_reduceDims
#SBATCH --output=run_02_terekhanova_reduceDims_%j.out
#SBATCH --error=run_02_terekhanova_reduceDims_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=60G
#SBATCH --time=04:00:00

# Load necessary modules (adjust as needed for your system)
module load R/4.2.1

# Run R script
Rscript - <<EOF

# load libraries
library(ArchR)
library(ggrastr)
library(Seurat)
set.seed(1)

# Set the number of threads for ArchR
addArchRThreads(threads = 18)

# Load the project
proj <- loadArchRProject(path = "Terekhanova")

# IterativeLSI used for dimensionality reduction
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")

# runs seurat's FindClusters function
proj <- addClusters(input = proj, 
                    reducedDims = "IterativeLSI" # name of the reducedDims object
                    )

print("harmony batch correction")

# run harmony batch correction 
proj <- addHarmony(
    ArchRProj = proj,
    reducedDims = "IterativeLSI",
    name = "Harmony",
    groupBy = "Sample" # the variable to correct for
)

proj <- addClusters(
    input = proj,
    reducedDims = "Harmony",
    method = "Seurat",
    name = "Clusters_harmony",
    resolution = 0.8
)

print("UMAP")

# UMAP embeddings with IterativeLSI
proj <- addUMAP(
    ArchRProj = proj, 
    reducedDims = "IterativeLSI", 
    name = "UMAP", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)

# adds harmony to embeddings
proj <- addUMAP(
    ArchRProj = proj, 
    reducedDims = "Harmony", 
    name = "UMAP_harmony", 
    nNeighbors = 30, 
    minDist = 0.5, 
    metric = "cosine"
)

# plot UMAP
print("plot umap")

# plot by sample (IterativeLSI vs Harmony embeddings)
u1 <- plotEmbedding(ArchRProj = proj, 
                    colorBy = "cellColData", 
                    name = "Sample", 
                    embedding = "UMAP"
                    )

u2 <- plotEmbedding(ArchRProj = proj, 
                    colorBy = "cellColData", 
                    name = "Sample", 
                    embedding = "UMAP_harmony"
                    )

# plot by clusters
u3 <- plotEmbedding(ArchRProj = proj, 
                    colorBy = "cellColData", 
                    name = "Clusters", 
                    embedding = "UMAP"
                    )

u4 <- plotEmbedding(ArchRProj = proj, 
                    colorBy = "cellColData", 
                    name = "Clusters_harmony", 
                    embedding = "UMAP_harmony"
                    )

# plot the LSI embedding with harmony cluster labels
u5 <- plotEmbedding(ArchRProj = proj, 
                    colorBy = "cellColData", 
                    name = "Clusters_harmony", 
                    embedding = "UMAP"
                    )
## save PDF of UMAPs

plotPDF(u1,u2,u3,u4,u5,
        name = "Terekhanova-all-Plot-UMAP-Sample-HarmonyClusters.pdf", 
        ArchRProj = proj, 
        addDOC = TRUE, # adds date of creation to end of filename
        width = 5, height = 5)

# Save the project
saveArchRProject(ArchRProj = proj, outputDirectory = "Terekhanova", load = TRUE)

EOF

echo "ArchR analysis completed"
