#!/bin/bash
#SBATCH --job-name=run_03_human_multiome_harmony
#SBATCH --output=run_03_human_multiome_harmony_%j.out
#SBATCH --error=run_03_human_multiome_harmony_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=60G
#SBATCH --time=01:00:00

# Load necessary modules (adjust as needed for your system)
module load R/4.4.1

# Run R script
Rscript - <<EOF

# load libraries
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Seurat)
library(here)
set.seed(1)

# Set the number of threads for ArchR
addArchRThreads(threads = 18)

# Load the project
proj <- loadArchRProject(path = "human_multiome")

print("Harmony batch correction")

print("Available LSI reducedDims:")
print(proj@reducedDims)
original_reducedDims <- names(proj@reducedDims)
original_embeddings <- names(proj@embeddings)

print("Correcting for groupBy variable: Sample")

# run harmony batch correction
for (i in seq_along(proj@reducedDims)) {
    print(paste("Adding harmony for reducedDims:", names(proj@reducedDims)[i]))
    proj <- addHarmony(
        ArchRProj = proj,
        reducedDims = names(proj@reducedDims)[i],
        name = paste0("Harmony_", names(proj@reducedDims)[i]),
        groupBy = "Sample" # the variable to correct for
    )
    # runs seurat's FindClusters function
    print(paste("Adding clusters for reducedDims:", names(proj@reducedDims)[i]))
    proj <- addClusters(
        input = proj,
        reducedDims = paste0("Harmony_", names(proj@reducedDims)[i]),
        method = "Seurat",
        name = paste0("Clusters_harmony_", names(proj@reducedDims)[i]),
        resolution = 0.8
    )
} 

print("UMAP")

# UMAP embeddings with each of the reducedDims, including harmony embeddings
for (i in seq_along(proj@reducedDims)) {
    print(paste("Adding UMAP for reducedDims:", names(proj@reducedDims)[i]))
    proj <- addUMAP(
        ArchRProj = proj, 
        reducedDims = names(proj@reducedDims)[i], 
        name = paste0("UMAP_", names(proj@reducedDims)[i]), 
        nNeighbors = 30, 
        minDist = 0.5, 
        metric = "cosine"
    )
}


# plot UMAP

new_reducedDims <- setdiff(names(proj@reducedDims), original_reducedDims)

print("Initiate UMAP plotting")
print("Previous available UMAP embeddings:")
print(original_embeddings)

print("plot umap for new embeddings only:")
print("Newly added Harmony UMAP embeddings are:")
new_embeddings <- setdiff(names(proj@embeddings), original_embeddings)
print(new_embeddings)

# plot by sample (IterativeLSI vs Harmony embeddings)
for (i in seq_along(new_embeddings)) {
    print(paste("Plotting UMAP for embedding:", names(new_embeddings)[i]))
    u <- plotEmbedding(ArchRProj = proj, 
                       colorBy = "cellColData", 
                       name = "Sample", 
                       embedding = names(new_embeddings)[i]
                       )
    # save each plot
    plotPDF(u, 
            name = paste0("Plot-UMAP-Sample-", names(new_embeddings)[i], ".pdf"), 
            ArchRProj = proj, 
            addDOC = TRUE, # adds date of creation to end of filename
            width = 5, height = 5)
}

# plot by harmony clusters for each reducedDims
for (i in seq_along(new_embeddings)) {
    print(paste("Plotting UMAP for clusters_harmony:", names(new_embeddings)[i]))
    u <- plotEmbedding(ArchRProj = proj, 
                       colorBy = "cellColData", 
                       name = paste0("Clusters_harmony_", names(new_reducedDims)[i]),
                       embedding = names(new_embeddings)[i]
                       )
    # save each plot
    plotPDF(u, 
            name = paste0("Harmony-Plot-UMAP-Clusters_harmony-", names(new_embeddings)[i], ".pdf"), 
            ArchRProj = proj, 
            addDOC = TRUE, # adds date of creation to end of filename
            width = 5, height = 5)
}

# Save the project
print("Saving Project")
# saveArchRProject(ArchRProj = proj, outputDirectory = "human_multiome", load = TRUE)
EOF

echo "ArchR analysis completed"
