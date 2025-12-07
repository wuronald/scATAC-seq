#!/bin/bash
#SBATCH --job-name=run_03_Gaiti_multiome_harmony
#SBATCH --output=run_03_Gaiti_multiome_harmony_%j.out
#SBATCH --error=run_03_Gaiti_multiome_harmony_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=60G
#SBATCH --time=04:00:00

# Load necessary modules (adjust as needed for your system)
module load R/4.4.1

# Run R scriptmv
Rscript - <<'EOF'

# load libraries
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Seurat)
library(here)
set.seed(1)

# Set the number of threads for ArchR
addArchRThreads(threads = 18)

# Load the project
print("Loading ArchR project")
proj <- loadArchRProject(path = "Gaiti_multiome")

# Check available reducedDims and embeddings:
print("Available Original LSI reducedDims:")
original_reducedDims <- names(proj@reducedDims)
original_embeddings <- names(proj@embeddings)
print(paste0("Available original reducedDims:", original_reducedDims))
print(paste0("Available original embeddings:", original_embeddings))

# add patient IDs as batch to cellColData:
# patient IDs: 6419, 6425, 6434, 6467, 6509, 6514 
cell_samples <- getCellColData(proj, select = "Sample", drop = TRUE)
batch_vector <- ifelse(grepl("6419", cell_samples), "batch1",
                ifelse(grepl("6425", cell_samples), "batch2",
                ifelse(grepl("6434", cell_samples), "batch3",
                ifelse(grepl("6467", cell_samples), "batch4",
                ifelse(grepl("6509", cell_samples), "batch5",   
                ifelse(grepl("6514", cell_samples), "batch6", NA))))))

proj <- addCellColData(
  ArchRProj = proj,
  data = batch_vector,
  name = "lib_batch",
  cells = getCellNames(proj),
  force = TRUE
)

print("Harmony batch correction")
# cellColData variable to correct for batch effects using Harmony:
groupByVar <- c("Sample","lib_batch")

# Check if the variable exists in cellColData
if (!all(groupByVar %in% names(proj@cellColData))) {
    stop(paste("The variable:", groupByVar, "does not exist in cellColData. Please check the variable name."))
} else {
    print(paste("The variable:", groupByVar, "exists in cellColData. Proceeding with Harmony batch correction."))
}       
print(paste0("Correcting for groupBy variable:", groupByVar))

# run harmony batch correction
for (i in seq_along(proj@reducedDims)) {
    print(paste("Adding harmony for reducedDims:", names(proj@reducedDims)[i]))
    proj <- addHarmony(
        ArchRProj = proj,
        reducedDims = names(proj@reducedDims)[i],
        name = paste0("Harmony_", names(proj@reducedDims)[i]),
        groupBy = groupByVar # the variable to correct for
    )
    # runs seurat's FindClusters function
    print(paste("Adding clusters for reducedDims:", names(proj@reducedDims)[i]))
    proj <- addClusters(
        input = proj,
        reducedDims = paste0("Harmony_", names(proj@reducedDims)[i]),
        method = "Seurat",
        name = paste0("Clusters_harmony_", names(proj@reducedDims)[i]),
        resolution = 0.4
    )
} 
print("Available LSI reducedDims after harmony:")
print(proj@reducedDims)
#original_reducedDims <- names(proj@reducedDims)

print("UMAP")

# UMAP embeddings with each of the reducedDims, including harmony embeddings
for (i in seq_along(proj@reducedDims)) {
    print(paste("Adding UMAP for reducedDims:", names(proj@reducedDims)[i]))
    proj <- addUMAP(
        ArchRProj = proj, 
        reducedDims = names(proj@reducedDims)[i], 
        name = paste0("UMAP_", names(proj@reducedDims)[i]), 
        nNeighbors = 40, 
        minDist = 0.5, 
        metric = "cosine"
    )
}


# plot UMAP

new_reducedDims <- setdiff(names(proj@reducedDims), original_reducedDims)
print("Newly added reducedDims:")
print(new_reducedDims)

print("Initiate UMAP plotting")
print("Previous available UMAP embeddings:")
print(original_embeddings)

print("plot umap for new embeddings only:")
print("Newly added Harmony UMAP embeddings are:")
new_embeddings <- setdiff(names(proj@embeddings), original_embeddings)
print(new_embeddings)

# plot by sample (IterativeLSI vs Harmony embeddings)
for (i in seq_along(new_embeddings)) {
    print(paste("Plotting UMAP for embedding:", new_embeddings[i], "with Sample as colorBy"))
    u <- plotEmbedding(ArchRProj = proj, 
                       colorBy = "cellColData", 
                       name = "Sample", 
                       embedding = new_embeddings[i]
                       )
    # save each plot
    plotPDF(u, 
            name = paste0("Plot-UMAP-Sample-", new_embeddings[i], ".pdf"), 
            ArchRProj = proj, 
            addDOC = TRUE, # adds date of creation to end of filename
            width = 5, height = 5)
}

# print cellColData to see clusters
print("Available cellColData:")
print(names(proj@cellColData))

# plot by harmony clusters for each reducedDims
for (i in seq_along(new_embeddings)) {
    for (j in seq_along(original_reducedDims)) {
    print(paste("Plotting UMAP for embedding:", new_embeddings[i], "with Clusters_harmony as colorBy ", original_reducedDims[j]))
    u <- plotEmbedding(ArchRProj = proj, 
                       colorBy = "cellColData", 
                       name = paste0("Clusters_harmony_", original_reducedDims[j]),
                       embedding = new_embeddings[i]
                       )
    # save each plot
    plotPDF(u, 
            name = paste0("Harmony-Plot-UMAP-Clusters_harmony-", new_embeddings[i],"-", original_reducedDims[j],".pdf"), 
            ArchRProj = proj, 
            addDOC = TRUE, # adds date of creation to end of filename
            width = 5, height = 5)
            }
}

# Save the project
print("Saving Project")
saveArchRProject(ArchRProj = proj, outputDirectory = "Gaiti_multiome_harmony", load = TRUE)
EOF

echo "ArchR analysis completed"
