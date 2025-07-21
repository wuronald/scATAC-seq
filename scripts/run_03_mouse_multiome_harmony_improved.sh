#!/bin/bash
#SBATCH --job-name=run_03_mouse_multiome_harmony_improved
#SBATCH --output=run_03_mouse_multiome_harmony_improved_%j.out
#SBATCH --error=run_03_mouse_multiome_harmony_improved_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=60G
#SBATCH --time=01:30:00

# Load necessary modules (adjust as needed for your system)
module load R/4.4.1

# Run R script
Rscript - <<'EOF'

# load libraries
library(ArchR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Seurat)
library(here)
set.seed(1)

# Set the number of threads for ArchR
addArchRThreads(threads = 18)

# Load the project
print("Loading ArchR project")
proj <- loadArchRProject(path = "mouse_multiome")

# Check available reducedDims and embeddings:
print("Available Original LSI reducedDims:")
original_reducedDims <- names(proj@reducedDims)
original_embeddings <- names(proj@embeddings)
print(paste0("Available original reducedDims:", original_reducedDims))
print(paste0("Available original embeddings:", original_embeddings))

# add libraries as batch to cellColData:
cell_samples <- getCellColData(proj, select = "Sample", drop = TRUE)
batch_vector <- ifelse(cell_samples == "SM222_GL1", "batch1",
                 ifelse(cell_samples %in% c("SM122_AR1", "SM122_BR4", "SM222_AR1"), "batch2",
                 ifelse(cell_samples %in% c("SM222_HL1", "SM922_AR4"), "batch3", 
                ifelse(cell_samples == "SM922_BR4", "batch4", NA))))

proj <- addCellColData(
  ArchRProj = proj,
  data = batch_vector,
  name = "lib_batch",
  cells = getCellNames(proj),
  force = TRUE
)

print("Harmony batch correction with parameter testing")
# cellColData variable to correct for batch effects using Harmony:
groupByVar <- c("Sample","lib_batch")

# Check if the variable exists in cellColData
if (!all(groupByVar %in% names(proj@cellColData))) {
    stop(paste("The variable:", groupByVar, "does not exist in cellColData. Please check the variable name."))
} else {
    print(paste("The variable:", groupByVar, "exists in cellColData. Proceeding with Harmony batch correction."))
}       
print(paste0("Correcting for groupBy variable:", groupByVar))

# Define parameter sets to test
# Each set is designed to reduce clusters and make them more cohesive
# theta and lambda need to be vectors with one value per groupBy variable (Sample, lib_batch)
param_sets <- list(
  default = list(sigma = 0.1, theta = NULL, lambda = NULL, max.iter.harmony = 10, resolution = 0.4, suffix = "default"),
  higher_sigma = list(sigma = 0.3, theta = NULL, lambda = NULL, max.iter.harmony = 10, resolution = 0.2, suffix = "sigma0p3_res0p2"),
  strong_correction = list(sigma = 0.2, theta = c(2, 2), lambda = c(1, 1), max.iter.harmony = 20, resolution = 0.15, suffix = "theta2_lambda1_res0p15"),
  very_cohesive = list(sigma = 0.5, theta = c(3, 3), lambda = c(1.5, 1.5), max.iter.harmony = 25, resolution = 0.1, suffix = "sigma0p5_theta3_lambda1p5_res0p1"),
  moderate_cohesive = list(sigma = 0.25, theta = c(1.5, 1.5), lambda = c(0.5, 0.5), max.iter.harmony = 15, resolution = 0.25, suffix = "sigma0p25_theta1p5_lambda0p5_res0p25")
)

# Store all new reducedDims and embeddings for plotting
all_new_reducedDims <- c()
all_new_embeddings <- c()

# run harmony batch correction with different parameter sets
for (param_name in names(param_sets)) {
    params <- param_sets[[param_name]]
    print(paste("Testing parameter set:", param_name))
    print(paste("Parameters: sigma =", params$sigma, ", theta =", params$theta, 
                ", lambda =", params$lambda, ", max.iter.harmony =", params$max.iter.harmony,
                ", resolution =", params$resolution))
    
    for (i in seq_along(proj@reducedDims)) {
        if (names(proj@reducedDims)[i] %in% original_reducedDims) {  # Only process original reducedDims
            reducedDim_name <- names(proj@reducedDims)[i]
            harmony_name <- paste0("Harmony_", reducedDim_name, "_", params$suffix)
            cluster_name <- paste0("Clusters_harmony_", reducedDim_name, "_", params$suffix)
            
            print(paste("Adding harmony for reducedDims:", reducedDim_name, "with suffix:", params$suffix))
            
            # Build harmony parameters dynamically
            harmony_args <- list(
                ArchRProj = proj,
                reducedDims = reducedDim_name,
                name = harmony_name,
                groupBy = groupByVar,
                sigma = params$sigma,
                max.iter.harmony = params$max.iter.harmony
            )
            
            # Add theta and lambda only if they're not NULL
            if (!is.null(params$theta)) {
                harmony_args$theta <- params$theta
            }
            if (!is.null(params$lambda)) {
                harmony_args$lambda <- params$lambda
            }
            
            # Run addHarmony with the specified parameters
            proj <- do.call(addHarmony, harmony_args)
            
            # runs seurat's FindClusters function with specified resolution
            print(paste("Adding clusters for reducedDims:", reducedDim_name, "with resolution:", params$resolution))
            proj <- addClusters(
                input = proj,
                reducedDims = harmony_name,
                method = "Seurat",
                name = cluster_name,
                resolution = params$resolution
            )
            
            # Track new reducedDims
            all_new_reducedDims <- c(all_new_reducedDims, harmony_name)
        }
    }
}

print("Available LSI reducedDims after harmony:")
print(names(proj@reducedDims))

print("UMAP")

# UMAP embeddings with each of the new harmony reducedDims
for (reducedDim in all_new_reducedDims) {
    umap_name <- paste0("UMAP_", reducedDim)
    print(paste("Adding UMAP for reducedDims:", reducedDim))
    proj <- addUMAP(
        ArchRProj = proj, 
        reducedDims = reducedDim, 
        name = umap_name, 
        nNeighbors = 40, 
        minDist = 0.5, 
        metric = "cosine"
    )
    all_new_embeddings <- c(all_new_embeddings, umap_name)
}

print("Initiate UMAP plotting")
print("Newly added Harmony UMAP embeddings:")
print(all_new_embeddings)

# plot by sample for all parameter sets
for (embedding in all_new_embeddings) {
    print(paste("Plotting UMAP for embedding:", embedding, "with Sample as colorBy"))
    u <- plotEmbedding(ArchRProj = proj, 
                       colorBy = "cellColData", 
                       name = "Sample", 
                       embedding = embedding
                       )
    # save each plot with parameter info in filename
    plotPDF(u, 
            name = paste0("Plot-UMAP-Sample-", embedding, "-colorBy-Sample.pdf"), 
            ArchRProj = proj, 
            addDOC = TRUE,
            width = 5, height = 5)
}

# print cellColData to see clusters
print("Available cellColData:")
cluster_columns <- names(proj@cellColData)[grepl("Clusters_harmony", names(proj@cellColData))]
print("Available cluster columns:")
print(cluster_columns)

# plot by harmony clusters for each parameter set
for (embedding in all_new_embeddings) {
    # Extract the parameter suffix from the embedding name
    embedding_parts <- strsplit(embedding, "_")[[1]]
    if (length(embedding_parts) >= 4) {
        # Find corresponding cluster column
        reducedDim_base <- paste(embedding_parts[3:(length(embedding_parts)-1)], collapse = "_")
        param_suffix <- embedding_parts[length(embedding_parts)]
        cluster_col <- paste0("Clusters_harmony_", reducedDim_base, "_", param_suffix)
        
        if (cluster_col %in% names(proj@cellColData)) {
            print(paste("Plotting UMAP for embedding:", embedding, "with clusters:", cluster_col))
            u <- plotEmbedding(ArchRProj = proj, 
                               colorBy = "cellColData", 
                               name = cluster_col,
                               embedding = embedding
                               )
            # save each plot
            plotPDF(u, 
                    name = paste0("Harmony-Plot-UMAP-", embedding, "-colorBy-", cluster_col, ".pdf"), 
                    ArchRProj = proj, 
                    addDOC = TRUE,
                    width = 5, height = 5)
        }
    }
}

# Generate summary of cluster numbers for each parameter set
print("=== CLUSTER SUMMARY ===")
for (cluster_col in cluster_columns) {
    cluster_data <- getCellColData(proj, select = cluster_col, drop = TRUE)
    n_clusters <- length(unique(cluster_data))
    print(paste("Parameter set", cluster_col, "produced", n_clusters, "clusters"))
}

# Save the project
# print("Saving Project")
# saveArchRProject(ArchRProj = proj, outputDirectory = "mouse_multiome_harmony_test", load = FALSE)

EOF

echo "ArchR Harmony parameter testing completed"