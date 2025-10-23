#!/bin/bash
#SBATCH --job-name=run_07c_human_multiome_projectBulkATACseq
#SBATCH --output=run_07c_human_multiome_projectBulkATACseq_%j.out
#SBATCH --error=run_07c_human_multiome_projectBulkATACseq_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=60G
#SBATCH --time=01:00:00

# Parse command line arguments
SEBULK_PATH="${1:-/cluster/projects/wouterslab/RW555/DiffBind/2023-03-27_seBulk_hg38_all_normoxia_vs_hypoxia_DiffBind.rds}"

# Export the parameter so R can read it
export SEBULK_PATH

echo "Using seBulk path: $SEBULK_PATH"

# Load necessary modules (adjust as needed for your system)
module load R/4.4.1

# Run R script with parameters
Rscript - <<'EOF'

# load libraries
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(here)
set.seed(1)

cat("\n=== Starting Bulk ATAC-seq Projection Analysis ===\n")

# Set the number of threads for ArchR
addArchRThreads(threads = 18)
cat("Set ArchR threads to 18\n")

# set genome to hg38
addArchRGenome("hg38")
cat("Set genome to hg38\n")

# Load the project
cat("\n=== Loading ArchR Project ===\n")
proj <- loadArchRProject(path = "human_multiome_harmony_merged_malig_peak")
cat("Successfully loaded ArchR project\n")
cat("Number of cells:", nCells(proj), "\n")

# Get seBulk path from environment variable
seBulk_path <- Sys.getenv("SEBULK_PATH")
cat("\n=== Loading Bulk ATAC-seq Data ===\n")
cat("Loading from:", seBulk_path, "\n")
seBulk <- readRDS(seBulk_path)
cat("Successfully loaded seBulk object\n")

# Print seBulk information for debugging
cat("\n=== seBulk Object Information ===\n")
cat("seBulk class:", class(seBulk), "\n")
cat("seBulk dimensions:", dim(seBulk), "\n")
cat("seBulk feature type:", class(rowRanges(seBulk)), "\n")
cat("Number of features in seBulk:", length(rowRanges(seBulk)), "\n")

# Get all available reductions and embeddings
cat("\n=== Getting Available Reductions and Embeddings ===\n")
available_reductions <- names(proj@reducedDims)
available_embeddings <- names(proj@embeddings)
cat("Available reductions:", paste(available_reductions, collapse = ", "), "\n")
cat("Available embeddings:", paste(available_embeddings, collapse = ", "), "\n")

# Initialize list to store plots
plot_list <- list()

# Loop over all reductions
for (reduction in available_reductions) {
  cat("\n=== Processing Reduction:", reduction, "===\n")
  
  # Only process Harmony reductions
  if (!is.null(reduction) && !is.na(reduction) && !grepl("^Harmony", reduction)) {
    cat("Skipping", reduction, "- only processing Harmony reductions\n")
    next
  }
  
  # Skip problematic combined reduction if needed
  if (!is.null(reduction) && !is.na(reduction) && reduction == "Harmony_LSI_Combined") {
    cat("Temporarily skipping", reduction, "- known to cause issues\n")
    next
  }
  
  # The corresponding UMAP should be "UMAP_" + reduction name
  umap_name <- paste0("UMAP_", reduction)
  
  # Check if UMAP exists for this reduction
  if (!umap_name %in% available_embeddings) {
    cat("Warning: UMAP embedding", umap_name, "not found. Skipping...\n")
    next
  }
  
  cat("Projecting bulk ATAC-seq using", reduction, "reduction with", umap_name, "embedding\n")
  
  # Add error handling around the projection
  tryCatch({
    # Print reduction dimensions for debugging
    if (reduction %in% names(proj@reducedDims)) {
      cat("Reduced dims dimensions:", dim(proj@reducedDims[[reduction]]), "\n")
    }
    
    # Validate embedding exists and has proper structure
    if (umap_name %in% names(proj@embeddings)) {
      cat("Embedding dimensions:", dim(proj@embeddings[[umap_name]]), "\n")
    } else {
      cat("WARNING: Embedding", umap_name, "not accessible. Skipping...\n")
      next
    }
    
    # Additional validation: check that reduction and embedding are compatible
    cat("Checking compatibility between reduction and embedding...\n")
    reduction_dim <- nrow(proj@reducedDims[[reduction]])
    embedding_dim <- nrow(proj@embeddings[[umap_name]])
    cat("Cells in reduction:", reduction_dim, "\n")
    cat("Cells in embedding:", embedding_dim, "\n")
    
    if (reduction_dim == 0 || embedding_dim == 0) {
      cat("WARNING: Empty reduction or embedding detected. Skipping...\n")
      next
    }
    
    if (reduction_dim != embedding_dim) {
      cat("WARNING: Mismatch between reduction and embedding dimensions. Skipping...\n")
      next
    }
    
    # project bulk ATAC-seq to scATAC-seq
    bulkPro <- projectBulkATAC(
      ArchRProj = proj,
      seATAC = seBulk,
      reducedDims = reduction,
      embedding = umap_name,
      n = 250,
      verbose = TRUE
    )
    
    cat("Projection complete\n")
    cat("Single cell UMAP dimensions:", dim(bulkPro$singleCellUMAP), "\n")
    cat("Simulated bulk UMAP dimensions:", dim(bulkPro$simulatedBulkUMAP), "\n")
    
    # Show first few rows
    cat("First few rows of simulatedBulkUMAP:\n")
    print(head(bulkPro$simulatedBulkUMAP))
    
    cat("\n=== Creating Plot for", reduction, "===\n")
    
    # Concatenate single-cell and pseudocell embedding positions
    pro_df <- rbind(bulkPro$singleCellUMAP, bulkPro$simulatedBulkUMAP)
    
    # Create a color palette and force the scATAC cells to be grey
    pal <- paletteDiscrete(values = unique(as.character(pro_df$Type)), set = "stallion")
    pal["scATAC"] <- "#BABABA"
    
    # Plot using ggPoint
    p <- ggPoint(
      x = pro_df[, 1],
      y = pro_df[, 2],
      discrete = TRUE,
      color = as.character(pro_df$Type),
      pal = pal,
      xlabel = paste0(colnames(pro_df)[1]),
      ylabel = paste0(colnames(pro_df)[2]),
      title = paste0("Bulk ATAC-seq Projection - ", reduction)
    )
    
    # Store plot
    plot_list[[reduction]] <- p
    cat("Plot created for", reduction, "\n")
    
    # Clean up to free memory
    rm(bulkPro, pro_df, p)
    gc()
    
  }, error = function(e) {
    cat("\n!!! ERROR occurred while processing", reduction, "!!!\n")
    cat("Error message:", conditionMessage(e), "\n")
    cat("Skipping this reduction and continuing...\n\n")
  })
}

# Check if we have any plots to save
if (length(plot_list) > 0) {
  # Save all plots to a single PDF using ArchR's plotPDF
  cat("\n=== Saving Plots to PDF ===\n")
  output_name <- "BulkATAC_projection_all_reductions"
  cat("Output file:", paste0(output_name, ".pdf"), "\n")
  cat("Number of plots to save:", length(plot_list), "\n")
  
  plotPDF(
    plotList = plot_list,
    name = output_name,
    ArchRProj = proj,
    addDOC = TRUE,
    width = 10,
    height = 8
  )
  
  cat("\n=== Analysis Complete ===\n")
  cat("PDF saved to:", paste0(getOutputDirectory(proj), "/", output_name, ".pdf"), "\n")
  cat("Successfully created plots for:", paste(names(plot_list), collapse = ", "), "\n")
} else {
  cat("\n=== WARNING: No plots were generated ===\n")
  cat("All reductions failed or were skipped.\n")
}

EOF
