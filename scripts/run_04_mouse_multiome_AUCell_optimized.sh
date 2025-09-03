#!/bin/bash
#SBATCH --job-name=run_04_mouse_multiome_AUCell_optimized
#SBATCH --output=run_04_mouse_multiome_AUCell_optimized_%j.out
#SBATCH --error=run_04_mouse_multiome_AUCell_optimized_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=60G
#SBATCH --time=02:00:00

# Set error handling
set -euo pipefail

# Configuration variables
PROJECT_DIR="mouse_multiome_harmony_merged_subset"
OUTPUT_DIR="AUCell_results"
THREADS=18
SEED=42
CURRENT_DATE=$(date +%Y-%m-%d)

# Load necessary modules (adjust as needed for your system)
module load R/4.4.1

echo "Starting AUCell analysis at $(date)"
echo "Using ${THREADS} threads"

# Run R script with improved error handling and logging
Rscript - <<'EOF'

# Load libraries with error handling
required_packages <- c("here", "ArchR", "AUCell", "ggrastr", "readr", "dplyr", "tools")
missing_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]

if(length(missing_packages) > 0) {
  stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
}


library(ArchR)
library(AUCell)
library(ggrastr)
library(readr)
library(dplyr)
library(tools)
library(here)

# Configuration
PROJECT_DIR <- "mouse_multiome_harmony_merged_subset"
OUTPUT_DIR <- "AUCell_results"
THREADS <- 18
SEED <- 42
CURRENT_DATE <- Sys.Date()

set.seed(SEED)
addArchRThreads(threads = THREADS)

cat("Starting AUCell analysis at", as.character(Sys.time()), "\n")

# Create output directories
output_path <- file.path(PROJECT_DIR, OUTPUT_DIR)
plots_path <- file.path(PROJECT_DIR, "Plots")
dir.create(output_path, showWarnings = FALSE, recursive = TRUE)
dir.create(plots_path, showWarnings = FALSE, recursive = TRUE)

# Load the project with error checking
cat("Loading ArchR project...\n")
if (!dir.exists(PROJECT_DIR)) {
  stop("Project directory not found: ", PROJECT_DIR)
}

proj <- loadArchRProject(path = PROJECT_DIR)
cat("Project loaded successfully. Number of cells:", nCells(proj), "\n")

# Function to load or create expression matrix
load_or_create_expr_matrix <- function(proj, output_path) {
  expr_file <- file.path(output_path, "mouse_multiome_ExprMatrix.RData")
  
  if (file.exists(expr_file)) {
    cat("Loading previously saved expression matrix...\n")
    load(expr_file, envir = .GlobalEnv)
    return(exprMatrix)
  } else {
    cat("Extracting expression matrix from project...\n")
    exprMatrix <- getMatrixFromProject(proj, useMatrix = "GeneExpressionMatrix")
    
    # Add gene names to rownames
    if (is.null(rownames(exprMatrix))) {
      rownames(exprMatrix) <- SummarizedExperiment::rowData(exprMatrix)$name
      cat("Added gene names to rownames\n")
    }
    
    cat("Saving expression matrix...\n")
    save(exprMatrix, file = expr_file)
    return(exprMatrix)
  }
}

# Load expression matrix
exprMatrix <- load_or_create_expr_matrix(proj, output_path)

# Validate expression matrix
if (!is(exprMatrix, "SummarizedExperiment")) {
  stop("Expression matrix is not a SummarizedExperiment object")
}

cat("Expression matrix dimensions:", nrow(exprMatrix), "x", ncol(exprMatrix), "\n")

# Function to load gene sets
load_gene_sets <- function() {
  cat("Loading gene sets...\n")
  
  # Load TSV gene sets
  signatures_dir <- here("data", "Signatures")
  if (!dir.exists(signatures_dir)) {
    stop("Signatures directory not found: ", signatures_dir)
  }
  
  txt_files <- list.files(path = signatures_dir, pattern = "\\.tsv$", full.names = TRUE)
  
  if (length(txt_files) == 0) {
    warning("No .tsv files found in signatures directory")
    geneSets <- list()
  } else {
    geneSets <- lapply(txt_files, readr::read_lines)
    names(geneSets) <- tools::file_path_sans_ext(basename(txt_files))
  }
  
  # Load mouse PIMO signature
  pimo_file <- here("data", "Signatures", "pimo_sig", "phoebe_mouse_pimo_sigs", "mouse_pimo_sig.csv")
  
  if (file.exists(pimo_file)) {
    PIMO_pos <- readr::read_csv(pimo_file, show_col_types = FALSE) %>%
      dplyr::filter(logfoldchanges >= 0) %>%
      dplyr::pull(names)
    
    geneSets <- append(geneSets, list(PIMO_pos = PIMO_pos))
    cat("Added PIMO signature with", length(PIMO_pos), "genes\n")
  } else {
    warning("PIMO signature file not found: ", pimo_file)
  }
  
  cat("Loaded", length(geneSets), "gene sets\n")
  return(geneSets)
}

# Function to calculate gene set overlap
calculate_overlap <- function(se, geneSets) {
  se_genes <- rownames(se)
  results <- data.frame(
    Gene_Set = character(),
    Set_Size = integer(),
    Overlap = integer(),
    Percentage_OL = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (set_name in names(geneSets)) {
    gene_set <- geneSets[[set_name]]
    overlap <- length(intersect(se_genes, gene_set))
    set_size <- length(gene_set)
    percentage <- (overlap / set_size) * 100
    
    results <- rbind(results, data.frame(
      Gene_Set = set_name,
      Set_Size = set_size,
      Overlap = overlap,
      Percentage_OL = percentage,
      stringsAsFactors = FALSE
    ))
  }
  return(results)
}

# Load gene sets and check overlap
geneSets <- load_gene_sets()
overlap_results <- calculate_overlap(exprMatrix, geneSets)
print(overlap_results)

# Save overlap results
write.csv(overlap_results, 
          file = file.path(output_path, paste0(CURRENT_DATE, "_gene_set_overlap.csv")), 
          row.names = FALSE)

# Function to calculate AUCell scores
calculate_aucell_scores <- function(exprMatrix, geneSets, output_path) {
  aucell_file <- file.path(output_path, "cells_AUC.RData")
  
  if (file.exists(aucell_file)) {
    cat("Loading previously calculated AUCell scores...\n")
    load(aucell_file, envir = .GlobalEnv)
    return(cells_AUC)
  } else {
    cat("Calculating cell rankings and AUC scores...\n")
    
    # Calculate cell rankings
    cell_rankings <- AUCell_buildRankings(exprMatrix, plotStats = FALSE)
    
    # Calculate AUC scores
    cells_AUC <- AUCell_calcAUC(geneSets, cell_rankings)
    
    # Save results
    save(cells_AUC, file = aucell_file)
    cat("AUCell scores saved to:", aucell_file, "\n")
    
    return(cells_AUC)
  }
}

# Calculate AUCell scores
cells_AUC <- calculate_aucell_scores(exprMatrix, geneSets, output_path)

# Explore thresholds and create assignments
cat("Exploring AUCell thresholds...\n")
set.seed(SEED)

cells_assignment <- AUCell_exploreThresholds(cells_AUC, 
                                             plotHist = TRUE,
                                             assign = TRUE)

# Check for warnings
warningMsg <- sapply(cells_assignment, function(x) x$aucThr$comment)
warnings_present <- warningMsg[which(warningMsg != "")]
if (length(warnings_present) > 0) {
  cat("Threshold warnings:\n")
  print(warnings_present)
}

# Save threshold plots
pdf_file <- file.path(plots_path, paste0(CURRENT_DATE, "_AUCell_thresholds.pdf"))
pdf(pdf_file, width = 12, height = 9)
  par(mfrow = c(3, 3))
  AUCell_exploreThresholds(cells_AUC, plotHist = TRUE, assign = TRUE)
dev.off()
cat("Threshold plots saved to:", pdf_file, "\n")

# Add AUCell assignments to ArchR project
cat("Adding AUCell assignments to ArchR project...\n")

# Add binary assignments
for (set_name in names(cells_assignment)) {
  assigned_cells <- cells_assignment[[set_name]]$assignment
  proj@cellColData[[set_name]] <- getCellNames(proj) %in% assigned_cells
  
  cat("Added binary assignment for", set_name, ":", 
      sum(proj@cellColData[[set_name]]), "cells assigned\n")
}

# Add continuous AUC scores
cat("Adding continuous AUC scores to ArchR project...\n")
auc_matrix <- assay(cells_AUC)

for (set_name in rownames(auc_matrix)) {
  score_name <- paste0("AUCell_", set_name)
  proj <- addCellColData(
    ArchRProj = proj,
    data = auc_matrix[set_name, ],
    name = score_name,
    cells = colnames(auc_matrix)
  )
  cat("Added continuous scores for", score_name, "\n")
}

# Save summary statistics
summary_stats <- data.frame(
  Gene_Set = names(cells_assignment),
  Threshold = sapply(cells_assignment, function(x) x$aucThr$selected),
  Assigned_Cells = sapply(cells_assignment, function(x) length(x$assignment)),
  Total_Cells = nCells(proj),
  Percentage_Assigned = sapply(cells_assignment, function(x) 
    round(length(x$assignment) / nCells(proj) * 100, 2))
)

write.csv(summary_stats, 
          file = file.path(output_path, paste0(CURRENT_DATE, "_aucell_summary.csv")), 
          row.names = FALSE)

cat("Summary statistics:\n")
print(summary_stats)

# Save the updated project
cat("Saving updated ArchR project...\n")
saveArchRProject(ArchRProj = proj, outputDirectory = PROJECT_DIR, load = FALSE)

cat("AUCell analysis completed successfully at", as.character(Sys.time()), "\n")
cat("Results saved in:", output_path, "\n")

EOF

# Check if R script completed successfully
if [ $? -eq 0 ]; then
    echo "AUCell analysis completed successfully at $(date)"
else
    echo "AUCell analysis failed at $(date)" >&2
    exit 1
fi