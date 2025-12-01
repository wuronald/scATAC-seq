#!/bin/bash
#SBATCH --job-name=run_01_Gaiti_multiome_createArrowFiles
#SBATCH --output=run_01_Gaiti_multiome_createArrowFiles_%j.out
#SBATCH --error=run_01_Gaiti_multiome_createArrowFiles_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=58G
#SBATCH --time=12:00:00

# Load necessary modules (adjust as needed for your system)
module load R/4.4.1

# Run R script
Rscript - <<'EOF'

# load libraries
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)  # human genome
library(Seurat)
library(here)
set.seed(1)

# Set the number of threads for ArchR
addArchRThreads(threads = 6)

# input files setup:
parent_dir <- here::here("data", "Gaiti_multiome")

# Check if directory exists
if (!dir.exists(parent_dir)) {
  stop(paste("Directory does not exist:", parent_dir))
}

# Get a list of subdirectories 
subdirs <- list.dirs(parent_dir, recursive = FALSE, full.names = FALSE)

if (length(subdirs) == 0) {
  stop("No subdirectories found in", parent_dir)
}

# Find ATAC fragments files
atacFiles <- list.files(
  path = parent_dir,
  pattern = "fragments.tsv.gz$",
  recursive = TRUE,
  full.names = TRUE
)
# Find RNA matrix files
rnaFiles <- list.files(
  path = parent_dir,
  pattern = "filtered_feature_bc_matrix.h5$",
  recursive = TRUE,
  full.names = TRUE
)

# Check if any atac files were found
if (length(atacFiles) == 0) {
  stop("No ATAC fragments files found in the specified directory.")
}
# Check if any RNA files were found
if (length(rnaFiles) == 0) {
  stop("No RNA matrix files found in the specified directory.")
}

# print("Finding matrix, barcode, and features files for RNA data")
# # Find matrix files
# matrix_file <- list.files(
#   path = parent_dir,
#   pattern = "matrix.mtx.gz$",
#   recursive = TRUE,
#   full.names = TRUE
# )
# # Find barcode files
# barcode_file <- list.files(
#   path = parent_dir,
#   pattern = "barcodes.tsv.gz$",
#   recursive = TRUE,
#   full.names = TRUE
# )
# # Find features files
# features_file <- list.files(
#   path = parent_dir,
#   pattern = "features.tsv.gz$",
#   recursive = TRUE,
#   full.names = TRUE
# )

# Extract sample names from file paths
sampleNames <- basename(dirname(atacFiles))

# Name the files with sample names
names(atacFiles) <- sampleNames
names(rnaFiles) <- sampleNames

# Set genome to hg (human genome)
addArchRGenome("hg38")

print("ATAC input files:")
print(atacFiles)
print("RNA input files:")
print(rnaFiles)

# Remove sample "Zadeh_Shelia__61" (deemed low quality in scRNA))
remove_sample <- "Zadeh_Shelia__61"
if (remove_sample %in% names(atacFiles)) {
  atacFiles <- atacFiles[names(atacFiles) != remove_sample]
}

print("Starting create arrowFiles")

# create arrow files
ArrowFiles <- createArrowFiles(
  inputFiles = atacFiles,
  sampleNames = names(atacFiles),
  minTSS = 4, # Dont set this too high because you can always increase later
  minFrags = 1000,
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

# vector of paths to all samples
print("Created ArrowFiles:")
print(ArrowFiles)
# Uncomment the following line if you want to use specific Arrow files
# ArrowFiles <- c("SM122_BR4.arrow","SM222_AR1.arrow","SM222_GL1.arrow")

# Create ArchR project
print("Creating ArchR Project")
projMulti1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "Gaiti_multiome",
  copyArrows = TRUE #This is recommended so that you maintain an unaltered copy for later usage.
)
# add doublet scores
print("Adding doublet scores to ArchR Project")
projMulti1 <- addDoubletScores(projMulti1, force = TRUE)

# show available matrices for the project
print("Available matrices in the ArchR Project:") 
getAvailableMatrices(projMulti1) #  "GeneScoreMatrix" "TileMatrix" 

# # Check if matrix_file, barcode_file, and features_file are of the same length
# if (length(matrix_file) != length(barcode_file) || length(matrix_file) != length(features_file)) {
#   stop("The number of matrix, barcode, and features files do not match!")
# }

# # Extract prefixes from ArrowFiles
# arrow_prefixes <- sub("\\.arrow$", "", basename(ArrowFiles))

# # Extract prefixes from matrix_file, barcode_file, and features_file
# matrix_prefixes <- sub("_matrix\\.mtx\\.gz$", "", basename(matrix_file))
# barcode_prefixes <- sub("_barcodes\\.tsv\\.gz$", "", basename(barcode_file))
# features_prefixes <- sub("_features\\.tsv\\.gz$", "", basename(features_file))

# # Check if all prefixes match
# if (!all(matrix_prefixes == barcode_prefixes & barcode_prefixes == features_prefixes)) {
#   stop("Prefixes of matrix, barcode, and features files do not match!")
# }

# # Check if prefixes match the ArrowFile prefixes
# if (!all(matrix_prefixes %in% arrow_prefixes)) {
#   stop("Prefixes of input files do not match the prefixes of ArrowFiles!")
# }

# print("All checks passed: matrix, barcode, and features files are consistent and match ArrowFiles.")

# # For each Arrow file, read the corresponding RNA data and add GeneExpressionMatrix

# for (i in seq_along(ArrowFiles)) {
#   message(paste0("Adding GeneExpression matrix to Arrow file for sample ", names(ArrowFiles)[i], ": ", ArrowFiles[i]))

# # Read the data
# message("Reading matrix...")
# print(matrix_file[i])
# mat <- readMM(matrix_file[i])
# message("Reading features and barcodes...")
# barcodes <- fread(barcode_file[i], header = FALSE)$V1
# features <- fread(features_file[i], header = FALSE)

# # Name the feature columns for clarity
# # 10x Multiome features.tsv has 3 columns: ID, Name, Type
# colnames(features) <- c("id", "name", "type","seqnames", "start", "end")

# message("Splitting the matrix...")

# # Find the row indices for "Gene Expression"
# gex_indices <- which(features$type == "Gene Expression")

# # Subset the main matrix to get *only* the GEX counts
# gex_mat <- mat[gex_indices, ]

# # Get the corresponding feature metadata for GEX
# gex_features <- features[gex_indices, ]

# # Set the row and column names for the GEX matrix
# # Use the gene *name* (e.g., "SOX2") as rownames, as this is standard.
# rownames(gex_mat) <- gex_features$name
# colnames(gex_mat) <- barcodes

# message("Creating SummarizedExperiment...")

# # Create the SummarizedExperiment (seGEX)
# seGEX <- SummarizedExperiment(
#   assays = list(counts = gex_mat),
#   rowData = gex_features,
#   colData = S4Vectors::DataFrame(row.names = colnames(gex_mat))
# )

# print("Verifying SummarizedExperiment object:")
# print(seGEX)
# print("Row data of SummarizedExperiment:")
# print(rowData(seGEX))

# Current directory
print(paste0("Current working directory:",getwd() ))

# Apply the patch
#print("Applying patch to import10xFeatureMatrix function")
#source("scripts/patch_function.r")
#patch_import10x_function("ArchR")

# import rna data
print("importing scRNA data")
seRNA <- import10xFeatureMatrix(
  input = rnaFiles,
  names = names(rnaFiles),
  #force = TRUE,
  strictMatch = FALSE
)

# Add GeneExpressionMatrix to ArrowFiles
print("Adding GeneExpressionMatrix to ArrowFile")

# addGeneExpressionMatrix(
#   input = projMulti1,
#   #input = paste0("./", ArrowFiles[i]),
#   seRNA = seGEX,
#   force = TRUE
# )
# message("Done! Your GEX data is now in the arrowfiles.")
# }


# Check if the cell names in the RNA data match the Arrow files:
# 1. Retrieve GeneExpressionMatrix from the ArchRProject
# print("Retrieving GeneExpressionMatrix from ArchRProject")
# seRNA <- getMatrixFromProject(ArchRProj = projMulti1, 
#                               useMatrix = "GeneExpressionMatrix")

print("Number of cell names in Arrow files not in RNA data:")
length(which(getCellNames(projMulti1) %ni% colnames(seRNA)))

# keep only the cells that are in the RNA data
print("Keeping cells that are in RNA data")
cellsToKeep <- which(getCellNames(projMulti1) %in% colnames(seRNA))

# Check if any cells to keep
if (length(cellsToKeep) == 0) {
  stop("No overlapping cells found between RNA data and Arrow files!")
}

# Subset the ArchR project to keep only these cells
print("Subsetting ArchR project to keep only overlapping cells")
projMulti2 <- subsetArchRProject(ArchRProj = projMulti1, cells = getCellNames(projMulti1)[cellsToKeep], outputDirectory = "Gaiti_multiome", force = TRUE)

# remove projMulti1 from memory
rm(projMulti1)

# add gene expression matrix to the ArchR project
print("Adding GeneExpressionMatrix to ArchR Project")
projMulti2 <- addGeneExpressionMatrix(input = projMulti2, seRNA = seRNA, strictMatch = TRUE, force = TRUE)

# Save the project
print("Saving the project")
saveArchRProject(ArchRProj = projMulti2, outputDirectory = "Gaiti_multiome", load = TRUE)

# filter doublets; note addDoubletScores must be run previously
## Default filterRatio = 1; this is a consistent filter applied on all samples
## Can be adjusted to filter more cells
print("filterDoublets")
# projMulti2 <- addDoubletScores(projMulti2, force = TRUE)
projMulti2 <- filterDoublets(ArchRProj = projMulti2)

# Save the project
print("Saving the project")
saveArchRProject(ArchRProj = projMulti2, outputDirectory = "Gaiti_multiome", load = TRUE)

EOF

echo "ArchR analysis completed"
