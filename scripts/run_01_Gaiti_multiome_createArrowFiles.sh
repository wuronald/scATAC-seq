#!/bin/bash
#SBATCH --job-name=run_01_Gaiti_multiome_createArrowFiles
#SBATCH --output=run_01_Gaiti_multiome_createArrowFiles_%j.out
#SBATCH --error=run_01_Gaiti_multiome_createArrowFiles_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=58G
#SBATCH --time=08:00:00

# ==============================================================================
# Purpose: SLURM batch submission script to create ArchR Arrow files and initialize
#          the integrated ArchR project for the Gaiti et al. 10x Multiome dataset (26 samples).
# Workflow:
#   1. Scans 'data/Gaiti_multiome/' for ATAC fragments and RNA H5 matrix pairs.
#   2. Generates Arrow files using the hg38 genome (minTSS = 4, minFrags = 1000).
#   3. Initializes ArchRProject and imports scRNA-seq GeneExpressionMatrix.
#   4. Subsets project to matching cells and adds GeneExpressionMatrix.
#   5. Computes doublet scores, exports predicted doublets to Gaiti_multiome_doublets.txt,
#      and filters doublets.
#   6. Saves processed ArchR project into 'Gaiti_multiome/'.
# Input:   'data/Gaiti_multiome/<Sample>/' containing:
#            - *fragments.tsv.gz & *fragments.tsv.gz.tbi
#            - filtered_feature_bc_matrix.h5 (from convert_10x_geo_files_to_h5.sh)
# Output:  Arrow files + 'Gaiti_multiome/' ArchR project directory + Gaiti_multiome_doublets.txt
#
# Example usage:
#   sbatch scripts/run_01_Gaiti_multiome_createArrowFiles.sh
# ==============================================================================

# Note on Sample Swap:
# Remember to rename and move swapped 6467_A_filtered_feature_bc_matrix.h5
# and 6467_B_filtered_feature_bc_matrix.h5 to their correct sample folders before running this script!
# mv /cluster/projects/wouterslab/ArchR103_4/data/Gaiti_multiome/GSM8820713_6467_A/6467_A_filtered_feature_bc_matrix.h5 /cluster/projects/wouterslab/ArchR103_4/data/Gaiti_multiome/GSM8820714_6467_B/temp.h5 && \
# mv /cluster/projects/wouterslab/ArchR103_4/data/Gaiti_multiome/GSM8820714_6467_B/6467_B_filtered_feature_bc_matrix.h5 /cluster/projects/wouterslab/ArchR103_4/data/Gaiti_multiome/GSM8820713_6467_A/6467_A_filtered_feature_bc_matrix.h5 && \
# mv /cluster/projects/wouterslab/ArchR103_4/data/Gaiti_multiome/GSM8820714_6467_B/temp.h5 /cluster/projects/wouterslab/ArchR103_4/data/Gaiti_multiome/GSM8820714_6467_B/6467_B_filtered_feature_bc_matrix.h5

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

# Extract sample names from file paths independently
names(atacFiles) <- basename(dirname(atacFiles))
names(rnaFiles)  <- basename(dirname(rnaFiles))

# Ensure symmetric sample matching
commonSamples <- intersect(names(atacFiles), names(rnaFiles))
if (length(commonSamples) == 0) {
  stop("No matching sample pairs found between ATAC fragments and RNA matrix files.")
}
atacFiles <- atacFiles[commonSamples]
rnaFiles  <- rnaFiles[commonSamples]

# Set genome to hg (human genome)
addArchRGenome("hg38")

print("ATAC input files:")
print(atacFiles)
print("RNA input files:")
print(rnaFiles)

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

# Create ArchR project
print("Creating ArchR Project")
projMulti1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "Gaiti_multiome",
  copyArrows = TRUE #This is recommended so that you maintain an unaltered copy for later usage.
)

# show available matrices for the project
print("Available matrices in the ArchR Project:") 
getAvailableMatrices(projMulti1) #  "GeneScoreMatrix" "TileMatrix" 

# Current directory
print(paste0("Current working directory:", getwd()))

# import rna data
print("importing scRNA data")
seRNA <- import10xFeatureMatrix(
  input = rnaFiles,
  names = names(rnaFiles),
  strictMatch = TRUE
)

print("Number of cells in ArchRProject:")
length(getCellNames(projMulti1))

print("Number of cell names in Arrow files not in RNA data:")
length(which(getCellNames(projMulti1) %ni% colnames(seRNA)))

# keep only the cells that are in the RNA data
print("Keeping cells that are in RNA data")
cellsToKeep <- which(getCellNames(projMulti1) %in% colnames(seRNA))
print(paste("Number of overlapping cells:", length(cellsToKeep)))
print("First few overlapping cell names:")
head(cellsToKeep)

# Check if any cells to keep
if (length(cellsToKeep) == 0) {
  stop("No overlapping cells found between RNA data and Arrow files!")
}

# Tabulate cell counts before subsetting
print("Cell counts before subsetting:")
table(projMulti1$Sample)

# Subset the ArchR project to keep only these cells
print("Subsetting ArchR project to keep only overlapping cells")
projMulti2 <- subsetArchRProject(ArchRProj = projMulti1, cells = getCellNames(projMulti1)[cellsToKeep], outputDirectory = "Gaiti_multiome", force = TRUE)

# Tabulate cell counts after subsetting
print("Cell counts after subsetting:")
table(projMulti2$Sample)

# remove projMulti1 from memory
rm(projMulti1)

# add gene expression matrix to the ArchR project
print("Adding GeneExpressionMatrix to ArchR Project")
projMulti2 <- addGeneExpressionMatrix(input = projMulti2, seRNA = seRNA, strictMatch = TRUE, force = TRUE)

# Tabulate cell counts after adding GeneExpressionMatrix
print("Cell counts after adding gene expression matrix:")
table(projMulti2$Sample)

# filter doublets; note addDoubletScores must be run previously
print("Adding doublet scores to ArchR Project")
projMulti2 <- addDoubletScores(projMulti2, force = TRUE)

# Function to preview and save doublets before filtering
saveDoubletBarcodes <- function(ArchRProj, cutEnrich = 1, cutScore = -Inf, 
                                 filterRatio = 1, outfile = "doublets.txt"){
  
  df <- getCellColData(ArchRProj, c("Sample", "DoubletEnrichment", "DoubletScore"))
  splitDF <- split(seq_len(nrow(df)), as.character(df$Sample))
  
  doublets <- lapply(splitDF, function(y){
    x <- df[y, , drop = FALSE]
    n <- nrow(x)
    x <- x[order(x$DoubletEnrichment, decreasing = TRUE), ]
    
    if(!is.null(cutEnrich)){
      x <- x[which(x$DoubletEnrichment >= cutEnrich), ]
    }
    if(!is.null(cutScore)){
      x <- x[which(x$DoubletScore >= cutScore), ]
    }
    
    if(nrow(x) > 0){
      head(rownames(x), filterRatio * n * (n / 100000))
    } else {
      NULL
    }
  }) %>% unlist(use.names = FALSE)
  
  write.table(doublets, file = outfile, quote = FALSE, 
              row.names = FALSE, col.names = FALSE)
  
  message("Saved ", length(doublets), " doublet barcodes to ", outfile)
  return(doublets)
}

# Use it:
print("Saving doublet barcodes before filtering")
doublet_list <- saveDoubletBarcodes(projMulti2, outfile = "Gaiti_multiome_doublets.txt")

print("filterDoublets")
projMulti2 <- filterDoublets(ArchRProj = projMulti2)

# Tabulate cell counts after doublet filtering
print("Cell counts after doublet filtering:")
table(projMulti2$Sample)

# Save the project
print("Saving the project")
saveArchRProject(ArchRProj = projMulti2, outputDirectory = "Gaiti_multiome", load = TRUE)

EOF

echo "ArchR analysis completed"
