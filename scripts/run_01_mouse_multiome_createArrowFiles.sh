#!/bin/bash
#SBATCH --job-name=run_01_mouse_multiome_createArrowFiles
#SBATCH --output=run_01_mouse_multiome_createArrowFiles_%j.out
#SBATCH --error=run_01_mouse_multiome_createArrowFiles_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=60G
#SBATCH --time=03:00:00

# Load necessary modules (adjust as needed for your system)
module load R/4.4.1

# Run R script
Rscript - <<'EOF'

# load libraries
library(ArchR)
library(BSgenome.Mmusculus.UCSC.mm10)  # Changed to mouse genome
library(Seurat)
library(here)
set.seed(1)

# Set the number of threads for ArchR
addArchRThreads(threads = 18)

# input files setup:
parent_dir <- here::here("data", "mouse_multiome")

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

# Extract sample names from file paths
sampleNames <- basename(dirname(atacFiles))
sampleNames2 <- basename(dirname(rnaFiles))

# Name the files with sample names
names(atacFiles) <- sampleNames
names(rnaFiles) <- sampleNames2

# Set genome to mm10 (mouse genome)
addArchRGenome("mm10")

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

# Uncomment the following line if you want to use specific Arrow files
# ArrowFiles <- c("SM122_BR4.arrow","SM222_AR1.arrow","SM222_GL1.arrow")

print("Creating ArchRProject Project with the following arrowFiles:")
print(ArrowFiles)

projMulti1 <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "mouse_multiome",
  copyArrows = TRUE #This is recommended so that you maintain an unaltered copy for later usage.
)

# show available matrices for the project 
getAvailableMatrices(projMulti1) #  "GeneScoreMatrix" "TileMatrix" 


# import rna data
print("importing scRNA data")
seRNA <- import10xFeatureMatrix(
  input = rnaFiles,
  names = names(rnaFiles),
  strictMatch = TRUE
)

# check if the cell names in the RNA data match the Arrow files

print("Number of cell names in Arrow files not in RNA data:")
length(which(getCellNames(projMulti1) %ni% colnames(seRNA)))

# keep only the cells that are in the RNA data
print("Keeping cells that are in RNA data")
cellsToKeep <- which(getCellNames(projMulti1) %in% colnames(seRNA))
projMulti2 <- subsetArchRProject(ArchRProj = projMulti1, cells = getCellNames(projMulti1)[cellsToKeep], outputDirectory = "mouse_multiome", force = TRUE)

# add gene expression matrix to the ArchR project
projMulti2 <- addGeneExpressionMatrix(input = projMulti2, seRNA = seRNA, strictMatch = TRUE, force = TRUE)

# filter doublets; note addDoubletScores must be run previously
## Default filterRatio = 1; this is a consistent filter applied on all samples
## Can be adjusted to filter more cells
print("Starting adddoubletscores and filterDoublets")
projMulti2 <- addDoubletScores(projMulti2, force = TRUE)
projMulti2 <- filterDoublets(ArchRProj = projMulti2)

# Save the project
print("Saving the project")
saveArchRProject(ArchRProj = projMulti2, outputDirectory = "mouse_multiome", load = TRUE)
EOF

echo "ArchR analysis completed"
