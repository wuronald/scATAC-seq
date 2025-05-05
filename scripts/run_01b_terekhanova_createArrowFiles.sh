#!/bin/bash
#SBATCH --job-name=run_01b_terekhanova_createArrowFiles
#SBATCH --output=run_01b_terekhanova_createArrowFiles_%j.out
#SBATCH --error=run_01b_terekhanova_createArrowFiles_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=60G
#SBATCH --time=02:00:00

# Load necessary modules (adjust as needed for your system)
#module load R/4.2.1
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

# input files setup:
parent_dir <- here::here("data","Terekhanova")

# Get a list of subdirectories 
subdirs <- list.dirs(parent_dir, recursive = FALSE, full.names = FALSE)
# Extract prefixes before the underscore to get sample names
prefixes <- sapply(subdirs, function(x) strsplit(x, "_")[[1]][1])

# Get list of relative path to fragments.tsv.gz for all Terekhanova samples
list_files <- function(dir_path) {
  full_paths <- list.files(dir_path, recursive = TRUE, full.names = FALSE, pattern = "fragments.tsv.gz$")
  return(full_paths)
}
matching_files <- list_files(dir_path = here::here("data","Terekhanova"))

# make character vector with sample name and full path to their respective fragments.tsv.gz
inputFiles <- paste0(here::here("data","Terekhanova",matching_files))
names(inputFiles) <- prefixes

# subset inputfiles for the 5 that failed previously
#inputFiles2 <- inputFiles1[c("C3L-03405", "C3L-03968", "C3N-01814", "C3N-02186", "C3N-01818")]
# set genome to hg38
addArchRGenome("hg38")

print("There are the inputfiles:")
print(inputFiles)

print("Starting create arrowFiles")

# create arrow files
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, # Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

# vector of paths to all 18 of Terekhanova samples
ArrowFiles <- c("C3L-02705.arrow","C3L-03405.arrow","C3L-03968.arrow","C3N-00662.arrow","C3N-00663.arrow","C3N-01334.arrow",
                "C3N-01518.arrow","C3N-01798.arrow","C3N-01814.arrow","C3N-01816.arrow","C3N-01818.arrow","C3N-02181.arrow",
                "C3N-02186.arrow","C3N-02188.arrow","C3N-02769.arrow","C3N-02783.arrow","C3N-02784.arrow","C3N-03186.arrow")
print("Starting add doubletscores")
print(ArrowFiles)

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "Terekhanova",
  copyArrows = TRUE #This is recommended so that you maintain an unaltered copy for later usage.
)

# show available matrices for the project 
getAvailableMatrices(proj) #  "GeneScoreMatrix" "TileMatrix" 

# filter doublets; note addDoubletScores must be run previously
## Default filterRatio = 1; this is a consistent filter applied on all samples
## Can be adjusted to filter more cells
proj <- filterDoublets(ArchRProj = proj) 

# Save the project
saveArchRProject(ArchRProj = proj, outputDirectory = "Terekhanova", load = TRUE)

EOF

echo "ArchR analysis completed"
