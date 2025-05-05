#!/bin/bash
#SBATCH --job-name=ArchR_analysis
#SBATCH --output=ArchR_analysis_%j.out
#SBATCH --error=ArchR_analysis_%j.err
#SBATCH --partition=all
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=20G
#SBATCH --time=02:00:00

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
addArchRThreads(threads = 12)

# input files setup:
inputFiles <- c("data/HemeFragments/scATAC_BMMC_R1.fragments.tsv.gz",
                "data/HemeFragments/scATAC_CD34_BMMC_R1.fragments.tsv.gz",
                "data/HemeFragments/scATAC_PBMC_R1.fragments.tsv.gz")
names(inputFiles) <- c("scATAC_BMMC_R1","scATAC_CD34_BMMC_R1", "scATAC_PBMC_R1")

# create arrow files

# set genome to hg19
addArchRGenome("hg19")

# create arrow files (subset, first sample only)
ArrowFiles <- createArrowFiles(
  inputFiles = inputFiles,
  sampleNames = names(inputFiles),
  minTSS = 4, # Dont set this too high because you can always increase later
  minFrags = 1000, 
  addTileMat = TRUE,
  addGeneScoreMat = TRUE
)

ArrowFiles

doubScores <- addDoubletScores(
  input = ArrowFiles,
  k = 10, #Refers to how many cells near a "pseudo-doublet" to count.
  knnMethod = "UMAP", #Refers to the embedding to use for nearest neighbor search.
  LSIMethod = 1
)

proj <- ArchRProject(
  ArrowFiles = ArrowFiles, 
  outputDirectory = "HemeTutorial",
  copyArrows = TRUE #This is recommended so that you maintain an unaltered copy for later usage.
)

# show available matrices for the project 
getAvailableMatrices(proj) #  "GeneScoreMatrix" "TileMatrix" 

# filter doublets; note addDoubletScores must be run previously
## Default filterRatio = 1; this is a consistent filter applied on all samples
## Can be adjusted to filter more cells
proj <- filterDoublets(ArchRProj = proj) 

# Save the project
saveArchRProject(ArchRProj = proj, outputDirectory = "HemeTutorial", load = TRUE)

# IterativeLSI used for dimensionality reduction
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")

# runs seurat's FindClusters function
proj <- addClusters(input = proj, 
                    reducedDims = "IterativeLSI" # name of the reducedDims object
                    )


# table showing # of cells per cluster
proj$Clusters %>% table

# total # of clusters
proj$Clusters %>% table %>% length # 16

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


# table of # of cells per cluster

proj$Clusters_harmony %>% table # harmony clusters

# total harmony clusters
proj$Clusters_harmony %>% table %>% length # 23

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
        name = "HemeTutorial-Plot-UMAP-Sample-HarmonyClusters.pdf", 
        ArchRProj = proj, 
        addDOC = TRUE, # adds date of creation to end of filename
        width = 5, height = 5)
EOF

echo "ArchR analysis completed"
