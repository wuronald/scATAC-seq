#!/bin/bash
#SBATCH --job-name=ArchR_analysis_clustering
#SBATCH --output=ArchR_analysis_clustering_%j.out
#SBATCH --error=ArchR_analysis_clustering_%j.err
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

# set genome to hg19
addArchRGenome("hg19")

# load ArchRProject
proj <-loadArchRProject(path = "HemeTutorial")


# show available matrices for the project 
getAvailableMatrices(proj) #  "GeneScoreMatrix" "TileMatrix" 

# IterativeLSI used for dimensionality reduction
proj <- addIterativeLSI(ArchRProj = proj, useMatrix = "TileMatrix", name = "IterativeLSI")

# runs seurat's FindClusters function
proj <- addClusters(input = proj, 
                    reducedDims = "IterativeLSI" # name of the reducedDims object
                    )
proj$Clusters

# table showing # of cells per cluster
# proj$Clusters %>% table

# total # of clusters
# proj$Clusters %>% table %>% length # 16

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

# save the project 
proj <- saveArchRProject(ArchRProj = proj)

# table of # of cells per cluster

# proj$Clusters_harmony %>% table # harmony clusters

# total harmony clusters
# proj$Clusters_harmony %>% table %>% length # 23

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
