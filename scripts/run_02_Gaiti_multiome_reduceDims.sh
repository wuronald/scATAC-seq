#!/bin/bash
#SBATCH --job-name=run_02_Gaiti_multiome_reduceDims
#SBATCH --output=run_02_Gaiti_multiome_reduceDims_%j.out
#SBATCH --error=run_02_Gaiti_multiome_reduceDims_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=58G
#SBATCH --time=08:00:00

# Load necessary modules (adjust as needed for your system)
module load R/4.4.1

# Run R script
Rscript - <<'EOF'

# load libraries
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Seurat)
library(here)
set.seed(1)

# Set the number of threads for ArchR
addArchRThreads(threads = 8)

# Load the project
projMulti2  <- loadArchRProject(path = "Gaiti_multiome")

# IterativeLSI used for dimensionality reduction

print("Adding Iterative LSI with TileMatrix")

projMulti2 <- addIterativeLSI(
  ArchRProj = projMulti2, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "TileMatrix", 
  depthCol = "nFrags",
  name = "LSI_ATAC",
  force = TRUE
)

print("Adding Iterative LSI with GeneExpressionMatrix")

projMulti2 <- addIterativeLSI(
  ArchRProj = projMulti2, 
  clusterParams = list(
    resolution = 0.2, 
    sampleCells = 10000,
    n.start = 10
  ),
  saveIterations = FALSE,
  useMatrix = "GeneExpressionMatrix", 
  depthCol = "Gex_nUMI",
  varFeatures = 2500,
  firstSelection = "variable",
  binarize = FALSE,
  name = "LSI_RNA",
  force = TRUE
)

# Check number of cells in each
print("Dimensions of reducedDims")
dim(projMulti2@reducedDims$LSI_ATAC$matSVD)
dim(projMulti2@reducedDims$LSI_RNA$matSVD)

# Check if cell names match
print("Check if cell names match between reducedDims")
identical(rownames(projMulti2@reducedDims$LSI_ATAC$matSVD), 
          rownames(projMulti2@reducedDims$LSI_RNA$matSVD))

# Check for any issues in the matrices
print("Check for NA values in reducedDims matrices")
sum(is.na(projMulti2@reducedDims$LSI_ATAC$matSVD))
sum(is.na(projMulti2@reducedDims$LSI_RNA$matSVD))

# Get cells present in both modalities
print("Get common cells present in both modalities")
commonCells <- intersect(
  rownames(projMulti2@reducedDims$LSI_ATAC$matSVD),
  rownames(projMulti2@reducedDims$LSI_RNA$matSVD)
)

print("number of getcellnames")
length(getCellNames(projMulti2))

print("Number of common cells:")
length(commonCells)
head(commonCells)

# Subset the ArchRProject to only include common cells
print("Subsetting ArchRProject to only include common cells")
# projMulti2 <- subsetArchRProject(ArchRProj = projMulti2, cells = getCellNames(projMulti2)[commonCells], outputDirectory = "Gaiti_multiome", force = TRUE)

projMulti2 <- subsetArchRProject(
  ArchRProj = projMulti2, 
  cells = commonCells,  # Just use commonCells directly, not as an index
  outputDirectory = "Gaiti_multiome", 
  force = TRUE
)

# reduce the number of dimensions with Both scATAC and scRNA data Combined
print("Reduced dimensions with Both scATAC and scRNA data Combined")
projMulti2 <- addCombinedDims(projMulti2, reducedDims = c("LSI_ATAC", "LSI_RNA"), name =  "LSI_Combined")

# create UMAP embeddings for each of the 3 reduced dimensions
print("Adding UMAP with Iterative LSI")
projMulti2 <- addUMAP(projMulti2, reducedDims = "LSI_ATAC", name = "UMAP_ATAC", minDist = 0.8, force = TRUE)
print("Adding UMAP with GeneExpressionMatrix")
projMulti2 <- addUMAP(projMulti2, reducedDims = "LSI_RNA", name = "UMAP_RNA", minDist = 0.8, force = TRUE)
print("Adding UMAP with combined reduced dimensions")
projMulti2 <- addUMAP(projMulti2, reducedDims = "LSI_Combined", name = "UMAP_Combined", minDist = 0.8, force = TRUE)

# call clusters
print("Adding clusters")
projMulti2 <- addClusters(projMulti2, reducedDims = "LSI_ATAC", name = "Clusters_ATAC", resolution = 0.4, force = TRUE)
projMulti2 <- addClusters(projMulti2, reducedDims = "LSI_RNA", name = "Clusters_RNA", resolution = 0.4, force = TRUE)
projMulti2 <- addClusters(projMulti2, reducedDims = "LSI_Combined", name = "Clusters_Combined", resolution = 0.4, force = TRUE)

# Save the project
print("Saving Project")
saveArchRProject(ArchRProj = projMulti2, outputDirectory = "Gaiti_multiome", load = TRUE)

# plotting clusters
p1 <- plotEmbedding(projMulti2, name = "Clusters_Combined", embedding = "UMAP_ATAC", size = 1, labelAsFactors=F, labelMeans=F)
p2 <- plotEmbedding(projMulti2, name = "Clusters_Combined", embedding = "UMAP_RNA", size = 1, labelAsFactors=F, labelMeans=F)
p3 <- plotEmbedding(projMulti2, name = "Clusters_Combined", embedding = "UMAP_Combined", size = 1, labelAsFactors=F, labelMeans=F)

p <- lapply(list(p1,p2,p3), function(x){
  x + guides(color = "none", fill = "none") + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    theme(
      axis.text.x=element_blank(), 
      axis.ticks.x=element_blank(), 
      axis.text.y=element_blank(), 
      axis.ticks.y=element_blank()
    )
})

# do.call(cowplot::plot_grid, c(list(ncol = 3),p))

## save PDF of UMAPs
print("save pdf of UMAP")
plotPDF(p1,p2,p3,p,
        name = "Gaiti_multiome-scATAC-scRNA-Combined.pdf", 
        ArchRProj = projMulti2, # PDF saved in the plots subfolder in the associated project directory
        addDOC = TRUE, # adds date of creation to end of filename
        width = 5, height = 5)

# Visualize Cluster Differences with Confusion Matrix
print("Visualize Cluster Differences with Confusion Matrix")
cM_atac_rna <- confusionMatrix(paste0(projMulti2$Clusters_ATAC), paste0(projMulti2$Clusters_RNA))
cM_atac_rna <- cM_atac_rna / Matrix::rowSums(cM_atac_rna)

library(pheatmap)
p_atac_rna <- pheatmap::pheatmap(
  mat = as.matrix(cM_atac_rna), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)
p_atac_rna


## save PDF of confusion matrix
print("save pdf of confusion matrix")
plotPDF(p_atac_rna,
        name = "Gaiti_multiome-scATAC-scRNA-Combined_ConfusionMatrix.pdf", 
        ArchRProj = projMulti2, # PDF saved in the plots subfolder in the associated project directory
        addDOC = TRUE, # adds date of creation to end of filename
        width = 5, height = 5)

# plots Sample for each of the 3 reduced dimensions
p4 <- plotEmbedding(projMulti2, embedding = "UMAP_ATAC", size = 1, labelAsFactors=F, labelMeans=F,
colorBy = "cellColData", name = "Sample")
p5 <- plotEmbedding(projMulti2, embedding = "UMAP_RNA", size = 1, labelAsFactors=F, labelMeans=F,
colorBy = "cellColData", name = "Sample")
p6 <- plotEmbedding(projMulti2, embedding = "UMAP_Combined", size = 1, labelAsFactors=F, labelMeans=F,
colorBy = "cellColData", name = "Sample")
plotPDF(p4,p5,p6,
        name = "Gaiti_multiome-scATAC-scRNA-colorby-Sample.pdf", 
        ArchRProj = projMulti2, # PDF saved in the plots subfolder in the associated project directory
        addDOC = TRUE, # adds date of creation to end of filename
        width = 5, height = 5)

EOF