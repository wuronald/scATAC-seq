#!/bin/bash
#SBATCH --job-name=run_02_mouse_multiome_reduceDims
#SBATCH --output=run_02_mouse_multiome_reduceDims_%j.out
#SBATCH --error=run_02_mouse_multiome_reduceDims_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=60G
#SBATCH --time=01:00:00

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
projMulti2  <- loadArchRProject(path = "mouse_multiome")

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
# check dimensions:
print("Dimensions of LSI_ATAC:")
print(dim(projMulti2@reducedDims$LSI_ATAC$matSVD))
print("Dimensions of LSI_RNA:")
print(dim(projMulti2@reducedDims$LSI_RNA$matSVD))

# Check which cells are present in each modality
atac_cells <- rownames(projMulti2@reducedDims$LSI_ATAC$matSVD)
rna_cells <- rownames(projMulti2@reducedDims$LSI_RNA$matSVD)

# Find missing cells
missing_in_rna <- setdiff(atac_cells, rna_cells)
missing_in_atac <- setdiff(rna_cells, atac_cells)

print(paste("Cells missing in RNA:", length(missing_in_rna)))
print(paste("Cells missing in ATAC:", length(missing_in_atac)))

# Find cells that are present in BOTH modalities
common_cells <- intersect(atac_cells, rna_cells)
print(paste("Cells present in both modalities:", length(common_cells)))

# Get all cell names from the ArchRProject
all_cells <- getCellNames(projMulti2)

# Create logical vector for subsetting - keep only cells present in both modalities
cells_to_keep <- all_cells %in% common_cells

# Subset the ArchRProject to keep only cells with both ATAC and RNA data
# projMulti2_filtered <- projMulti2[cells_to_keep]
projMulti2 <- subsetArchRProject(ArchRProj = projMulti2, cells = getCellNames(projMulti2)[cells_to_keep], outputDirectory = "mouse_multiome", force = TRUE)

# Verify the dimensions after filtering
print("Dimensions after filtering:")
print("LSI_ATAC dimensions:")
print(dim(projMulti2@reducedDims$LSI_ATAC$matSVD))
print("LSI_RNA dimensions:")
print(dim(projMulti2@reducedDims$LSI_RNA$matSVD))


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
saveArchRProject(ArchRProj = projMulti2, outputDirectory = "mouse_multiome", load = TRUE)

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
        name = "mouse_multiome-scATAC-scRNA-Combined.pdf", 
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
        name = "mouse_multiome-scATAC-scRNA-Combined_ConfusionMatrix.pdf", 
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
        name = "mouse_multiome-scATAC-scRNA-colorby-Sample.pdf", 
        ArchRProj = projMulti2, # PDF saved in the plots subfolder in the associated project directory
        addDOC = TRUE, # adds date of creation to end of filename
        width = 5, height = 5)

EOF