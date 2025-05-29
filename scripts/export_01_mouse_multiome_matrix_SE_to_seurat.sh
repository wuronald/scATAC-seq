#!/bin/bash
#SBATCH --job-name=export_01_mouse_multiome_matrix_SE_to_seurat
#SBATCH --output=export_01_mouse_multiome_matrix_SE_to_seurat_%j.out
#SBATCH --error=export_01_mouse_multiome_matrix_SE_to_seurat_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=30G
#SBATCH --time=01:00:00

# Load necessary modules (adjust as needed for your system)
module load R/4.4.1

# Run R script
Rscript - <<'EOF'

# load libraries
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(Seurat)
library(SingleCellExperiment)
library(here)
set.seed(1)

# Set the number of threads for ArchR
addArchRThreads(threads = 18)

# Load the project
projMulti2  <- loadArchRProject(path = "mouse_multiome")

# Exporting Matrix Level data
print("Exporting Matrix Level data")
GSM_se <- getMatrixFromProject(
    ArchRProj = projMulti2,
    useMatrix <- "GeneExpressionMatrix"
)

GSM_se

print("Converting to SingleCellExperiment")
GSM_sce <- as(GSM_se, "SingleCellExperiment")
counts <- as.matrix(assay(GSM_sce, "GeneExpressionMatrix"))
assays(GSM_sce)$counts <- counts
libSizes <- colSums(counts)
sizeFactors <- libSizes/mean(libSizes)
assays(GSM_sce)$logcounts <- log2(t(t(counts)/sizeFactors) + 1)
rownames(GSM_sce) <- rowData(GSM_sce)$name

print("Converting to Seurat object")
seuratObj <- as.Seurat(GSM_sce, counts = "counts", data = "logcounts")

seuratObj



# Add the ArchR metadata to the Seurat object
print("Adding ArchR metadata to Seurat object")
seuratObj[["LSI_ATAC"]] <- CreateDimReducObject(embeddings = getReducedDims(projMulti2, reducedDims = "LSI_ATAC"), key = "LSIatac_", assay = DefaultAssay(seuratObj))
seuratObj[["LSI_RNA"]] <- CreateDimReducObject(embeddings = getReducedDims(projMulti2, reducedDims = "LSI_RNA"), key = "LSIrna_", assay = DefaultAssay(seuratObj))
seuratObj[["LSI_Combined"]] <- CreateDimReducObject(embeddings = getReducedDims(projMulti2, reducedDims = "LSI_Combined"), key = "LSIcombined_", assay = DefaultAssay(seuratObj))

print("Getting UMAP embeddings for the Seurat object")
UMAP_ATAC_matrix <- getEmbedding(
  ArchRProj = projMulti2,
  embedding = "UMAP_ATAC", # or "UMAPHarmony" if you used harmony integration
  returnDF = TRUE
)
UMAP_ATAC_matrix <- as.matrix(UMAP_ATAC_matrix)

UMAP_RNA_matrix <- getEmbedding(
  ArchRProj = projMulti2,
  embedding = "UMAP_RNA", # or "UMAPHarmony" if you used harmony integration
  returnDF = TRUE
)
UMAP_RNA_matrix <- as.matrix(UMAP_RNA_matrix)

UMAP_Combined_matrix <- getEmbedding(
  ArchRProj = projMulti2,
  embedding = "UMAP_Combined", # or "UMAPHarmony" if you used harmony integration
  returnDF = TRUE
)
UMAP_Combined_matrix <- as.matrix(UMAP_Combined_matrix)

# Add the UMAP embeddings to the Seurat object
print("Adding UMAP embeddings to Seurat object")
seuratObj[["UMAP_ATAC"]] <- CreateDimReducObject(embeddings = UMAP_ATAC_matrix, key = "UMAPatac_", assay = DefaultAssay(seuratObj))
seuratObj[["UMAP_RNA"]] <- CreateDimReducObject(embeddings = UMAP_RNA_matrix, key = "UMAPrna_", assay = DefaultAssay(seuratObj))
seuratObj[["UMAP_Combined"]] <- CreateDimReducObject(embeddings = UMAP_Combined_matrix, key = "UMAPCombined_", assay = DefaultAssay(seuratObj))


# save the Seurat object
print("Saving Seurat object")
saveRDS(seuratObj, file = here("mouse_multiome","seurat_object.rds"))


# readr::read_tsv(here("data","Signatures","Neftel_MES_ms.tsv"), col_names = FALSE)

# seuratObj <- readRDS(here("mouse_multiome","seurat_object.rds"))
EOF
