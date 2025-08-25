#!/bin/bash
#SBATCH --job-name=export_02_merged_metadata_seurat_archR_objects
#SBATCH --output=export_02_merged_metadata_seurat_archR_objects_%j.out
#SBATCH --error=export_02_merged_metadata_seurat_archR_objects_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=20G
#SBATCH --time=00:30:00

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
print("Loading ArchR project")
proj <- loadArchRProject(path = "mouse_multiome_harmony_test")

# Load Seurat object metadata
seurat_metadata <- readr::read_csv(here("mouse_multiome_harmony_test","merged_seurat_metadata.csv"))

# get barcodes from Seurat object
print("Getting barcodes from Seurat metadata")
seurat_atac_barcodes <- seurat_metadata$seurat_atac_barcode
seurat_rna_barcodes <- seurat_metadata$seurat_gex_barcode
print(paste("Number of cells in the seurat metadata:", length(seurat_rna_barcodes))) # 24644

# get barcodes from ArchR object
print("Getting barcodes from ArchR project")
cells = getCellNames(proj)
print(paste("Number of cells in the ArchR project:", length(cells))) # 26622

# filter cells based on Seurat rna barcodes
# atac barcodes seems to be missing from the ArchR object, so we will filter based on RNA barcodes
#cellsToKeep <- which(cells %in% seurat_atac_barcodes)
cellsToKeep <- which(cells %in% seurat_rna_barcodes)

print(paste("Number of cells to keep:", length(cellsToKeep))) # 21559

# Subset ArchR project to keep only the cells present in the Seurat object
print("Subsetting ArchR project to keep only the cells present in the Seurat object")

projSubset <- subsetArchRProject(
  ArchRProj = proj,
  cells = proj$cellNames[cellsToKeep],
  outputDirectory = "mouse_multiome_harmony_merged_subset",
  dropCells = TRUE,
  force = TRUE
  )
# projSubset <- proj[cellsToKeep, ] # improper way to subset ArchR project

# subset seurat metadata to keep only the cells present in the ArchR project
seurat_metadata_subset <- seurat_metadata[seurat_metadata$seurat_gex_barcode %in% projSubset$cellNames, ] 

# re-order the Seurat metadata to match the order of cells in the ArchR project
seurat_metadata_subset <- seurat_metadata_subset[match(projSubset$cellNames, seurat_metadata_subset$seurat_gex_barcode), ]

# check if order of cells is the same
are_identical <- all(
projSubset$cellNames == seurat_metadata_subset$seurat_gex_barcode
)
print(paste("Are the cell names in the same order?", are_identical))

# Dynamically add columns from seurat_metadata_subset to ArchR metadata
cols_to_add <- c("hybrid_pair","Azimuth_class", "Azimuth_subclass") # add more column names as needed

for (col in cols_to_add) {
    projSubset <- addCellColData(
        projSubset,
        data = seurat_metadata_subset[[col]],
        name = col,
        cells = getCellNames(projSubset)
    )
}

# Save the project
print("Saving the project")
saveArchRProject(ArchRProj = projSubset, outputDirectory = "mouse_multiome_harmony_merged_subset", load = TRUE)

print("Script completed!")
EOF