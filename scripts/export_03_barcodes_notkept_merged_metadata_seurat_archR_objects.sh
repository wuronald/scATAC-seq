#!/bin/bash
#SBATCH --job-name=export_03_barcodes_notkept_merged_metadata_seurat_archR_objects
#SBATCH --output=export_03_barcodes_notkept_merged_metadata_seurat_archR_objects_%j.out
#SBATCH --error=export_03_barcodes_notkept_merged_metadata_seurat_archR_objects_%j.err
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
library(Seurat)
library(here)
set.seed(1)

# Load the project
print(paste("Loading ArchR project from", archr_dir))
proj <- loadArchRProject(path = archr_dir)

# get barcodes from ArchR object
print("Getting barcodes from ArchR project")
cells = getCellNames(proj)

# Load Seurat object metadata
metadata_path <- here(archr_dir, "merged_seurat_metadata.csv")
if (!file.exists(metadata_path)) {
    stop(paste("Seurat Metadata file does not exist:", metadata_path))
}
seurat_metadata <- readr::read_csv(metadata_path)
print(paste("Seurat metadata loaded from", metadata_path))

# get barcodes from Seurat object
print("Getting barcodes from Seurat metadata")
seurat_atac_barcodes <- seurat_metadata$seurat_atac_barcode
seurat_rna_barcodes <- seurat_metadata$seurat_gex_barcode
print(paste("Number of cells in the seurat metadata:", length(seurat_rna_barcodes))) # 24644

# Export ArchR barcodes not kept
cellsNotKept <- cells[!(cells %in% seurat_rna_barcodes)]
readr::write_csv(
  data.frame(archr_barcode = cellsNotKept),
  file = file.path(output_dir, "archr_barcodes_not_in_seurat.csv")
)

# Export Seurat metadata for barcodes not in ArchR
seuratBarcodesNotInArchR <- seurat_rna_barcodes[!(seurat_rna_barcodes %in% cells)]
seurat_metadata_not_in_archr <- seurat_metadata[seurat_metadata$seurat_gex_barcode %in% seuratBarcodesNotInArchR, ]
readr::write_csv(
  seurat_metadata_not_in_archr,
  file = file.path(output_dir, "seurat_metadata_not_in_archr.csv")
)

# Import Seurat unfiltered barcodes
seurat_unfiltered_metadata <- readr::read_csv("raw_metadata.csv")

# Add prefixed barcode column (barcode format eg.: Zadeh__C0736__5117#AAACAGCCAGGACACA-1)
seurat_unfiltered_metadata$prefixed_barcode <- paste0(
  seurat_unfiltered_metadata$orig.ident, "#", seurat_unfiltered_metadata$gex_barcode
)

# Compare prefixed_barcode with barcodes from cellsNotKept and seuratBarcodesNotInArchR
print(paste("number of cells in original unfiltered seurat:", nrow(seurat_unfiltered_metadata)))
print(paste("number of cells in ArchR not kept:", length(cellsNotKept)))
print(paste("number of cells in Seurat filtered not kept:", length(seuratBarcodesNotInArchR)))

# 1. Find prefixed_barcodes in seurat_unfiltered_metadata that are in cellsNotKept

prefixed_in_cellsNotKept <- seurat_unfiltered_metadata$prefixed_barcode[
  seurat_unfiltered_metadata$prefixed_barcode %in% cellsNotKept
]
cat("Number of prefixed_barcodes in cellsNotKept:", length(prefixed_in_cellsNotKept), "\n")
if (length(prefixed_in_cellsNotKept) > 0) {
  print(head(prefixed_in_cellsNotKept))
}

# 2. Find prefixed_barcodes in seurat_unfiltered_metadata that are in seuratBarcodesNotInArchR
prefixed_in_seuratBarcodesNotInArchR <- seurat_unfiltered_metadata$prefixed_barcode[
  seurat_unfiltered_metadata$prefixed_barcode %in% seuratBarcodesNotInArchR
]
cat("Number of prefixed_barcodes in seuratBarcodesNotInArchR:", length(prefixed_in_seuratBarcodesNotInArchR), "\n")
if (length(prefixed_in_seuratBarcodesNotInArchR) > 0) {
  print(head(prefixed_in_seuratBarcodesNotInArchR))
}
EOF