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

# Usage: sbatch scripts/export_02_merged_metadata_seurat_archR_objects.sh human_multiome_harmony
# Usage: sbatch scripts/export_02_merged_metadata_seurat_archR_objects.sh mouse_multiome_harmony_test
ARCHR_DIR="$1"
export ARCHR_DIR
echo "ARCHR_DIR is set to '$ARCHR_DIR'"
if [ -z "$ARCHR_DIR" ]; then
    echo "Usage: $0 <ARCHR_DIR>"
    exit 1
fi

# Run R script
Rscript - <<'EOF'

# load libraries
library(ArchR)
library(Seurat)
library(here)
set.seed(1)

archr_dir <- Sys.getenv("ARCHR_DIR")
print(paste("ARCHR_DIR is:", archr_dir))
output_dir <- paste0(archr_dir, "_merged")

# Detect species from folder name
if (grepl("human", archr_dir, ignore.case=TRUE)) {
    library(BSgenome.Hsapiens.UCSC.hg38)
    genome <- "hg38"
} else if (grepl("mouse", archr_dir, ignore.case=TRUE)) {
    library(BSgenome.Mmusculus.UCSC.mm10)
    genome <- "mm10"
} else {
    stop("Could not determine species from ARCHR_DIR name.")
}


# Set the number of threads for ArchR
addArchRThreads(threads = 18)

# Load the project
print(paste("Loading ArchR project from", archr_dir))
proj <- loadArchRProject(path = archr_dir)

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
  outputDirectory = output_dir,
  dropCells = TRUE,
  force = TRUE
  )

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
cols_to_add <- c("hybrid_pair","neftel_4_state",
"Azimuth_class", "Azimuth_subclass",
"Ambiguous", "PIMO_status", "PIMO_up_status") # add more column names as needed

# Only keep columns that exist in seurat_metadata_subset
cols_to_add <- cols_to_add[cols_to_add %in% colnames(seurat_metadata_subset)]
print(paste("Columns to add to ArchR metadata:", paste(cols_to_add, collapse = ", ")))

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
saveArchRProject(ArchRProj = projSubset, outputDirectory = output_dir, load = TRUE)

print("Script completed!")
EOF