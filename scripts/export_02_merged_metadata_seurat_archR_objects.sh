#!/bin/bash
#SBATCH --job-name=export_02_merged_metadata_seurat_archR_objects
#SBATCH --output=export_02_merged_metadata_seurat_archR_objects_%j.out
#SBATCH --error=export_02_merged_metadata_seurat_archR_objects_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=9
#SBATCH --mem=30G
#SBATCH --time=01:00:00

# Load necessary modules (adjust as needed for your system)
module load R/4.4.1

# Note: Run R markdown notebook to generate/pre-process seurat object metadata before running this script.
# eg. human_multiome_merged_seurat_archr.Rmd for human_multiome
# eg. mouse_multiome_merged_seurat_archr.Rmd for mouse_multiome
# eg. Gaiti_multiome_merged_seurat_archr.Rmd for Gaiti_multiome

# Usage:
#   sbatch scripts/export_02_merged_metadata_seurat_archR_objects.sh <ARCHR_DIR> [OUTPUT_DIR] [METADATA_PATH] [COLUMNS]
#
# Arguments (all positional):
#   $1  ARCHR_DIR      Path to the ArchR project directory (required)
#   $2  OUTPUT_DIR     Output directory for the subsetted ArchR project
#                      (default: <ARCHR_DIR>_merged)
#   $3  METADATA_PATH  Path to the Seurat metadata CSV file
#                      (default: <ARCHR_DIR>/merged_seurat_metadata.csv)
#   $4  COLUMNS        Comma-separated list of metadata columns to add to ArchR object
#                      (default: hybrid_pair,neftel_4_state,Region.annotation,
#                       mislabelled,Azimuth_class,Azimuth_subclass,
#                       Ambiguous,PIMO_status,PIMO_up_status)
#
# Examples:
#   sbatch scripts/export_02_merged_metadata_seurat_archR_objects.sh human_multiome_harmony_merged_malig_peak
#   sbatch scripts/export_02_merged_metadata_seurat_archR_objects.sh human_multiome_harmony_merged_malig_peak human_multiome_hmmp_nomura
#   sbatch scripts/export_02_merged_metadata_seurat_archR_objects.sh human_multiome_harmony_merged_malig_peak "" data/metadata/human_multiome/merged_seurat_nomura_metadata_subset.csv
#   sbatch scripts/export_02_merged_metadata_seurat_archR_objects.sh human_multiome_harmony_merged_malig_peak human_multiome_hmmp_nomura data/metadata/human_multiome/merged_seurat_nomura_metadata_subset.csv "meta_module"

# --- Parse positional arguments ---
ARCHR_DIR="$1"
OUTPUT_DIR="$2"
METADATA_PATH="$3"
COLUMNS="$4"

if [ -z "$ARCHR_DIR" ]; then
    echo "Error: ARCHR_DIR is required."
    echo "Usage: $0 <ARCHR_DIR> [OUTPUT_DIR] [METADATA_PATH] [COLUMNS]"
    exit 1
fi

# --- Apply defaults if not set ---
if [ -z "$OUTPUT_DIR" ]; then
    OUTPUT_DIR="${ARCHR_DIR}_merged"
fi

if [ -z "$METADATA_PATH" ]; then
    METADATA_PATH="${ARCHR_DIR}/merged_seurat_metadata.csv"
fi

# --- Export variables for R ---
export ARCHR_DIR OUTPUT_DIR METADATA_PATH COLUMNS

echo "ARCHR_DIR     : $ARCHR_DIR"
echo "OUTPUT_DIR    : $OUTPUT_DIR"
echo "METADATA_PATH : $METADATA_PATH"
echo "COLUMNS       : ${COLUMNS:-<using defaults>}"

# Run R script
Rscript - <<'EOF'

# load libraries
library(ArchR)
library(Seurat)
library(here)
set.seed(1)

archr_dir     <- Sys.getenv("ARCHR_DIR")
output_dir    <- Sys.getenv("OUTPUT_DIR")
metadata_path <- Sys.getenv("METADATA_PATH")
columns_env   <- Sys.getenv("COLUMNS")

print(paste("ARCHR_DIR     :", archr_dir))
print(paste("OUTPUT_DIR    :", output_dir))
print(paste("METADATA_PATH :", metadata_path))

# Helper: print a tidy summary block for an ArchR project
print_archr_summary <- function(proj, label) {
    cd <- as.data.frame(proj@cellColData)
    cat("\n", strrep("=", 60), "\n", sep = "")
    cat(" SUMMARY:", label, "\n")
    cat(strrep("=", 60), "\n", sep = "")
    cat(sprintf("  %-30s %d\n",   "Total cells:",          nrow(cd)))
    cat(sprintf("  %-30s %d\n",   "Number of samples:",    length(unique(cd$Sample))))
    cat("\n  Cells per sample:\n")
    sample_tbl <- sort(table(cd$Sample), decreasing = TRUE)
    for (nm in names(sample_tbl)) {
        cat(sprintf("    %-28s %d\n", nm, sample_tbl[[nm]]))
    }
    # Fragment size / TSS stats if available
    if ("TSSEnrichment" %in% colnames(cd)) {
        cat(sprintf("\n  %-30s mean=%.2f, median=%.2f, min=%.2f, max=%.2f\n",
            "TSS Enrichment:",
            mean(cd$TSSEnrichment,   na.rm=TRUE),
            median(cd$TSSEnrichment, na.rm=TRUE),
            min(cd$TSSEnrichment,    na.rm=TRUE),
            max(cd$TSSEnrichment,    na.rm=TRUE)))
    }
    if ("nFrags" %in% colnames(cd)) {
        cat(sprintf("  %-30s mean=%.0f, median=%.0f, min=%.0f, max=%.0f\n",
            "nFrags:",
            mean(cd$nFrags,   na.rm=TRUE),
            median(cd$nFrags, na.rm=TRUE),
            min(cd$nFrags,    na.rm=TRUE),
            max(cd$nFrags,    na.rm=TRUE)))
    }
    if ("DoubletScore" %in% colnames(cd)) {
        cat(sprintf("  %-30s mean=%.4f, median=%.4f\n",
            "Doublet Score:",
            mean(cd$DoubletScore,   na.rm=TRUE),
            median(cd$DoubletScore, na.rm=TRUE)))
    }
    cat(strrep("=", 60), "\n\n", sep = "")
}

# Detect species from folder name
if (grepl("human", archr_dir, ignore.case=TRUE)) {
    library(BSgenome.Hsapiens.UCSC.hg38)
    genome <- "hg38"
} else if (grepl("mouse", archr_dir, ignore.case=TRUE)) {
    library(BSgenome.Mmusculus.UCSC.mm10)
    genome <- "mm10"
} else {
    print("Could not determine species from ARCHR_DIR name. Default to human")
    library(BSgenome.Hsapiens.UCSC.hg38)
    genome <- "hg38"
}


# Set the number of threads for ArchR
addArchRThreads(threads = 9)

# Load the project
print(paste("Loading ArchR project from", archr_dir))
proj <- loadArchRProject(path = archr_dir)

# Load Seurat object metadata
if (!file.exists(metadata_path)) {
    stop(paste("Seurat Metadata file does not exist:", metadata_path))
}
seurat_metadata <- readr::read_csv(metadata_path)
print(paste("Seurat metadata loaded from", metadata_path))

# --- BEFORE summary ---
print_archr_summary(proj, paste("BEFORE filtering —", archr_dir))
cat(sprintf("  %-30s %d\n", "Cells in Seurat metadata:", nrow(seurat_metadata)))
cat(sprintf("  %-30s %d\n\n", "Cells in ArchR project:",  length(getCellNames(proj))))

# check if columns exist in seurat_metadata
required_columns <- c("seurat_atac_barcode", "seurat_gex_barcode")
missing_columns <- setdiff(required_columns, colnames(seurat_metadata))
if (length(missing_columns) > 0) {
    stop(paste("Missing required columns in Seurat metadata:", paste(missing_columns, collapse = ", ")))
} else {
    print("All required columns are present in Seurat metadata.")
}

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

n_archr      <- length(cells)
n_seurat     <- length(seurat_rna_barcodes)
n_overlap    <- length(cellsToKeep)
n_archr_only <- n_archr - n_overlap
n_seurat_only <- n_seurat - n_overlap
cat(sprintf("  %-30s %d\n",   "Cells in ArchR:",           n_archr))
cat(sprintf("  %-30s %d\n",   "Cells in Seurat metadata:", n_seurat))
cat(sprintf("  %-30s %d\n",   "Overlapping cells (kept):", n_overlap))
cat(sprintf("  %-30s %d\n",   "ArchR-only (dropped):",     n_archr_only))
cat(sprintf("  %-30s %d\n\n", "Seurat-only (not in ArchR):", n_seurat_only))

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

# Resolve which columns to add to ArchR metadata:
# Use columns passed via --columns if provided, otherwise fall back to defaults
if (nchar(columns_env) > 0) {
    cols_to_add <- trimws(strsplit(columns_env, ",")[[1]])
    print("Using user-supplied columns.")
} else {
    cols_to_add <- c("hybrid_pair","neftel_4_state",
    "Region.annotation", "mislabelled",
    "Azimuth_class", "Azimuth_subclass",
    "Ambiguous", "PIMO_status", "PIMO_up_status") # add more column names as needed
    print("No columns specified, using defaults.")
}

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

# --- AFTER summary ---
print_archr_summary(projSubset, paste("AFTER filtering & metadata transfer —", output_dir))

# Per-column value breakdowns for all added metadata columns
cat(strrep("-", 60), "\n", sep = "")
cat(" Metadata column distributions (added columns):\n")
cat(strrep("-", 60), "\n", sep = "")
cd_after <- as.data.frame(projSubset@cellColData)
for (col in cols_to_add) {
    if (col %in% colnames(cd_after)) {
        col_vals <- cd_after[[col]]
        cat(sprintf("\n  %s:\n", col))
        if (is.numeric(col_vals)) {
            cat(sprintf("    mean=%.4f, median=%.4f, min=%.4f, max=%.4f, NAs=%d\n",
                mean(col_vals, na.rm=TRUE), median(col_vals, na.rm=TRUE),
                min(col_vals,  na.rm=TRUE), max(col_vals,   na.rm=TRUE),
                sum(is.na(col_vals))))
        } else {
            tbl <- sort(table(col_vals, useNA = "ifany"), decreasing = TRUE)
            for (v in names(tbl)) {
                pct <- 100 * tbl[[v]] / length(col_vals)
                cat(sprintf("    %-30s %5d  (%5.1f%%)\n", v, tbl[[v]], pct))
            }
        }
    }
}
cat(strrep("-", 60), "\n\n", sep = "")

# Save the project
print("Saving the project")
saveArchRProject(ArchRProj = projSubset, outputDirectory = output_dir, load = TRUE)

print("Script completed!")
EOF