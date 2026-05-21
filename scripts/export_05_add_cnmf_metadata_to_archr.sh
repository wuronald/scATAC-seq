#!/bin/bash
#SBATCH --job-name=export_05_add_cnmf_metadata_to_archr
#SBATCH --output=export_05_add_cnmf_metadata_to_archr_%j.out
#SBATCH --error=export_05_add_cnmf_metadata_to_archr_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=60G
#SBATCH --time=01:00:00

# Add cNMF metadata columns to an ArchR project and save the result to a new
# directory. No subsetting is performed — all cells in the original project are
# retained. The cNMF_program column is renamed to match the supplied groupBy
# name (e.g. "k3_program") and values are prefixed with "program_" so that ArchR
# does not receive purely numeric group names.
#
# Arguments
# ---------
#   1  groupBy          Column name to use for the cNMF program (e.g. "k3_program")
#   2  base_archr_dir   Path to the existing (source) ArchR project directory
#   3  output_archr_dir Path to the new ArchR project directory that will be created
#   4  metadata_csv     Path to the cNMF metadata CSV (must contain a "cell_id" column)
#
# Example usage
# -------------
  # sbatch scripts/export_05_add_cnmf_metadata_to_archr.sh \
  #     k3_program \
  #     human_multiome_harmony_merged_malig_peak \
  #     human_multiome_hmmp_k3_program_all_cells \
  #     data/metadata/human_multiome/cNMF_cell_program_assignment_with_usages_K3.csv

# ---------------------------------------------------------------------------
# Parse arguments
# ---------------------------------------------------------------------------
GROUP_BY="${1}"
BASE_ARCHR_PROJECT="${2}"
OUTPUT_ARCHR_PROJECT="${3}"
METADATA_CSV="${4}"

# ---------------------------------------------------------------------------
# Validate required arguments
# ---------------------------------------------------------------------------
if [ -z "${GROUP_BY}" ] || [ -z "${BASE_ARCHR_PROJECT}" ] || \
   [ -z "${OUTPUT_ARCHR_PROJECT}" ] || [ -z "${METADATA_CSV}" ]; then
    echo "ERROR: All four arguments are required."
    echo ""
    echo "Usage:"
    echo "  sbatch $0 <groupBy> <base_archr_dir> <output_archr_dir> <metadata_csv>"
    echo ""
    echo "Example:"
    echo "  sbatch $0 k3_program \\"
    echo "      human_multiome_harmony_merged_malig_peak \\"
    echo "      human_multiome_hmmp_k3_program_all_cells \\"
    echo "      data/metadata/human_multiome/cNMF_cell_program_assignment_with_usages_K3.csv"
    exit 1
fi

# ---------------------------------------------------------------------------
# Export parameters so R (inline heredoc) can read them
# ---------------------------------------------------------------------------
export GROUP_BY
export BASE_ARCHR_PROJECT
export OUTPUT_ARCHR_PROJECT
export METADATA_CSV

echo "=========================================="
echo "Adding cNMF metadata to ArchR project"
echo "  (no cell subsetting)"
echo "=========================================="
echo "  groupBy:              ${GROUP_BY}"
echo "  Base ArchR project:   ${BASE_ARCHR_PROJECT}"
echo "  Output ArchR project: ${OUTPUT_ARCHR_PROJECT}"
echo "  Metadata CSV:         ${METADATA_CSV}"
echo "=========================================="
echo ""

# ---------------------------------------------------------------------------
# Load R module
# ---------------------------------------------------------------------------
module load R/4.4.1

# ---------------------------------------------------------------------------
# Run R inline
# ---------------------------------------------------------------------------
Rscript - <<'EOF'

library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(here)
set.seed(1)

# ---- Read parameters from environment ----
groupBy           <- Sys.getenv("GROUP_BY")
baseArchRProject  <- Sys.getenv("BASE_ARCHR_PROJECT")
outputArchRProject <- Sys.getenv("OUTPUT_ARCHR_PROJECT")
metadataCsv       <- Sys.getenv("METADATA_CSV")

cat("groupBy:              ", groupBy, "\n")
cat("Base ArchR project:   ", baseArchRProject, "\n")
cat("Output ArchR project: ", outputArchRProject, "\n")
cat("Metadata CSV:         ", metadataCsv, "\n\n")

# ---- ArchR setup ----
addArchRThreads(threads = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "8")))
addArchRGenome("hg38")

# ---- Load the base ArchR project ----
cat("Loading base ArchR project...\n")
proj <- loadArchRProject(path = baseArchRProject)
cat(sprintf("Loaded project with %d cells.\n\n", nCells(proj)))

# ---- Copy project to output directory before modifying ----
# Saving first ensures a clean, isolated copy and avoids mutating the original.
cat("Saving copy of project to output directory...\n")
saveArchRProject(proj, outputDirectory = outputArchRProject, load = FALSE)
cat("Reloading project from output directory...\n")
proj <- loadArchRProject(path = outputArchRProject)
cat(sprintf("Reloaded project with %d cells.\n\n", nCells(proj)))

# ---- Validate and read metadata CSV ----
cat("========================================\n")
cat("Adding external metadata from CSV\n")
cat("========================================\n")
cat("CSV path:", metadataCsv, "\n")

if (!file.exists(metadataCsv)) {
  stop(paste("Metadata CSV not found:", metadataCsv))
}

meta <- read.csv(metadataCsv, stringsAsFactors = FALSE)
cat(sprintf("Loaded metadata: %d rows, %d columns\n", nrow(meta), ncol(meta)))
cat("Columns:", paste(colnames(meta), collapse = ", "), "\n\n")

if (!"cell_id" %in% colnames(meta)) {
  stop("Metadata CSV must contain a 'cell_id' column matching ArchR cell names.")
}

# ---- Reformat cell_id: CSV uses '____', ArchR uses '#' ----
# CSV format:   "Zadeh__C0736__5117____BARCODE-1"
# ArchR format: "Zadeh__C0736__5117#BARCODE-1"
meta$cell_id <- sub("____", "#", meta$cell_id)
cat("Reformatted cell_id: replaced '____' with '#' to match ArchR cell names.\n")
cat("Example reformatted cell_id:", head(meta$cell_id, 1), "\n\n")

# ---- Rename cNMF_program -> groupBy and prefix values with "program_" ----
# Purely numeric group names cause problems in ArchR (e.g. for pseudobulking).
if ("cNMF_program" %in% colnames(meta)) {
  meta[[groupBy]] <- paste0("program_", meta$cNMF_program)
  cat(sprintf("Renamed 'cNMF_program' -> '%s' (values prefixed with 'program_').\n\n", groupBy))
}

# ---- Match metadata rows to ArchR cells ----
archrCells <- rownames(proj@cellColData)
matchIdx   <- match(archrCells, meta$cell_id)

nMatched   <- sum(!is.na(matchIdx))
nTotal     <- length(archrCells)
cat(sprintf("Matched %d / %d ArchR cells to metadata rows.\n", nMatched, nTotal))
if (nMatched == 0) {
  stop("No ArchR cells matched the metadata cell_id values. Check cell_id formatting.")
}
cat(sprintf("  (%d cells will have NA for metadata columns — they are kept in the project.)\n\n",
            nTotal - nMatched))

# ---- Add all metadata columns (except cell_id) to the ArchR project ----
colsToAdd <- setdiff(colnames(meta), "cell_id")
cat(sprintf("Adding %d metadata column(s): %s\n\n", length(colsToAdd),
            paste(colsToAdd, collapse = ", ")))

for (col in colsToAdd) {
  values <- meta[[col]][matchIdx]  # NA for unmatched cells (preserved in project)
  proj <- addCellColData(
    ArchRProj = proj,
    data      = values,
    name      = col,
    cells     = archrCells,
    force     = TRUE
  )
  cat(sprintf("  ✓ Added column '%s'\n", col))
}
cat("\n")

# ---- Summary of program assignments ----
if (groupBy %in% colnames(proj@cellColData)) {
  progTable <- table(proj@cellColData[[groupBy]], useNA = "ifany")
  cat(sprintf("'%s' distribution across all cells (no subsetting):\n", groupBy))
  for (i in seq_along(progTable)) {
    label <- if (is.na(names(progTable)[i])) "<NA>" else names(progTable)[i]
    cat(sprintf("  %-20s : %d cells\n", label, progTable[[i]]))
  }
  cat("\n")
}

# ---- Build combined_program column (coalesce: k3_program over PIMO_up_status) ----
# Logic: start with PIMO_up_status; overwrite with groupBy (e.g. k3_program) wherever
# groupBy is non-NA. This means cNMF-assigned cells are labelled by their program,
# and all remaining cells retain their PIMO_up_status label.
pimoCol    <- "PIMO_up_status"
combinedCol <- "combined_program"

if (pimoCol %in% colnames(proj@cellColData) && groupBy %in% colnames(proj@cellColData)) {
  cat("========================================\n")
  cat(sprintf("Building '%s' column\n", combinedCol))
  cat(sprintf("  Base layer : '%s'\n", pimoCol))
  cat(sprintf("  Override   : '%s' (where non-NA)\n", groupBy))
  cat("========================================\n")

  pimoVals    <- as.character(proj@cellColData[[pimoCol]])
  programVals <- as.character(proj@cellColData[[groupBy]])

  # Coalesce: use programVals where it is not NA, otherwise fall back to pimoVals
  combinedVals <- ifelse(!is.na(programVals), programVals, pimoVals)

  proj <- addCellColData(
    ArchRProj = proj,
    data      = combinedVals,
    name      = combinedCol,
    cells     = archrCells,
    force     = TRUE
  )
  cat(sprintf("  ✓ Added column '%s'\n\n", combinedCol))

  # Summary
  combTable <- table(proj@cellColData[[combinedCol]], useNA = "ifany")
  cat(sprintf("'%s' distribution:\n", combinedCol))
  for (i in seq_along(combTable)) {
    label <- if (is.na(names(combTable)[i])) "<NA>" else names(combTable)[i]
    cat(sprintf("  %-20s : %d cells\n", label, combTable[[i]]))
  }
  cat("\n")
} else {
  missing <- c(
    if (!pimoCol  %in% colnames(proj@cellColData)) pimoCol,
    if (!groupBy  %in% colnames(proj@cellColData)) groupBy
  )
  cat(sprintf(
    "⚠ Skipping '%s': required column(s) not present in project: %s\n\n",
    combinedCol, paste(missing, collapse = ", ")
  ))
}

# ---- Save the final project in-place ----
cat("Saving final ArchR project to:", outputArchRProject, "\n")
saveArchRProject(proj, outputDirectory = outputArchRProject, load = FALSE)

cat("========================================\n")
cat(sprintf("Done. Project saved to: %s\n", outputArchRProject))
cat(sprintf("Total cells retained: %d (no subsetting performed)\n", nCells(proj)))
cat("========================================\n")

EOF

echo ""
echo "=========================================="
echo "COMPLETE"
echo "=========================================="
echo "  Output ArchR project: ${OUTPUT_ARCHR_PROJECT}"
echo "=========================================="