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

# ==============================================================================
# Purpose: SLURM batch script to append cNMF program assignment and usage metadata
#          to an ArchR project across all cells (no subsetting) and construct a
#          coalesced 'combined_program' annotation column for downstream analyses
#          (e.g., chromVAR motif deviations, group comparisons, browser tracks).
#
# Workflow:
#   1. Loads the source/base ArchR project without modifying original files.
#   2. Detects reference genome (hg38 for human, mm10 for mouse).
#   3. Clones/saves the project into a dedicated output directory.
#   4. Ingests the cNMF metadata CSV, reformats cell IDs ('____' -> '#').
#   5. Safely prefixes numeric program identifiers with 'program_' while strictly
#      preserving NA values (avoiding 'program_NA' string conversion).
#   6. Matches ArchR cell barcodes and appends program assignment and usage columns
#      via addCellColData().
#   7. Constructs the 'combined_program' column by coalescing:
#        - Base layer: PIMO_up_status ('PIMOdown', 'PIMOinter', 'PIMOup')
#        - Overwrite: cNMF program ('program_1', 'program_2', etc.) where non-NA.
#   8. Saves the final metadata-augmented ArchR project.
#
# Usage:
#   sbatch scripts/export_05_add_cnmf_metadata_to_archr.sh \
#       [GROUP_BY] [BASE_ARCHR_PROJECT] [OUTPUT_ARCHR_PROJECT] [METADATA_CSV]
#
# Arguments (all positional, optional with sensible defaults):
#   $1  GROUP_BY              Column name for the cNMF program (default: k3_program)
#   $2  BASE_ARCHR_PROJECT    Source ArchR project directory
#                             (default: human_multiome_harmony_merged_malig_peak)
#   $3  OUTPUT_ARCHR_PROJECT  Destination directory for new ArchR project
#                             (default: human_multiome_hmmp_${GROUP_BY}_all_cells)
#   $4  METADATA_CSV          Path to cNMF cell metadata CSV (must have cell_id)
#                             (default: data/metadata/human_multiome/cNMF_cell_program_assignment_with_usages_K3.csv)
#
# Examples:
#   # Standard K=3 cNMF metadata transfer:
#   sbatch scripts/export_05_add_cnmf_metadata_to_archr.sh
#
#   # Explicit K=3 invocation:
#   sbatch scripts/export_05_add_cnmf_metadata_to_archr.sh \
#       k3_program \
#       human_multiome_harmony_merged_malig_peak \
#       human_multiome_hmmp_k3_program_all_cells \
#       data/metadata/human_multiome/cNMF_cell_program_assignment_with_usages_K3.csv
#
#   # K=9 cNMF metadata transfer:
#   sbatch scripts/export_05_add_cnmf_metadata_to_archr.sh \
#       k9_program \
#       human_multiome_harmony_merged_malig_peak \
#       human_multiome_hmmp_k9_program_all_cells \
#       data/metadata/human_multiome/cNMF_cell_program_assignment_with_usages_K9.csv
# ==============================================================================

# ---------------------------------------------------------------------------
# Parse arguments with default values
# ---------------------------------------------------------------------------
GROUP_BY="${1:-k3_program}"
BASE_ARCHR_PROJECT="${2:-human_multiome_harmony_merged_malig_peak}"
OUTPUT_ARCHR_PROJECT="${3:-human_multiome_hmmp_${GROUP_BY}_all_cells}"
[ -z "${OUTPUT_ARCHR_PROJECT}" ] && OUTPUT_ARCHR_PROJECT="human_multiome_hmmp_${GROUP_BY}_all_cells"

METADATA_CSV="${4:-data/metadata/human_multiome/cNMF_cell_program_assignment_with_usages_K3.csv}"
[ -z "${METADATA_CSV}" ] && METADATA_CSV="data/metadata/human_multiome/cNMF_cell_program_assignment_with_usages_K3.csv"

# ---------------------------------------------------------------------------
# Export parameters so R (inline heredoc) can read them
# ---------------------------------------------------------------------------
export GROUP_BY
export BASE_ARCHR_PROJECT
export OUTPUT_ARCHR_PROJECT
export METADATA_CSV

echo "=========================================="
echo "Adding cNMF metadata to ArchR project"
echo "  (all cells retained, no subsetting)"
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
library(here)
set.seed(1)

# ---- Read parameters from environment ----
groupBy            <- Sys.getenv("GROUP_BY",            unset = "k3_program")
baseArchRProject   <- Sys.getenv("BASE_ARCHR_PROJECT",   unset = "human_multiome_harmony_merged_malig_peak")
outputArchRProject <- Sys.getenv("OUTPUT_ARCHR_PROJECT", unset = "human_multiome_hmmp_k3_program_all_cells")
metadataCsv        <- Sys.getenv("METADATA_CSV",         unset = "data/metadata/human_multiome/cNMF_cell_program_assignment_with_usages_K3.csv")

cat("groupBy:              ", groupBy, "\n")
cat("Base ArchR project:   ", baseArchRProject, "\n")
cat("Output ArchR project: ", outputArchRProject, "\n")
cat("Metadata CSV:         ", metadataCsv, "\n\n")

# ---- Resolve paths ----
baseProjectPath   <- if (grepl("^/", baseArchRProject))   baseArchRProject   else here(baseArchRProject)
outputProjectPath <- if (grepl("^/", outputArchRProject)) outputArchRProject else here(outputArchRProject)
metadataCsvPath   <- if (grepl("^/", metadataCsv))        metadataCsv        else here(metadataCsv)

# ---- ArchR setup ----
addArchRThreads(threads = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "8")))

# ---- Dynamic species / genome detection ----
if (grepl("mouse", baseArchRProject, ignore.case = TRUE)) {
  library(BSgenome.Mmusculus.UCSC.mm10)
  genome <- "mm10"
} else {
  library(BSgenome.Hsapiens.UCSC.hg38)
  genome <- "hg38"
}
cat(sprintf("Configuring ArchR reference genome: %s\n", genome))
addArchRGenome(genome)

# ---- Validate input existence ----
if (!file.exists(baseProjectPath)) {
  stop(paste("Base ArchR project directory not found:", baseProjectPath))
}
if (!file.exists(metadataCsvPath)) {
  stop(paste("Metadata CSV not found:", metadataCsvPath))
}

# ---- Load the base ArchR project ----
cat("Loading base ArchR project from:", baseProjectPath, "\n")
proj <- loadArchRProject(path = baseProjectPath)
cat(sprintf("Loaded project with %d cells.\n\n", nCells(proj)))

# ---- Copy project to output directory before modifying ----
# Saving first ensures a clean, isolated copy and avoids mutating the original.
cat("Saving copy of project to output directory:", outputProjectPath, "\n")
saveArchRProject(proj, outputDirectory = outputProjectPath, load = FALSE)
cat("Reloading project from output directory...\n")
proj <- loadArchRProject(path = outputProjectPath)
cat(sprintf("Reloaded project with %d cells.\n\n", nCells(proj)))

# ---- Validate and read metadata CSV ----
cat("========================================\n")
cat("Adding external metadata from CSV\n")
cat("========================================\n")
cat("CSV path:", metadataCsvPath, "\n")

meta <- read.csv(metadataCsvPath, stringsAsFactors = FALSE)
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

# ---- Rename cNMF_program -> groupBy and prefix values with 'program_' ----
# Safely handle NA values so that paste0 does not turn NA into "program_NA" string.
if ("cNMF_program" %in% colnames(meta)) {
  progVals <- as.character(meta$cNMF_program)
  progValsClean <- ifelse(is.na(progVals) | progVals %in% c("NA", ""), NA_character_,
                          ifelse(grepl("^program_", progVals), progVals, paste0("program_", progVals)))
  meta[[groupBy]] <- progValsClean
  cat(sprintf("Renamed 'cNMF_program' -> '%s' (values prefixed with 'program_', preserving NA).\n\n", groupBy))
} else if (groupBy %in% colnames(meta)) {
  progVals <- as.character(meta[[groupBy]])
  progValsClean <- ifelse(is.na(progVals) | progVals %in% c("NA", ""), NA_character_,
                          ifelse(grepl("^program_", progVals), progVals, paste0("program_", progVals)))
  meta[[groupBy]] <- progValsClean
  cat(sprintf("Ensured values in '%s' are prefixed with 'program_' (preserving NA).\n\n", groupBy))
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

# ---- Build combined_program column (coalesce: groupBy over PIMO_up_status) ----
# Logic: start with PIMO_up_status; overwrite with groupBy (e.g. k3_program) wherever
# groupBy is non-NA and not 'program_NA'. This means cNMF-assigned cells are labelled
# by their program, and all remaining cells retain their PIMO_up_status label.
pimoCol     <- "PIMO_up_status"
combinedCol <- "combined_program"

if (pimoCol %in% colnames(proj@cellColData) && groupBy %in% colnames(proj@cellColData)) {
  cat("========================================\n")
  cat(sprintf("Building '%s' column\n", combinedCol))
  cat(sprintf("  Base layer : '%s'\n", pimoCol))
  cat(sprintf("  Override   : '%s' (where non-NA)\n", groupBy))
  cat("========================================\n")

  pimoVals    <- as.character(proj@cellColData[[pimoCol]])
  programVals <- as.character(proj@cellColData[[groupBy]])

  # Coalesce: use programVals where it is not NA/program_NA, otherwise fall back to pimoVals
  isProgramAssigned <- !is.na(programVals) & !(programVals %in% c("NA", "program_NA", ""))
  combinedVals      <- ifelse(isProgramAssigned, programVals, pimoVals)

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
cat("Saving final ArchR project to:", outputProjectPath, "\n")
saveArchRProject(proj, outputDirectory = outputProjectPath, load = FALSE)

cat("========================================\n")
cat(sprintf("Done. Project saved to: %s\n", outputProjectPath))
cat(sprintf("Total cells retained: %d (no subsetting performed)\n", nCells(proj)))
cat("========================================\n")

EOF

echo ""
echo "=========================================="
echo "COMPLETE"
echo "=========================================="
echo "  Output ArchR project: ${OUTPUT_ARCHR_PROJECT}"
echo "=========================================="