#!/bin/bash
#SBATCH --job-name=export_04_pseudobulk_bigwig_heatmaps
#SBATCH --output=export_04_pseudobulk_bigwig_heatmaps_%j.out
#SBATCH --error=export_04_pseudobulk_bigwig_heatmaps_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=60G
#SBATCH --time=01:00:00

# Parse command line arguments
# First argument: groupBy variable (e.g., "PIMO_up_status", "hybrid_pair", "k3_program")
# Second argument: path to directory containing BED files
# Third argument (optional): pattern to match BED files (e.g., "PIMO_MP", "Nomura_")
# Fourth argument (optional): "force" to re-export bigWig files even if they exist
# Fifth argument (optional): path to a metadata CSV to add to the ArchR object before export.
#                            The CSV must contain a "cell_id" column matching ArchR cell names.
# Sixth argument (optional): BASE_ARCHR_PROJECT path (defaults to "human_multiome_harmony_merged_malig_peak")
# Seventh argument (optional): regular expression to exclude certain bigwig groups (e.g., "(PIMOinterTC|PIMOinterEB|PIMOinterPT)")

# Example usage:
# sbatch scripts/export_04_pseudobulk_bigwig_heatmaps_v2.sh PIMO_up_status data/Signatures/PIMO_metaprogram/ PIMO_MP
# sbatch scripts/export_04_pseudobulk_bigwig_heatmaps_v2.sh k3_program data/Signatures/PIMO_metaprogram PIMO_MP no data/metadata/human_multiome/cNMF_cell_program_assignment_with_usages_K3.csv human_multiome_harmony_merged_malig_peak "(PIMOinterTC|PIMOinterEB|PIMOinterPT)"

GROUP_BY="${1:-PIMO_up_status}"  # Default to "PIMO_up_status" if no argument provided
BED_DIR="${2:-data/Signatures}"  # Default BED directory
BED_PATTERN="${3:-PIMO_MP}"  # Default to PIMO Metaprograms
FORCE_EXPORT="${4:-no}"  # Default to not forcing re-export
METADATA_CSV="${5:-}"  # Optional path to metadata CSV
BASE_ARCHR_PROJECT="${6:-human_multiome_harmony_merged_malig_peak}"  # Optional base project path
EXCLUDE_PATTERN="${7:-}"  # Optional regex pattern to exclude certain bigwig groups

# If a metadata CSV is provided but no groupBy was given, warn and exit
if [ -n "${METADATA_CSV}" ] && [ -z "${1}" ]; then
    echo "ERROR: When providing a metadata CSV (arg 5), you must also specify a groupBy"
    echo "       value as the first argument (e.g. k3_program, k9_program)."
    echo "       This should match the K value of the cNMF solution in the CSV."
    exit 1
fi

# Determine working ArchR project directory
if [ -n "${METADATA_CSV}" ]; then
    BASE_NAME=$(basename "${BASE_ARCHR_PROJECT}")
    ARCHR_PROJECT_DIR="${BASE_NAME}_subset_${GROUP_BY}"
else
    ARCHR_PROJECT_DIR="${BASE_ARCHR_PROJECT}"
fi

# Export the parameters so R can access them
export GROUP_BY
export BED_DIR
export BED_PATTERN
export FORCE_EXPORT
export METADATA_CSV
export BASE_ARCHR_PROJECT
export ARCHR_PROJECT_DIR

echo "=========================================="
echo "Running pseudo-bulk bigWig export and heatmap generation (v2)"
echo "=========================================="
echo "  groupBy: ${GROUP_BY}"
echo "  BED directory: ${BED_DIR}"
echo "  BED file pattern: ${BED_PATTERN}"
echo "  Force re-export: ${FORCE_EXPORT}"
echo "  Metadata CSV: ${METADATA_CSV:-<none>}"
echo "  Base ArchR project: ${BASE_ARCHR_PROJECT}"
echo "  Working ArchR project: ${ARCHR_PROJECT_DIR}"
echo "  Exclude pattern: ${EXCLUDE_PATTERN:-<none>}"
echo "=========================================="
echo ""

# Load necessary modules
module load R/4.4.1

########################################################
# PART 1: Export pseudo-bulked data to bigWig files
########################################################

echo "PART 1: Exporting pseudo-bulked bigWig files using ArchR"
echo "=========================================="

# Run R script to export bigWig files
Rscript - <<'EOF'

# Load libraries
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(here)
set.seed(1)

# Get parameters from environment variables
groupBy          <- Sys.getenv("GROUP_BY",           unset = "PIMO_up_status")
forceExport      <- Sys.getenv("FORCE_EXPORT",       unset = "no")
baseArchRProject <- Sys.getenv("BASE_ARCHR_PROJECT", unset = "human_multiome_harmony_merged_malig_peak")
archrProjectDir  <- Sys.getenv("ARCHR_PROJECT_DIR",  unset = "human_multiome_harmony_merged_malig_peak")

cat("Using groupBy:", groupBy, "\n")
cat("Force re-export:", forceExport, "\n")
cat("Base ArchR project:", baseArchRProject, "\n")
cat("Working ArchR project:", archrProjectDir, "\n")

# Set the number of threads for ArchR
addArchRThreads(threads = 18)

# Set genome to hg38
addArchRGenome("hg38")

# Load the original (unmodified) ArchR project
cat("Loading ArchR project...\n")
proj <- loadArchRProject(path = baseArchRProject)

print("InitiaL SAVE before subsetting...")
saveArchRProject(proj, outputDirectory = archrProjectDir, load = FALSE)
print("InitiaL load before subsetting...")
proj <- loadArchRProject(path = archrProjectDir)
print("Passed load before subsetting...")

########################################################
# Optional: Add external metadata and subset cells
########################################################

metadataCsv <- Sys.getenv("METADATA_CSV", unset = "")

if (nchar(metadataCsv) > 0) {
  cat("\n========================================\n")
  cat("Adding external metadata from CSV\n")
  cat("========================================\n")
  cat("CSV path:", metadataCsv, "\n")
  
  if (!file.exists(metadataCsv)) {
    stop(paste("Metadata CSV not found:", metadataCsv))
  }
  
  # Read the metadata CSV
  meta <- read.csv(metadataCsv, stringsAsFactors = FALSE)
  cat(sprintf("Loaded metadata: %d rows, %d columns\n", nrow(meta), ncol(meta)))
  cat("Columns:", paste(colnames(meta), collapse = ", "), "\n")
  
  if (!"cell_id" %in% colnames(meta)) {
    stop("Metadata CSV must contain a 'cell_id' column matching ArchR cell names.")
  }
  
  # Reformat cell_id from CSV format to ArchR format
  meta$cell_id <- sub("____", "#", meta$cell_id)
  cat("Reformatted cell_id: replaced '____' separator with '#' to match ArchR cell names\n")
  cat("Example reformatted cell_id:", head(meta$cell_id, 1), "\n")
  
  # Rename cNMF_program -> groupBy column name
  if ("cNMF_program" %in% colnames(meta)) {
    meta[[groupBy]] <- paste0("program_", meta$cNMF_program)
    cat(sprintf("Renamed 'cNMF_program' -> '%s' (values prefixed with 'program_')\n", groupBy))
  }
  
  # Match metadata to ArchR cells
  archrCells <- rownames(proj@cellColData)
  matchIdx   <- match(archrCells, meta$cell_id)
  
  nMatched <- sum(!is.na(matchIdx))
  cat(sprintf("Matched %d / %d ArchR cells to metadata\n", nMatched, length(archrCells)))
  
  # Add all metadata columns (except cell_id) to the ArchR project
  colsToAdd <- setdiff(colnames(meta), "cell_id")

  for (col in colsToAdd) {
    values <- meta[[col]][matchIdx]
    proj <- addCellColData(
      ArchRProj = proj,
      data      = values,
      name      = col,
      cells     = archrCells,
      force     = TRUE
    )
    cat(sprintf("  ✓ Added column '%s'\n", col))
  }

  # Subset to cells that have a non-NA cNMF program assignment
  if (groupBy %in% colsToAdd) {
    cellsBefore <- nCells(proj)
    cellsWithProgram <- archrCells[!is.na(meta[[groupBy]][matchIdx])]
    print(paste("Number of cells with program:", length(cellsWithProgram)))
    print(head(cellsWithProgram))
    cells = getCellNames(proj)
    print("Example ArchR cell names:")
    print(head(cells))

    if (length(cellsWithProgram) == 0) {
      stop(sprintf(
        "No cells with a non-NA '%s' found after metadata join. Check that cell_id values match ArchR cell names.",
        groupBy
      ))
    }
    
    tryCatch({
      proj <- subsetArchRProject(
        ArchRProj       = proj,
        cells           = cellsWithProgram,
        dropCells       = TRUE,
        outputDirectory = archrProjectDir,
        force           = TRUE
      )
    }, error = function(e) {
      arrowFiles <- list.files(file.path(archrProjectDir, "ArrowFiles"), pattern = "\\.arrow$")
      if (length(arrowFiles) == 0) {
        stop(sprintf("subsetArchRProject failed and no Arrow files found in %s/ArrowFiles: %s",
                     archrProjectDir, e$message))
      }
      cat(sprintf("⚠ Caught expected ArchR reload bug after subsetting (Arrow files present on disk). Reloading manually.\n"))
    })
    
    cat(sprintf("\n✓ Subsetted ArchR project: %d → %d cells (kept cells with non-NA %s)\n",
                cellsBefore, nCells(proj), groupBy))

    # Reload the subset project from disk after pause
    Sys.sleep(2)
    print("Reloading subsetted ArchR project from disk after pause...")
    proj <- loadArchRProject(path = archrProjectDir)
       
    # Summary of program assignments
    progTable <- table(proj@cellColData[[groupBy]])
    cat(sprintf("\n%s distribution after subsetting:\n", groupBy))
    for (prog in names(progTable)) {
      cat(sprintf("  %s: %d cells\n", prog, progTable[[prog]]))
    }
  }
  
  cat("========================================\n\n")
} else {
  cat("No metadata CSV provided. Skipping metadata addition and cell subsetting.\n\n")
}

# Create output directory for bigWig files
filePrefix <- gsub("_", "", groupBy)
outDir <- here(paste0(archrProjectDir, "/pseudobulk_bigwig_", groupBy, "/"))
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

cat("Output directory:", outDir, "\n")

# Get the unique groups
groups <- unique(proj@cellColData[[groupBy]])
groups <- groups[!is.na(groups)]
cat("Groups found:", paste(groups, collapse = ", "), "\n")

########################################################
# Export pseudo-bulk bigWig files for each group
########################################################

cat("\n========================================\n")
cat("Exporting pseudo-bulk bigWig files\n")
cat("========================================\n\n")

# Check if bigWig files already exist
bigwigDir <- file.path(proj@projectMetadata$outputDirectory, "GroupBigWigs", groupBy)
cat("BigWig output directory:", bigwigDir, "\n")

# Force export if requested
if (tolower(forceExport) == "force" || tolower(forceExport) == "yes") {
  cat("\n⚠ Force re-export enabled. Will regenerate all bigWig files.\n\n")
  allGroupsHaveBigWigs <- FALSE
} else {
  # Check for existing bigWig files
  existingBigWigs <- list.files(bigwigDir, pattern = "\\.bw$", full.names = FALSE)
  
  if (length(existingBigWigs) > 0) {
    cat("\nFound existing bigWig files:\n")
    for (bw in existingBigWigs) {
      cat(sprintf("  - %s\n", bw))
    }
    cat("\nChecking if all groups have bigWig files...\n")
    
    # Check if we have a bigWig for each group
    allGroupsHaveBigWigs <- TRUE
    for (group in groups) {
      expectedFile <- paste0(group, "-TileSize-100-normMethod-ReadsInTSS-ArchR.bw")
      if (!file.exists(file.path(bigwigDir, expectedFile))) {
        cat(sprintf("  Missing bigWig for group: %s\n", group))
        allGroupsHaveBigWigs <- FALSE
      } else {
        cat(sprintf("  ✓ Found bigWig for group: %s\n", group))
      }
    }
    
    if (allGroupsHaveBigWigs) {
      cat("\n✓ All bigWig files already exist. Skipping export.\n")
    } else {
      cat("\nSome bigWig files are missing. Exporting all groups...\n\n")
      allGroupsHaveBigWigs <- FALSE
    }
  } else {
    cat("\nNo existing bigWig files found. Will export all groups.\n\n")
    allGroupsHaveBigWigs <- FALSE
  }
}

# Only export if we need to
if (!allGroupsHaveBigWigs) {
  cat("Exporting bigWig files for all groups in", groupBy, "...\n")
  cat("Groups to be exported:\n")
  for (group in groups) {
    cellsInGroup <- sum(proj@cellColData[[groupBy]] == group, na.rm = TRUE)
    cat(sprintf("  - %s (%d cells)\n", group, cellsInGroup))
  }
  cat("\n")
  
  tryCatch({
    getGroupBW(
      ArchRProj = proj,
      groupBy = groupBy,
      normMethod = "ReadsInTSS",
      tileSize = 100,
      maxCells = 500,
      ceiling = 4,
      verbose = TRUE,
      threads = getArchRThreads()
    )
    
    cat("\n✓ Successfully exported bigWig files for all groups\n")
  }, error = function(e) {
    cat(sprintf("\n✗ Error exporting bigWig files: %s\n", e$message))
    stop("BigWig export failed")
  })
  
} else {
  cat("\nSkipping bigWig export (files already exist)\n")
}

cat("========================================\n")
if (!allGroupsHaveBigWigs) {
  cat("BigWig export complete!\n")
} else {
  cat("BigWig files already exist (skipped export)\n")
}
cat("========================================\n\n")

bigwigDir <- file.path(proj@projectMetadata$outputDirectory, "GroupBigWigs", groupBy)
cat("BigWig files location:", bigwigDir, "\n")

# List final bigWig files
finalBigWigs <- list.files(bigwigDir, pattern = "\\.bw$", full.names = FALSE)
cat(sprintf("\nTotal bigWig files available: %d\n", length(finalBigWigs)))
for (bw in finalBigWigs) {
  bwPath <- file.path(bigwigDir, bw)
  bwSize <- file.info(bwPath)$size / (1024^3)
  cat(sprintf("  - %s (%.2f GB)\n", bw, bwSize))
}

# Save the bigwig directory path for the bash script
cat(bigwigDir, file = file.path(outDir, "bigwig_dir.txt"))

EOF

echo ""
echo "PART 1 COMPLETE: BigWig files exported"
echo "=========================================="
echo ""

########################################################
# PART 2: Run deeptools computeMatrix and plotHeatmap
########################################################

echo "PART 2: Running deeptools to create heatmaps"
echo "=========================================="

# Load deeptools module
module load deeptools/3.2.1

# Create timestamp for output files
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")
echo "Timestamp: ${TIMESTAMP}"

# Read the bigwig directory path
BIGWIG_DIR=$(cat ${ARCHR_PROJECT_DIR}/pseudobulk_bigwig_${GROUP_BY}/bigwig_dir.txt)

echo "BigWig directory: ${BIGWIG_DIR}"
echo "BED directory: ${BED_DIR}"
echo "BED pattern: ${BED_PATTERN}"

# Check if bigwig directory exists
if [ ! -d "${BIGWIG_DIR}" ]; then
    echo "ERROR: BigWig directory not found: ${BIGWIG_DIR}"
    exit 1
fi

# Check if BED directory exists
if [ ! -d "${BED_DIR}" ]; then
    echo "ERROR: BED directory not found: ${BED_DIR}"
    exit 1
fi

# Find all bigWig files, applying exclusion pattern if provided
if [ -n "${EXCLUDE_PATTERN}" ]; then
    echo "Applying exclude pattern: ${EXCLUDE_PATTERN}"
    BIGWIG_FILES=$(find ${BIGWIG_DIR} -name "*.bw" | grep -vE "${EXCLUDE_PATTERN}" | sort)
else
    BIGWIG_FILES=$(find ${BIGWIG_DIR} -name "*.bw" | sort)
fi

if [ -z "${BIGWIG_FILES}" ]; then
    echo "ERROR: No bigWig files found in ${BIGWIG_DIR} (or all were excluded)"
    exit 1
fi

echo ""
echo "Found bigWig files for heatmap generation:"
echo "${BIGWIG_FILES}"
echo ""

# Find all BED files matching the pattern
BED_FILES=$(find ${BED_DIR} -name "${BED_PATTERN}*.bed" | sort)

if [ -z "${BED_FILES}" ]; then
    echo "ERROR: No BED files found matching pattern ${BED_PATTERN}*.bed in ${BED_DIR}"
    exit 1
fi

echo "Found BED files:"
echo "${BED_FILES}"
echo ""

# Create output directory for deeptools results
DEEPTOOLS_OUT="${ARCHR_PROJECT_DIR}/deeptools_${GROUP_BY}_${BED_PATTERN}"
mkdir -p ${DEEPTOOLS_OUT}

echo "Output directory: ${DEEPTOOLS_OUT}"
echo ""

########################################################
# Run computeMatrix
########################################################

echo "Running computeMatrix..."
echo "=========================================="

# Convert file lists to space-separated strings
BIGWIG_LIST=$(echo ${BIGWIG_FILES} | tr '\n' ' ')
BED_LIST=$(echo ${BED_FILES} | tr '\n' ' ')

# Extract sample labels from bigwig filenames
SAMPLE_LABELS=$(for f in ${BIGWIG_FILES}; do basename $f .bw | cut -d'-' -f1; done | tr '\n' ' ')

# Extract region labels from BED filenames
REGION_LABELS=$(for f in ${BED_FILES}; do basename $f _genes.bed; done | tr '\n' ' ')

echo "Sample labels: ${SAMPLE_LABELS}"
echo "Region labels: ${REGION_LABELS}"
echo ""

# Run computeMatrix with reference-point mode (TSS)
computeMatrix reference-point \
    --referencePoint TSS \
    -b 2000 -a 2000 \
    -R ${BED_LIST} \
    -S ${BIGWIG_LIST} \
    -o ${DEEPTOOLS_OUT}/matrix_refpt_${GROUP_BY}_${BED_PATTERN}.gz \
    --outFileNameMatrix ${DEEPTOOLS_OUT}/matrix_refpt_${GROUP_BY}_${BED_PATTERN}.tab \
    --outFileSortedRegions ${DEEPTOOLS_OUT}/regions_refpt_${GROUP_BY}_${BED_PATTERN}.bed \
    -p ${SLURM_CPUS_PER_TASK:-4} \
    --verbose

if [ $? -eq 0 ]; then
    echo "✓ computeMatrix reference-point completed successfully"
else
    echo "✗ computeMatrix reference-point failed"
    exit 1
fi

echo ""

# Run computeMatrix with scale-regions mode (gene bodies)
computeMatrix scale-regions \
    -R ${BED_LIST} \
    -S ${BIGWIG_LIST} \
    -o ${DEEPTOOLS_OUT}/matrix_scale_region_${GROUP_BY}_${BED_PATTERN}.gz \
    --outFileNameMatrix ${DEEPTOOLS_OUT}/matrix_scale_region_${GROUP_BY}_${BED_PATTERN}.tab \
    --outFileSortedRegions ${DEEPTOOLS_OUT}/matrix_scale_region_${GROUP_BY}_${BED_PATTERN}.bed \
    -p ${SLURM_CPUS_PER_TASK:-4} \
    --beforeRegionStartLength 2000 \
    --regionBodyLength 5000 \
    --afterRegionStartLength 2000 \
    --verbose 

if [ $? -eq 0 ]; then
    echo "✓ computeMatrix scale-regions completed successfully"
else    
    echo "✗ computeMatrix scale-regions failed"
    exit 1
fi
echo ""

########################################################
# Run plotHeatmap
########################################################

echo "Running plotHeatmap for reference-point mode..."
echo "=========================================="

plotHeatmap \
    -m ${DEEPTOOLS_OUT}/matrix_refpt_${GROUP_BY}_${BED_PATTERN}.gz \
    -out ${DEEPTOOLS_OUT}/heatmap_refpt_${GROUP_BY}_${BED_PATTERN}_${TIMESTAMP}.pdf \
    --colorMap RdYlBu_r \
    --heatmapHeight 15 \
    --refPointLabel "TSS" \
    --samplesLabel ${SAMPLE_LABELS} \
    --regionsLabel ${REGION_LABELS} \
    --plotTitle "ATAC-seq signal at gene signatures (${BED_PATTERN})" \
    --sortRegions descend \
    --sortUsing mean \
    --averageTypeSummaryPlot mean

if [ $? -eq 0 ]; then
    echo "✓ plotHeatmap completed successfully"
else
    echo "✗ plotHeatmap failed"
    exit 1
fi

echo ""

# create heatmap for scale-regions mode
echo "Running plotHeatmap for scale-regions mode..."
plotHeatmap \
    -m ${DEEPTOOLS_OUT}/matrix_scale_region_${GROUP_BY}_${BED_PATTERN}.gz \
    -out ${DEEPTOOLS_OUT}/heatmap_scale_region_${GROUP_BY}_${BED_PATTERN}_${TIMESTAMP}.pdf \
    --colorMap RdYlBu_r \
    --heatmapHeight 15 \
    --samplesLabel ${SAMPLE_LABELS} \
    --regionsLabel ${REGION_LABELS} \
    --plotTitle "ATAC-seq signal across gene bodies (${BED_PATTERN})" \
    --sortRegions descend \
    --sortUsing mean \
    --averageTypeSummaryPlot mean

if [ $? -eq 0 ]; then
    echo "✓ plotHeatmap for scale-regions completed successfully"
else
    echo "✗ plotHeatmap for scale-regions failed"
    exit 1
fi  
echo ""

# Profile plot for reference-point mode
echo "Creating profile plot for reference-point mode..."
plotProfile \
    -m ${DEEPTOOLS_OUT}/matrix_refpt_${GROUP_BY}_${BED_PATTERN}.gz \
    -out ${DEEPTOOLS_OUT}/profile_refpt_${GROUP_BY}_${BED_PATTERN}_${TIMESTAMP}.pdf \
    --perGroup \
    --colors blue yellow red \
    --refPointLabel "TSS" \
    --samplesLabel ${SAMPLE_LABELS} \
    --regionsLabel ${REGION_LABELS} \
    --plotTitle "Average ATAC-seq signal (${BED_PATTERN})" \
    --yAxisLabel "Mean ATAC signal"

if [ $? -eq 0 ]; then
    echo "✓ plotProfile completed successfully"
else
    echo "✗ plotProfile failed (non-critical)"
fi

# Profile plot for scale-regions mode
echo "Creating profile plot for scale-regions mode..."
plotProfile \
    -m ${DEEPTOOLS_OUT}/matrix_scale_region_${GROUP_BY}_${BED_PATTERN}.gz \
    -out ${DEEPTOOLS_OUT}/profile_scale_region_${GROUP_BY}_${BED_PATTERN}_${TIMESTAMP}.pdf \
    --perGroup \
    --colors blue yellow red \
    --samplesLabel ${SAMPLE_LABELS} \
    --regionsLabel ${REGION_LABELS} \
    --plotTitle "Average ATAC-seq signal across gene bodies (${BED_PATTERN})" \
    --yAxisLabel "Mean ATAC signal"

if [ $? -eq 0 ]; then
    echo "✓ plotProfile for scale-regions completed successfully"
else
    echo "✗ plotProfile for scale-regions failed (non-critical)"
fi

echo ""
echo "=========================================="
echo "ANALYSIS COMPLETE!"
echo "=========================================="
echo ""
echo "Output files:"
echo "  BigWig files: ${BIGWIG_DIR}"
echo "  Heatmap (refpt): ${DEEPTOOLS_OUT}/heatmap_refpt_${GROUP_BY}_${BED_PATTERN}_${TIMESTAMP}.pdf"
echo "  Heatmap (scale): ${DEEPTOOLS_OUT}/heatmap_scale_region_${GROUP_BY}_${BED_PATTERN}_${TIMESTAMP}.pdf"
echo "  Profile (refpt): ${DEEPTOOLS_OUT}/profile_refpt_${GROUP_BY}_${BED_PATTERN}_${TIMESTAMP}.pdf"
echo "  Profile (scale): ${DEEPTOOLS_OUT}/profile_scale_region_${GROUP_BY}_${BED_PATTERN}_${TIMESTAMP}.pdf"
echo ""
