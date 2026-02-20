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
# First argument: groupBy variable (e.g., "PIMO_up_status", "hybrid_pair")
# Second argument: path to directory containing BED files
# Third argument (optional): pattern to match BED files (e.g., "PIMO_MP", "Nomura_")
# Fourth argument (optional): "force" to re-export bigWig files even if they exist

# Example usage:
# sbatch scripts/export_04_pseudobulk_bigwig_heatmaps.sh PIMO_up_status data/Signatures/PIMO_metaprogram/ PIMO_MP
# sbatch scripts/export_04_pseudobulk_bigwig_heatmaps.sh PIMO_up_status data/Signatures/nomura Nomura_

GROUP_BY="${1:-PIMO_up_status}"  # Default to "PIMO_up_status" if no argument provided
BED_DIR="${2:-data/Signatures}"  # Default BED directory
BED_PATTERN="${3:-PIMO_MP}"  # Default to PIMO Metaprograms
FORCE_EXPORT="${4:-no}"  # Default to not forcing re-export

# Export the parameters so R can access them
export GROUP_BY
export BED_DIR
export BED_PATTERN
export FORCE_EXPORT

echo "=========================================="
echo "Running pseudo-bulk bigWig export and heatmap generation"
echo "=========================================="
echo "  groupBy: ${GROUP_BY}"
echo "  BED directory: ${BED_DIR}"
echo "  BED file pattern: ${BED_PATTERN}"
echo "  Force re-export: ${FORCE_EXPORT}"
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
groupBy <- Sys.getenv("GROUP_BY", unset = "PIMO_up_status")
forceExport <- Sys.getenv("FORCE_EXPORT", unset = "no")

cat("Using groupBy:", groupBy, "\n")
cat("Force re-export:", forceExport, "\n")

# Set the number of threads for ArchR
addArchRThreads(threads = 18)

# Set genome to hg38
addArchRGenome("hg38")

# Load the ArchR project
cat("Loading ArchR project...\n")
proj <- loadArchRProject(path = "human_multiome_harmony_merged_malig_peak")

# Create output directory for bigWig files
filePrefix <- gsub("_", "", groupBy)
outDir <- here(paste0("human_multiome_harmony_merged_malig_peak/pseudobulk_bigwig_", groupBy, "/"))
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
      # ArchR creates filenames like: "GroupName-TileSize-100-normMethod-ReadsInTSS-ArchR.bw"
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
      cat("To force re-export, use the 'force' argument:\n")
      cat("  sbatch script.sh PIMO_up_status /path/to/beds PIMO_MP force\n")
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
  
  # getGroupBW creates bigWig files for ALL groups at once
  # It does not have a useGroups parameter - it processes all groups in groupBy
  tryCatch({
    getGroupBW(
      ArchRProj = proj,
      groupBy = groupBy,
      normMethod = "ReadsInTSS",  # Normalize by reads in TSS
      tileSize = 100,  # 100bp bins
      maxCells = 500,  # Subsample if more than 500 cells per group
      ceiling = 4,  # Set ceiling for coverage
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

# The bigWig files are saved in: proj@projectMetadata$outputDirectory/GroupBigWigs/
bigwigDir <- file.path(proj@projectMetadata$outputDirectory, "GroupBigWigs", groupBy)
cat("BigWig files location:", bigwigDir, "\n")

# List final bigWig files
finalBigWigs <- list.files(bigwigDir, pattern = "\\.bw$", full.names = FALSE)
cat(sprintf("\nTotal bigWig files available: %d\n", length(finalBigWigs)))
for (bw in finalBigWigs) {
  bwPath <- file.path(bigwigDir, bw)
  bwSize <- file.info(bwPath)$size / (1024^3)  # Convert to GB
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
BIGWIG_DIR=$(cat human_multiome_harmony_merged_malig_peak/pseudobulk_bigwig_${GROUP_BY}/bigwig_dir.txt)

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

# Find all bigWig files
BIGWIG_FILES=$(find ${BIGWIG_DIR} -name "*.bw" | sort)

if [ -z "${BIGWIG_FILES}" ]; then
    echo "ERROR: No bigWig files found in ${BIGWIG_DIR}"
    exit 1
fi

echo ""
echo "Found bigWig files:"
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
DEEPTOOLS_OUT="human_multiome_harmony_merged_malig_peak/deeptools_${GROUP_BY}_${BED_PATTERN}"
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
    -p ${SLURM_CPUS_PER_TASK} \
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
    -p ${SLURM_CPUS_PER_TASK} \
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

# Create heatmap for reference-point mode
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
echo "  Heatmap: ${DEEPTOOLS_OUT}/heatmap_${GROUP_BY}_${BED_PATTERN}_${TIMESTAMP}.pdf"
echo "  Profile: ${DEEPTOOLS_OUT}/profile_${GROUP_BY}_${BED_PATTERN}_${TIMESTAMP}.pdf"
echo "  Matrix: ${DEEPTOOLS_OUT}/matrix_${GROUP_BY}_${BED_PATTERN}_${TIMESTAMP}.gz"
echo ""