#!/bin/bash
#SBATCH --job-name=run_05_human_multiome_peakCalling
#SBATCH --output=run_05_human_multiome_peakCalling_%j.out
#SBATCH --error=run_05_human_multiome_peakCalling_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=60G
#SBATCH --time=03:00:00

# ─── Usage ────────────────────────────────────────────────────────────────────
# Parameters are passed via sbatch --export. All are optional; defaults are used
# for any variable not supplied.
#
# Variables:
#   ARCHR_PROJECT   Path to load the ArchR project
#                   (default: human_multiome_harmony_merged)
#   SKIP_SUBSET     Skip the malignant cell subsetting step; true or false
#                   (default: false)
#   SUBSET_OUTDIR   Output directory name for the subsetted project
#                   (default: human_multiome_harmony_merged_malig_peak)
#   GROUP_BY_COLS   Comma-separated metadata columns for peak calling
#                   (default: hybrid_pair,PIMO_up_status)
#
# Examples:
#   sbatch run_05_human_multiome_peakCalling.sh
#   sbatch --export=ALL,SKIP_SUBSET=true,ARCHR_PROJECT=human_multiome_hmmp_k3_program scripts/run_05_human_multiome_peakCalling.sh
#   sbatch --export=ALL,ARCHR_PROJECT=human_multiome_hmmp_k3_program,SKIP_SUBSET=true,SUBSET_OUTDIR=human_multiome_hmmp_k3_program_peak,GROUP_BY_COLS=k3_program scripts/run_05_human_multiome_peakCalling.sh
#   sbatch --export=ALL,ARCHR_PROJECT=human_multiome_hmmp_nomura,SKIP_SUBSET=true,SUBSET_OUTDIR=human_multiome_hmmp_nomura_peak,GROUP_BY_COLS=meta_module scripts/run_05_human_multiome_peakCalling.sh
# ──────────────────────────────────────────────────────────────────────────────

# ─── Parameters (set via --export at sbatch submission) ───────────────────────
# Defaults are used for any variable not supplied via --export.
ARCHR_PROJECT="${ARCHR_PROJECT:-human_multiome_harmony_merged}"
SKIP_SUBSET="${SKIP_SUBSET:-false}"
SUBSET_OUTDIR="${SUBSET_OUTDIR:-human_multiome_harmony_merged_malig_peak}"
GROUP_BY_COLS="${GROUP_BY_COLS:-hybrid_pair,PIMO_up_status}"

echo "=========================================="
echo "  ArchR Peak Calling Pipeline"
echo "=========================================="
echo "  ArchR project path : $ARCHR_PROJECT"
echo "  Skip subset step   : $SKIP_SUBSET"
echo "  Subset output dir  : $SUBSET_OUTDIR"
echo "  GroupBy columns    : $GROUP_BY_COLS"
echo "=========================================="

# Load necessary modules (adjust as needed for your system)
module load R/4.4.1
module load python3/3.7.2
module load MACS/2.2.7.1

# Export shell variables so R can read them via Sys.getenv()
export ARCHR_PROJECT SKIP_SUBSET SUBSET_OUTDIR GROUP_BY_COLS

# Run R script
Rscript - <<'EOF'

# load libraries
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(here)
set.seed(1)

# ─── Parameters passed in from shell ─────────────────────────────────────────
archr_project_path <- Sys.getenv("ARCHR_PROJECT")
skip_subset        <- as.logical(Sys.getenv("SKIP_SUBSET"))
subset_outdir      <- Sys.getenv("SUBSET_OUTDIR")
groupBy_list       <- strsplit(Sys.getenv("GROUP_BY_COLS"), ",")[[1]]
groupBy_list       <- trimws(groupBy_list)  # remove any accidental whitespace

# Set the number of threads for ArchR
addArchRThreads(threads = 18)

# set genome to hg38
addArchRGenome("hg38")

# ─── Load the project ────────────────────────────────────────────────────────
print(paste("Loading ArchR project from:", archr_project_path))
proj_hyp <- loadArchRProject(path = archr_project_path)

# ─── Optional: subset to malignant cells ─────────────────────────────────────
if (!skip_subset) {
  print("Subsetting ArchR object to only malignant cells and non-ambiguous: ")
  malignant_cells <- proj_hyp$cellNames[which(proj_hyp$Azimuth_class == "Malignant" & proj_hyp$Ambiguous == FALSE)]
  print(paste("Number of malignant non-ambiguous cells:", length(malignant_cells)))

  proj_hyp <- subsetArchRProject(
    ArchRProj = proj_hyp,
    cells = malignant_cells,
    outputDirectory = subset_outdir,
    dropCells = TRUE,
    force = TRUE
  )
} else {
  print("Skipping subset step. Using project as-is.")
  # When skipping subset, treat the loaded project directory as the output dir
  # (used for saving outputs below)
  subset_outdir <- archr_project_path
}

# ─── Validate groupBy columns ─────────────────────────────────────────────────
print(paste("GroupBy list:", paste(groupBy_list, collapse = ", ")))
available_cols <- names(proj_hyp@cellColData)
missing_cols <- groupBy_list[!groupBy_list %in% available_cols]
if (length(missing_cols) > 0) {
  stop(paste(
    "The following groupBy column(s) were not found in the ArchR project metadata:",
    paste(missing_cols, collapse = ", "),
    "\nAvailable columns are:",
    paste(available_cols, collapse = ", ")
  ))
}

# ─── Find MACS2 path ─────────────────────────────────────────────────────────
if (exists("findMacs2")) {
  pathToMacs2 <- findMacs2()
  if (!is.null(pathToMacs2)) {
    print(paste("MACS2 found at:", pathToMacs2))
    pathToMacs2 <- "/cluster/tools/software/rocky9/MACS/2.2.7.1/bin/macs2"
  } else {
    stop("MACS2 not found or path is invalid. Please check your MACS2 installation.")
  }
} else {
  stop("findMacs2() function not found. Please ensure ArchR is installed correctly.")
}

# ─── Peak calling loop ────────────────────────────────────────────────────────
for (groupBy in groupBy_list) {
  print(paste("Processing groupBy:", groupBy))

  # Add group coverages
  proj_hyp <- addGroupCoverages(
    ArchRProj = proj_hyp,
    groupBy = groupBy
  )

  # Peak calling
  print(paste("Call Peaks with MACS2 for", groupBy))
  proj_hyp2 <- addReproduciblePeakSet(
    ArchRProj = proj_hyp,
    groupBy = groupBy,
    pathToMacs2 = pathToMacs2
  )

  # Get peak set
  print("get Peaks set")
  myPeakSet <- getPeakSet(proj_hyp2)

  # Tabulate type of peaks
  print("tabulate type of peaks")
  print(table(myPeakSet$peakType))

  # Save the peakSet gr object
  peakset_path <- file.path(subset_outdir, "PeakCalls", paste0("PeakSet_gr_", groupBy, ".rds"))
  print(paste("Saved peakSet as .RData for", groupBy))
  saveRDS(myPeakSet, file = peakset_path)

  # Check available matrices
  print("check available matrices")
  print(getAvailableMatrices(proj_hyp2))

  # Add peak matrix
  print("adding peak matrix")
  proj_hyp2 <- addPeakMatrix(proj_hyp2)

  # Check available matrices after adding peak matrix
  print("check available matrices after adding peak matrix")
  print(getAvailableMatrices(proj_hyp2))

  # Fix incompatible dimensions error (change factor back to logical if binary)
  if (length(levels(proj_hyp2@cellColData[[groupBy]])) == 2) {
    print("fix incompatible dimensions error")
    proj_hyp2@cellColData[[groupBy]] <- as.logical(proj_hyp2@cellColData[[groupBy]])
  }

  # Check dimensions of PeakMatrix and cellColData
  print("Check dimensions of PeakMatrix and cellColData: ")
  peakMat <- getMatrixFromProject(ArchRProj = proj_hyp2, useMatrix = "PeakMatrix")
  print(dim(assay(peakMat)))
  print(length(proj_hyp2@cellColData[[groupBy]]))

  # Get marker peaks for groupBy
  print(paste("Get marker peaks for", groupBy, "groups"))
  markersPeaks <- getMarkerFeatures(
    ArchRProj = proj_hyp2,
    useMatrix = "PeakMatrix",
    groupBy = groupBy,
    maxCells = 1000,
    bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
    testMethod = "wilcoxon"
  )

  # Save the markersPeaks SE object
  markers_path <- file.path(subset_outdir, "PeakCalls", paste0("markersPeaks_", groupBy, ".rds"))
  print(paste("Saved markersPeaks as .rds for", groupBy))
  saveRDS(markersPeaks, file = markers_path)

  # Save the markersPeaks GR object
  markers_gr_path <- file.path(subset_outdir, "PeakCalls", paste0("markersPeaks_GR_", groupBy, ".rds"))
  print(paste("Saved markersPeaks GR as .rds for", groupBy))
  markersPeaks_GR <- getMarkers(markersPeaks, cutOff = "FDR <=  1 & abs(Log2FC) >= 0", returnGR = TRUE)
  saveRDS(markersPeaks_GR, file = markers_gr_path)
}

# ─── Pairwise test: PIMOup vs PIMOdown (only if PIMO_up_status was in groupBy) ─
if ("PIMO_up_status" %in% groupBy_list) {
  print("Pairwise test between PIMO_up_status groups: PIMOup vs PIMOdown")
  markerTest <- getMarkerFeatures(
    ArchRProj = proj_hyp2,
    useMatrix = "PeakMatrix",
    groupBy = "PIMO_up_status",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
    useGroups = "PIMOup",
    bgdGroups = "PIMOdown"
  )

  # Save markerTest SE object
  markertest_path <- file.path(subset_outdir, "PeakCalls", "markersTest_PIMO_up_status_PIMOup_vs_PIMOdown.rds")
  print("Saved markerTest as .rds for PIMO_up_status PIMOup vs PIMOdown")
  saveRDS(markerTest, file = markertest_path)

  # Extract and save markerTest GR object
  print("Extracting markerTestGR object for PIMO_up_status")
  markerTest_GR <- getMarkers(markerTest, cutOff = "FDR <= 1 & abs(Log2FC) >= 0", returnGR = TRUE)

  print("Number of markerTest peaks identified per group:")
  print(sapply(markerTest_GR, length))

  markertest_gr_path <- file.path(subset_outdir, "PeakCalls", "markersTest_GR_PIMO_up_status_PIMOup_vs_PIMOdown.rds")
  saveRDS(markerTest_GR, file = markertest_gr_path)
} else {
  print("PIMO_up_status not in groupBy list -- skipping pairwise PIMOup vs PIMOdown test.")
}

# ─── Save the project ─────────────────────────────────────────────────────────
proj_hyp2 <- saveArchRProject(ArchRProj = proj_hyp2, outputDirectory = subset_outdir, load = TRUE)

EOF

echo "Peak Calling analysis completed"