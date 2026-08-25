#!/bin/bash
#SBATCH --job-name=run_05_Gaiti_multiome_peakCalling
#SBATCH --output=run_05_Gaiti_multiome_peakCalling_%j.out
#SBATCH --error=run_05_Gaiti_multiome_peakCalling_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=60G
#SBATCH --time=06:00:00

# ─── Usage ────────────────────────────────────────────────────────────────────
# Parameters are passed via sbatch --export. All are optional; defaults are used
# for any variable not supplied.
#
# Variables:
#   ARCHR_PROJECT   Path to load the ArchR project
#                   (default: Gaiti_multiome_harmony_merged)
#   SKIP_SUBSET     Skip the malignant cell subsetting step; true or false
#                   (default: false)
#   SUBSET_OUTDIR   Output directory name for the subsetted project
#                   (default: Gaiti_multiome_harmony_merged_malig_peak)
#   GROUP_BY_COLS   Comma-separated metadata columns for peak calling
#                   (default: PIMO_up_status,PIMO_Region,Region.annotation)
#
# Examples:
#   sbatch scripts/run_05_Gaiti_multiome_peakCalling.sh
#   sbatch --export=ALL,ARCHR_PROJECT=Gaiti_multiome_hmmp_k3_program_all_cells,SKIP_SUBSET=true,SUBSET_OUTDIR=Gaiti_multiome_hmmp_k3_program_all_cells,GROUP_BY_COLS=combined_program scripts/run_05_Gaiti_multiome_peakCalling.sh
#   sbatch --export=ALL,ARCHR_PROJECT=Gaiti_multiome_harmony_merged_malig_peak,SKIP_SUBSET=true,SUBSET_OUTDIR=Gaiti_multiome_harmony_merged_malig_peak,GROUP_BY_COLS=PIMO_Region scripts/run_05_Gaiti_multiome_peakCalling.sh
# ──────────────────────────────────────────────────────────────────────────────

ARCHR_PROJECT="${ARCHR_PROJECT:-Gaiti_multiome_harmony_merged}"
SKIP_SUBSET="${SKIP_SUBSET:-false}"
SUBSET_OUTDIR="${SUBSET_OUTDIR:-Gaiti_multiome_harmony_merged_malig_peak}"
GROUP_BY_COLS="${GROUP_BY_COLS:-PIMO_up_status,PIMO_Region,Region.annotation}"

echo "=========================================="
echo "  ArchR Peak Calling Pipeline (Gaiti)"
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
archr_project_path <- Sys.getenv("ARCHR_PROJECT", unset = "Gaiti_multiome_harmony_merged")
skip_subset        <- as.logical(Sys.getenv("SKIP_SUBSET", unset = "false"))
subset_outdir      <- Sys.getenv("SUBSET_OUTDIR", unset = "Gaiti_multiome_harmony_merged_malig_peak")
groupBy_list       <- strsplit(Sys.getenv("GROUP_BY_COLS", unset = "PIMO_up_status,PIMO_Region,Region.annotation"), ",")[[1]]
groupBy_list       <- trimws(groupBy_list)

# Set the number of threads for ArchR
addArchRThreads(threads = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "8")))

# set genome to hg38
addArchRGenome("hg38")

# ─── Load the project ────────────────────────────────────────────────────────
print(paste("Loading ArchR project from:", archr_project_path))
proj_hyp <- loadArchRProject(path = archr_project_path)

# ─── Optional: subset to malignant cells ─────────────────────────────────────
if (!skip_subset) {
  print("Subsetting ArchR object to only malignant cells and not mislabelled cells")
  malignant_cells <- proj_hyp$cellNames[which(proj_hyp$Azimuth_class == "Malignant" & proj_hyp$mislabelled == "Correct")]
  print(paste("Number of malignant non-ambiguous cells:", length(malignant_cells)))

  # add additional metadata columns if not already present
  if (!"PIMO_Region" %in% colnames(proj_hyp@cellColData)) {
    print("Adding PIMO_Region metadata column based on PIMO_up_status and Region.annotation")
    proj_hyp@cellColData$PIMO_Region <- paste(
      proj_hyp$PIMO_up_status,
      proj_hyp$Region.annotation,
      sep = "_"
    )
  } else {
    print("PIMO_Region metadata column already exists")
  }

  proj_hyp <- subsetArchRProject(
    ArchRProj = proj_hyp,
    cells = malignant_cells,
    outputDirectory = subset_outdir,
    dropCells = TRUE,
    force = TRUE
  )
} else {
  print("Skipping subset step. Using project as-is.")
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

# Initialise list to store one proj_hyp2 per groupBy for use in pairwise tests below
proj_hyp2_list <- list()

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
  peak_calls_dir <- file.path(subset_outdir, "PeakCalls")
  dir.create(peak_calls_dir, showWarnings = FALSE, recursive = TRUE)
  
  peakset_path <- file.path(peak_calls_dir, paste0("PeakSet_gr_", groupBy, ".rds"))
  print(paste("Saved peakSet as .rds for", groupBy, "to", peakset_path))
  saveRDS(myPeakSet, file = peakset_path)
  
  # Check available matrices
  print("check available matrices")
  print(getAvailableMatrices(proj_hyp2))
  
  # Add peak matrix
  print("adding peak matrix")
  proj_hyp2 <- addPeakMatrix(proj_hyp2, force = TRUE)
  
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
  markers_path <- file.path(peak_calls_dir, paste0("markersPeaks_", groupBy, ".rds"))
  print(paste("Saved markersPeaks as .rds for", groupBy))
  saveRDS(markersPeaks, file = markers_path)

  # Save the markersPeaks GR object
  markers_gr_path <- file.path(peak_calls_dir, paste0("markersPeaks_GR_", groupBy, ".rds"))
  print(paste("Saved markersPeaks GR as .rds for", groupBy))
  markersPeaks_GR <- getMarkers(markersPeaks, cutOff = "FDR <= 1 & abs(Log2FC) >= 0", returnGR = TRUE)
  saveRDS(markersPeaks_GR, file = markers_gr_path)

  # Store proj_hyp2 for pairwise tests below
  proj_hyp2_list[[groupBy]] <- proj_hyp2
}

# ─── Pairwise tests for PIMO_Region (if PIMO_Region in groupBy_list) ─────────
if ("PIMO_Region" %in% groupBy_list) {
  print("Loading PIMO_Region proj_hyp2 for pairwise tests")
  proj_hyp2_pimo_reg <- proj_hyp2_list[["PIMO_Region"]]

  # 1. PIMOup_EB vs PIMOup_TC
  print("Pairwise test between PIMO_Region groups: PIMOup_EB vs PIMOup_TC")
  markerTest <- getMarkerFeatures(
    ArchRProj = proj_hyp2_pimo_reg, 
    useMatrix = "PeakMatrix",
    groupBy = "PIMO_Region",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
    useGroups = "PIMOup_EB",
    bgdGroups = "PIMOup_TC"
  )
  saveRDS(markerTest, file = file.path(subset_outdir, "PeakCalls", "markersTest_PIMO_Region_PIMOup_EB_vs_PIMOup_TC.rds"))
  markerTest_GR <- getMarkers(markerTest, cutOff = "FDR <= 1 & abs(Log2FC) >= 0", returnGR = TRUE)
  saveRDS(markerTest_GR, file = file.path(subset_outdir, "PeakCalls", "markersTest_GR_PIMO_Region_PIMOup_EB_vs_PIMOup_TC.rds"))

  # 2. PIMOup_EB vs PIMOdown_EB
  print("Pairwise test between PIMO_Region groups: PIMOup_EB vs PIMOdown_EB")
  markerTest <- getMarkerFeatures(
    ArchRProj = proj_hyp2_pimo_reg, 
    useMatrix = "PeakMatrix",
    groupBy = "PIMO_Region",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
    useGroups = "PIMOup_EB",
    bgdGroups = "PIMOdown_EB"
  )
  saveRDS(markerTest, file = file.path(subset_outdir, "PeakCalls", "markersTest_PIMO_Region_PIMOup_EB_vs_PIMOdown_EB.rds"))
  markerTest_GR <- getMarkers(markerTest, cutOff = "FDR <= 1 & abs(Log2FC) >= 0", returnGR = TRUE)
  saveRDS(markerTest_GR, file = file.path(subset_outdir, "PeakCalls", "markersTest_GR_PIMO_Region_PIMOup_EB_vs_PIMOdown_EB.rds"))

  # 3. PIMOup_TC vs PIMOdown_TC
  print("Pairwise test between PIMO_Region groups: PIMOup_TC vs PIMOdown_TC")
  markerTest <- getMarkerFeatures(
    ArchRProj = proj_hyp2_pimo_reg, 
    useMatrix = "PeakMatrix",
    groupBy = "PIMO_Region",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
    useGroups = "PIMOup_TC",
    bgdGroups = "PIMOdown_TC"
  )
  saveRDS(markerTest, file = file.path(subset_outdir, "PeakCalls", "markersTest_PIMO_Region_PIMOup_TC_vs_PIMOdown_TC.rds"))
  markerTest_GR <- getMarkers(markerTest, cutOff = "FDR <= 1 & abs(Log2FC) >= 0", returnGR = TRUE)
  saveRDS(markerTest_GR, file = file.path(subset_outdir, "PeakCalls", "markersTest_GR_PIMO_Region_PIMOup_TC_vs_PIMOdown_TC.rds"))

  # 4. PIMOup_PT vs PIMOdown_PT
  print("Pairwise test between PIMO_Region groups: PIMOup_PT vs PIMOdown_PT")
  markerTest <- getMarkerFeatures(
    ArchRProj = proj_hyp2_pimo_reg, 
    useMatrix = "PeakMatrix",
    groupBy = "PIMO_Region",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
    useGroups = "PIMOup_PT",
    bgdGroups = "PIMOdown_PT"
  )
  saveRDS(markerTest, file = file.path(subset_outdir, "PeakCalls", "markersTest_PIMO_Region_PIMOup_PT_vs_PIMOdown_PT.rds"))
  markerTest_GR <- getMarkers(markerTest, cutOff = "FDR <= 1 & abs(Log2FC) >= 0", returnGR = TRUE)
  saveRDS(markerTest_GR, file = file.path(subset_outdir, "PeakCalls", "markersTest_GR_PIMO_Region_PIMOup_PT_vs_PIMOdown_PT.rds"))
}

# ─── Pairwise tests for Region.annotation (if Region.annotation in groupBy_list) ─
if ("Region.annotation" %in% groupBy_list) {
  print("Loading Region.annotation proj_hyp2 for pairwise tests")
  proj_hyp2_reg_anno <- proj_hyp2_list[["Region.annotation"]]

  # 1. EB vs TC
  print("Pairwise test between Region.annotation groups: EB vs TC") 
  markerTest <- getMarkerFeatures(
    ArchRProj = proj_hyp2_reg_anno, 
    useMatrix = "PeakMatrix",
    groupBy = "Region.annotation",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
    useGroups = "EB",
    bgdGroups = "TC"
  )
  saveRDS(markerTest, file = file.path(subset_outdir, "PeakCalls", "markersTest_Region.annotation_EB_vs_TC.rds"))
  markerTest_GR <- getMarkers(markerTest, cutOff = "FDR <= 1 & abs(Log2FC) >= 0", returnGR = TRUE)
  saveRDS(markerTest_GR, file = file.path(subset_outdir, "PeakCalls", "markersTest_GR_Region.annotation_EB_vs_TC.rds"))

  # 2. EB vs PT
  print("Pairwise test between Region.annotation groups: EB vs PT")
  markerTest <- getMarkerFeatures(
    ArchRProj = proj_hyp2_reg_anno, 
    useMatrix = "PeakMatrix",
    groupBy = "Region.annotation",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
    useGroups = "EB",
    bgdGroups = "PT"
  )
  saveRDS(markerTest, file = file.path(subset_outdir, "PeakCalls", "markersTest_Region.annotation_EB_vs_PT.rds"))
  markerTest_GR <- getMarkers(markerTest, cutOff = "FDR <= 1 & abs(Log2FC) >= 0", returnGR = TRUE)
  saveRDS(markerTest_GR, file = file.path(subset_outdir, "PeakCalls", "markersTest_GR_Region.annotation_EB_vs_PT.rds"))

  # 3. TC vs PT
  print("Pairwise test between Region.annotation groups: TC vs PT")
  markerTest <- getMarkerFeatures(
    ArchRProj = proj_hyp2_reg_anno, 
    useMatrix = "PeakMatrix",
    groupBy = "Region.annotation",
    testMethod = "wilcoxon",
    bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
    useGroups = "TC",
    bgdGroups = "PT"
  )
  saveRDS(markerTest, file = file.path(subset_outdir, "PeakCalls", "markersTest_Region.annotation_TC_vs_PT.rds"))
  markerTest_GR <- getMarkers(markerTest, cutOff = "FDR <= 1 & abs(Log2FC) >= 0", returnGR = TRUE)
  saveRDS(markerTest_GR, file = file.path(subset_outdir, "PeakCalls", "markersTest_GR_Region.annotation_TC_vs_PT.rds"))
}

# ─── Save the project ─────────────────────────────────────────────────────────
if (exists("proj_hyp2")) {
  proj_hyp2 <- saveArchRProject(ArchRProj = proj_hyp2, outputDirectory = subset_outdir, load = TRUE)
}

EOF

echo "Peak Calling analysis completed"