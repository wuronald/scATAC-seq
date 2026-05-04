#!/bin/bash
#SBATCH --job-name=run_05b_Gaiti_multiome_pairwiseTesting
#SBATCH --output=run_05b_Gaiti_multiome_pairwiseTesting_%j.out
#SBATCH --error=run_05b_Gaiti_multiome_pairwiseTesting_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=58G
#SBATCH --time=06:00:00

# Load necessary modules
module load R/4.4.1

# Run R script
Rscript - <<'EOF'

# load libraries
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(here)
set.seed(1)

# Set the number of threads for ArchR
addArchRThreads(threads = 6)

# Set genome to hg38
addArchRGenome("hg38")

# ─────────────────────────────────────────────────────────────────────────────
# Load the already-subsetted malignant-cell project (output of peak calling run)
# ─────────────────────────────────────────────────────────────────────────────
print("Loading ArchR project from Gaiti_multiome_harmony_merged_malig_peak")
proj_hyp2 <- loadArchRProject(path = "Gaiti_multiome_harmony_merged_malig_peak")
print(proj_hyp2)

# Convenience path to the pre-computed peak call outputs
peakCallDir <- "Gaiti_multiome_harmony_merged_malig_peak/PeakCalls"

# ─────────────────────────────────────────────────────────────────────────────
# Helper: load a saved peakSet GR, inject it into the project via addPeakSet(),
# and rebuild the PeakMatrix so the project is ready for getMarkerFeatures().
# ─────────────────────────────────────────────────────────────────────────────
loadPeakSetAndMatrix <- function(proj, groupBy) {
  peakSetPath <- file.path(peakCallDir, paste0("PeakSet_gr_", groupBy, ".rds"))
  print(paste("Loading peakSet for groupBy:", groupBy))
  print(paste("  Path:", peakSetPath))
  peakSet <- readRDS(peakSetPath)
  print(paste("  Peaks loaded:", length(peakSet)))
  print(paste("  Peak types:"))
  print(table(peakSet$peakType))
  proj <- addPeakSet(ArchRProj = proj, peakSet = peakSet, force = TRUE)
  proj <- addPeakMatrix(proj)
  print(paste("  PeakMatrix rebuilt for groupBy:", groupBy))
  print(paste("  Available matrices:", paste(getAvailableMatrices(proj), collapse = ", ")))
  return(proj)
}

# ─────────────────────────────────────────────────────────────────────────────
# Pairwise tests — PIMO_Region peak set
# ─────────────────────────────────────────────────────────────────────────────
print("=== Loading PIMO_Region peak set for pairwise tests ===")
proj_hyp2 <- loadPeakSetAndMatrix(proj_hyp2, "PIMO_Region")

# 1. PIMOup_EB vs PIMOup_TC
print("Pairwise test: PIMOup_EB vs PIMOup_TC")
markerTest <- getMarkerFeatures(
  ArchRProj  = proj_hyp2,
  useMatrix  = "PeakMatrix",
  groupBy    = "PIMO_Region",
  testMethod = "wilcoxon",
  bias       = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
  useGroups  = "PIMOup_EB",
  bgdGroups  = "PIMOup_TC"
)
saveRDS(markerTest,
  file = file.path(peakCallDir, "markersTest_PIMO_Region_PIMOup_EB_vs_PIMOup_TC.rds"))
markerTest_GR <- getMarkers(markerTest, cutOff = "FDR <= 1 & abs(Log2FC) >= 0", returnGR = TRUE)
saveRDS(markerTest_GR,
  file = file.path(peakCallDir, "markersTest_GR_PIMO_Region_PIMOup_EB_vs_PIMOup_TC.rds"))
print("Number of markerTest peaks identified per group:")
print(sapply(markerTest_GR, length))

# 2. PIMOup_EB vs PIMOdown_EB
print("Pairwise test: PIMOup_EB vs PIMOdown_EB")
markerTest <- getMarkerFeatures(
  ArchRProj  = proj_hyp2,
  useMatrix  = "PeakMatrix",
  groupBy    = "PIMO_Region",
  testMethod = "wilcoxon",
  bias       = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
  useGroups  = "PIMOup_EB",
  bgdGroups  = "PIMOdown_EB"
)
saveRDS(markerTest,
  file = file.path(peakCallDir, "markersTest_PIMO_Region_PIMOup_EB_vs_PIMOdown_EB.rds"))
markerTest_GR <- getMarkers(markerTest, cutOff = "FDR <= 1 & abs(Log2FC) >= 0", returnGR = TRUE)
saveRDS(markerTest_GR,
  file = file.path(peakCallDir, "markersTest_GR_PIMO_Region_PIMOup_EB_vs_PIMOdown_EB.rds"))
print("Number of markerTest peaks identified per group:")
print(sapply(markerTest_GR, length))

# 3. PIMOup_TC vs PIMOdown_TC
print("Pairwise test: PIMOup_TC vs PIMOdown_TC")
markerTest <- getMarkerFeatures(
  ArchRProj  = proj_hyp2,
  useMatrix  = "PeakMatrix",
  groupBy    = "PIMO_Region",
  testMethod = "wilcoxon",
  bias       = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
  useGroups  = "PIMOup_TC",
  bgdGroups  = "PIMOdown_TC"
)
saveRDS(markerTest,
  file = file.path(peakCallDir, "markersTest_PIMO_Region_PIMOup_TC_vs_PIMOdown_TC.rds"))
markerTest_GR <- getMarkers(markerTest, cutOff = "FDR <= 1 & abs(Log2FC) >= 0", returnGR = TRUE)
saveRDS(markerTest_GR,
  file = file.path(peakCallDir, "markersTest_GR_PIMO_Region_PIMOup_TC_vs_PIMOdown_TC.rds"))
print("Number of markerTest peaks identified per group:")
print(sapply(markerTest_GR, length))

# 4. PIMOup_PT vs PIMOdown_PT
print("Pairwise test: PIMOup_PT vs PIMOdown_PT")
markerTest <- getMarkerFeatures(
  ArchRProj  = proj_hyp2,
  useMatrix  = "PeakMatrix",
  groupBy    = "PIMO_Region",
  testMethod = "wilcoxon",
  bias       = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
  useGroups  = "PIMOup_PT",
  bgdGroups  = "PIMOdown_PT"
)
saveRDS(markerTest,
  file = file.path(peakCallDir, "markersTest_PIMO_Region_PIMOup_PT_vs_PIMOdown_PT.rds"))
markerTest_GR <- getMarkers(markerTest, cutOff = "FDR <= 1 & abs(Log2FC) >= 0", returnGR = TRUE)
saveRDS(markerTest_GR,
  file = file.path(peakCallDir, "markersTest_GR_PIMO_Region_PIMOup_PT_vs_PIMOdown_PT.rds"))
print("Number of markerTest peaks identified per group:")
print(sapply(markerTest_GR, length))

# ─────────────────────────────────────────────────────────────────────────────
# Pairwise tests — Region.annotation peak set
# ─────────────────────────────────────────────────────────────────────────────
print("=== Loading Region.annotation peak set for pairwise tests ===")
proj_hyp2 <- loadPeakSetAndMatrix(proj_hyp2, "Region.annotation")

# 1. EB vs TC
print("Pairwise test: EB vs TC")
markerTest <- getMarkerFeatures(
  ArchRProj  = proj_hyp2,
  useMatrix  = "PeakMatrix",
  groupBy    = "Region.annotation",
  testMethod = "wilcoxon",
  bias       = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
  useGroups  = "EB",
  bgdGroups  = "TC"
)
saveRDS(markerTest,
  file = file.path(peakCallDir, "markersTest_Region.annotation_EB_vs_TC.rds"))
markerTest_GR <- getMarkers(markerTest, cutOff = "FDR <= 1 & abs(Log2FC) >= 0", returnGR = TRUE)
saveRDS(markerTest_GR,
  file = file.path(peakCallDir, "markersTest_GR_Region.annotation_EB_vs_TC.rds"))
print("Number of markerTest peaks identified per group:")
print(sapply(markerTest_GR, length))

# 2. EB vs PT
print("Pairwise test: EB vs PT")
markerTest <- getMarkerFeatures(
  ArchRProj  = proj_hyp2,
  useMatrix  = "PeakMatrix",
  groupBy    = "Region.annotation",
  testMethod = "wilcoxon",
  bias       = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
  useGroups  = "EB",
  bgdGroups  = "PT"
)
saveRDS(markerTest,
  file = file.path(peakCallDir, "markersTest_Region.annotation_EB_vs_PT.rds"))
markerTest_GR <- getMarkers(markerTest, cutOff = "FDR <= 1 & abs(Log2FC) >= 0", returnGR = TRUE)
saveRDS(markerTest_GR,
  file = file.path(peakCallDir, "markersTest_GR_Region.annotation_EB_vs_PT.rds"))
print("Number of markerTest peaks identified per group:")
print(sapply(markerTest_GR, length))

# 3. TC vs PT
print("Pairwise test: TC vs PT")
markerTest <- getMarkerFeatures(
  ArchRProj  = proj_hyp2,
  useMatrix  = "PeakMatrix",
  groupBy    = "Region.annotation",
  testMethod = "wilcoxon",
  bias       = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
  useGroups  = "TC",
  bgdGroups  = "PT"
)
saveRDS(markerTest,
  file = file.path(peakCallDir, "markersTest_Region.annotation_TC_vs_PT.rds"))
markerTest_GR <- getMarkers(markerTest, cutOff = "FDR <= 1 & abs(Log2FC) >= 0", returnGR = TRUE)
saveRDS(markerTest_GR,
  file = file.path(peakCallDir, "markersTest_GR_Region.annotation_TC_vs_PT.rds"))
print("Number of markerTest peaks identified per group:")
print(sapply(markerTest_GR, length))

# ─────────────────────────────────────────────────────────────────────────────
# Save the project with the last-loaded peak set (Region.annotation) in place
# ─────────────────────────────────────────────────────────────────────────────
# print("Saving ArchR project")
# proj_hyp2 <- saveArchRProject(
#   ArchRProj       = proj_hyp2,
#   outputDirectory = "Gaiti_multiome_harmony_merged_malig_peak",
#   load            = TRUE
# )

EOF

echo "Pairwise testing analysis completed"