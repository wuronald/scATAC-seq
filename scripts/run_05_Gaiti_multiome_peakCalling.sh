#!/bin/bash
#SBATCH --job-name=run_05_Gaiti_multiome_peakCalling
#SBATCH --output=run_05_Gaiti_multiome_peakCalling_%j.out
#SBATCH --error=run_05_Gaiti_multiome_peakCalling_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=58G
#SBATCH --time=06:00:00

# Load necessary modules (adjust as needed for your system)
module load R/4.4.1
module load python3/3.7.2
module load MACS/2.2.7.1

# Run R script
Rscript - <<'EOF'

# load libraries
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(here)
set.seed(1)

# Set the number of threads for ArchR
addArchRThreads(threads = 6)

# set genome to hg38
addArchRGenome("hg38")

# Load the project
proj_hyp <- loadArchRProject(path = "Gaiti_multiome_harmony_merged")

# subset to only malignant cells
print("Subsetting ArchR object to only malignant cells and not mislabelled cells")
malignant_cells <- proj_hyp$cellNames[which(proj_hyp$Azimuth_class == "Malignant" & proj_hyp$mislabelled == "Correct")]
print(paste("Number of malignant non-ambiguous cells:", length(malignant_cells)))

# add additional metadata columns if not already present
# 1. PIMO_Region: is a combination of PIMO_up_status and Region.annotation columns
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
  outputDirectory = "Gaiti_multiome_harmony_merged_malig_peak",
  dropCells = TRUE,
  force = TRUE
  )

# Find MACS2 path and check if it exists
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

# Set of ArchR object metadata columns to iterate over for peak calling
# groupBy_list <- c("hybrid_pair", "PIMO_up_status","PIMO_Region") # hybrid_pair has NAs causing issues atm
groupBy_list <- c("PIMO_up_status","PIMO_Region","Region.annotation")

print(paste("GroupBy list:", paste(groupBy_list, collapse = ", ")))

for (groupBy in groupBy_list) {
  print(paste("Processing groupBy:", groupBy))

  # Ensure grouping column is a factor
  # proj_hyp@cellColData[[groupBy]] <- as.factor(proj_hyp@cellColData[[groupBy]])

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

  # Save the peakSet gr object for easier load in the future
  print(paste("Saved peakSet as .RData for", groupBy))
  saveRDS(myPeakSet, file = paste0("Gaiti_multiome_harmony_merged_malig_peak/PeakCalls/PeakSet_gr_", groupBy, ".rds"))
  
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
# Get PeakMatrix as a matrix object
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
  
# Save the markersPeaks SE object for easier load in the future
print(paste("Saved markersPeaks as .rds for", groupBy))
saveRDS(markersPeaks, file = paste0("Gaiti_multiome_harmony_merged_malig_peak/PeakCalls/markersPeaks_", groupBy, ".rds"))

# Save the markersPeaks GR object for easier load in the future
print(paste("Saved markersPeaks GR as .rds for", groupBy))
markersPeaks_GR <- getMarkers(markersPeaks, cutOff = "FDR <=  1 & abs(Log2FC) >= 0", returnGR = TRUE)
saveRDS(markersPeaks_GR, file = paste0("Gaiti_multiome_harmony_merged_malig_peak/PeakCalls/markersPeaks_GR_", groupBy, ".rds"))

}

# pairwise test between PIMO_Region groups: 
# 1. Within PIMOup EB vs TC
print("Pairwise test between PIMO_Region groups: PIMOup_EB vs PIMOup_TC")
markerTest <- getMarkerFeatures(
  ArchRProj = proj_hyp2, 
  useMatrix = "PeakMatrix",
  groupBy = "PIMO_Region",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
  useGroups = "PIMOup_EB",
  bgdGroups = "PIMOup_TC"
)
# Save the markerTest SE object for easier load in the future
print(paste("Saved markerTest as .rds for PIMO_Region PIMOup_EB vs PIMOup_TC" ))
saveRDS(markerTest, file = paste0("Gaiti_multiome_harmony_merged_malig_peak/PeakCalls/markersTest_PIMO_Region_", "PIMOup_EB_vs_PIMOup_TC", ".rds"))

# Extract markerTest GR object for PIMO_Region groups
print("Extracting markerTestGR object for PIMO_Region")
markerTest_GR <- getMarkers(markerTest, cutOff = "FDR <= 1 & abs(Log2FC) >= 0", returnGR = TRUE)

# Save the markerTest GR object for easier load in the future
print(paste("Saved markerTest GR as .rds for PIMO_Region PIMOup_EB vs PIMOup_TC" ))
saveRDS(markerTest_GR, file = "Gaiti_multiome_harmony_merged_malig_peak/PeakCalls/markersTest_GR_PIMO_Region_PIMOup_EB_vs_PIMOup_TC.rds")

print("Number of markerTest peaks identified per group:")
print(sapply(markerTest_GR, length))

# 2. PIMOup vs PIMOdown for EB regions only
print("Pairwise test between PIMO_Region groups: PIMOup_EB vs PIMOdown_EB")
markerTest <- getMarkerFeatures(
  ArchRProj = proj_hyp2, 
  useMatrix = "PeakMatrix",
  groupBy = "PIMO_Region",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
  useGroups = "PIMOup_EB",
  bgdGroups = "PIMOdown_EB"
)
# Save the markerTest SE object for easier load in the future
print(paste("Saved markerTest as .rds for PIMO_Region :", "PIMOup_EB vs PIMOdown_EB" ))
saveRDS(markerTest, file = paste0("Gaiti_multiome_harmony_merged_malig_peak/PeakCalls/markersTest_PIMO_Region_", "PIMOup_EB_vs_PIMOdown_EB", ".rds"))

# Extract markerTest GR object for PIMO_Region groups
print("Extracting markerTestGR object for PIMO_Region")
markerTest_GR <- getMarkers(markerTest, cutOff = "FDR <= 1 & abs(Log2FC) >= 0", returnGR = TRUE)

# Save the markerTest SE object for easier load in the future
print(paste("Saved markerTest as .rds for PIMO_Region PIMOup_EB vs PIMOdown_EB" ))
saveRDS(markerTest_GR, file = "Gaiti_multiome_harmony_merged_malig_peak/PeakCalls/markersTest_GR_PIMO_Region_PIMOup_EB_vs_PIMOdown_EB.rds")

print("Number of markerTest peaks identified per group:")
print(sapply(markerTest_GR, length))

# 3. PIMOup vs PIMOdown for TC regions only
print("Pairwise test between PIMO_Region groups: PIMOup_TC vs PIMOdown_TC")
markerTest <- getMarkerFeatures(
  ArchRProj = proj_hyp2, 
  useMatrix = "PeakMatrix",
  groupBy = "PIMO_Region",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
  useGroups = "PIMOup_TC",
  bgdGroups = "PIMOdown_TC"
)
# Save the markerTest SE object for easier load in the future
print(paste("Saved markerTest as .rds for PIMO_Region :", "PIMOup_TC vs PIMOdown_TC" ))
saveRDS(markerTest, file = paste0("Gaiti_multiome_harmony_merged_malig_peak/PeakCalls/markersTest_PIMO_Region_", "PIMOup_TC_vs_PIMOdown_TC", ".rds"))

# extract markerTest GR object for PIMO_Region groups
print("Extracting markerTestGR object for PIMO_Region")
markerTest_GR <- getMarkers(markerTest, cutOff = "FDR <= 1 & abs(Log2FC) >= 0", returnGR = TRUE)

# Save the markerTest SE object for easier load in the future
print(paste("Saved markerTest as .rds for PIMO_Region PIMOup_TC vs PIMOdown_TC" ))
saveRDS(markerTest_GR, file = "Gaiti_multiome_harmony_merged_malig_peak/PeakCalls/markersTest_GR_PIMO_Region_PIMOup_TC_vs_PIMOdown_TC.rds")

print("Number of markerTest peaks identified per group:")
print(sapply(markerTest_GR, length))

# 4. PIMOup vs PIMOdown for PT regions only
print("Pairwise test between PIMO_Region groups: PIMOup_PT vs PIMOdown_PT")
markerTest <- getMarkerFeatures(
  ArchRProj = proj_hyp2, 
  useMatrix = "PeakMatrix",
  groupBy = "PIMO_Region",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
  useGroups = "PIMOup_PT",
  bgdGroups = "PIMOdown_PT"
)
# Save the markerTest SE object for easier load in the future
print(paste("Saved markerTest as .rds for PIMO_Region :", "PIMOup_PT vs PIMOdown_PT" ))
saveRDS(markerTest, file = paste0("Gaiti_multiome_harmony_merged_malig_peak/PeakCalls/markersTest_PIMO_Region_", "PIMOup_PT_vs_PIMOdown_PT", ".rds"))

# extract markerTest GR object for PIMO_Region groups
print("Extracting markerTestGR object for PIMO_Region")
markerTest_GR <- getMarkers(markerTest, cutOff = "FDR <= 1 & abs(Log2FC) >= 0", returnGR = TRUE)

# Save the markerTest SE object for easier load in the future
print(paste("Saved markerTest as .rds for PIMO_Region PIMOup_PT vs PIMOdown_PT" ))
saveRDS(markerTest_GR, file = "Gaiti_multiome_harmony_merged_malig_peak/PeakCalls/markersTest_GR_PIMO_Region_PIMOup_PT_vs_PIMOdown_PT.rds")

print("Number of markerTest peaks identified per group:")
print(sapply(markerTest_GR, length))

# pairwise test between Region.annotation groups:
# 1. EB vs TC
print("Pairwise test between Region.annotation groups: EB vs TC") 
markerTest <- getMarkerFeatures(
  ArchRProj = proj_hyp2, 
  useMatrix = "PeakMatrix",
  groupBy = "Region.annotation",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
  useGroups = "EB",
  bgdGroups = "TC"
)
# Save the markerTest SE object for easier load in the future
print(paste("Saved markerTest as .rds for Region.annotation :", "EB vs TC" ))
saveRDS(markerTest, file = paste0("Gaiti_multiome_harmony_merged_malig_peak/PeakCalls/markersTest_Region.annotation_", "EB_vs_TC", ".rds"))
# extract markerTest GR object for Region.annotation groups
print("Extracting markerTestGR object for Region.annotation")
markerTest_GR <- getMarkers(markerTest, cutOff = "FDR <=  1 & abs(Log2FC) >= 0", returnGR = TRUE)
# Save the markerTest SE object for easier load in the future
print(paste("Saved markerTest as .rds for Region.annotation EB vs TC" ))
saveRDS(markerTest_GR, file = "Gaiti_multiome_harmony_merged_malig_peak/PeakCalls/markersTest_GR_Region.annotation_EB_vs_TC.rds")
print("Number of markerTest peaks identified per group:")
print(sapply(markerTest_GR, length))

# 2. EB vs PT
print("Pairwise test between Region.annotation groups: EB vs PT")
markerTest <- getMarkerFeatures(
  ArchRProj = proj_hyp2, 
  useMatrix = "PeakMatrix",
  groupBy = "Region.annotation",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
  useGroups = "EB",
  bgdGroups = "PT"
)
# Save the markerTest SE object for easier load in the future
print(paste("Saved markerTest as .rds for Region.annotation :", "EB vs PT"))
saveRDS(markerTest, file = paste0("Gaiti_multiome_harmony_merged_malig_peak/PeakCalls/markersTest_Region.annotation_EB_vs_PT.rds"))

# extract markerTest GR object for Region.annotation groups
print("Extracting markerTestGR object for Region.annotation")
markerTest_GR <- getMarkers(markerTest, cutOff = "FDR <= 1 & abs(Log2FC) >= 0", returnGR = TRUE)

# Save the markerTest SE object for easier load in the future
print(paste("Saved markerTest as .rds for Region.annotation EB vs PT"))
saveRDS(markerTest_GR, file = "Gaiti_multiome_harmony_merged_malig_peak/PeakCalls/markersTest_GR_Region.annotation_EB_vs_PT.rds")

print("Number of markerTest peaks identified per group:")
print(sapply(markerTest_GR, length))

# 3. TC vs PT
print("Pairwise test between Region.annotation groups: TC vs PT")
markerTest <- getMarkerFeatures(
  ArchRProj = proj_hyp2, 
  useMatrix = "PeakMatrix",
  groupBy = "Region.annotation",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
  useGroups = "TC",
  bgdGroups = "PT"
)
# Save the markerTest SE object for easier load in the future
print(paste("Saved markerTest as .rds for Region.annotation :", "TC vs PT"))
saveRDS(markerTest, file = paste0("Gaiti_multiome_harmony_merged_malig_peak/PeakCalls/markersTest_Region.annotation_TC_vs_PT.rds"))
# extract markerTest GR object for Region.annotation groups
print("Extracting markerTestGR object for Region.annotation")
markerTest_GR <- getMarkers(markerTest, cutOff = "FDR <= 1 & abs(Log2FC) >= 0", returnGR = TRUE)
# Save the markerTest SE object for easier load in the future
print(paste("Saved markerTest as .rds for Region.annotation TC vs PT"))
saveRDS(markerTest_GR, file = "Gaiti_multiome_harmony_merged_malig_peak/PeakCalls/markersTest_GR_Region.annotation_TC_vs_PT.rds")
print("Number of markerTest peaks identified per group:")
print(sapply(markerTest_GR, length))

# Save the project
proj_hyp2 <- saveArchRProject(ArchRProj = proj_hyp2, outputDirectory = "Gaiti_multiome_harmony_merged_malig_peak", load = TRUE)
EOF

echo "Peak Calling analysis completed"