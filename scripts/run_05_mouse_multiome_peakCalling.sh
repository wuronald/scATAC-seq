#!/bin/bash
#SBATCH --job-name=run_05_mouse_multiome_peakCalling
#SBATCH --output=run_05_mouse_multiome_peakCalling_%j.out
#SBATCH --error=run_05_mouse_multiome_peakCalling_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=60G
#SBATCH --time=02:00:00

# Load necessary modules (adjust as needed for your system)
module load R/4.4.1
module load python3/3.7.2
module load MACS/2.2.7.1

# Run R script
Rscript - <<'EOF'

# load libraries
library(ArchR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(here)
set.seed(1)

# Set the number of threads for ArchR
addArchRThreads(threads = 18)

# set genome to mm10
addArchRGenome("mm10")

# Load the project
proj_hyp <- loadArchRProject(path = "mouse_multiome_harmony_merged_subset")

# subset to only malignant cells
print("Subsetting ArchR object to only malignant cells: ")
malignant_cells <- proj_hyp$cellNames[which(proj_hyp$Azimuth_class == "Malignant")]
print(paste("Number of malignant cells:", length(malignant_cells)))

proj_hyp <- subsetArchRProject(
  ArchRProj = proj_hyp,
  cells = malignant_cells,
  outputDirectory = "mouse_multiome_harmony_merged_malig_peak_subset",
  dropCells = TRUE,
  force = TRUE
  )

# Save the project to a new directory
# proj_hyp <- saveArchRProject(ArchRProj = proj_hyp, outputDirectory = "mouse_multiome_harmony_merged_malig_peak_subset", load = TRUE)

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
groupBy_list <- c("PIMO_pos", "hybrid_pair")
# groupBy_list <- c("hybrid_pair")
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
  
# TEMPORARY FOR DEBUGGING: Save the project prior to where error occurs
#print("Saving ArchR project before getMarkerFeatures step")
#proj_hyp2 <- saveArchRProject(ArchRProj = proj_hyp2, outputDirectory = "mouse_multiome_harmony_merged_malig_peak_subset", load = TRUE)
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
    bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
    # bias = c("TSSEnrichment", "nFrags", "Gex_nUMI"),  
    testMethod = "wilcoxon"
  )
  
  # Save the markersPeaks SE object for easier load in the future
  print(paste("Saved markersPeaks as .RData for", groupBy))
  save(markersPeaks, file = paste0("mouse_multiome_harmony_merged_malig_peak_subset/PeakCalls/markersPeaks_", groupBy, ".RData"))

}
# Save the project after each groupBy
proj_hyp2 <- saveArchRProject(ArchRProj = proj_hyp2, outputDirectory = "mouse_multiome_harmony_merged_malig_peak_subset", load = TRUE)
EOF

echo "Peak Calling analysis completed"