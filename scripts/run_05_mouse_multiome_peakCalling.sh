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
library(here)
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
set.seed(1)

# Set the number of threads for ArchR
addArchRThreads(threads = 18)

# set genome to hg38
addArchRGenome("hg38")

# Load the project
proj_hyp <- loadArchRProject(path = "mouse_multiome_harmony_merged_subset")

# Save the project to a new directory
proj_hyp <- saveArchRProject(ArchRProj = proj_hyp, outputDirectory = "mouse_multiome_harmony_merged_subset_peak", load = TRUE)

# Save pseudobulk info as new ArchRProject object
print("Creating pseudobulk based on: PIMO_pos")

proj_hyp$PIMO_pos <- as.factor(proj_hyp$PIMO_pos)
proj_hyp <- addGroupCoverages(ArchRProj = proj_hyp,
                           groupBy = "PIMO_pos" # group for which the pseudobulk rep is made
                           )

# MACS2 Peak Calling

# # PATH to MACS2
# #pathToMacs2 <- "/cluster/tools/software/centos7/MACS/2.2.5/bin/macs2"
# #pathToMacs2 <- "/cluster/home/rwu/.local/lib/python3.7/site-packages/MACS2"
# pathToMacs2 <- "/cluster/home/rwu/.local/bin/macs2"

pathToMacs2 <- findMacs2()

# # call peaks
print("Call Peaks with MACS2")
proj_hyp2 <- addReproduciblePeakSet(ArchRProj = proj_hyp,
                        groupBy = "PIMO_pos",
                        pathToMacs2 = pathToMacs2)


print("get Peaks set")
# get the peak set
myPeakSet <- getPeakSet(proj_hyp2)

# tabulate type of peaks
print("tabulate type of peaks")
table(myPeakSet$peakType)

# check available matrices
print("check available matrices")
getAvailableMatrices(proj_hyp2)

# add peak Matrix
print("adding peak matrix")
proj_hyp2 <- addPeakMatrix(proj_hyp2)

# check available matrices
print("check available matrices after adding peak matrix")
getAvailableMatrices(proj_hyp2)

# fix incompatible dimensions error (change factor back to logical)
print("fix incompatible dimensions error")
proj_hyp2$proj_hyp2 <- as.logical(proj_hyp2$proj_hyp2)

# get marker peaks for PIMO_pos groups
print("Get marker peaks for PIMO_pos groups")

markersPeaks_PIMO_pos <- getMarkerFeatures(
  ArchRProj = proj_hyp2,
  useMatrix = "PeakMatrix",
  groupBy = "PIMO_pos",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

# save the markersPeaks SE object for easier load in the future
save(markersPeaks_PIMO_pos, file = "mouse_multiome_harmony_merged_subset_peak/PeakCalls/markersPeaks_PIMO_pos.RData")

# correct groupCoverages paths
print("output groupCoverages paths:")
print(proj_hyp2@projectMetadata$GroupCoverages[[1]]$coverageMetadata$File)

# my_dir <- "/cluster/projects/wouterslab/ArchR103_4/Terekhanova_hypoxia"
# old_project_dir <- "/cluster/projects/wouterslab/ArchR103_4/Terekhanova"

# old_paths <- proj_hyp2@projectMetadata$GroupCoverages[[1]]$coverageMetadata$File
# new_paths <- gsub(old_project_dir, my_dir, old_paths)
# proj_hyp2@projectMetadata$GroupCoverages[[1]]$coverageMetadata$File <- new_paths

print("new groupCoverages paths:")
#print(new_paths)
print("confirm assignment of new paths:")
print(proj_hyp2@projectMetadata$GroupCoverages[[1]]$coverageMetadata$File)

# Save the project
print("Saving the project")

proj_hyp2 <- saveArchRProject(ArchRProj = proj_hyp2, outputDirectory = "mouse_multiome_harmony_merged_subset_peak", load = TRUE)
print(paste0("current output dir:",getOutputDirectory(proj_hyp2)))

EOF

echo "Peak Calling analysis completed"