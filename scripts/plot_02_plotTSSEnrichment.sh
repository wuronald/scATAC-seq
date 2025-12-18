#!/bin/bash
#SBATCH --job-name=plot_02_plotTSSEnrichment
#SBATCH --output=plot_02_plotTSSEnrichment_%j.out
#SBATCH --error=plot_02_plotTSSEnrichment_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=6
#SBATCH --mem=58G
#SBATCH --time=04:00:00

# Parse command line arguments
# First argument: ARCHR_PROJECT_PATH - path to the ArchR project
# Second argument: GROUP_BY - groupBy variable for TSS enrichment plotting (e.g., "PIMO_up_status", "hybrid_pair", "PIMO_Region")
# Example usage:
# sbatch scripts/plot_02_plotTSSEnrichment.sh "human_multiome_harmony_merged_malig_peak" "PIMO_up_status"
# sbatch scripts/plot_02_plotTSSEnrichment.sh "Gaiti_multiome_harmony_merged_malig_peak" "PIMO_Region"

ARCHR_PROJECT_PATH="${1:-human_multiome_harmony_merged_malig_peak}" # Default to "human_multiome_harmony_merged_malig_peak" if no argument provided
GROUP_BY="${2:-PIMO_up_status}"  # Default to "PIMO_up_status" if no argument provided

# Export the parameters so R can access them
export ARCHR_PROJECT_PATH
export GROUP_BY
echo "Running TSS enrichment plotting with:"
echo "  ArchR project path: ${ARCHR_PROJECT_PATH}"
echo "  groupBy: ${GROUP_BY}"

# Load necessary modules (adjust as needed for your system)
module load R/4.4.1

# Run R script
Rscript - <<'EOF'

# load libraries
library(ArchR)
library(here)
set.seed(1)

# Get parameters from environment variables
archrProjectPath <- Sys.getenv("ARCHR_PROJECT_PATH", unset = "human_multiome_harmony_merged_malig_peak")
groupBy <- Sys.getenv("GROUP_BY", unset = "PIMO_up_status")
cat("Using ArchR project path:", archrProjectPath, "\n")
cat("Using groupBy:", groupBy, "\n")

# Set the number of threads for ArchR
addArchRThreads(threads = 6)

# Load the project
proj <- loadArchRProject(path = archrProjectPath)

# Plot fragment size distribution
fragsizes <- plotFragmentSizes(ArchRProj = proj)

# Plot TSS enrichment by group
tssPlot <- plotTSSEnrichment(ArchRProj = proj, groupBy = groupBy)

# Save the plot
print("Saving TSS enrichment plot and fragment size distribution plot as PDF")
plotPDF(fragsizes,tssPlot, 
        name = paste0(groupBy,"-QC-Sample-FragSizes-TSS_Enrichment_Plot.pdf"), 
        ArchRProj = proj, 
        addDOC = TRUE, 
        width = 5, 
        height = 5
        )
EOF