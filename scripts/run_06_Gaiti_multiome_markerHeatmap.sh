#!/bin/bash
#SBATCH --job-name=run_06_Gaiti_multiome_markerHeatmap
#SBATCH --output=run_06_Gaiti_multiome_markerHeatmap_%j.out
#SBATCH --error=run_06_Gaiti_multiome_markerHeatmap_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=16G
#SBATCH --time=01:00:00

# Parse command line arguments
# First argument: GROUPBY - which groupBy condition to analyze (PIMO_up_status or PIMO_Region)
# Second argument: ANALYSIS_TYPE - either "markers" for general markers or "pairwise" for pairwise test
# Third argument: COMPARISON - for pairwise tests, specify the comparison (optional, only used if ANALYSIS_TYPE="pairwise")
# Fourth argument: GENES - genes for plotting browser tracks
# 
# Example usage for general markers:
# sbatch scripts/run_06_Gaiti_multiome_markerHeatmap.sh "PIMO_up_status" "markers"
# sbatch scripts/run_06_Gaiti_multiome_markerHeatmap.sh "PIMO_Region" "markers"
#
# Example usage for pairwise tests:
# sbatch scripts/run_06_Gaiti_multiome_markerHeatmap.sh "PIMO_Region" "pairwise" "PIMOup_EB_vs_PIMOup_TC"
# sbatch scripts/run_06_Gaiti_multiome_markerHeatmap.sh "PIMO_Region" "pairwise" "PIMOup_EB_vs_PIMOdown_EB"
# sbatch scripts/run_06_Gaiti_multiome_markerHeatmap.sh "PIMO_Region" "pairwise" "PIMOup_TC_vs_PIMOdown_TC"
# sbatch scripts/run_06_Gaiti_multiome_markerHeatmap.sh "PIMO_Region" "pairwise" "PIMOup_PT_vs_PIMOdown_PT"
#
# Example with custom genes:
# sbatch scripts/run_06_Gaiti_multiome_markerHeatmap.sh "PIMO_Region" "pairwise" "PIMOup_EB_vs_PIMOup_TC" "HIF1A,CA9,VEGFA"

GROUPBY="${1:-PIMO_up_status}"
ANALYSIS_TYPE="${2:-markers}"
COMPARISON="${3:-}"
GENES="${4:-CTSB,OLIG1,OLIG2,SOX2,CD109,CD44,RND3,STMN2,NGFR,SOX10,ID1,ID2,ID3,CA9,VEGFA,SLC2A1,EMX2,HES1,NR4A2,NR4A3,KLF9,NRN1,DNER,DPYSL4,RND3,CADM3,PLPPR3,SLC2A3,SLC5A3,PCP4,SOX12,SV2A,VAMP1,PTPRN,PPFIA3,PPP1R1A,HOMER2,NMB,ADM,ENO2,CA11,HMOX1,CRYAB,TIPARP,SOD2,XBP1,ATF3,DDIT4,DDIT4L}"

# Export the parameters so R can access them
export GROUPBY
export ANALYSIS_TYPE
export COMPARISON
export GENES

echo "Running marker heatmap analysis with:"
echo "  groupBy: ${GROUPBY}"
echo "  analysis_type: ${ANALYSIS_TYPE}"
echo "  comparison: ${COMPARISON}"
echo "  genes: ${GENES}"

# Load necessary modules (adjust as needed for your system)
module load R/4.4.1

# Run R script
Rscript - <<'EOF'

# load libraries
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(here)
set.seed(1)

# Get parameters from environment variables
groupBy <- Sys.getenv("GROUPBY", unset = "PIMO_up_status")
analysis_type <- Sys.getenv("ANALYSIS_TYPE", unset = "markers")
comparison <- Sys.getenv("COMPARISON", unset = "")
genes <- Sys.getenv("GENES", unset = "CTSB,OLIG1,OLIG2,SOX2,CD109,CD44,RND3,STMN2,NGFR,SOX10,ID1,ID2,ID3,CA9,VEGFA,SLC2A1,EMX2,HES1,NR4A2,NR4A3,KLF9,NRN1,DNER,DPYSL4,RND3,CADM3,PLPPR3,SLC2A3,SLC5A3,PCP4,SOX12,SV2A,VAMP1,PTPRN,PPFIA3,PPP1R1A,HOMER2,NMB,ADM,ENO2,CA11,HMOX1,CRYAB,TIPARP,SOD2,XBP1,ATF3,DDIT4,DDIT4L")

# Convert comma-separated string to vector
genes <- unlist(strsplit(genes, ","))

print(paste("Analyzing groupBy:", groupBy))
print(paste("Analysis type:", analysis_type))
if (comparison != "") {
    print(paste("Comparison:", comparison))
}
print(paste("Number of genes:", length(genes)))

# Set the number of threads for ArchR
addArchRThreads(threads = 18)

# set genome to hg38
addArchRGenome("hg38")

# Load the project
proj_hyp <- loadArchRProject(path = "Gaiti_multiome_harmony_merged_malig_peak")

# Determine which markersPeaks RDS file to load based on analysis type
if (analysis_type == "pairwise") {
    # For pairwise tests
    if (comparison == "") {
        stop("COMPARISON argument is required when ANALYSIS_TYPE='pairwise'.\n",
             "Available comparisons:\n",
             "  - PIMOup_EB_vs_PIMOup_TC\n",
             "  - PIMOup_EB_vs_PIMOdown_EB\n",
             "  - PIMOup_TC_vs_PIMOdown_TC\n",
             "  - PIMOup_PT_vs_PIMOdown_PT")
    }
    
    markersPeaks_rds <- paste0("Gaiti_multiome_harmony_merged_malig_peak/PeakCalls/markersTest_", 
                                groupBy, "_", comparison, ".rds")
    output_suffix <- paste0(groupBy, "_", comparison)
} else {
    # For general markers across all groups
    markersPeaks_rds <- paste0("Gaiti_multiome_harmony_merged_malig_peak/PeakCalls/markersPeaks_", 
                                groupBy, ".rds")
    output_suffix <- paste0(groupBy, "_markers")
}

print(paste("Loading markers from:", markersPeaks_rds))

# Load markersPeaks
if (file.exists(markersPeaks_rds)) {
    print("Loading previously saved markersPeaks")
    markersPeaks <- readRDS(file = markersPeaks_rds)
} else {
    stop(paste("File not found:", markersPeaks_rds, 
               "\nPlease run the peak calling script first to generate this file."))
}

# Check if markersPeaks is a SummarizedExperiment object
print("Checking if markersPeaks is a SummarizedExperiment object")
if (!is(markersPeaks, "SummarizedExperiment")) {
    stop("markersPeaks is not a SummarizedExperiment object. Please check your markersPeaks data.")
} else {
    print("markersPeaks is a valid SummarizedExperiment object.")
}

# Plot marker heatmap
print(paste("Plotting marker heatmap for", output_suffix))
cutOff <- "FDR <= 0.1 & abs(Log2FC) >= 0.5"
print(paste("Using cutOff:", cutOff))

heatmap <- plotMarkerHeatmap(
    seMarker = markersPeaks,
    cutOff = cutOff,
    limits = c(-2, 2),
    plotLog2FC = TRUE,
    transpose = TRUE,
    returnMatrix = FALSE
)

# Save the heatmap as a PDF
print("Saving heatmap as PDF")
plotPDF(heatmap, 
        name = paste0("markerPeaks-Heatmap_", output_suffix), 
        width = 8, 
        height = 6, 
        ArchRProj = proj_hyp, 
        addDOC = TRUE)

# Plot MA and Volcano plots for marker peaks
print("Plotting MA and Volcano plots for marker peaks")

# Iterate over each group in the markersPeaks object
for (group in colnames(markersPeaks)) {
    print(paste("Creating MA plot for group:", group))
    pma <- plotMarkers(
        seMarker = markersPeaks,
        name = group,
        cutOff = cutOff,
        plotAs = "MA"
    )
    
    print(paste("Creating Volcano plot for group:", group))
    pvolcano <- plotMarkers(
        seMarker = markersPeaks,
        name = group,
        cutOff = cutOff,
        plotAs = "Volcano"
    )
    
    # Save each MA and Volcano plot as a PDF
    print(paste("Saving plots for group:", group))
    plotPDF(pma, pvolcano, 
            name = paste0("markerPeaks-MAplot_", output_suffix, "_", group), 
            width = 5, 
            height = 5, 
            ArchRProj = proj_hyp, 
            addDOC = TRUE)
}
print("All MA and Volcano plots created and saved.")

# Plot marker peaks in browser tracks for genes of interest
print("Plotting marker peaks in browser tracks for genes of interest")
print(paste("Genes of interest:", paste(genes, collapse = ", ")))

# For pairwise comparisons, we need to use the appropriate groupBy for browser tracks
# Extract the groups being compared for proper visualization
if (analysis_type == "pairwise") {
    # Use the groupBy variable for browser tracks (e.g., "PIMO_Region")
    browser_groupBy <- groupBy
} else {
    browser_groupBy <- groupBy
}

p <- plotBrowserTrack(
    ArchRProj = proj_hyp, 
    groupBy = browser_groupBy, 
    geneSymbol = genes,
    features = getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", returnGR = TRUE),
    upstream = 50000,
    downstream = 50000
)

plotPDF(plotList = p, 
        name = paste0("markerPeaks_browserTrack_", output_suffix), 
        ArchRProj = proj_hyp, 
        addDOC = TRUE, 
        width = 5, 
        height = 5)

print(paste("Analysis complete for", output_suffix))

EOF

echo "Marker heatmap analysis completed"