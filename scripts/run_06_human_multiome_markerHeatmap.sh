#!/bin/bash
#SBATCH --job-name=run_06_human_multiome_markerHeatmap
#SBATCH --output=run_06_human_multiome_markerHeatmap_%j.out
#SBATCH --error=run_06_human_multiome_markerHeatmap_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=16G
#SBATCH --time=00:30:00

# Load necessary modules (adjust as needed for your system)
module load R/4.4.1

# Run R script
Rscript - <<'EOF'

# load libraries
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(here)
set.seed(1)

# Set the number of threads for ArchR
addArchRThreads(threads = 18)

# set genome to hg38
addArchRGenome("hg38")

# Load the project
proj_hyp <- loadArchRProject(path = "human_multiome_harmony_merged_malig_peak")

# path to saved markersPeaks RDS file
#markersPeaks_rds <- "human_multiome_harmony_merged_malig_peak/PeakCalls/markersPeaks_PIMO_up_status.rds"
markersPeaks_rds <- "human_multiome_harmony_merged_malig_peak/PeakCalls/markersTest_PIMO_up_status_PIMOup_vs_PIMOdown.rds"

if (!exists("markersPeaks")) {
    print("loading markersPeaks")
    if (file.exists(markersPeaks_rds)) {
        print("previously saved markersPeaks loaded")
        markersPeaks <- readRDS(file = markersPeaks_rds)
    } else {
        print("extracting markersPeaks")
        markersPeaks <- getMarkerFeatures(
            ArchRProj = proj_hyp,
            useMatrix = "PeakMatrix",
            groupBy = "PIMO_up_status",
            bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
            testMethod = "wilcoxon"
        )
    }
}

# check if markersPeaks is a SummarizedExperiment object; needed for plotting
print("Checking if markersPeaks is a SummarizedExperiment object")
if (!is(markersPeaks, "SummarizedExperiment")) {
    stop("markersPeaks is not a SummarizedExperiment object. Please check your markersPeaks data.")
} else {
    print("markersPeaks is a valid SummarizedExperiment object.")
}

# Extract GR object from SE object for PIMO_up_status groups (used later for MA plots in rstudio)
print("Extracting marker peaks GR object for PIMO_up_status groups")

markersPeaks_GR_rds <- "human_multiome_harmony_merged_malig_peak/PeakCalls/markersPeaks_GR_PIMO_up_status.rds"
if (file.exists(markersPeaks_GR_rds)) {
    print(paste("File already exists:", markersPeaks_GR_rds))
    markersPeaks_GR <- readRDS(file = markersPeaks_GR_rds)
} else {
    markersPeaks_GR <- getMarkers(markersPeaks, cutOff = "FDR <= 1 & abs(Log2FC) >= 0", returnGR = TRUE)
    saveRDS(markersPeaks_GR, file = markersPeaks_GR_rds)
    print(paste("Saved markersPeaks_GR to:", markersPeaks_GR_rds))
}

print(paste("Number of marker peaks identified:", sapply(markersPeaks_GR, length)))


# Plot marker heatmap for Azimuth_class
print("Plotting marker heatmap for PIMO_up_status groups")
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
plotPDF(heatmap, name = "markerPeaks-Heatmap", width = 8, height = 6, ArchRProj = proj_hyp, addDOC = TRUE)

# Plot MA plots for marker peaks
print("Plotting MA plots for marker peaks")

# Iterate over each group in the markersPeaks object and create MA plots
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
    # Save each MA plot as a PDF
    print(paste("Saving MA plot for group:", group))
    plotPDF(pma, pvolcano, 
    name = paste0("markerPeaks-MAplot_", group), width = 5, height = 5, ArchRProj = proj_hyp, addDOC = TRUE)
}
print("All MA and Volcano plots created and saved.")

# get gene annotation for hg38
geneAnno <- getGeneAnnotation(proj_hyp)

# Plot marker peaks in browser tracks for genes of interest:
print("Plotting marker peaks in browser tracks for genes of interest")
genes <- c("CTSB","OLIG1", "OLIG2",
"SOX2", "CD109", "CD44", "RND3", "STMN2", "NGFR", "SOX10", 
"ID2", "ID3", "CA9", "VEGFA", "SLC2A1") 
print(paste("Genes of interest:", paste(genes, collapse = ", ")))

p <- plotBrowserTrack(
    ArchRProj = proj_hyp, 
    groupBy = "PIMO_up_status", 
    geneSymbol = genes,
    features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", returnGR = TRUE),
    upstream = 50000,
    downstream = 50000
)
plotPDF(plotList = p, 
    name = paste0("markerPeaks_browserTrack-PIMO_up_status"), 
    ArchRProj = projHeme5, 
    addDOC = TRUE, 
    width = 5, 
    height = 5
    )

# for (gene in genes) {
# print(paste("Plotting marker peaks in browser tracks for: ", gene))

# p <- plotBrowserTrack(
#     ArchRProj = proj_hyp, 
#     groupBy = "PIMO_up_status", 
#     geneSymbol = gene,
#     features =  getMarkers(markersPeaks, cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", returnGR = TRUE),
#     upstream = 50000,
#     downstream = 50000
# )

# plotPDF(p, name = paste0("markerPeaks_browserTrack-PIMO_up_status_", gene), width = 5, height = 5, ArchRProj = proj_hyp, addDOC = TRUE)
# }

EOF