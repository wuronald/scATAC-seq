#!/bin/bash
#SBATCH --job-name=run_07_Gaiti_multiome_motifEnrichment
#SBATCH --output=run_07_Gaiti_multiome_motifEnrichment_%j.out
#SBATCH --error=run_07_Gaiti_multiome_motifEnrichment_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=58G
#SBATCH --time=01:30:00

# Parse command line arguments
# First argument: motif set (homer, encode, JASPAR2020, JASPAR2018, JASPAR2016, cisbp)
# Second argument: groupBy variable (e.g., "PIMO_up_status", "PIMO_Region")
# Third argument: comparison name (e.g., "PIMOup_vs_PIMOdown", "PIMOup_EB_vs_PIMOup_TC", "PIMOup_EB_vs_PIMOdown_EB")
#                 If not provided, defaults based on GROUP_BY:
#                   - PIMO_up_status -> PIMOup_vs_PIMOdown
#                   - PIMO_Region -> PIMOup_EB_vs_PIMOup_TC
# Example usage:
# sbatch scripts/run_07_Gaiti_multiome_motifEnrichment.sh homer PIMO_up_status PIMOup_vs_PIMOdown
# sbatch scripts/run_07_Gaiti_multiome_motifEnrichment.sh homer PIMO_Region PIMOup_EB_vs_PIMOup_TC
# sbatch scripts/run_07_Gaiti_multiome_motifEnrichment.sh cisbp PIMO_Region PIMOup_EB_vs_PIMOdown_EB

MOTIF_SET="${1:-homer}"  # Default to "homer" if no argument provided
GROUP_BY="${2:-PIMO_Region}"  # Default to "PIMO_Region" if no argument provided

# Set default comparison based on GROUP_BY if not provided
if [ -z "$3" ]; then
    if [ "$GROUP_BY" == "PIMO_up_status" ]; then
        COMPARISON="PIMOup_vs_PIMOdown"
    else
        COMPARISON="PIMOup_EB_vs_PIMOup_TC"
    fi
else
    COMPARISON="$3"
fi

# Export the parameters so R can access them
export MOTIF_SET
export GROUP_BY
export COMPARISON

echo "Running motif enrichment analysis with:"
echo "  motifSet: ${MOTIF_SET}"
echo "  groupBy: ${GROUP_BY}"
echo "  comparison: ${COMPARISON}"

# Load necessary modules (adjust as needed for your system)
module load R/4.4.1

# Run R script with parameters
Rscript - <<'EOF'

# load libraries
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(here)
set.seed(1)

# Get parameters from environment variables
motifSet <- Sys.getenv("MOTIF_SET", unset = "homer")
groupBy <- Sys.getenv("GROUP_BY", unset = "PIMO_Region")
comparison <- Sys.getenv("COMPARISON", unset = "PIMOup_EB_vs_PIMOup_TC")
cat("Using motifSet:", motifSet, "\n")
cat("Using groupBy:", groupBy, "\n")
cat("Using comparison:", comparison, "\n")

# Set the number of threads for ArchR
addArchRThreads(threads = 8)

# set genome to hg38
addArchRGenome("hg38")

# Load the project
proj <- loadArchRProject(path = "Gaiti_multiome_harmony_merged_malig_peak")

# Create file naming prefix based on groupBy variable and comparison
# Replace underscores and special chars to create clean file names
filePrefix <- paste0(gsub("_", "", groupBy), "_", gsub("_", "", comparison))

# Load and assign appropriate PeakSet based on groupBy
peakSetPath <- paste0("/cluster/projects/wouterslab/ArchR103_4/Gaiti_multiome_harmony_merged_malig_peak/PeakCalls/PeakSet_gr_", groupBy, ".rds")
print(paste("Loading PeakSet from:", peakSetPath))

if (file.exists(peakSetPath)) {
    peakSet_gr <- readRDS(peakSetPath)
    print(paste("PeakSet loaded successfully for", groupBy))
    
    # Add the PeakSet to the ArchR project
    proj <- addPeakSet(
        ArchRProj = proj,
        peakSet = peakSet_gr,
        force = TRUE
    )
    print("PeakSet added to ArchR project")
} else {
    stop(paste("PeakSet file not found:", peakSetPath))
}

# Load markerTest from pairwise differential peaks getMarkerFeatures()
# Construct path dynamically based on groupBy and comparison
markersTest_rds <- paste0("Gaiti_multiome_harmony_merged_malig_peak/PeakCalls/markersTest_", groupBy, "_", comparison, ".rds")

print(paste("Attempting to load markerTest from:", markersTest_rds))

if (file.exists(markersTest_rds)) {
    print(paste("Loading markerTest for comparison:", comparison))
    markerTest <- readRDS(file = markersTest_rds)
    print("markerTest loaded successfully")
} else {
    # Provide helpful error message based on groupBy
    if (groupBy == "PIMO_up_status") {
        stop(paste("markerTest file not found:", markersTest_rds, 
                   "\nFor PIMO_up_status, expected comparison:",
                   "\n  - PIMOup_vs_PIMOdown"))
    } else if (groupBy == "PIMO_Region") {
        stop(paste("markerTest file not found:", markersTest_rds, 
                   "\nFor PIMO_Region, available comparisons should be:",
                   "\n  - PIMOup_EB_vs_PIMOup_TC",
                   "\n  - PIMOup_EB_vs_PIMOdown_EB", 
                   "\n  - PIMOup_TC_vs_PIMOdown_TC",
                   "\n  - PIMOup_PT_vs_PIMOdown_PT"))
    } else {
        stop(paste("markerTest file not found:", markersTest_rds))
    }
}

# Add motif annotations
print(paste("Adding motif annotations to ArchR project with motifSet:", motifSet))

# Create annotation name based on motifSet
annoName <- paste0("Motif_", motifSet)

print(paste("Checking if motif annotations already exist for", annoName))
if (annoName %in% names(proj@peakAnnotation)) {
    print(paste("Motif annotations", annoName, "already exist in the ArchR project. Skipping addition."))
} else {
    print(paste("Motif annotations do not exist. Proceeding to add them with motifSet:", motifSet))
    proj <- addMotifAnnotations(ArchRProj = proj, motifSet = motifSet, annoName = annoName)
}   

# Determine motifs in a given gene locus
# Example for CEBPA

# get peakset
# note: addPeakSet() can be used to assign the peakset if multiple peaksets exist
# only one peakset can be active at a time

print("Getting peakset from ArchR project")
pSet <- getPeakSet(ArchRProj = proj)
pSet$name <- paste(seqnames(pSet), start(pSet), end(pSet), sep = "_")

# get motif matches for the peakset
print(paste("Getting motif matches for the peakset using", annoName))
matches <- getMatches(ArchRProj = proj, name = annoName)
rownames(matches) <- paste(seqnames(matches), start(matches), end(matches), sep = "_")
matches <- matches[pSet$name]

# GR object corresponding to CEBPA promoter
print("Finding motifs in CEBPA promoter region")
gr <- GRanges(seqnames = c("chr19"), ranges = IRanges(start = c(33792929), end = c(33794030)))
queryHits <- queryHits(findOverlaps(query = pSet, subject = gr, type = "within"))

# Motif Enrichment in differentially Accessible Peaks

print("Performing motif enrichment for upregulated peaks")
motifsUp <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = proj,
    peakAnnotation = annoName,
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

print("Performing motif enrichment for downregulated peaks")
motifsDown <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = proj,
    peakAnnotation = annoName,
    cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
  )

# Export SummarizedExperiment Motif enrichment results
print("Exporting motif enrichment results")
outDir <- here(paste0("Gaiti_multiome_harmony_merged_malig_peak/motifEnrichment_", motifSet, "_", groupBy, "_", comparison, "/"))
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
saveRDS(motifsUp, file = file.path(outDir, paste0("motifsUp_", filePrefix, "_", motifSet, ".rds")))
saveRDS(motifsDown, file = file.path(outDir, paste0("motifsDown_", filePrefix, "_", motifSet, ".rds")))

# Function: motif enrichment SE objects for plotting
print("Plotting ggplot of motif enrichment results")

plotMotifEnrichments <- function(motifs) {
    # Convert to data frame for ggplot
    df <- data.frame(TF = rownames(motifs), mlog10Padj = assay(motifs)[,1])
    df <- df[order(df$mlog10Padj, decreasing = TRUE),]
    df$rank <- seq_len(nrow(df))
    
    gg <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
      geom_point(size = 1) +
      ggrepel::geom_label_repel(
            data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
            size = 1.5,
            nudge_x = 150,
            force = 5,
            direction = "both",
            hjust = 0.5,
            color = "black"
      ) + theme_ArchR() + 
      ylab("-log10(P-adj) Motif Enrichment") + 
      xlab("Rank Sorted TFs Enriched") +
      scale_colour_gradientn(colors = paletteContinuous(set = "comet"))
    return(gg)
}

ggUp <- plotMotifEnrichments(motifsUp)
ggDown <- plotMotifEnrichments(motifsDown)

plotPDF(ggUp, ggDown, name = paste0(filePrefix, "-Markers-Motifs-Enriched_", motifSet), width = 5, height = 5, ArchRProj = proj, addDOC = TRUE)

print("Plotting heatmap of motif enrichment results")
heatmapEM_up <- plotEnrichHeatmap(motifsUp, n = 10, transpose = TRUE)
heatmapEM_down <- plotEnrichHeatmap(motifsDown, n = 10, transpose = TRUE)
plotPDF(heatmapEM_up, heatmapEM_down, name = paste0(groupBy, "_", comparison, "-MotifsUp_MotifsDown-Enriched-Marker-Heatmap_", motifSet), width = 8, height = 6, ArchRProj = proj, addDOC = TRUE)

print("Motif enrichment analysis completed successfully!")

EOF