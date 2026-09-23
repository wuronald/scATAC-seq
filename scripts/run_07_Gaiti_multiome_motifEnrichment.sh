#!/bin/bash
#SBATCH --job-name=run_07_Gaiti_multiome_motifEnrichment
#SBATCH --output=run_07_Gaiti_multiome_motifEnrichment_%j.out
#SBATCH --error=run_07_Gaiti_multiome_motifEnrichment_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=58G
#SBATCH --time=00:30:00

## ==============================================================================
# Purpose: SLURM batch submission script to perform transcription factor (TF)
#          motif enrichment analysis on pairwise differential accessible peaks
#          for the Gaiti multiome ArchR project.
# Workflow:
#   1. Loads peak-annotated ArchR project ('Gaiti_multiome_harmony_merged_malig_peak/').
#   2. Loads and attaches the corresponding PeakSet ('PeakSet_gr_<groupBy>.rds').
#   3. Loads pre-computed pairwise marker features ('markersTest_<groupBy>_<comparison>.rds').
#   4. Adds motif annotations for the specified motif set (e.g. homer, cisbp, encode, JASPAR).
#   5. Computes hypergeometric motif enrichment (peakAnnoEnrichment) for upregulated
#      (FDR <= 0.1 & Log2FC >= 0.5) and downregulated (FDR <= 0.1 & Log2FC <= -0.5) peaks.
#   6. Exports motif enrichment SummarizedExperiment objects (.rds).
#   7. Generates rank-ordered TF enrichment ggplot PDFs with ggrepel labels and
#      heatmaps (plotEnrichHeatmap).
#
# Usage:
#   sbatch scripts/run_07_Gaiti_multiome_motifEnrichment.sh \
#       [MOTIF_SET] [GROUP_BY] [COMPARISON]
#
# Arguments (all positional, optional):
#   $1  MOTIF_SET   Motif database (homer, cisbp, encode, JASPAR2020, JASPAR2018, JASPAR2016)
#                   (default: homer)
#   $2  GROUP_BY    Metadata column used for peak calling (default: PIMO_Region)
#   $3  COMPARISON  Pairwise comparison name (default: auto-detected based on GROUP_BY)
#
# Examples:
#   # PIMO_Region pairwise tests:
#   sbatch scripts/run_07_Gaiti_multiome_motifEnrichment.sh homer PIMO_Region PIMOup_EB_vs_PIMOup_TC
#   sbatch scripts/run_07_Gaiti_multiome_motifEnrichment.sh homer PIMO_Region PIMOup_EB_vs_PIMOdown_EB
#   sbatch scripts/run_07_Gaiti_multiome_motifEnrichment.sh cisbp PIMO_Region PIMOup_EB_vs_PIMOdown_EB
#
#   # Region.annotation pairwise tests:
#   sbatch scripts/run_07_Gaiti_multiome_motifEnrichment.sh homer Region.annotation EB_vs_TC
#   sbatch scripts/run_07_Gaiti_multiome_motifEnrichment.sh homer Region.annotation EB_vs_PT
#   sbatch scripts/run_07_Gaiti_multiome_motifEnrichment.sh homer Region.annotation TC_vs_PT
#
#   # PIMO_up_status pairwise test:
#   sbatch scripts/run_07_Gaiti_multiome_motifEnrichment.sh homer PIMO_up_status PIMOup_vs_PIMOdown
# ==============================================================================

MOTIF_SET="${1:-homer}"  # Default to "homer" if no argument provided
GROUP_BY="${2:-PIMO_Region}"  # Default to "PIMO_Region" if no argument provided

# Set default comparison based on GROUP_BY if not provided
if [ -z "$3" ]; then
    if [ "$GROUP_BY" == "PIMO_up_status" ]; then
        COMPARISON="PIMOup_vs_PIMOdown"
    elif [ "$GROUP_BY" == "Region.annotation" ]; then
        COMPARISON="EB_vs_TC"
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
addArchRThreads(threads = as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "8")))

# set genome to hg38
addArchRGenome("hg38")

# Load the project
archrPath <- "Gaiti_multiome_harmony_merged_malig_peak"
proj <- loadArchRProject(path = archrPath)

# Create file naming prefix based on groupBy variable and comparison
# Replace underscores and special chars to create clean file names
filePrefix <- paste0(gsub("_|\\.", "", groupBy), "_", gsub("_", "", comparison))

# Load and assign appropriate PeakSet based on groupBy
localPeakSetPath <- file.path(archrPath, "PeakCalls", paste0("PeakSet_gr_", groupBy, ".rds"))
clusterPeakSetPath <- file.path("/cluster/projects/wouterslab/ArchR103_4", archrPath, "PeakCalls", paste0("PeakSet_gr_", groupBy, ".rds"))

if (file.exists(localPeakSetPath)) {
    peakSetPath <- localPeakSetPath
} else if (file.exists(clusterPeakSetPath)) {
    peakSetPath <- clusterPeakSetPath
} else {
    stop(paste("PeakSet file not found in local path:", localPeakSetPath, "or cluster path:", clusterPeakSetPath))
}

print(paste("Loading PeakSet from:", peakSetPath))
peakSet_gr <- readRDS(peakSetPath)
print(paste("PeakSet loaded successfully for", groupBy))

# Add the PeakSet to the ArchR project
proj <- addPeakSet(
    ArchRProj = proj,
    peakSet = peakSet_gr,
    force = TRUE
)
print("PeakSet added to ArchR project")

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
    } else if (groupBy == "Region.annotation") {
        stop(paste("markerTest file not found:", markersTest_rds, 
                   "\nFor Region.annotation, available comparisons should be:",
                   "\n  - EB_vs_TC",
                   "\n  - EB_vs_PT",
                   "\n  - TC_vs_PT"))
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

# ── Motif Enrichment in differentially Accessible Peaks ───────────────────────

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

# Motif Enrichment heatmaps with tryCatch
print("Plotting heatmap of motif enrichment results")
heatmapEM_up <- tryCatch(
    plotEnrichHeatmap(motifsUp, n = 10, transpose = TRUE),
    error = function(e) {
        message("WARNING: plotEnrichHeatmap for motifsUp failed — skipping. Reason: ", conditionMessage(e))
        NULL
    }
)
heatmapEM_down <- tryCatch(
    plotEnrichHeatmap(motifsDown, n = 10, transpose = TRUE),
    error = function(e) {
        message("WARNING: plotEnrichHeatmap for motifsDown failed — skipping. Reason: ", conditionMessage(e))
        NULL
    }
)

heatmapList <- Filter(Negate(is.null), list(heatmapEM_up, heatmapEM_down))
if (length(heatmapList) > 0) {
    do.call(plotPDF, c(heatmapList, list(
        name = paste0(groupBy, "_", comparison, "-MotifsUp_MotifsDown-Enriched-Marker-Heatmap_", motifSet),
        width = 8,
        height = 6,
        ArchRProj = proj,
        addDOC = TRUE
    )))
} else {
    message("WARNING: Both heatmaps failed — no heatmap PDF will be written.")
}

print("Motif enrichment analysis completed successfully!")

EOF