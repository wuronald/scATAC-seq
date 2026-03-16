#!/bin/bash
#SBATCH --job-name=run_07_human_multiome_motifEnrichment
#SBATCH --output=run_07_human_multiome_motifEnrichment_%j.out
#SBATCH --error=run_07_human_multiome_motifEnrichment_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=32G
#SBATCH --time=01:30:00

# Parse command line arguments
# First argument:  motif set      (homer, encode, JASPAR2020, JASPAR2018, JASPAR2016, cisbp)
# Second argument: groupBy        (e.g., "PIMO_up_status", "hybrid_pair", "k3_program")
# Third argument:  archr_path     path to ArchR project directory
# Fourth argument: pairwise       "true" (default when groupBy=PIMO_up_status) or "false"
#
# Example usage:
#   sbatch scripts/run_07_human_multiome_motifEnrichment.sh homer PIMO_up_status
#   sbatch scripts/run_07_human_multiome_motifEnrichment.sh homer k3_program human_multiome_hmmp_k3_program false

DEFAULT_ARCHR_PATH="human_multiome_harmony_merged_malig_peak"

MOTIF_SET="${1:-cisbp}"               # Default: cisbp
GROUP_BY="${2:-hybrid_pair}"          # Default: hybrid_pair
ARCHR_PATH="${3:-${DEFAULT_ARCHR_PATH}}"  # Default: human_multiome_harmony_merged_malig_peak

# Pairwise default: true when groupBy is PIMO_up_status, false otherwise
if [ -z "${4}" ]; then
    if [ "${GROUP_BY}" = "PIMO_up_status" ]; then
        PAIRWISE="true"
    else
        PAIRWISE="false"
    fi
else
    PAIRWISE="${4}"
fi

# Export the parameters so R can access them
export MOTIF_SET
export GROUP_BY
export ARCHR_PATH
export PAIRWISE

echo "Running motif enrichment analysis with:"
echo "  motifSet:   ${MOTIF_SET}"
echo "  groupBy:    ${GROUP_BY}"
echo "  archrPath:  ${ARCHR_PATH}"
echo "  pairwise:   ${PAIRWISE}"

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
motifSet  <- Sys.getenv("MOTIF_SET",   unset = "cisbp")
groupBy   <- Sys.getenv("GROUP_BY",    unset = "hybrid_pair")
archrPath <- Sys.getenv("ARCHR_PATH",  unset = "human_multiome_harmony_merged_malig_peak")
pairwise  <- tolower(Sys.getenv("PAIRWISE", unset = "false")) == "true"

cat("Using motifSet:  ", motifSet,  "\n")
cat("Using groupBy:   ", groupBy,   "\n")
cat("Using archrPath: ", archrPath, "\n")
cat("Using pairwise:  ", pairwise,  "\n")

# Set the number of threads for ArchR
addArchRThreads(threads = 18)

# set genome to hg38
addArchRGenome("hg38")

# Load the project
proj <- loadArchRProject(path = archrPath)

# Create file naming prefix based on groupBy variable
# Replace underscores and special chars to create clean file names
filePrefix <- gsub("_", "", groupBy)

# Load and assign appropriate PeakSet based on groupBy
peakSetPath <- paste0("/cluster/projects/wouterslab/ArchR103_4/", archrPath, "/PeakCalls/PeakSet_gr_", groupBy, ".rds")
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

# ── Load marker features ───────────────────────────────────────────────────────
# Pairwise mode  → markersTest_{groupBy}PIMOup_vs_PIMOdown.rds
# Non-pairwise   → markersPeaks_{groupBy}.rds

if (pairwise) {
    markers_rds <- paste0(archrPath, "/PeakCalls/markersTest_", groupBy, "PIMOup_vs_PIMOdown", ".rds")
    print(paste("Pairwise mode: loading markerTest from", markers_rds))

    if (!exists("markerTest")) {
        if (file.exists(markers_rds)) {
            print("Previously saved markerTest loaded.")
            markerTest <- readRDS(file = markers_rds)
        } else {
            print("RDS not found — extracting pairwise markerTest for PIMO_up_status PIMOup vs PIMOdown")
            markerTest <- getMarkerFeatures(
                ArchRProj  = proj,
                useMatrix  = "PeakMatrix",
                groupBy    = "PIMO_up_status",
                testMethod = "wilcoxon",
                bias       = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
                useGroups  = "PIMOup",
                bgdGroups  = "PIMOdown"
            )
        }
    }
} else {
    markers_rds <- paste0(archrPath, "/PeakCalls/markersPeaks_", groupBy, ".rds")
    print(paste("Non-pairwise mode: loading markersPeaks from", markers_rds))

    if (!exists("markerTest")) {
        if (file.exists(markers_rds)) {
            print("Previously saved markersPeaks loaded.")
            markerTest <- readRDS(file = markers_rds)
        } else {
            print("RDS not found — extracting non-pairwise markersPeaks")
            markerTest <- getMarkerFeatures(
                ArchRProj  = proj,
                useMatrix  = "PeakMatrix",
                groupBy    = groupBy,
                bias       = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
                testMethod = "wilcoxon"
            )
        }
    }
}

# ── Add motif annotations ──────────────────────────────────────────────────────
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

# ── Motif Enrichment in differentially Accessible Peaks ───────────────────────

motifsUp <- peakAnnoEnrichment(
    seMarker      = markerTest,
    ArchRProj     = proj,
    peakAnnotation = annoName,
    cutOff        = "FDR <= 0.1 & Log2FC >= 0.5"
)

motifsDown <- peakAnnoEnrichment(
    seMarker      = markerTest,
    ArchRProj     = proj,
    peakAnnotation = annoName,
    cutOff        = "FDR <= 0.1 & Log2FC <= -0.5"
)

# Export SummarizedExperiment Motif enrichment results
print("Exporting motif enrichment results")
outDir <- here(paste0(archrPath, "/motifEnrichment_", motifSet, "_", groupBy, "/"))
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
saveRDS(motifsUp,   file = file.path(outDir, paste0("motifsUp_",   filePrefix, "_", motifSet, ".rds")))
saveRDS(motifsDown, file = file.path(outDir, paste0("motifsDown_", filePrefix, "_", motifSet, ".rds")))

# Function: motif enrichment SE objects for plotting
print("Plotting ggplot of motif enrichment results")

# plotMotifEnrichments <- function(motifs) {
#     # Convert to data frame for ggplot
#     df <- data.frame(TF = rownames(motifs), mlog10Padj = assay(motifs)[,1])
#     df <- df[order(df$mlog10Padj, decreasing = TRUE),]
#     df$rank <- seq_len(nrow(df))

#     gg <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) +
#       geom_point(size = 1) +
#       ggrepel::geom_label_repel(
#             data    = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF),
#             size    = 1.5,
#             nudge_x = 150,
#             force   = 5,
#             direction = "both",
#             hjust   = 0.5,
#             color   = "black"
#       ) + theme_ArchR() +
#       ylab("-log10(P-adj) Motif Enrichment") +
#       xlab("Rank Sorted TFs Enriched") +
#       scale_colour_gradientn(colors = paletteContinuous(set = "comet"))
#     return(gg)
# }

# ggUp   <- plotMotifEnrichments(motifsUp)
# ggDown <- plotMotifEnrichments(motifsDown)

# plotPDF(ggUp, ggDown, name = paste0(filePrefix, "-Markers-Motifs-Enriched_", motifSet), width = 5, height = 5, ArchRProj = proj, addDOC = TRUE)
# Function: motif enrichment SE objects for plotting
print("Plotting ggplot of motif enrichment results")

plotMotifEnrichments <- function(motifs, groupName = NULL) {
    df <- data.frame(TF = rownames(motifs), mlog10Padj = assay(motifs)[, groupName])
    df <- df[order(df$mlog10Padj, decreasing = TRUE),]
    df$rank <- seq_len(nrow(df))

    title <- if (!is.null(groupName)) groupName else "Enrichment"

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
      ggtitle(title) +
      ylab("-log10(P-adj) Motif Enrichment") +
      xlab("Rank Sorted TFs Enriched") +
      scale_colour_gradientn(colors = paletteContinuous(set = "comet"))
    return(gg)
}

# Generate one plot per group for motifsUp and motifsDown
groups <- colnames(assay(motifsUp))
ggUp   <- lapply(groups, function(g) plotMotifEnrichments(motifsUp,   groupName = g))
ggDown <- lapply(groups, function(g) plotMotifEnrichments(motifsDown, groupName = g))

# Flatten into a single list: [up_prog1, up_prog2, ..., down_prog1, down_prog2, ...]
allPlots <- c(ggUp, ggDown)

do.call(plotPDF, c(allPlots, list(
    name      = paste0(filePrefix, "-Markers-Motifs-Enriched_", motifSet),
    width     = 5,
    height    = 5,
    ArchRProj = proj,
    addDOC    = TRUE
)))

# Motif Enrichment heatmaps
print("Plotting heatmap of motif enrichment results")
heatmapEM_up <- tryCatch(
    plotEnrichHeatmap(motifsUp,   n = 10, transpose = TRUE),
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

# Motif Enrichment heatmaps (relaxed thresholds)
heatmapEM_up_relaxed <- tryCatch(
    plotEnrichHeatmap(motifsUp,   n = 10, transpose = TRUE, cutOff = 5),
    error = function(e) {
        message("WARNING: plotEnrichHeatmap for motifsUp (relaxed) failed — skipping. Reason: ", conditionMessage(e))
        NULL
    }
)
heatmapEM_down_relaxed <- tryCatch(
    plotEnrichHeatmap(motifsDown, n = 10, transpose = TRUE, cutOff = 5),
    error = function(e) {
        message("WARNING: plotEnrichHeatmap for motifsDown (relaxed) failed — skipping. Reason: ", conditionMessage(e))
        NULL
    }
)

# Only plot whichever heatmaps succeeded
heatmapList <- Filter(Negate(is.null), list(heatmapEM_up, heatmapEM_down, heatmapEM_up_relaxed, heatmapEM_down_relaxed))
if (length(heatmapList) > 0) {
    do.call(plotPDF, c(heatmapList, list(
        name      = paste0(groupBy, "-_MotifsUp_MotifsDown-Enriched-Marker-Heatmap_", motifSet),
        width     = 8,
        height    = 6,
        ArchRProj = proj,
        addDOC    = TRUE
    )))
} else {
    message("WARNING: Both heatmaps failed — no heatmap PDF will be written.")
}
EOF