#!/bin/bash
#SBATCH --job-name=run_06_mouse_multiome_markerHeatmap
#SBATCH --output=run_06_mouse_multiome_markerHeatmap_%j.out
#SBATCH --error=run_06_mouse_multiome_markerHeatmap_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=60G
#SBATCH --time=04:00:00

# =============================================================================
# Parse command line arguments
# Usage:
#   sbatch scripts/run_06_mouse_multiome_markerHeatmap.sh \
#       [ARCHR_PATH] [GROUP_BY] [GENES] [PAIRWISE]
#
# Arguments (all optional, positional):
#   $1  ARCHR_PATH  Path to the ArchRProject directory
#                   Default: mouse_multiome_harmony_merged_malig_peak_subset
#
#   $2  GROUP_BY    Metadata column to group cells by for marker analysis
#                   Default: hybrid_pair
#
#   $3  GENES       Comma-separated gene symbols for browser tracks
#                   Default: Olig1,Olig2,... (see below)
#
#   $4  PAIRWISE    Whether to load a pairwise markersTest RDS (true/false)
#                   Default: true when GROUP_BY=PIMO_up_status, false otherwise
#                   true  → loads markersTest_{GROUP_BY}_PIMOup_vs_PIMOdown.rds
#                   false → loads markersPeaks_{GROUP_BY}.rds
#
# Examples:
#   # Run with all defaults (hybrid_pair, pairwise=false):
#   sbatch scripts/run_06_mouse_multiome_markerHeatmap.sh
#
#   # Custom groupBy, pairwise auto-defaults to false:
#   sbatch scripts/run_06_mouse_multiome_markerHeatmap.sh \
#       "mouse_multiome_harmony_merged_malig_peak_subset" "neftel_4_state"
#
#   # PIMO_up_status groupBy (pairwise auto-defaults to true):
#   sbatch scripts/run_06_mouse_multiome_markerHeatmap.sh \
#       "mouse_multiome_harmony_merged_malig_peak_subset" "PIMO_up_status"
#
#   # Force pairwise off for PIMO_up_status:
#   sbatch scripts/run_06_mouse_multiome_markerHeatmap.sh \
#       "mouse_multiome_harmony_merged_malig_peak_subset" "PIMO_up_status" "" "false"
#
#   # All arguments specified:
#   sbatch scripts/run_06_mouse_multiome_markerHeatmap.sh \
#       "mouse_multiome_harmony_merged_malig_peak_subset" "hybrid_pair" \
#       "Olig1,Olig2,Sox2,Sox10,Id2,Id3,Vegfa,Ctsb,Ngfr,Ca9" "false"
# =============================================================================

DEFAULT_GENES="Olig1,Olig2,Sox2,Sox10,Id1,Id2,Id3,Vegfa,Ctsb,Ngfr,Ca9,Slc2a1,\
Hes1,Nr4a2,Nr4a3,Klf9,Nrn1,Dner,Dpysl4,Rnd3,Cadm3,Plppr3,Slc2a3,Pcp4,Sox12,\
Sv2a,Vamp1,Ptprn,Homer2,Nmb,Adm,Eno2,Hmox1,Cryab,Tiparp,Sod2,Xbp1,\
Atf3,Ddit4"

ARCHR_PATH="${1:-mouse_multiome_harmony_merged_malig_peak_subset}"
GROUP_BY="${2:-hybrid_pair}"
GENES="${3:-${DEFAULT_GENES}}"

# PAIRWISE: default true only when GROUP_BY is PIMO_up_status
if [ -n "$4" ]; then
    PAIRWISE="$4"
elif [ "${GROUP_BY}" = "PIMO_up_status" ]; then
    PAIRWISE="true"
else
    PAIRWISE="false"
fi

# Validate that the ArchRProject directory exists before launching R
if [ ! -d "${ARCHR_PATH}" ]; then
    echo "ERROR: ArchRProject directory not found: ${ARCHR_PATH}"
    exit 1
fi

# Validate PAIRWISE value
if [ "${PAIRWISE}" != "true" ] && [ "${PAIRWISE}" != "false" ]; then
    echo "ERROR: PAIRWISE must be 'true' or 'false', got: ${PAIRWISE}"
    exit 1
fi

# Export parameters so the R heredoc can read them via Sys.getenv()
export ARCHR_PATH
export GROUP_BY
export GENES
export PAIRWISE

echo "============================================="
echo "Running markerHeatmap analysis with:"
echo "  ArchRProject path : ${ARCHR_PATH}"
echo "  Group-by variable : ${GROUP_BY}"
echo "  Pairwise mode     : ${PAIRWISE}"
echo "  Genes             : ${GENES}"
echo "============================================="

module load R/4.4.1

Rscript - <<'EOF'

library(ArchR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(here)
set.seed(1)

# ── Read parameters from environment ──────────────────────────────────────────
archr_path <- Sys.getenv("ARCHR_PATH",
    unset = "mouse_multiome_harmony_merged_malig_peak_subset")
group_by   <- Sys.getenv("GROUP_BY",
    unset = "hybrid_pair")
pairwise   <- tolower(Sys.getenv("PAIRWISE", unset = "false")) == "true"
genes_raw  <- Sys.getenv("GENES",
    unset = "Olig1,Olig2,Sox2,Sox10,Id2,Id3,Vegfa")

genes <- unlist(strsplit(genes_raw, ","))

# ── Derived file paths ────────────────────────────────────────────────────────
peak_calls_dir <- file.path(archr_path, "PeakCalls")

# markersPeaks RDS: pairwise test file vs. standard markers file
if (pairwise) {
    markersPeaks_rds <- file.path(peak_calls_dir,
        paste0("markersTest_", group_by, "_PIMOup_vs_PIMOdown.rds"))
} else {
    markersPeaks_rds <- file.path(peak_calls_dir,
        paste0("markersPeaks_", group_by, ".rds"))
}

markersPeaks_GR_rds <- file.path(peak_calls_dir,
    paste0("markersPeaks_GR_", group_by, ".rds"))

cat("Resolved parameters:\n")
cat("  ArchRProject    :", archr_path, "\n")
cat("  PeakCalls dir   :", peak_calls_dir, "\n")
cat("  Pairwise mode   :", pairwise, "\n")
cat("  markersPeaks RDS:", markersPeaks_rds, "\n")
cat("  markersPeaks GR :", markersPeaks_GR_rds, "\n")
cat("  group_by        :", group_by, "\n")
cat("  genes           :", paste(genes, collapse = ", "), "\n\n")

# ── ArchR setup ───────────────────────────────────────────────────────────────
addArchRThreads(threads = 18)
addArchRGenome("mm10")

# ── Load ArchRProject ─────────────────────────────────────────────────────────
if (!dir.exists(archr_path)) stop("ArchRProject directory not found: ", archr_path)
proj_hyp <- loadArchRProject(path = archr_path)

# ── Load or compute markersPeaks ──────────────────────────────────────────────
if (!exists("markersPeaks")) {
    print("Loading markersPeaks ...")
    if (file.exists(markersPeaks_rds)) {
        print(paste("Previously saved markersPeaks found; loading from:", markersPeaks_rds))
        markersPeaks <- readRDS(file = markersPeaks_rds)
    } else {
        print(paste("No saved RDS at:", markersPeaks_rds, "- computing now."))
        print(paste("Running getMarkerFeatures() with groupBy =", group_by))
        markersPeaks <- getMarkerFeatures(
            ArchRProj  = proj_hyp,
            useMatrix  = "PeakMatrix",
            groupBy    = group_by,
            bias       = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
            testMethod = "wilcoxon"
        )
        dir.create(peak_calls_dir, showWarnings = FALSE, recursive = TRUE)
        saveRDS(markersPeaks, file = markersPeaks_rds)
        print(paste("Saved markersPeaks to:", markersPeaks_rds))
    }
}

print("Checking if markersPeaks is a SummarizedExperiment object ...")
if (!is(markersPeaks, "SummarizedExperiment")) {
    stop("markersPeaks is not a SummarizedExperiment object. Please check your data.")
} else {
    print("markersPeaks is a valid SummarizedExperiment object.")
}

# ── Load or compute markersPeaks GR ──────────────────────────────────────────
print(paste("Extracting marker peaks GR object for", group_by, "groups"))
if (file.exists(markersPeaks_GR_rds)) {
    print(paste("Previously saved GR found; loading from:", markersPeaks_GR_rds))
    markersPeaks_GR <- readRDS(file = markersPeaks_GR_rds)
} else {
    print(paste("No saved GR RDS at:", markersPeaks_GR_rds, "- computing now."))
    markersPeaks_GR <- getMarkers(markersPeaks,
        cutOff = "FDR <= 1 & abs(Log2FC) >= 0", returnGR = TRUE)
    saveRDS(markersPeaks_GR, file = markersPeaks_GR_rds)
    print(paste("Saved markersPeaks_GR to:", markersPeaks_GR_rds))
}

print(paste("Number of marker peaks identified:", sapply(markersPeaks_GR, length)))

# ── Marker heatmap ────────────────────────────────────────────────────────────
print(paste("Plotting marker heatmap for", group_by, "groups"))
cutOff <- "FDR <= 0.1 & abs(Log2FC) >= 0.5"
print(paste("Using cutOff:", cutOff))

heatmap <- plotMarkerHeatmap(
    seMarker     = markersPeaks,
    cutOff       = cutOff,
    limits       = c(-2, 2),
    plotLog2FC   = TRUE,
    transpose    = TRUE,
    returnMatrix = FALSE
)
plotPDF(heatmap,
    name      = paste0("markerPeaks-Heatmap_", group_by),
    width     = 8,
    height    = 6,
    ArchRProj = proj_hyp,
    addDOC    = TRUE)

# ── MA & Volcano plots ────────────────────────────────────────────────────────
print("Plotting MA and Volcano plots for marker peaks ...")
for (group in colnames(markersPeaks)) {
    print(paste("Creating MA + Volcano plots for group:", group))
    pma <- plotMarkers(seMarker = markersPeaks, name = group,
        cutOff = cutOff, plotAs = "MA")
    pvolcano <- plotMarkers(seMarker = markersPeaks, name = group,
        cutOff = cutOff, plotAs = "Volcano")
    plotPDF(pma, pvolcano,
        name      = paste0("markerPeaks-MAplot_", group_by, "_", group),
        width     = 5,
        height    = 5,
        ArchRProj = proj_hyp,
        addDOC    = TRUE)
}
print("All MA and Volcano plots created and saved.")

# ── Browser tracks ────────────────────────────────────────────────────────────
print("Plotting browser tracks for genes of interest ...")
print(paste("Genes:", paste(genes, collapse = ", ")))

# PIMO_up_status-specific palette; NULL (ArchR default) for any other groupBy
PIMO_up_status_colors <- c("PIMOdown" = "blue", "PIMOinter" = "gold", "PIMOup" = "red")
track_pal <- if (group_by == "PIMO_up_status") PIMO_up_status_colors else NULL

p <- plotBrowserTrack(
    ArchRProj  = proj_hyp,
    groupBy    = group_by,
    geneSymbol = genes,
    pal        = track_pal,
    features   = getMarkers(markersPeaks,
        cutOff = "FDR <= 0.1 & abs(Log2FC) >= 1", returnGR = TRUE),
    upstream   = 50000,
    downstream = 50000
)
plotPDF(plotList = p,
    name      = paste0("markerPeaks_browserTrack-", group_by),
    ArchRProj = proj_hyp,
    addDOC    = TRUE,
    width     = 5,
    height    = 5)

print("Done! All outputs saved.")

EOF