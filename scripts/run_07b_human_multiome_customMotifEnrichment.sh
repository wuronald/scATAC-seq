#!/bin/bash
#SBATCH --job-name=run_07b_human_multiome_customMotifEnrichment
#SBATCH --output=run_07b_human_multiome_customMotifEnrichment_%j.out
#SBATCH --error=run_07b_human_multiome_customMotifEnrichment_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=16G
#SBATCH --time=00:30:00

# ==============================================================================
# Purpose: SLURM batch submission script to perform custom peak annotation
#          enrichment analysis comparing single-cell marker peaks in the Human
#          Multiome ArchR project against custom Differentially Accessible (DA)
#          peaks from bulk ATAC-seq (G361 and G411 GSC lines under hypoxia/anoxia).
# Workflow:
#   1. Loads peak-annotated ArchR project (default: 'human_multiome_harmony_merged_malig_peak/').
#   2. Loads and attaches the corresponding PeakSet ('PeakSet_gr_<groupBy>.rds').
#   3. Loads pre-computed marker features (SE object) for the specified groupBy.
#   4. Registers custom bulk ATAC-seq DiffBind BED peak sets (addPeakAnnotations).
#   5. Computes hypergeometric peak annotation enrichment (peakAnnoEnrichment)
#      for upregulated (FDR <= 0.1 & Log2FC >= 0.5) and downregulated
#      (FDR <= 0.1 & Log2FC <= -0.5) peaks.
#   6. Exports custom enrichment SummarizedExperiment objects (.rds).
#   7. Generates heatmap PDFs (plotEnrichHeatmap).
#
# Usage:
#   sbatch scripts/run_07b_human_multiome_customMotifEnrichment.sh \
#       [GROUP_BY] [ARCHR_PATH]
#
# Arguments (all positional, optional):
#   $1  GROUP_BY    Metadata column used for grouping/markers (default: PIMO_up_status)
#   $2  ARCHR_PATH  Path to ArchR project directory
#                   (default: human_multiome_harmony_merged_malig_peak)
#
# Examples:
#   # Standard run with PIMO_up_status:
#   sbatch scripts/run_07b_human_multiome_customMotifEnrichment.sh PIMO_up_status
#
#   # Run with hybrid_pair grouping:
#   sbatch scripts/run_07b_human_multiome_customMotifEnrichment.sh hybrid_pair
# ==============================================================================

# Parse command line arguments
GROUP_BY="${1:-PIMO_up_status}"
ARCHR_PATH="${2:-human_multiome_harmony_merged_malig_peak}"

# Export the parameters so R can access them
export GROUP_BY
export ARCHR_PATH

echo "Running custom peak annotation enrichment analysis with:"
echo "  groupBy:   ${GROUP_BY}"
echo "  archrPath: ${ARCHR_PATH}"

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
groupBy   <- Sys.getenv("GROUP_BY",   unset = "PIMO_up_status")
archrPath <- Sys.getenv("ARCHR_PATH", unset = "human_multiome_harmony_merged_malig_peak")

cat("Using groupBy:  ", groupBy,   "\n")
cat("Using archrPath:", archrPath, "\n")

# Set the number of threads for ArchR
addArchRThreads(threads = 18)

# set genome to hg38
addArchRGenome("hg38")

# Load the project
proj <- loadArchRProject(path = archrPath)

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

# Load marker features (SE object) obtained previously by getMarkerFeatures(): 
localMarkersPath <- file.path(archrPath, "PeakCalls", paste0("markersPeaks_", groupBy, ".rds"))
clusterMarkersPath <- file.path("/cluster/projects/wouterslab/ArchR103_4", archrPath, "PeakCalls", paste0("markersPeaks_", groupBy, ".rds"))

if (file.exists(localMarkersPath)) {
    markersPeaks_rds <- localMarkersPath
} else if (file.exists(clusterMarkersPath)) {
    markersPeaks_rds <- clusterMarkersPath
} else {
    markersPeaks_rds <- localMarkersPath
}

if (!exists("markersPeaks")) {
    print("loading markersPeaks")
    if (file.exists(markersPeaks_rds)) {
        print(paste("previously saved markersPeaks loaded:", markersPeaks_rds))
        markersPeaks <- readRDS(file = markersPeaks_rds)
    } else {
        print("extracting markersPeaks")
        markersPeaks <- getMarkerFeatures(
            ArchRProj = proj,
            useMatrix = "PeakMatrix",
            groupBy = groupBy,
            bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
            testMethod = "wilcoxon"
        )
        saveRDS(markersPeaks, file = markersPeaks_rds) # save for future use
    }
}

# Load Custom Peaks (bed files)
print("Loading custom peaks from bed files")

# Resolve custom BED files directory (check local first, then cluster)
bed_dir <- "data/GSC_diffbind_bed"
if (!dir.exists(bed_dir)) {
    cluster_bed_dir <- "/cluster/projects/wouterslab/Data/GSC_diffbind_bed"
    if (dir.exists(cluster_bed_dir)) {
        bed_dir <- cluster_bed_dir
    }
}

customPeaks <- c(
    "G361_N_72_vs_A_24_hypoxia_enrich"  = file.path(bed_dir, "2023-01-26_G361_N_72_vs_A_24_hypoxia_enrich_hg38.bed"),
    "G361_N_72_vs_A_24_normoxia_enrich" = file.path(bed_dir, "2023-01-26_G361_N_72_vs_A_24_normoxia_enrich_hg38.bed"),
    "G361_N_72_vs_H_24_hypoxia_enrich"  = file.path(bed_dir, "2023-01-26_G361_N_72_vs_H_24_hypoxia_enrich_hg38.bed"),
    "G361_N_72_vs_H_24_normoxia_enrich" = file.path(bed_dir, "2023-01-26_G361_N_72_vs_H_24_normoxia_enrich_hg38.bed"),
    "G361_N_72_vs_H_72_hypoxia_enrich"  = file.path(bed_dir, "2023-01-26_G361_N_72_vs_H_72_hypoxia_enrich_hg38.bed"),
    "G361_N_72_vs_H_72_normoxia_enrich" = file.path(bed_dir, "2023-01-26_G361_N_72_vs_H_72_normoxia_enrich_hg38.bed"),
    "G411_N_72_vs_A_24_hypoxia_enrich"  = file.path(bed_dir, "2023-01-26_G411_N_72_vs_A_24_hypoxia_enrich_hg38.bed"),
    "G411_N_72_vs_A_24_normoxia_enrich" = file.path(bed_dir, "2023-01-26_G411_N_72_vs_A_24_normoxia_enrich_hg38.bed"),
    "G411_N_72_vs_H_24_hypoxia_enrich"  = file.path(bed_dir, "2023-01-26_G411_N_72_vs_H_24_hypoxia_enrich_hg38.bed"),
    "G411_N_72_vs_H_24_normoxia_enrich" = file.path(bed_dir, "2023-01-26_G411_N_72_vs_H_24_normoxia_enrich_hg38.bed"),
    "G411_N_72_vs_H_72_hypoxia_enrich"  = file.path(bed_dir, "2023-01-26_G411_N_72_vs_H_72_hypoxia_enrich_hg38.bed"),
    "G411_N_72_vs_H_72_normoxia_enrich" = file.path(bed_dir, "2023-01-26_G411_N_72_vs_H_72_normoxia_enrich_hg38.bed")
)

missing_beds <- customPeaks[!file.exists(customPeaks)]
if (length(missing_beds) > 0) {
    warning(paste("The following BED files were not found:\n",
                  paste(missing_beds, collapse = "\n")))
}

# Add custom peak annotations to the ArchR project
print("Adding custom peak annotations to ArchR project")
proj <- addPeakAnnotations(
    ArchRProj = proj,
    regions = customPeaks,
    name = "bulk_ATAC_DA_peaks",
    force = TRUE
)

# perform peak annotation enrichment for custom peaks
print("Performing peak annotation enrichment for custom peaks")
enrichRegions <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj,
    peakAnnotation = "bulk_ATAC_DA_peaks",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
)
enrichRegions

depletedRegions <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj,
    peakAnnotation = "bulk_ATAC_DA_peaks",
    cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
)
depletedRegions

# Export custom enrichment results
print("Exporting custom enrichment results")
motifSet <- "bulk_ATAC_DA_peaks"
outDir <- file.path(archrPath, paste0("customMotifEnrichment_", motifSet))
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
saveRDS(enrichRegions, file = file.path(outDir, paste0("enrichRegions_", motifSet, "_vs_", groupBy, ".rds")))
saveRDS(depletedRegions, file = file.path(outDir, paste0("depletedRegions_",  motifSet, "_vs_", groupBy, ".rds")))

# plot heatmaps for enriched and depleted regions
print("Plotting heatmaps for enriched and depleted regions")
enrichHeatmapRegions <- tryCatch({
    plotEnrichHeatmap(enrichRegions, n = 7, transpose = TRUE)
}, error = function(e) {
    cat("Warning: plotEnrichHeatmap failed for enriched regions:", conditionMessage(e), "\n")
    NULL
})

depletedHeatmapRegions <- tryCatch({
    plotEnrichHeatmap(depletedRegions, n = 7, transpose = TRUE)
}, error = function(e) {
    cat("Warning: plotEnrichHeatmap failed for depleted regions:", conditionMessage(e), "\n")
    NULL
})

# save heatmaps to pdf
plot_list <- list(enrichHeatmapRegions, depletedHeatmapRegions)
plot_list <- plot_list[!sapply(plot_list, is.null)]

if (length(plot_list) > 0) {
    print("Saving heatmaps to PDF")
    plotPDF(
        plotList = plot_list,
        name = paste0("custom-Motif-Enriched-Marker-Heatmap-", motifSet, "_vs_", groupBy),
        width = 8,
        height = 6,
        ArchRProj = proj,
        addDOC = TRUE
    )
} else {
    cat("No heatmaps generated to save.\n")
}

EOF