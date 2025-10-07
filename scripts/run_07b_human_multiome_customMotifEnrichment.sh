#!/bin/bash
#SBATCH --job-name=run_07b_human_multiome_CustomMotifEnrichment
#SBATCH --output=run_07b_human_multiome_CustomMotifEnrichment_%j.out
#SBATCH --error=run_07b_human_multiome_CustomMotifEnrichment_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=16G
#SBATCH --time=00:30:00

# Parse command line arguments

# First argument: groupBy variable (e.g., "PIMO_up_status", "hybrid_pair")
# Example usage:
# sbatch scripts/run_07b_human_multiome_motifEnrichment.sh PIMO_up_status

GROUP_BY="${1:-PIMO_up_status}"  # Default to "PIMO_up_status" if no argument provided

# Export the parameters so R can access them
export GROUP_BY

echo "Running custom motif enrichment analysis with:"
echo "  groupBy: ${GROUP_BY}"

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
groupBy <- Sys.getenv("GROUP_BY", unset = "PIMO_up_status")
cat("Using groupBy:", groupBy, "\n")

# Set the number of threads for ArchR
addArchRThreads(threads = 18)

# set genome to hg38
addArchRGenome("hg38")

# Load the project
proj <- loadArchRProject(path = "human_multiome_harmony_merged_malig_peak")

# Load and assign appropriate PeakSet based on groupBy
peakSetPath <- paste0("human_multiome_harmony_merged_malig_peak/PeakCalls/PeakSet_gr_", groupBy, ".rds")
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

# Load marker features (SE object) obtained previously by getMarkerFeatures(): 
# Construct path dynamically based on groupBy
markersPeaks_rds <- paste0("human_multiome_harmony_merged_malig_peak/PeakCalls/markersPeaks_", groupBy, ".rds")

if (!exists("markersPeaks")) {
    print("loading markersPeaks")
    if (file.exists(markersPeaks_rds)) {
        print("previously saved markersPeaks loaded:")
        print(markersPeaks_rds)
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

customPeaks <- c(
    "G361_N_72_vs_A_24_hypoxia_enrich" = "data/GSC_diffbind_bed/2023-01-26_G361_N_72_vs_A_24_hypoxia_enrich_hg38.bed",
    "G361_N_72_vs_A_24_normoxia_enrich" = "data/GSC_diffbind_bed/2023-01-26_G361_N_72_vs_A_24_normoxia_enrich_hg38.bed",
    "G361_N_72_vs_H_24_hypoxia_enrich" = "data/GSC_diffbind_bed/2023-01-26_G361_N_72_vs_H_24_hypoxia_enrich_hg38.bed",
    "G361_N_72_vs_H_24_normoxia_enrich" = "data/GSC_diffbind_bed/2023-01-26_G361_N_72_vs_H_24_normoxia_enrich_hg38.bed",
    "G361_N_72_vs_H_72_hypoxia_enrich" = "data/GSC_diffbind_bed/2023-01-26_G361_N_72_vs_H_72_hypoxia_enrich_hg38.bed",
    "G361_N_72_vs_H_72_normoxia_enrich" = "data/GSC_diffbind_bed/2023-01-26_G361_N_72_vs_H_72_normoxia_enrich_hg38.bed",
    "G411_N_72_vs_A_24_hypoxia_enrich" = "data/GSC_diffbind_bed/2023-01-26_G411_N_72_vs_A_24_hypoxia_enrich_hg38.bed",
    "G411_N_72_vs_A_24_normoxia_enrich" = "data/GSC_diffbind_bed/2023-01-26_G411_N_72_vs_A_24_normoxia_enrich_hg38.bed",
    "G411_N_72_vs_H_24_hypoxia_enrich" = "data/GSC_diffbind_bed/2023-01-26_G411_N_72_vs_H_24_hypoxia_enrich_hg38.bed",
    "G411_N_72_vs_H_24_normoxia_enrich" = "data/GSC_diffbind_bed/2023-01-26_G411_N_72_vs_H_24_normoxia_enrich_hg38.bed",
    "G411_N_72_vs_H_72_hypoxia_enrich" = "data/GSC_diffbind_bed/2023-01-26_G411_N_72_vs_H_72_hypoxia_enrich_hg38.bed",
    "G411_N_72_vs_H_72_normoxia_enrich" = "data/GSC_diffbind_bed/2023-01-26_G411_N_72_vs_H_72_normoxia_enrich_hg38.bed"
)

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
outDir <- here(paste0("human_multiome_harmony_merged_malig_peak/customMotifEnrichment_", motifSet,"/"))
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
saveRDS(enrichRegions, file = file.path(outDir, paste0("enrichRegions_", motifSet,"_vs_", groupBy,".rds")))
saveRDS(depletedRegions, file = file.path(outDir, paste0("depletedRegions_",  motifSet,"_vs_", groupBy, ".rds")))


# plot heatmaps for enriched and depleted regions
print("Plotting heatmaps for enriched and depleted regions")
enrichHeatmapRegions <- plotEnrichHeatmap(enrichRegions, n = 7, transpose = TRUE)
depletedHeatmapRegions <- plotEnrichHeatmap(depletedRegions, n = 7, transpose = TRUE)

# save heatmaps to pdf
print("Saving heatmaps to PDF")
plotPDF(enrichHeatmapRegions, depletedHeatmapRegions,
name = paste0("custom-Motif-Enriched-Marker-Heatmap-", motifSet, "_vs_", groupBy),  width = 8, height = 6, ArchRProj = proj, addDOC = TRUE)

EOF