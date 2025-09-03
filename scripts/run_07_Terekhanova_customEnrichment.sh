#!/bin/bash
#SBATCH --job-name=run_07_Terekhanova_customEnrichment
#SBATCH --output=run_07_Terekhanova_customEnrichment_%j.out
#SBATCH --error=run_07_Terekhanova_customEnrichment_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=60G
#SBATCH --time=01:00:00

# Load necessary modules (adjust as needed for your system)
module load R/4.2.1

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
proj <- loadArchRProject(path = "Terekhanova_hypoxia")

# Load bulk-ATAC-seq DiffBind peaks
# note: we previously converted the bed file from hg19 to hg38
print("Load bed files from previous DiffBind bulk ATAC-seq analysis:")
diffPeaks <- c(
			   G361_N_72_vs_A_24_hypoxia_enrich = "data/bed/2023-01-26_G361_N_72_vs_A_24_hypoxia_enrich_hg38.bed",
			   G361_N_72_vs_A_24_normoxia_enrich = "data/bed/2023-01-26_G361_N_72_vs_A_24_normoxia_enrich_hg38.bed",
			   G361_N_72_vs_H_24_hypoxia_enrich = "data/bed/2023-01-26_G361_N_72_vs_H_24_hypoxia_enrich_hg38.bed",
			   G361_N_72_vs_H_24_normoxia_enrich = "data/bed/2023-01-26_G361_N_72_vs_H_24_normoxia_enrich_hg38.bed",
			   G361_N_72_vs_H_72_hypoxia_enrich = "data/bed/2023-01-26_G361_N_72_vs_H_72_hypoxia_enrich_hg38.bed",
			   G361_N_72_vs_H_72_normoxia_enrich = "data/bed/2023-01-26_G361_N_72_vs_H_72_normoxia_enrich_hg38.bed",
			   G411_N_72_vs_A_24_hypoxia_enrich = "data/bed/2023-01-26_G411_N_72_vs_A_24_hypoxia_enrich_hg38.bed",
			   G411_N_72_vs_A_24_normoxia_enrich = "data/bed/2023-01-26_G411_N_72_vs_A_24_normoxia_enrich_hg38.bed",
			   G411_N_72_vs_H_24_hypoxia_enrich = "data/bed/2023-01-26_G411_N_72_vs_H_24_hypoxia_enrich_hg38.bed",
			   G411_N_72_vs_H_24_normoxia_enrich = "data/bed/2023-01-26_G411_N_72_vs_H_24_normoxia_enrich_hg38.bed",
			   G411_N_72_vs_H_72_hypoxia_enrich = "data/bed/2023-01-26_G411_N_72_vs_H_72_hypoxia_enrich_hg38.bed",
			   G411_N_72_vs_H_72_normoxia_enrich = "data/bed/2023-01-26_G411_N_72_vs_H_72_normoxia_enrich_hg38.bed"
			   )
# subset diffPeaks as G361 beds only   
diffPeaks <- diffPeaks[1]
# Add Custom Annotations from bed files
print("add Custom annotations from DiffBind bed files")
proj <- addPeakAnnotations(ArchRProj = proj, regions = diffPeaks, name = "bulkATAC")
						
# Load Marker Peaks SE object(previously computed and saved)
## eg. markersPeaks_DAEG_UP_24 from run_04_terekhanova_peakCalling.sh
print("load markers peaks")
load(file = "Terekhanova_hypoxia/PeakCalls/markersPeaks_DAEG_UP_24.RData")

# Calculate Enrichment of EncodeTFBS with marker Peak set
print("calculate enrichment of custom peaks via peakAnnoEnrichment")
enrichRegions <- peakAnnoEnrichment(
    seMarker = markersPeaks_DAEG_UP_24,
    ArchRProj = proj,
    peakAnnotation = "bulkATAC",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

# Examine results:
df <- data.frame(TF = rownames(enrichRegions), mlog10Padj = assay(enrichRegions)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
head(df)


# Save SE object results of peakAnnoEnrichment()
print("save results of peakAnnoEnrichment")
save(enrichRegions, file = "Terekhanova_hypoxia/motifs/enrichRegions.RData")
						
EOF

echo "Custom Peak Enrichment analysis completed"