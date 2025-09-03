#!/bin/bash
#SBATCH --job-name=run_05_Terekhanova_motifEnrichment
#SBATCH --output=run_05_Terekhanova_motifEnrichment_%j.out
#SBATCH --error=run_05_Terekhanova_motifEnrichment_%j.err
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

# check peakAnnotation
print("checking peak annotation slot")
print(proj@peakAnnotation)

# Add Motif Annotations from homer
print("add motif annotations")
proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "homer", name = "Motif")

# check peakAnnotation
print("checking peak annotation slot after adding motifs")
print(proj@peakAnnotation)

# Perform pairwise testing between group (hypoxic vs non-hypoxic)
proj$DAEG_UP_24 <- as.logical(proj$DAEG_UP_24)

print("pairwise testing b/w AUCell assigned group: DAEG_UP_24")
markerTest <- getMarkerFeatures(
  ArchRProj = proj,
  useMatrix = "PeakMatrix",
  groupBy = "DAEG_UP_24",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  maxCells = 5000,
  useGroups = "TRUE",
  bgdGroups = "FALSE"
)

# Compute enrichment
print("Compute upregulated motifs")
motifsUp <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = proj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )
  
print("Compute downregulated motifs")  

motifsDown <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = proj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
  )
 
# Examine Motif results
motifsUp
df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
head(df)

motifsDown
df <- data.frame(TF = rownames(motifsDown), mlog10Padj = assay(motifsDown)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
head(df)

# create a folder for motif results
dir.create("Terekhanova_hypoxia/motifs/")

# Save SE object of motifsUP and MotifsDown
save(motifsUp, file = "Terekhanova_hypoxia/motifs/motifsUp.RData")
save(motifsDown, file = "Terekhanova_hypoxia/motifs/motifsDown.RData")
 
EOF

echo "Motif Enrichment analysis completed"