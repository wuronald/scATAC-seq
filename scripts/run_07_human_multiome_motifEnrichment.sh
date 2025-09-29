#!/bin/bash
#SBATCH --job-name=run_06_human_multiome_motifEnrichment
#SBATCH --output=run_06_human_multiome_motifEnrichment_%j.out
#SBATCH --error=run_06_human_multiome_motifEnrichment_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=16G
#SBATCH --time=01:30:00

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
proj <- loadArchRProject(path = "human_multiome_harmony_merged_malig_peak")

# Load markerTest from getMarkerFeatures(): 
# path to saved markersPeaks RDS file
markersPeaks_rds <- "human_multiome_harmony_merged_malig_peak/PeakCalls/markersTest_PIMO_up_status_PIMOup_vs_PIMOdown.rds"

if (!exists("markerTest")) {
    print("loading markerTest")
    if (file.exists(markersPeaks_rds)) {
        print("previously saved markersPeaks loaded")
        markerTest <- readRDS(file = markersPeaks_rds)
    } else {
        print("extracting markersPeaks")
        markerTest <- getMarkerFeatures(
            ArchRProj = proj,
            useMatrix = "PeakMatrix",
            groupBy = "PIMO_up_status",
            bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
            testMethod = "wilcoxon"
        )
    }
}

# Add motif annotations
print("Adding motif annotations to ArchR project")

print("Checking if motif annotations already exist")
if ("Motif" %in% names(proj@peakAnnotation)) {
    print("Motif annotations already exist in the ArchR project. Skipping addition.")
} else {
    print("Motif annotations do not exist. Proceeding to add them.")
    proj <- addMotifAnnotations(ArchRProj = proj, motifSet = "cisbp", annoName = "Motif")
}   

# Determine motifs in a given gene locus
# Example for CEBPA

# get peakset
# note: addPeakSet() can be used to assign the peakset if multiple peaksets exist
# only one peakset can be active at a time

print("Getting peakset from ArchR project")
pSet <- getPeakSet(ArchRProj = proj)
pSet$name <- paste(seqnames(pSet), start(pSet), end(pSet), sep = "_")
print(paste("Retrieved the peakset:", pSet))

# get motif matches for the peakset
print("Getting motif matches for the peakset")
matches <- getMatches(ArchRProj = proj, name = "Motif")
rownames(matches) <- paste(seqnames(matches), start(matches), end(matches), sep = "_")
matches <- matches[pSet$name]

# GR object corresponding to CEBPA promoter
print("Finding motifs in CEBPA promoter region")
gr <- GRanges(seqnames = c("chr19"), ranges = IRanges(start = c(33792929), end = c(33794030)))
queryHits <- queryHits(findOverlaps(query = pSet, subject = gr, type = "within"))


# Motif Enrichment in differentially Accessible Peaks

motifsUp <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = proj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

motifsDown <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = proj,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
  )

# Prepare motif enrichment SE objects for plotting
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
            nudge_x = 2,
            color = "black"
      ) + theme_ArchR() + 
      ylab("-log10(P-adj) Motif Enrichment") + 
      xlab("Rank Sorted TFs Enriched") +
      scale_color_gradientn(colors = paletteContinuous(set = "comet"))
    return(gg)
}

ggUp <- plotMotifEnrichments(motifsUp)
ggDoown <- plotMotifEnrichments(motifsDown)

plotPDF(ggUp, ggDown, name = "PIMOup-vs-PIMOdown-Markers-Motifs-Enriched", width = 5, height = 5, ArchRProj = proj, addDOC = TRUE)


EOF