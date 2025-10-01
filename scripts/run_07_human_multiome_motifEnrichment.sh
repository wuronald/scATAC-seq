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
# options include: homer, encode, JASPAR2020, JASPAR2018, and JASPAR2016
MOTIF_SET="${1:-cisbp}"  # Default to "cisbp" if no argument provided

# Export the parameter so R can access it
export MOTIF_SET

echo "Running motif enrichment analysis with motifSet: ${MOTIF_SET}"

# Load necessary modules (adjust as needed for your system)
module load R/4.4.1

# Run R script with motifSet parameter
Rscript - <<'EOF'

# load libraries
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(here)
set.seed(1)

# Get motifSet from environment variable
motifSet <- Sys.getenv("MOTIF_SET", unset = "cisbp")
cat("Using motifSet:", motifSet, "\n")

# Set the number of threads for ArchR
addArchRThreads(threads = 18)

# set genome to hg38
addArchRGenome("hg38")

# Load the project
proj <- loadArchRProject(path = "human_multiome_harmony_merged_malig_peak")

# Load markerTest from pairwise differential peaks getMarkerFeatures(): 
# path to saved markersPeaks RDS file
markersTest_rds <- "human_multiome_harmony_merged_malig_peak/PeakCalls/markersTest_PIMO_up_status_PIMOup_vs_PIMOdown.rds"

if (!exists("markerTest")) {
    print("loading markerTest")
    if (file.exists(markersTest_rds)) {
        print("previously saved markersPeaks loaded:")
        print(markersTest_rds)
        markerTest <- readRDS(file = markersTest_rds)
    } else {
        print("extracting markersPeaks")
        markerTest <- getMarkerFeatures(
            ArchRProj = proj,
            useMatrix = "PeakMatrix",
            groupBy = "PIMO_up_status", # this needs to be changed based on the comparison of interest
            bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
            testMethod = "wilcoxon"
        )
    }
}

# load markerPeaks from non-pairwise getMarkerFeatures():
# path to saved markersPeaks RDS file
markersPeaks_rds <- "human_multiome_harmony_merged_malig_peak/PeakCalls/markersPeaks_PIMO_up_status.rds"
if (!exists("markersPeak")) {
    print("loading markersPeaks")
    if (file.exists(markersPeaks_rds)) {
        print("previously saved markersPeaks loaded")
        markersPeaks <- readRDS(file = markersPeaks_rds)
    } else {
        print("extracting markersPeaks")
        markersPeaks <- getMarkerFeatures(
            ArchRProj = proj,
            useMatrix = "PeakMatrix",
            groupBy = "PIMO_up_status",       # this needs to be changed based on the comparison of interest
            bias = c("TSSEnrichment", "log10(nFrags)", "log10(Gex_nUMI)"),
            testMethod = "wilcoxon"
        )
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

motifsUp <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = proj,
    peakAnnotation = annoName,
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

motifsDown <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = proj,
    peakAnnotation = annoName,
    cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
  )

# Export SummarizedExperiment Motif enrichment results
print("Exporting motif enrichment results")
outDir <- here(paste0("human_multiome_harmony_merged_malig_peak/motifEnrichment_", motifSet, "/"))
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)
saveRDS(motifsUp, file = file.path(outDir, paste0("motifsUp_PIMOup_vs_PIMOdown_", motifSet, ".rds")))
saveRDS(motifsDown, file = file.path(outDir, paste0("motifsDown_PIMOup_vs_PIMOdown_", motifSet, ".rds")))

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

plotPDF(ggUp, ggDown, name = paste0("PIMOup-vs-PIMOdown-Markers-Motifs-Enriched_", motifSet), width = 5, height = 5, ArchRProj = proj, addDOC = TRUE)


# Motif Enrichment in marker peaks (non-pairwise))
print("Plotting ggplot of motif enrichment results in non-pairwise marker peaks")
enrichMotifs <- peakAnnoEnrichment(
    seMarker = markersPeaks,
    ArchRProj = proj,
    peakAnnotation = annoName,
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5" # only upregulated motifs
  )

heatmapEM <- plotEnrichHeatmap(enrichMotifs, n = 7, transpose = TRUE)
plotPDF(heatmapEM, name = paste0("PIMO_up_status-upregulated_Motifs-Enriched-Marker-Heatmap_", motifSet), width = 8, height = 6, ArchRProj = proj, addDOC = TRUE)

########################################################
# Compute chromeVar Deviations
print("Computing chromVar deviations")
proj <- addBgdPeaks(proj, force = TRUE)
proj <- addDeviationsMatrix(ArchRProj = proj, peakAnnotation = annoName,
        matrixName = paste0("MotifMatrix_", motifSet), # name of the deviations matrix
        force = TRUE
        )

# Plot Variability of Motif Deviations
print("Plotting variability of motif deviations")
plotVarDev <- getVarDeviations(proj, name = paste0("MotifMatrix_", motifSet), plot = FALSE)

print("Saving variability of motif deviations plot and data")
saveRDS(plotVarDev, file = file.path(outDir, paste0("chromVarDeviations_PIMOup_vs_PIMOdown_", motifSet, ".rds")))

plotVarDev <- getVarDeviations(proj, name = paste0("MotifMatrix_", motifSet),
        n = 25, # label the top 25 most variable motifs
        plot = TRUE
        ) # set plot = TRUE to get the ggplot object
plotPDF(plotVarDev, name = paste0("PIMOup-vs-PIMOdown-Variable-Motif-Deviation-Scores_", motifSet), width = 5, height = 5, ArchRProj = proj, addDOC = TRUE)
########################################################
moi <- c("SOX","HIF","ARNT","NF1","NFI")

# add impute weights
proj <- addImputeWeights(proj,
reducedDims = "Harmony_LSI_Combined"
)

print(paste("Finding motif deviations for motifs of interest:", paste(moi, collapse=", ")))
markerMotifs <- getFeatures(proj, select = paste(moi, collapse="|"), useMatrix = paste0("MotifMatrix_", motifSet))
markerMotifs <- grep("z:", markerMotifs, value = TRUE)
print(paste("Marker motifs found:", paste(markerMotifs, collapse = ", ")))

# Plot motif deviations for motifs of interest
print("Plotting motif deviations for motifs of interest")
p <- plotGroups(ArchRProj = proj, 
  groupBy = "PIMO_up_status", 
  colorBy = paste0("MotifMatrix_", motifSet), 
  name = markerMotifs,
  imputeWeights = getImputeWeights(proj) # might be missing
)
# customize the plots
p2 <- lapply(seq_along(p), function(x){
  if(x != 1){
    p[[x]] + guides(color = "none", fill = "none") + 
    theme_ArchR(baseSize = 6) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    theme(
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
    ) + ylab("")
  }else{
    p[[x]] + guides(color = "none", fill = "none") + 
    theme_ArchR(baseSize = 6) +
    theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
    theme(
        axis.ticks.y=element_blank(),
        axis.title.y=element_blank()
    ) + ylab("")
  }
})
# do.call(cowplot::plot_grid, c(list(nrow = 1, rel_widths = c(2, rep(1, length(p2) - 1))),p2))

plotPDF(p, name = paste0("PIMOup-vs-PIMOdown-Groups-Deviations-w-Imputation_", motifSet), width = 5, height = 5, ArchRProj = proj, addDOC = TRUE)

p <- plotEmbedding(
    ArchRProj = proj, 
    colorBy = paste0("MotifMatrix_", motifSet), 
    name = sort(markerMotifs), 
    embedding = "UMAP_Harmony_LSI_Combined"
)
p2 <- lapply(p, function(x){
    x + guides(color = "none", fill = "none") + 
    theme_ArchR(baseSize = 6.5) +
    theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
    theme(
        axis.text.x=element_blank(), 
        axis.ticks.x=element_blank(), 
        axis.text.y=element_blank(), 
        axis.ticks.y=element_blank()
    )
})

EOF