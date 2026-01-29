#!/bin/bash
#SBATCH --job-name=run_08_Gaiti_multiome_chromVarDeviations
#SBATCH --output=run_08_Gaiti_multiome_chromVarDeviations_%j.out
#SBATCH --error=run_08_Gaiti_multiome_chromVarDeviations_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=58G
#SBATCH --time=02:00:00

# Parse command line arguments
# First argument: motif set (homer, encode, JASPAR2020, JASPAR2018, JASPAR2016, cisbp)
# Second argument: groupBy variable (e.g., "PIMO_up_status", "PIMO_Region")
# Third argument (optional): comma-separated motifs of interest (e.g., "SOX,HIF,ARNT,NF1,NFI")
# Example usage:
# sbatch scripts/run_08_Gaiti_multiome_chromVarDeviations.sh cisbp PIMO_Region "SOX,HIF,ARNT,NF1,NFI"
# sbatch scripts/run_08_Gaiti_multiome_chromVarDeviations.sh homer PIMO_up_status "HIF,ATF,FOS,JUN,AP-1,AP1,Bach"

MOTIF_SET="${1:-homer}"  # Default to "homer" if no argument provided
GROUP_BY="${2:-PIMO_Region}"  # Default to "PIMO_Region" if no argument provided
MOTIFS_OF_INTEREST="${3:-HIF,ATF,FOS,FRA,JUN,AP-1,AP1,Bach}"  # Default motifs if not provided

# Export the parameters so R can access them
export MOTIF_SET
export GROUP_BY
export MOTIFS_OF_INTEREST

echo "Running chromVAR deviations analysis with:"
echo "  motifSet: ${MOTIF_SET}"
echo "  groupBy: ${GROUP_BY}"
echo "  motifs of interest: ${MOTIFS_OF_INTEREST}"

# Load necessary modules (adjust as needed for your system)
module load R/4.4.1

# Run R script with parameters
Rscript - <<'EOF'

# load libraries
library(ArchR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(here)
library(ggplot2)
set.seed(1)

# Get parameters from environment variables
motifSet <- Sys.getenv("MOTIF_SET", unset = "cisbp")
groupBy <- Sys.getenv("GROUP_BY", unset = "PIMO_Region")
motifsOfInterest <- Sys.getenv("MOTIFS_OF_INTEREST", unset = "SOX,HIF,ARNT,NF1,NFI")

# Convert comma-separated string to vector
moi <- trimws(unlist(strsplit(motifsOfInterest, ",")))

cat("Using motifSet:", motifSet, "\n")
cat("Using groupBy:", groupBy, "\n")
cat("Motifs of interest:", paste(moi, collapse = ", "), "\n")

# Set the number of threads for ArchR
addArchRThreads(threads = 8)

# set genome to hg38
addArchRGenome("hg38")

# Load the project
proj <- loadArchRProject(path = "Gaiti_multiome_harmony_merged_malig_peak")

# Create file naming prefix based on groupBy variable
filePrefix <- gsub("_", "", groupBy)

# Create output directory
annoName <- paste0("Motif_", motifSet)
outDir <- here(paste0("Gaiti_multiome_harmony_merged_malig_peak/chromVarDeviations_", motifSet, "_", groupBy, "/"))
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

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
     # Add peak matrix
    print("adding peak matrix")
    proj <- addPeakMatrix(proj, force = TRUE)
    print("Peak matrix added to ArchR project")
    
} else {
    stop(paste("PeakSet file not found:", peakSetPath))
}

# Check if motif annotations exist, if not add them
print(paste("Checking if motif annotations exist for", annoName))
if (annoName %in% names(proj@peakAnnotation)) {
    print(paste("Motif annotations", annoName, "already exist in the ArchR project."))
} else {
    print(paste("Motif annotations do not exist. Adding them with motifSet:", motifSet))
    proj <- addMotifAnnotations(ArchRProj = proj, motifSet = motifSet, annoName = annoName)
}

########################################################
# Compute chromVAR Deviations
########################################################

print("Adding background peaks for chromVAR")
proj <- addBgdPeaks(proj, force = TRUE)

print(paste("Computing chromVAR deviations for motif set:", motifSet))
deviationsMatrixName <- paste0("MotifMatrix_", motifSet)
proj <- addDeviationsMatrix(
    ArchRProj = proj, 
    peakAnnotation = annoName,
    matrixName = deviationsMatrixName,
    force = TRUE
)

########################################################
# Plot Variability of Motif Deviations
########################################################
print("Computing variability of motif deviations")
plotVarDev <- getVarDeviations(proj, name = deviationsMatrixName, plot = FALSE)

print("Saving variability of motif deviations data")
saveRDS(plotVarDev, file = file.path(outDir, paste0("chromVarDeviations_", filePrefix, "_", motifSet, ".rds")))

print("Plotting variability of motif deviations (top 25)")
plotVarDev_plot <- getVarDeviations(
    proj, 
    name = deviationsMatrixName,
    n = 25, # label the top 25 most variable motifs
    plot = TRUE
)

plotPDF(
    plotVarDev_plot, 
    name = paste0(filePrefix, "-Variable-Motif-Deviation-Scores_", motifSet), 
    width = 5, 
    height = 5, 
    ArchRProj = proj, 
    addDOC = TRUE
)

########################################################
# Add Imputation Weights
########################################################
print("Adding imputation weights to help smooth plotting of motif deviations")

# Check if the reducedDims "Harmony_LSI_Combined" exists
if ("Harmony_LSI_Combined" %in% names(proj@reducedDims)) {
    print("ReducedDims 'Harmony_LSI_Combined' exists.")
} else {
    stop("ReducedDims 'Harmony_LSI_Combined' does not exist in the ArchR project.")
}

proj <- addImputeWeights(
    proj,
    reducedDims = "Harmony_LSI_Combined"
)

########################################################
# Plot Motif Deviations for Motifs of Interest
########################################################
print(paste("Finding motif deviations for motifs of interest:", paste(moi, collapse = ", ")))
markerMotifs <- getFeatures(proj, select = paste(moi, collapse = "|"), useMatrix = deviationsMatrixName)
markerMotifs <- grep("z:", markerMotifs, value = TRUE)
print(paste("Marker motifs found:", paste(markerMotifs, collapse = ", ")))

if (length(markerMotifs) == 0) {
    warning("No marker motifs found matching the motifs of interest!")
} else {
    # Plot motif deviations by group
    print("Plotting motif deviations by group")
    p <- plotGroups(
        ArchRProj = proj, 
        groupBy = groupBy, 
        colorBy = deviationsMatrixName, 
        name = markerMotifs,
        imputeWeights = getImputeWeights(proj)
    )
    
    # Customize the plots
    p2 <- lapply(seq_along(p), function(x){
        if(x != 1){
            p[[x]] + guides(color = "none", fill = "none") + 
            theme_ArchR(baseSize = 6) +
            theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
            theme(
                axis.text.y = element_blank(), 
                axis.ticks.y = element_blank(),
                axis.title.y = element_blank()
            ) + ylab("")
        } else {
            p[[x]] + guides(color = "none", fill = "none") + 
            theme_ArchR(baseSize = 6) +
            theme(plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "cm")) +
            theme(
                axis.ticks.y = element_blank(),
                axis.title.y = element_blank()
            ) + ylab("")
        }
    })
    
    plotPDF(
        p, 
        name = paste0(filePrefix, "-Groups-Deviations-w-Imputation_", motifSet), 
        width = 5, 
        height = 5, 
        ArchRProj = proj, 
        addDOC = TRUE
    )
    
    # Plot motif deviations on UMAP embedding
    print("Plotting motif deviations on UMAP embedding")
    p_umap <- plotEmbedding(
        ArchRProj = proj,
        colorBy = deviationsMatrixName, 
        name = sort(markerMotifs), 
        embedding = "UMAP_Harmony_LSI_Combined",
        imputeWeights = getImputeWeights(proj)
    )
    
    # Customize UMAP plots
    p_umap2 <- lapply(p_umap, function(x){
        x + 
        theme_ArchR(baseSize = 6.5) +
        theme(plot.margin = unit(c(0, 0, 0, 0), "cm")) +
        theme(
            axis.text.x = element_blank(), 
            axis.ticks.x = element_blank(), 
            axis.text.y = element_blank(), 
            axis.ticks.y = element_blank()
        )
    })
    print("Saving UMAP motif deviation plots")
    plotPDF(
        p_umap2, 
        name = paste0(filePrefix, "-UMAP-Deviations-w-Imputation_", motifSet), 
        width = 5, 
        height = 5, 
        ArchRProj = proj, 
        addDOC = TRUE
    )

}
print("chromVAR deviations analysis complete!")

########################################################
# Identification of Positive TF-regulators
########################################################

# Step 1. Identify deviant TF motifs
print("Step 1: Identifying deviant TF motifs")

seGroupMotif <- getGroupSE(
    ArchRProj = proj, 
    groupBy = groupBy, 
    useMatrix = deviationsMatrixName
)

seGroupMotif

# subset to only deviation z-scores
print("Subsetting to deviation z-scores only")
seZ <- seGroupMotif[rowData(seGroupMotif)$seqnames=="z",]

# Calculate max delta deviation for each motif
print("Calculating max delta deviation for each motif")
rowData(seZ)$maxDelta <- lapply(seq_len(ncol(seZ)), function(x){
  rowMaxs(assay(seZ) - assay(seZ)[,x])
}) %>% Reduce("cbind", .) %>% rowMaxs


# Step 2. Identify Correlated TF Motifs and TF Gene Score/Expression
print("Step 2: Identifying correlated TF motifs and TF gene expression")

corGSM_MM <- correlateMatrices(
    ArchRProj = proj,
    useMatrix1 = "GeneExpressionMatrix",
    useMatrix2 = deviationsMatrixName,
    removeFromName2 = "dot", # keeps motif to the left of the dot (for homer motifs) 
    reducedDims = "Harmony_LSI_Combined"
)

corGSM_MM

# Step 3. Add Maximum Delta Deviation to the Correlation Data Frame
# note: MotifMatrix_name might not exist? 
print("Step 3: Adding maximum delta deviation to the correlation data frame")
corGSM_MM$maxDelta <- rowData(seZ)[match(corGSM_MM[[paste0(deviationsMatrixName, "_name")]], rowData(seZ)$name), "maxDelta"]


# Step 4. Identify Positive TF Regulators
# threshold requirements:
# 1. TFs whose correlation between motif and gene score (or gene expression) is greater than 0.5
# 2. adjusted p-value less than 0.01
# 3. maximum inter-cluster difference in deviation z-score that is in the top quartile.
print("Step 4: Identifying positive TF regulators based on thresholds")

corGSM_MM <- corGSM_MM[order(abs(corGSM_MM$cor), decreasing = TRUE), ]
corGSM_MM <- corGSM_MM[which(!duplicated(gsub("\\-.*","",corGSM_MM[,paste0(deviationsMatrixName, "_name")]))), ]
corGSM_MM$TFRegulator <- "NO"
corGSM_MM$TFRegulator[which(corGSM_MM$cor > 0.5 & corGSM_MM$padj < 0.01 & corGSM_MM$maxDelta > quantile(corGSM_MM$maxDelta, 0.75))] <- "YES"
sort(corGSM_MM[corGSM_MM$TFRegulator=="YES",1])

# export rds of DataFrame
print("Saving positive TF regulators DataFrame")
saveRDS(corGSM_MM, file = file.path(outDir, paste0("Positive-TF-Regulators_", filePrefix, "_", motifSet, ".rds")))

# Plot positive TF regulators
print("Plotting positive TF regulators")

p <- ggplot(data.frame(corGSM_MM), aes(cor, maxDelta, color = TFRegulator)) +
  geom_point() + 
  theme_ArchR() +
  geom_vline(xintercept = 0, lty = "dashed") + 
  scale_color_manual(values = c("NO"="darkgrey", "YES"="firebrick3")) +
  xlab("Correlation To Gene Score") +
  ylab("Max TF Motif Delta") +
  scale_y_continuous(
    expand = c(0,0), 
    limits = c(0, max(corGSM_MM$maxDelta)*1.05)
  )

plotPDF(p, name = paste0(filePrefix, "-Positive-TF-regulators-Motifs-Enriched_", motifSet), width = 5, height = 5, ArchRProj = proj, addDOC = TRUE)

# Plot UMAPS of Example Gene Expression and Groupings
print("Plotting UMAPs of example gene expression and groupings")
g1 <- plotEmbedding(proj, 
colorBy = "GeneExpressionMatrix", 
name = "VEGFA", 
embedding = "UMAP_Harmony_LSI_Combined")

# set custom discrete color palette
PIMO_up_status_colors <- c("PIMOdown" = "blue", "PIMOinter" = "gold", "PIMOup" = "red")

g2 <- plotEmbedding(proj, 
colorBy = "cellColData", 
name = "PIMO_up_status",
pal = PIMO_up_status_colors, 
embedding = "UMAP_Harmony_LSI_Combined")

g3 <- plotEmbedding(proj, 
colorBy = "cellColData", 
name = "Sample", 
embedding = "UMAP_Harmony_LSI_Combined")

plotPDF(list(g1, g2, g3), name = paste0(filePrefix, "-Example-GeneExpression-and-Groupings_", motifSet), width = 5, height = 5, ArchRProj = proj, addDOC = TRUE)

EOF