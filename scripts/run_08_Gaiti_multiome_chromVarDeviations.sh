#!/bin/bash
#SBATCH --job-name=run_08_Gaiti_multiome_chromVarDeviations
#SBATCH --output=run_08_Gaiti_multiome_chromVarDeviations_%j.out
#SBATCH --error=run_08_Gaiti_multiome_chromVarDeviations_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=60G
#SBATCH --time=02:00:00

# Parse command line arguments
# First argument:  motif set (homer, encode, JASPAR2020, JASPAR2018, JASPAR2016, cisbp)
# Second argument: groupBy variable (e.g., "PIMO_up_status", "PIMO_Region", "Region.annotation")
# Third argument:  ArchR project name (e.g., "Gaiti_multiome_harmony_merged_malig_peak")
# Fourth argument (optional): comma-separated motifs of interest (e.g., "SOX,HIF,ARNT,NF1,NFI")
# Fifth argument  (optional): comma-separated groupBy levels to exclude from plotGroups (e.g., "PIMOup,PIMOinter")
# Sixth argument  (optional): comma-separated metadata column names to combine (e.g., "combined_program,Region.annotation")
# Example usage:
# sbatch scripts/run_08_Gaiti_multiome_chromVarDeviations.sh homer PIMO_up_status
# sbatch scripts/run_08_Gaiti_multiome_chromVarDeviations.sh cisbp PIMO_Region Gaiti_multiome_harmony_merged_malig_peak "SOX,HIF,ARNT,NF1,NFI"
# sbatch scripts/run_08_Gaiti_multiome_chromVarDeviations.sh homer combined_program Gaiti_multiome_hmmp_k3_program_all_cells "" "PIMOup"
# sbatch scripts/run_08_Gaiti_multiome_chromVarDeviations.sh cisbp combined_program Gaiti_multiome_hmmp_k3_program_all_cells "SOX,HIF" "PIMOup,PIMOinter"
# sbatch scripts/run_08_Gaiti_multiome_chromVarDeviations.sh homer "" Gaiti_multiome_hmmp_k3_program_all_cells "SOX,HIF" "" "combined_program,Region.annotation"
# sbatch scripts/run_08_Gaiti_multiome_chromVarDeviations.sh homer "" Gaiti_multiome_hmmp_k3_program_all_cells "HIF,ATF,FOS,JUN,SOX,OLIG" "program_1_PT,program_2_PT,program_3_PT,PIMOup_PT,PIMOinter_PT,PIMOdown_PT" "combined_program,Region.annotation"

MOTIF_SET="${1:-homer}"  # Default to "homer" if no argument provided
GROUP_BY="${2:-PIMO_Region}"  # Default to "PIMO_Region" if no argument provided
ARCHR_PROJECT="${3:-Gaiti_multiome_harmony_merged_malig_peak}"  # Default project name if not provided
MOTIFS_OF_INTEREST="${4:-HIF,ATF,FOS,FRA,JUN,AP-1,AP1,Bach,NRF,MAF,RFX,NF,SOX,OLIG,Neuro,ASCL}"  # Default motifs if not provided
EXCLUDE_GROUPS="${5:-}"  # Default to empty (no exclusions) if not provided
COMBINE_COLUMNS="${6:-}" # Default to empty (no combined columns) if not provided

# Export the parameters so R can access them
export MOTIF_SET
export GROUP_BY
export ARCHR_PROJECT
export MOTIFS_OF_INTEREST
export EXCLUDE_GROUPS
export COMBINE_COLUMNS

echo "Running chromVAR deviations analysis with:"
echo "  motifSet: ${MOTIF_SET}"
echo "  groupBy: ${GROUP_BY}"
echo "  ArchR project: ${ARCHR_PROJECT}"
echo "  motifs of interest: ${MOTIFS_OF_INTEREST}"
echo "  exclude groups: ${EXCLUDE_GROUPS:-none}"
echo "  combine columns: ${COMBINE_COLUMNS:-none}"

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
motifSet          <- Sys.getenv("MOTIF_SET",          unset = "homer")
groupBy           <- Sys.getenv("GROUP_BY",           unset = "PIMO_Region")
archrProject      <- Sys.getenv("ARCHR_PROJECT",      unset = "Gaiti_multiome_harmony_merged_malig_peak")
motifsOfInterest  <- Sys.getenv("MOTIFS_OF_INTEREST", unset = "SOX,HIF,ARNT,NF1,NFI")
excludeGroupsEnv  <- Sys.getenv("EXCLUDE_GROUPS",     unset = "")
combineColumnsEnv <- Sys.getenv("COMBINE_COLUMNS",    unset = "")

# Convert comma-separated string to vector
moi <- trimws(unlist(strsplit(motifsOfInterest, ",")))

# Convert comma-separated exclude groups to vector (empty vector if unset)
excludeGroups <- if (nchar(excludeGroupsEnv) > 0) trimws(unlist(strsplit(excludeGroupsEnv, ","))) else character(0)

# Convert comma-separated combine columns to vector (empty vector if unset)
combineCols <- if (nchar(combineColumnsEnv) > 0) trimws(unlist(strsplit(combineColumnsEnv, ","))) else character(0)

cat("Using motifSet:", motifSet, "\n")
cat("Using initial groupBy:", groupBy, "\n")
cat("Using ArchR project:", archrProject, "\n")
cat("Motifs of interest:", paste(moi, collapse = ", "), "\n")
cat("Exclude groups from plotGroups:", if (length(excludeGroups) > 0) paste(excludeGroups, collapse = ", ") else "none", "\n")
cat("Combine columns:", if (length(combineCols) > 0) paste(combineCols, collapse = ", ") else "none", "\n")

# Set the number of threads for ArchR
addArchRThreads(threads = 18)

# set genome to hg38
addArchRGenome("hg38")

# Load the project
proj <- loadArchRProject(path = archrProject)

# If combine columns are specified, check existence and create new metadata column
if (length(combineCols) > 0) {
    cat("Checking columns to combine:", paste(combineCols, collapse = ", "), "\n")
    availableCols <- names(proj@cellColData)
    missingCols <- combineCols[!combineCols %in% availableCols]
    
    if (length(missingCols) > 0) {
        stop(paste0(
            "The following column(s) specified in COMBINE_COLUMNS do not exist in project metadata: ",
            paste(missingCols, collapse = ", "),
            "\nAvailable columns are: ",
            paste(availableCols, collapse = ", ")
        ))
    }
    
    # Create the new combined column name and values
    newColName <- paste(combineCols, collapse = "_")
    cat("Creating new metadata column:", newColName, "\n")
    
    # Combine the column contents with an underscore
    combinedVals <- do.call(paste, c(lapply(combineCols, function(col) proj@cellColData[[col]]), sep = "_"))
    
    # Add new column to cellColData
    proj <- addCellColData(
        ArchRProj = proj,
        data      = combinedVals,
        name      = newColName,
        cells     = getCellNames(proj),
        force     = TRUE
    )
    
    # Override groupBy with the newly created combined column
    groupBy <- newColName
    cat("groupBy has been updated to:", groupBy, "\n")
}

# Create file naming prefix based on groupBy variable
filePrefix <- gsub("_", "", groupBy)

# Create output directory
annoName <- paste0("Motif_", motifSet)
outDir <- here(paste0(archrProject, "/chromVarDeviations_", motifSet, "_", groupBy, "/"))
dir.create(outDir, showWarnings = FALSE, recursive = TRUE)

# Load and assign appropriate PeakSet based on groupBy
peakSetPath <- paste0("/cluster/projects/wouterslab/ArchR103_4/", archrProject, "/PeakCalls/PeakSet_gr_", groupBy, ".rds")
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
    if (!is.null(proj@peakSet) && length(proj@peakSet) > 0) {
        print(paste("Warning: PeakSet file not found at", peakSetPath, "- proceeding with existing peakSet in ArchR project (", length(proj@peakSet), "peaks)"))
    } else {
        stop(paste("PeakSet file not found:", peakSetPath, "and no peakSet exists in ArchR project."))
    }
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
    # set custom discrete color palettes
    PIMO_up_status_colors <- c("PIMOdown" = "blue", "PIMOinter" = "gold", "PIMOup" = "red")
    combined_program_colors <- c(
        "PIMOdown"  = "blue",
        "PIMOinter" = "gold",
        "PIMOup"    = "red",
        "program_1" = "#FA8072",  # salmon
        "program_2" = "#01796F",  # pine green
        "program_3" = "#796A5B"   # greyish brown
    )
    group_pal <- if (groupBy == "PIMO_up_status") PIMO_up_status_colors else
                 if (groupBy == "combined_program") combined_program_colors else NULL

    # Subset the project to exclude unwanted groups from plotGroups, if any
    # FIX: Give the folder an completely unique prefix so ArchR's internal gsub doesn't break
    subsetTmpDir <- here(paste0("TMP_SUBSET_FOR_PLOTTING_", filePrefix))
    if (length(excludeGroups) > 0) {
        print(paste("Excluding groups from plotGroups:", paste(excludeGroups, collapse = ", ")))
        cellsToKeep <- proj$cellNames[!(proj@cellColData[[groupBy]] %in% excludeGroups)]
        print(paste("Cells before exclusion:", length(proj$cellNames), "| After:", length(cellsToKeep)))
        
        projSubset <- subsetArchRProject(
            ArchRProj = proj,
            cells = cellsToKeep,
            outputDirectory = subsetTmpDir,
            dropCells = TRUE,
            force = TRUE
        )
        # Re-add imputation weights on the subsetted project
        print("Re-adding imputation weights on subsetted project")
        projSubset <- addImputeWeights(
            projSubset,
            reducedDims = "Harmony_LSI_Combined"
        )
        # Drop excluded levels from the palette so the legend stays clean
        group_pal_subset <- group_pal[!names(group_pal) %in% excludeGroups]
    } else {
        print("No groups excluded; using full project for plotGroups")
        projSubset <- proj
        group_pal_subset <- group_pal
    }

    p <- plotGroups(
        ArchRProj = projSubset,
        pal = group_pal_subset,
        groupBy = groupBy, 
        colorBy = deviationsMatrixName, 
        name = markerMotifs,
        imputeWeights = getImputeWeights(projSubset)
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

    # Clean up the temporary subset project directory if it was created
    if (length(excludeGroups) > 0 && dir.exists(subsetTmpDir)) {
        print(paste("Removing temporary subset directory:", subsetTmpDir))
        unlink(subsetTmpDir, recursive = TRUE)
        print("Temporary subset directory removed")
    }
    
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

# set custom discrete color palettes
PIMO_up_status_colors <- c("PIMOdown" = "blue", "PIMOinter" = "gold", "PIMOup" = "red")
combined_program_colors <- c(
    "PIMOdown"  = "blue",
    "PIMOinter" = "gold",
    "PIMOup"    = "red",
    "program_1" = "#FA8072",  # salmon
    "program_2" = "#01796F",  # pine green
    "program_3" = "#796A5B"   # greyish brown
)
group_pal <- if (groupBy == "PIMO_up_status") PIMO_up_status_colors else
             if (groupBy == "combined_program") combined_program_colors else NULL

g2 <- plotEmbedding(proj, 
colorBy = "cellColData", 
name = groupBy,
pal = group_pal, 
embedding = "UMAP_Harmony_LSI_Combined")

g3 <- plotEmbedding(proj, 
colorBy = "cellColData", 
name = "Sample", 
embedding = "UMAP_Harmony_LSI_Combined")

plotPDF(list(g1, g2, g3), name = paste0(filePrefix, "-Example-GeneExpression-and-Groupings_", motifSet), width = 5, height = 5, ArchRProj = proj, addDOC = TRUE)

EOF