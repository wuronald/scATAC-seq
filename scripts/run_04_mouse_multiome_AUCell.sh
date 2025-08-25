#!/bin/bash
#SBATCH --job-name=run_04_mouse_multiome_AUCell
#SBATCH --output=run_04_mouse_multiome_AUCell_%j.out
#SBATCH --error=run_04_mouse_multiome_AUCell_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=60G
#SBATCH --time=02:00:00

# Load necessary modules (adjust as needed for your system)
module load R/4.4.1

# Run R script
Rscript - <<'EOF'

# load libraries
library(here)
library(ArchR)
library(AUCell)
#library(NMF)
#library(tidyverse)
library(ggrastr)
#library(Seurat)
set.seed(1)

# Set the number of threads for ArchR
addArchRThreads(threads = 18)

# Load the project
proj <- loadArchRProject(path = "mouse_multiome_harmony_merged_subset")

# create a folder for AUCell results
dir.create("mouse_multiome_harmony_merged_subset/AUCell/")

# load exprMatrix if not already in memory; otherwise create it

if (!exists("exprMatrix")) {
	print("loading exprMatrix")
  if (file.exists("mouse_multiome_harmony_merged_subset/AUCell/mouse_multiome_ExprMatrix.RData")) {
  	print("previously saved exprMatrix loaded")
    load(file = "mouse_multiome_harmony_merged_subset/AUCell/mouse_multiome_ExprMatrix.RData")
  }
else{
  	print("extracting exprMatrix")
  # load expression matrix (from ArchRProject object)
  ## Access the GeneScoreMatrix as SummarizedExperiment
  exprMatrix <- getMatrixFromProject(proj,
                                     useMatrix = "GeneExpressionMatrix")
  	print("saving exprMatrix")  
  # save the matrix for easier load in the future
  save(exprMatrix, file = "mouse_multiome_harmony_merged_subset/AUCell/mouse_multiome_ExprMatrix.RData") # xx gb
  
}
}

# check if exprMatrix is a SummarizedExperiment object; needed for AUCell
print("Checking if exprMatrix is a SummarizedExperiment object")

class(exprMatrix)
rownames(exprMatrix) # null

exprMatrix@elementMetadata$name
SummarizedExperiment::rowData(exprMatrix)
SummarizedExperiment::rowData(exprMatrix)$name

rownames(rowData(exprMatrix))

# # add gene names to the rownames slot
print("add gene names to the rownames slot") 
# #rownames(assay(exprMatrix)) <- as.character(rowData(exprMatrix)$name)
rownames(exprMatrix) <- SummarizedExperiment::rowData(exprMatrix)$name

## Load gene sets

print("load gene sets for AUCell")
# Load Neftel mouse gene sets from .tsv files in the Signatures directory
txt_files <- list.files(path = here::here("data","Signatures"), pattern = "\\.tsv$", full.names = TRUE)

## read the .txt files into a list
geneSets <- lapply(txt_files, readr::read_lines)

## name the elements of the list based on the original file name (without path and extension)
names(geneSets) <- tools::file_path_sans_ext(basename(txt_files))

# *** Load additional gene sets (not txt files) *** #

PIMO_pos <- readr::read_csv(
  here::here(
    "data",
    "Signatures",
    "pimo_sig",
    "ST_PIMO_FFPE_spot_gene_signature_HGNC.csv"
  )
) %>%
  dplyr::filter(logFC >=0) %>%
  dplyr::pull(Gene)
  
PIMO_pos %>% length # 327
  

# Append additional gene sets to list (if needed)
geneSets <- append(geneSets, list(PIMO_pos = PIMO_pos))

# Function: calculate_overlap 
## Checks if genes in the geneSets are in the GeneScoreMatrix
## Note: Consider other gene sets if many genes are missing from the GeneScoreMatrix

calculate_overlap <- function(se, geneSets) {
  # extract gene symbols from the SummarizedExperiment
  se_genes <- rownames(se)
  # make data frame to store info
  results <- data.frame(
    "Gene_Set" = character(),
    "Set_Size" = integer(),
    "Overlap" = integer(),
    "Percentage_OL" = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Loop through each gene set
  for (set_name in names(geneSets)) {
    gene_set <- geneSets[[set_name]]
    
    ## Calculate overlap
    overlap <- length(intersect(se_genes, gene_set))
    set_size <- length(gene_set)
    percentage <- (overlap / set_size) * 100
    
    ## Add results to the data frame
    results <- rbind(
      results,
      data.frame(
        "Gene_Set" = set_name,
        "Set_Size" = set_size,
        "Overlap" = overlap,
        "Percentage_OL" = percentage,
        stringsAsFactors = FALSE
      )
    )
  }
  return(results)
}

calculate_overlap(exprMatrix, geneSets) # check if geneset genes are in the matrix


## score gene signatures
print("Scoring gene signatures via AUCell")

if (!exists("cells_AUC")) {
  if (file.exists("Terekhanova/AUCell/cells_AUC.RData")) {
    load(file = "Terekhanova/AUCell/cells_AUC.RData") # 2.2 mb
  }
  else{
    print("file not found. creating cells rankings and cell AUC")
    
    # calculates both cell ranking and enrichment
    # cells_AUC <- AUCell_run(exprMatrix, geneSets)
    
    # alternative: calculate cell ranking and enrichment seperately
    cell_rankings <- AUCell_buildRankings(exprMatrix, plotStats = FALSE)
    cells_AUC <- AUCell_calcAUC(geneSets, cell_rankings)
    
    # save cells_AUC
    save(cells_AUC, file = "Terekhanova/AUCell/cells_AUC.RData")
  }
}

## Explore AUCell Thresholds

set.seed(42)


# plot all AUC histograms (for all geneSets) and calculates likely thresholds

cells_assignment <- AUCell_exploreThresholds(cells_AUC, 
                                             plotHist = TRUE,
                                             assign = TRUE)
# display warning messages
warningMsg <- sapply(cells_assignment, function(x) x$aucThr$comment)
warningMsg[which(warningMsg!="")]

# look at possible thresholds for (DAEG_UP_24)
cells_assignment$DAEG_UP_24$aucThr$thresholds
cells_assignment$Neftel_MES2$aucThr$thresholds

# look at auto-selected thresholds for (DAEG_UP_24)
cells_assignment$DAEG_UP_24$aucThr$selected # minimumDens
# of cells assigned 
cells_assignment$DAEG_UP_24$assignment %>% length() # 23353

# plot all AUC histograms (for all geneSets) and calculates likely thresholds
	par(mfrow=c(3,3))

	# Get the current date
	current_date <- Sys.Date()

	# Format the date as "YYYY-MM-DD"
	formatted_date <- format(current_date, "%Y-%m-%d")

	# Print the formatted date
	print(formatted_date)

pdf(here::here("Terekhanova","Plots",paste0(formatted_date,"_Terekhanova_all_AUCell_cells_assignment.pdf")))
  AUCell_exploreThresholds(cells_AUC,
                           plotHist = TRUE,
                           assign = TRUE)
dev.off()

##############
# Add AUCell assignments to ArchRProject object
print("Adding AUCell assignments to ArchRProject Object")

# name of gene sets to add as new columns
new_cols <- cells_AUC@NAMES

# determine membership of each barcode for each assigned gene Set
new_cols_list <- list() # initialize list
  # loop through each gene set
for (i in names(cells_AUC)) {
  new_cols_list[[i]] <- cells_assignment[[i]]$assignment # extract assigned barcodes
  # creates new column in cellColData
  proj@cellColData[[i]] <-
    ifelse(getCellNames(proj) %in% new_cols_list[[i]], TRUE, FALSE) # determine membership
}

## add numeric AUC scores to the cellColData
print("Adding AUCell scores to ArchRProject Object")

auc_matrix <- assay(cells_AUC)

for (set_name in rownames(auc_matrix)) {
  proj <- addCellColData(
    ArchRProj = proj,
    data = auc_matrix[set_name, ],
    name = paste0("AUCell_", set_name),
    cells = colnames(auc_matrix)
  )
}


# Save the project 
print("Saving ArchRProject object")
saveArchRProject(ArchRProj = proj, outputDirectory = "Terekhanova", load = TRUE)

EOF

echo "AUCell analysis completed"