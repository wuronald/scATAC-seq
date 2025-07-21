#!/bin/bash
#SBATCH --job-name=plot_01_markerGenes_modules_optimized
#SBATCH --output=plot_01_markerGenes_modules_optimized_%j.out
#SBATCH --error=plot_01_markerGenes_modules_optimized_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=20G
#SBATCH --time=01:00:00

# Load necessary modules (adjust as needed for your system)
module load R/4.4.1

# Run R script
Rscript - <<'EOF'

# load libraries
library(ArchR)
library(BSgenome.Mmusculus.UCSC.mm10)
library(Seurat)
library(here)
set.seed(1)

# Set the number of threads for ArchR
addArchRThreads(threads = 18)

# Load the project
print("Loading ArchR project")
proj <- loadArchRProject(path = "mouse_multiome")

# extract gene symbols from the project with more robust approach
print("Extracting gene information from project")
genes_gr <- getGenes(proj)
gene_symbols <- mcols(genes_gr)$symbol

# Remove any NA or empty gene symbols
gene_symbols <- gene_symbols[!is.na(gene_symbols) & gene_symbols != ""]
gene_symbols <- unique(gene_symbols)

print(paste("Total valid gene symbols found:", length(gene_symbols)))

# Select marker genes for different cell types
# Below are Mouse Motor Cortex cell type makers from Hubmap/Azimuth
features <- list(
  # From Hamed et al 2025
  Choroid_plexus = c("Ttr", "Enpp2", "2900040C04Rik", "1500015O10Rik", "Calml4", "Igfbp2", "Chchd10", "Rbp1", "Folr1", "Ppp1r1b"),
  Ependymal = c("Rsph1", "Tmem212", "Rarres2", "Ccdc153", "Dbi", "Fam183b", "Pltp", "1110017D15Rik", "Mns1", "Dynlrb2"),
  Immature_OL =  c("Gpr17", "Fyn", "Sirt2", "Opcml", "Nfasc", "Bcas1", "Mpzl1", "Tnr", "Enpp6", "Col9a3"),
  Macrophages = c("Cd74", "H2-Aa", "Lyz2", "H2-Eb1", "H2-Ab1", "Cxcl2", "Ifitm3", "Ifitm2", "Crip1", "Lgals3"),
  MDSCs = c("S100a9", "S100a8", "Retnlg", "Wfdc17", "Ifitm2", "S100a11", "Hp", "Fxyd5", "Wfdc21", "Ifitm1"),
  Microglia = c("Hexb", "C1qa", "Ctss", "C1qb", "Csf1r", "P2ry12", "C1qc", "Cx3cr1", "Ccl4", "Selplg"),
  NBs = c("Meg3", "Sox11", "Stmn2", "Tubb3", "Tmsb10", "Tubb2b", "Igfbpl1", "Meis2", "Dlx6os1", "Ccnd2"),
  OLs = c("Plp1", "Mal", "Ptgds", "Cldn11", "Mag", "Cryab", "Qdpr", "Mog", "Car2", "Mbp"),
  T_cells = c("Trbc2", "Nkg7", "AW112010", "Ms4a4b", "Ets1", "Cd3g", "Emb", "Trbc1", "Crip1", "Fxyd5"),
  TAPs = c("Hmgb2", "Top2a", "Cenpf", "Ube2c", "Hist1h2ap", "Hist1h2ae", "Tubb5", "Hmgn2", "Hist1h1b", "H2afz")
)

print("Original features to be added as module scores:")
print(lapply(features, length))

# Robust gene filtering function
filter_genes <- function(gene_list, available_genes) {
  # Remove any NA or empty values
  gene_list <- gene_list[!is.na(gene_list) & gene_list != ""]
  # Filter to only include genes present in the dataset
  valid_genes <- gene_list[gene_list %in% available_genes]
  return(valid_genes)
}

# Check if all features are in gene_symbols and report missing genes
print("Checking gene availability and filtering...")
missing_genes <- list()
features_filtered <- list()

for(feature_name in names(features)) {
  original_genes <- features[[feature_name]]
  valid_genes <- filter_genes(original_genes, gene_symbols)
  missing <- original_genes[!original_genes %in% gene_symbols]
  
  # Store results
  features_filtered[[feature_name]] <- valid_genes
  missing_genes[[feature_name]] <- missing
  
  # Report results
  cat("\n", feature_name, ":\n")
  cat("  Original genes:", length(original_genes), "\n")
  cat("  Valid genes:", length(valid_genes), "\n")
  cat("  Missing genes:", length(missing), "\n")
  if(length(missing) > 0) {
    cat("  Missing:", paste(missing, collapse=", "), "\n")
  }
}

# Remove any feature sets that have fewer than 3 genes (too few for meaningful module score)
min_genes <- 3
features_filtered <- features_filtered[sapply(features_filtered, length) >= min_genes]

print(paste("\nFeature sets with at least", min_genes, "genes:"))
print(sapply(features_filtered, length))

# Check if we have any valid feature sets left
if(length(features_filtered) == 0) {
  stop("No valid feature sets found after filtering. Check your gene names and minimum gene threshold.")
}

# Verify gene expression matrix is available
print("Checking available matrices...")
print(getAvailableMatrices(proj))

# Add module scores one by one with error handling
print("Adding module scores individually...")
successful_modules <- c()
failed_modules <- c()

for(i in seq_along(features_filtered)) {
  feature_name <- names(features_filtered)[i]
  feature_genes <- features_filtered[[i]]
  
  print(paste("Adding module score for:", feature_name, "with", length(feature_genes), "genes"))
  print(paste("Genes:", paste(feature_genes, collapse=", ")))
  
  # Create a temporary feature list with just this one feature
  temp_features <- list()
  temp_features[[feature_name]] <- feature_genes
  
  # Add module score with error handling
  tryCatch({
    proj <- addModuleScore(proj,
        useMatrix = "GeneExpressionMatrix",
        name = paste0("Module_", feature_name),
        features = temp_features,
        nBin = 25,
        nBgd = 100,
        seed = 1
        )
    
    successful_modules <- c(successful_modules, feature_name)
    print(paste("✓ Successfully added module score for:", feature_name))
    
  }, error = function(e) {
    failed_modules <- c(failed_modules, feature_name)
    print(paste("✗ Failed to add module score for:", feature_name))
    print(paste("Error:", e$message))
  })
}

# Report results
print(paste("Successfully added", length(successful_modules), "module scores"))
print(paste("Failed to add", length(failed_modules), "module scores"))

if(length(failed_modules) > 0) {
  print("Failed modules:")
  print(failed_modules)
}

# Verify that the module scores were added
print("Checking cellColData for module scores:")
cell_col_names <- names(proj@cellColData)
module_cols <- cell_col_names[grepl("Module_", cell_col_names)]
print(module_cols)

# Initialize plot_list outside the conditional block
plot_list <- list()

# Only proceed with plotting if we have successful modules
if(length(successful_modules) > 0) {
  
  # Create plots for successful module scores
  print("Creating plots for successful module scores")
  
  for(feature_name in successful_modules) {
    module_name <- paste0("Module_", feature_name, "_", feature_name)
    
    print(paste("Creating plot for:", module_name))
    
    tryCatch({
      p <- plotEmbedding(proj,
        embedding = "UMAP_Harmony_LSI_Combined",
        colorBy = "cellColData",
        name = module_name,
        plotAs = "points",
        size = 0.5,
        baseSize = 10,
        legendSize = 8)
      
      plot_list[[feature_name]] <- p
      
    }, error = function(e) {
      print(paste("Failed to create plot for:", module_name))
      print(paste("Error:", e$message))
    })
  }
  
  # Save plots to PDF if we have any successful plots
  if(length(plot_list) > 0) {
    print("Saving plots to PDF")
    do.call(plotPDF, c(plot_list, 
                       list(name = "mouse_multiome-scATAC-scRNA-Harmony_Combined_colorby-moduleScore_ALL.pdf", 
                            ArchRProj = proj,
                            addDOC = TRUE,
                            width = 5, height = 5)))
  }
  
} else {
  print("No successful module scores to plot")
}

# Create a comprehensive summary
print("Creating comprehensive summary")
summary_data <- data.frame(
  Module = names(features),
  Original_Gene_Count = sapply(features, length),
  Valid_Gene_Count = sapply(names(features), function(x) {
    if(x %in% names(features_filtered)) length(features_filtered[[x]]) else 0
  }),
  Missing_Gene_Count = sapply(missing_genes, length),
  Module_Added = names(features) %in% successful_modules,
  Plot_Created = names(features) %in% names(plot_list)
)

print(summary_data)

# Save summary to file
write.csv(summary_data, "module_gene_summary.csv", row.names = FALSE)

# Save detailed gene information
detailed_summary <- list()
for(feature_name in names(features)) {
  detailed_summary[[feature_name]] <- list(
    original_genes = features[[feature_name]],
    valid_genes = if(feature_name %in% names(features_filtered)) features_filtered[[feature_name]] else character(0),
    missing_genes = missing_genes[[feature_name]]
  )
}

# Save detailed summary as RDS for future reference
saveRDS(detailed_summary, "detailed_gene_summary.rds")

print("Script completed!")
print(paste("Total successful module scores:", length(successful_modules)))
print(paste("Total plots created:", length(plot_list)))

EOF