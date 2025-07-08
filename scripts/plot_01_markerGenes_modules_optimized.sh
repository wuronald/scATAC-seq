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

# extract gene symbols from the project
genes_gr <- getGenes(proj)
gene_symbols <- mcols(genes_gr)$symbol

# Select marker genes for different cell types
features <- list(
  Glutamatergic_score = c("Celf2", "Ptprd", "Arpp21", "Sv2b", "Pcsk2", "Phactr1", "Kalrn", "Nrg1", "Ano3", "Satb2"),
  GABAergic_score  = c("Grip1", "Galntl6", "Gad1", "Gad2", "Dlx6os1", "Cntnap2", "Gm26905", "Kcnmb2", "Kcnip1", "Erbb4"),
  Endo_score  = c("Flt1", "Slco1a4", "Adgrl4", "Ly6c1", "Slc2a1", "Klf2", "Mecom", "Bsg", "Ly6a", "Pltp"),
  Oligo_score  = c("Plp1", "Mbp", "St18", "Prr5l", "Mobp", "Mal", "Mog", "Cldn11", "Pde4b", "Mag"),
  Micro_PVM_score = c("Inpp5d", "Hexb", "Tgfbr1", "C1qa", "Ctss", "C1qb", "Zfhx3", "C1qc", "Selplg", "Cx3cr1"),
  OPC_score  = c("Lhfpl3", "Vcan", "Tnr", "Ptprz1", "Gm4876", "Xylt1", "Pdgfra", "Epn2", "Cacng4", "Megf11"),
  Vlmc_score  =c("Ptgds", "Bnc2", "Cped1", "Slc7a11", "Bmp6", "Apod", "Mgp", "Eya2", "Ranbp3l", "Adam12")
)

print("Features to be added as module scores:")
print(features)

# Check if all features are in gene_symbols and report missing genes
print("Checking if all features are in gene symbols")
missing_genes <- lapply(names(features), function(name) {
  missing <- features[[name]][!features[[name]] %in% gene_symbols]
  if(length(missing) > 0) {
    cat("Missing genes in", name, ":", paste(missing, collapse=", "), "\n")
  }
  return(missing)
})
names(missing_genes) <- names(features)

# Filter features to only include genes present in the dataset
features_filtered <- lapply(features, function(x) x[x %in% gene_symbols])

# Remove any feature sets that have no genes left
features_filtered <- features_filtered[sapply(features_filtered, length) > 0]

print("Filtered features (genes present in dataset):")
print(sapply(features_filtered, length))

# Check if we have any valid feature sets left
if(length(features_filtered) == 0) {
  stop("No valid feature sets found after filtering. Check your gene names.")
}

# Add module scores one by one to avoid batch processing issues
print("Adding module scores individually...")
for(i in seq_along(features_filtered)) {
  feature_name <- names(features_filtered)[i]
  feature_genes <- features_filtered[[i]]
  
  print(paste("Adding module score for:", feature_name, "with", length(feature_genes), "genes"))
  
  # Create a temporary feature list with just this one feature
  temp_features <- list()
  temp_features[[feature_name]] <- feature_genes
  
  # Add module score
  proj <- addModuleScore(proj,
      useMatrix = "GeneExpressionMatrix",
      name = "Module",
      features = temp_features)
}
# Verify that the module scores were added
print("Check if Module scores added successfully to cellColData")
print(names(proj@cellColData))

# Create plots for ALL module scores efficiently
print("Creating plots for all module scores")

# Get all module score names (updated naming scheme)
module_names <- paste0("Module_", names(features_filtered), "_", names(features_filtered))

# Create all plots at once using lapply
plot_list <- lapply(seq_along(features_filtered), function(i) {
  feature_name <- names(features_filtered)[i]
  module_name <- paste0("Module.", feature_name)
  
  print(paste("Creating plot for:", module_name))
  
  plotEmbedding(proj,
    embedding = "UMAP_Harmony_LSI_Combined",
    colorBy = "cellColData",
    name = module_name
  )
})

# Name the plots for better organization
names(plot_list) <- names(features_filtered)

# Save all plots to PDF
print("Saving all plots to PDF")
do.call(plotPDF, c(plot_list, 
                   list(name = "mouse_multiome-scATAC-scRNA-Harmony_Combined_colorby-moduleScore_ALL.pdf", 
                        ArchRProj = proj,
                        addDOC = TRUE,
                        width = 5, height = 5)))

# Optional: Create a summary plot showing the number of genes per module
print("Creating summary of module composition")
module_summary <- data.frame(
  Module = names(features_filtered),
  Genes_Count = sapply(features_filtered, length),
  Missing_Count = sapply(missing_genes, length)
)
print(module_summary)

# Save summary to file
write.csv(module_summary, "module_gene_summary.csv", row.names = FALSE)

print("Script completed successfully!")
print(paste("Total plots created:", length(plot_list)))

EOF