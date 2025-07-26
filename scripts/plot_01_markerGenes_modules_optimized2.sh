#!/bin/bash
#SBATCH --job-name=plot_01_markerGenes_modules_optimized2
#SBATCH --output=plot_01_markerGenes_modules_optimized2_%j.out
#SBATCH --error=plot_01_markerGenes_modules_optimized2_%j.err
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
proj <- loadArchRProject(path = "mouse_multiome_harmony_test")

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
  Glutamatergic_score = c("Celf2", "Ptprd", "Arpp21", "Sv2b", "Pcsk2", "Phactr1", "Kalrn", "Nrg1", "Ano3", "Satb2"),
  GABAergic_score  = c("Grip1", "Galntl6", "Gad1", "Gad2", "Dlx6os1", "Cntnap2", "Gm26905", "Kcnmb2", "Kcnip1", "Erbb4"),
  Astro_score = c("Gpc5", "Slc1a2", "Slc1a3", "Apoe", "Wdr17", "Plpp3", "Rorb", "Rmst", "Slc4a4", "Htra1"),
  Endo_score  = c("Flt1", "Slco1a4", "Adgrl4", "Ly6c1", "Slc2a1", "Klf2", "Mecom", "Bsg", "Ly6a", "Pltp"),
  Micro_PVM_score = c("Inpp5d", "Hexb", "Tgfbr1", "C1qa", "Ctss", "C1qb", "Zfhx3", "C1qc", "Selplg", "Cx3cr1"),
  Oligo_score  = c("Plp1", "Mbp", "St18", "Prr5l", "Mobp", "Mal", "Mog", "Cldn11", "Pde4b", "Mag"),
  OPC_score  = c("Lhfpl3", "Vcan", "Tnr", "Ptprz1", "Gm4876", "Xylt1", "Pdgfra", "Epn2", "Cacng4", "Megf11"),
  Peri_score = c("Atp13a5", "Vtn", "Cald1", "Ebf1", "Abcc9", "Dlc1", "Pdgfrb", "Tbx3os1", "Pde8b", "Slc38a11"),
  SMC_score = c("Acta2", "Map3k7cl", "Gm6249", "Myh11", "Tagln", "Pdlim3", "Ephx3", "Olfr558", "Crip1", "Tbx2"),
  Vlmc_score  =c("Ptgds", "Bnc2", "Cped1", "Slc7a11", "Bmp6", "Apod", "Mgp", "Eya2", "Ranbp3l", "Adam12"),
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

# Additional diagnostics before adding module scores
print("Performing additional diagnostics...")

# Check gene expression matrix dimensions and gene names
gene_exp_matrix <- getMatrixFromProject(proj, useMatrix = "GeneExpressionMatrix")
print(paste("Gene expression matrix dimensions:", paste(dim(gene_exp_matrix), collapse=" x ")))

# Get the actual gene names from the matrix
matrix_genes <- rowData(gene_exp_matrix)$name
print(paste("Number of genes in matrix:", length(matrix_genes)))

# Check if matrix genes match our gene symbols
print("Checking gene name format in matrix...")
print("First 10 genes from matrix:")
print(head(matrix_genes, 10))
print("First 10 genes from getGenes:")
print(head(gene_symbols, 10))

# Add module scores one by one with error handling
print("Adding module scores individually...")
successful_modules <- c()
failed_modules <- c()

for(i in seq_along(features_filtered)) {
  feature_name <- names(features_filtered)[i]
  feature_genes <- features_filtered[[i]]
  
  print(paste("Adding module score for:", feature_name, "with", length(feature_genes), "genes"))
  print(paste("Genes:", paste(feature_genes, collapse=", ")))
  
  # Double-check that genes exist in the matrix
  genes_in_matrix <- feature_genes %in% matrix_genes
  print(paste("Genes found in matrix:", sum(genes_in_matrix), "out of", length(feature_genes)))
  
  # Show which specific genes are missing
  missing_from_matrix <- feature_genes[!genes_in_matrix]
  if(length(missing_from_matrix) > 0) {
    print(paste("Genes missing from matrix:", paste(missing_from_matrix, collapse=", ")))
  }
  
  if(sum(genes_in_matrix) < 3) {
    print(paste("Skipping", feature_name, "- insufficient genes in matrix"))
    failed_modules <- c(failed_modules, feature_name)
    next
  }
  
  # Use only genes that are actually in the matrix
  valid_genes <- feature_genes[genes_in_matrix]
  print(paste("Using genes:", paste(valid_genes, collapse=", ")))
  
  # Additional check for problematic gene names (especially for Choroid_plexus)
  if(feature_name == "Choroid_plexus") {
    print("=== SPECIAL DEBUGGING FOR CHOROID_PLEXUS ===")
    
    # Check each gene individually with safer approach
    for(gene in valid_genes) {
      gene_idx <- which(matrix_genes == gene)
      print(paste("Gene:", gene, "- Found in matrix:", length(gene_idx) > 0))
      
      if(length(gene_idx) > 0) {
        print(paste("  Gene index:", gene_idx[1]))
      }
    }
    
    # Try with a smaller subset of genes first
    print("Trying with subset of genes...")
    safe_genes <- valid_genes[1:min(3, length(valid_genes))]  # Even smaller subset
    print(paste("Safe genes:", paste(safe_genes, collapse=", ")))
    
    temp_features_safe <- list()
    temp_features_safe[[paste0(feature_name, "_safe")]] <- safe_genes
    
    tryCatch({
      proj <- addModuleScore(proj,
          useMatrix = "GeneExpressionMatrix",
          name = "Module",
          features = temp_features_safe,
          nBin = 5,   # Very small bin size
          nBgd = 25,  # Small background
          seed = 1
      )
      print("✓ Safe subset worked! Issue might be with specific genes or parameters.")
      successful_modules <- c(successful_modules, paste0(feature_name, "_safe"))
      
    }, error = function(e) {
      print(paste("Even safe subset failed:", e$message))
      
      # Try with just one gene
      if(length(valid_genes) > 0) {
        print("Trying with single gene...")
        single_gene_features <- list()
        single_gene_features[[paste0(feature_name, "_single")]] <- valid_genes[1]
        
        tryCatch({
          proj <- addModuleScore(proj,
              useMatrix = "GeneExpressionMatrix",
              name = paste0("Module_", feature_name, "_single"),
              features = single_gene_features,
              nBin = 5,
              nBgd = 25,
              seed = 1
              )
          
          print("✓ Single gene worked!")
          successful_modules <- c(successful_modules, paste0(feature_name, "_single"))
          
        }, error = function(e2) {
          print(paste("Single gene also failed:", e2$message))
        })
      }
    })
    
    print("=== END SPECIAL DEBUGGING ===")
    next  # Skip the normal processing for Choroid_plexus
  }
  
  # Create a temporary feature list with just this one feature
  temp_features <- list()
  temp_features[[feature_name]] <- valid_genes
  
  # Add module score with error handling and more detailed error info
  tryCatch({
    proj <- addModuleScore(proj,
        useMatrix = "GeneExpressionMatrix",
        name = "Module",
        features = temp_features,
        nBin = 25,
        nBgd = 100,
        seed = 1)
    
    successful_modules <- c(successful_modules, feature_name)
    print(paste("✓ Successfully added module score for:", feature_name))
    
  }, error = function(e) {
    failed_modules <- c(failed_modules, feature_name)
    print(paste("✗ Failed to add module score for:", feature_name))
    print(paste("Full error message:", e$message))
    print(paste("Error call:", deparse(e$call)))
    
    # Try alternative approach with different parameters
    print("Trying alternative approach with reduced parameters...")
    tryCatch({
      proj <- addModuleScore(proj,
          useMatrix = "GeneExpressionMatrix",
          name = "Module",
          features = temp_features,
          nBin = 10,
          nBgd = 50,
          seed = 1
          )
      
      successful_modules <- c(successful_modules, feature_name)
      print(paste("✓ Successfully added module score for:", feature_name, "(with alternative parameters)"))
      
    }, error = function(e2) {
      print(paste("Alternative approach also failed:", e2$message))
    })
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
module_cols <- cell_col_names[grepl("Module", cell_col_names)]
print(module_cols)

# Initialize plot_list outside the conditional block
plot_list <- list()
# Select Embedding for plotting
# embedding_name <- "UMAP_Harmony_LSI_Combined"
embedding_name <- "UMAP_Harmony_LSI_Combined_default"
print(paste("Using embedding for plot:", embedding_name))

# Only proceed with plotting if we have successful modules
if(length(successful_modules) > 0) {
  
  # Create plots for successful module scores
  print("Creating plots for successful module scores")
  
  for(feature_name in successful_modules) {
    module_name <- paste0("Module.", feature_name)
    
    print(paste("Creating plot for:", module_name))
    
    tryCatch({
      p <- plotEmbedding(proj,
        embedding = embedding_name,
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