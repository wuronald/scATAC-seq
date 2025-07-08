#!/bin/bash
#SBATCH --job-name=plot_01_markerGenes_modules
#SBATCH --output=plot_01_markerGenes_modules_%j.out
#SBATCH --error=plot_01_markerGenes_modules_%j.err
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

# Select marker genes for B and T cells
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

# check if all features are in gene_symbols
print("Checking if all features are in gene symbols")
lapply(features, function(x) x %in% gene_symbols)

# add module scores for each feature set
proj <- addModuleScore(proj,
    useMatrix = "GeneExpressionMatrix",
    name = "Module",
    features = features)

# plot
print("Plotting module scores")
p1 <- plotEmbedding(proj,
    embedding = "UMAP_Harmony_LSI_Combined",
    colorBy = "cellColData",
    name="Module.Glutamatergic_score"
)

p2 <- plotEmbedding(proj,
    embedding = "UMAP_Harmony_LSI_Combined",
    colorBy = "cellColData",
    name="Module.GABAergic_score"
)

p3 <- plotEmbedding(proj,
    embedding = "UMAP_Harmony_LSI_Combined",
    colorBy = "cellColData",
    name="Module.Micro_PVM_score"
)

plotPDF(p1,p2,p3,
        name = "mouse_multiome-scATAC-scRNA-Harmony_Combined_colorby-moduleScore.pdf", 
        ArchRProj = proj, # PDF saved in the plots subfolder in the associated project directory
        addDOC = TRUE, # adds date of creation to end of filename
        width = 5, height = 5)

EOF