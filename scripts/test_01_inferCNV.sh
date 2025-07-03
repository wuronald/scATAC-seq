#!/bin/bash
#SBATCH --job-name=test_01_inferCNV
#SBATCH --output=test_01_inferCNV_%j.out
#SBATCH --error=test_01_inferCNV_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=12G
#SBATCH --time=01:00:00

# Load necessary modules (adjust as needed for your system)
module load R/4.1.0 # this version of R works with infercnv on HPC; other versions may not work
module load JAGS/4.3.0 # required for infercnv

mkdir -p infercnv_test/

# Run R script
Rscript - <<'EOF'

# load libraries
library(infercnv)
set.seed(1)

out_dir <- "infercnv_test/"
out_dir2 <- "infercnv_test2/"

print(paste("Output directory:", out_dir))
# Create the infercnv object using example data
# Note: The example data is provided in the infercnv package, but you can replace
# it with your own data by specifying the paths to your raw counts matrix,
# annotations file, and gene order file.
# The example data is a downsampled version of the oligodendroglioma dataset.
# For more information, see the infercnv documentation

print("Creating infercnv object...")
infercnv_obj = CreateInfercnvObject(raw_counts_matrix=system.file("extdata", "oligodendroglioma_expression_downsampled.counts.matrix.gz", package = "infercnv"),
                                    annotations_file=system.file("extdata", "oligodendroglioma_annotations_downsampled.txt", package = "infercnv"),
                                    delim="\t",
                                    gene_order_file=system.file("extdata", "gencode_downsampled.EXAMPLE_ONLY_DONT_REUSE.txt", package = "infercnv"),
                                    ref_group_names=c("Microglia/Macrophage","Oligodendrocytes (non-malignant)")) 
print("InferCNV object created successfully.")
# # Run inferCNV analysis
# # Note: The parameters used here are for demonstration purposes. You may need to adjust them based on your specific dataset and analysis requirements.
# # The cutoff parameter is set to 1 for Smart-seq2 data and 0.1 for 10x Genomics data.
# # The cluster_by_groups parameter is set to TRUE to cluster the cells by their annotations.
# # The denoise parameter is set to TRUE to apply denoising to the data.
# # The HMM parameter is set to TRUE to use a Hidden Markov Model for inference.


# print("Running inferCNV analysis...")
# start_time <- Sys.time()
# infercnv_obj_run = infercnv::run(infercnv_obj,
#                              cutoff=1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
#                              out_dir=out_dir, 
#                              cluster_by_groups=TRUE, 
#                              denoise=TRUE,
#                              HMM=TRUE)
# print("InferCNV analysis completed successfully.")

#time_taken <- Sys.time() - start_time
#print(paste("Total time taken for inferCNV analysis:", time_taken))

print("Running inferCNV analysis with parallel processing...")
start_time <- Sys.time()
infercnv_obj_parallel = infercnv::run(infercnv_obj,
                             cutoff=1, # cutoff=1 works well for Smart-seq2, and cutoff=0.1 works well for 10x Genomics
                             out_dir=out_dir2, 
                             cluster_by_groups=TRUE, 
                             denoise=TRUE,
                             num_threads=18,  # Set the number of threads to use for parallel processing
                             HMM=TRUE)
time_taken <- Sys.time() - start_time
print(paste("Total time taken for inferCNV analysis:", time_taken))

EOF
