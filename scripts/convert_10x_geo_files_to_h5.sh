#!/bin/bash
#SBATCH --job-name=convert_10x_geo_files_to_h5
#SBATCH --output=convert_10x_geo_files_to_h5_%j.out
#SBATCH --error=convert_10x_geo_files_to_h5_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=18
#SBATCH --mem=30G
#SBATCH --time=03:00:00

# Load necessary modules (adjust as needed for your system)
module load R/4.1.0

# Run R script
Rscript - <<'EOF'
library(DropletUtils)
library(Seurat)

# Define the base directory containing sample subdirectories
base_dir <- "geo_downloads"

# List of sample names
samples <- c(
    "GSM8820702_6419_A",
    "GSM8820703_6419_B",
    "GSM8820704_6419_C",
    "GSM8820705_6419_D",
    "GSM8820706_6425_A",
    "GSM8820707_6425_C",
    "GSM8820708_6425_D",
    "GSM8820709_6425_E",
    "GSM8820710_6425_F",
    "GSM8820711_6425_G",
    "GSM8820712_6434_A",
    "GSM8820713_6467_A",
    "GSM8820714_6467_B",
    "GSM8820715_6467_C",
    "GSM8820716_6467_D",
    "GSM8820717_6467_E",
    "GSM8820718_6467_F",
    "GSM8820719_6509_A",
    "GSM8820720_6509_B",
    "GSM8820721_6509_C",
    "GSM8820722_6509_D",
    "GSM8820723_6509_E",
    "GSM8820724_6514_A",
    "GSM8820725_6514_B",
    "GSM8820726_6514_C",
    "GSM8820727_6514_D"
)

# Counter for progress
total_samples <- length(samples)
current <- 0

cat("Starting processing of", total_samples, "samples...\n")
cat("Base directory:", base_dir, "\n\n")

# Loop through each sample
for (sample in samples) {
    current <- current + 1
    
    cat(sprintf("[%d/%d] Processing %s...\n", current, total_samples, sample))
    
    # Define input directory path
    input_dir <- file.path(base_dir, sample)
    
    # Check if directory exists
    if (!dir.exists(input_dir)) {
        cat("  ✗ Directory not found:", input_dir, "\n\n")
        next
    }
    
    # Check if required files exist with sample prefix
    original_files <- c(
        barcodes = file.path(input_dir, paste0(sample, "_barcodes.tsv.gz")),
        features = file.path(input_dir, paste0(sample, "_features.tsv.gz")),
        matrix = file.path(input_dir, paste0(sample, "_matrix.mtx.gz"))
    )
    
    missing_files <- !file.exists(original_files)
    if (any(missing_files)) {
        cat("  ✗ Missing files:\n")
        for (f in original_files[missing_files]) {
            cat("    -", basename(f), "\n")
        }
        cat("\n")
        next
    }
    
    # Create symbolic links with standard 10X names
    standard_files <- c(
        barcodes = file.path(input_dir, "barcodes.tsv.gz"),
        features = file.path(input_dir, "features.tsv.gz"),
        matrix = file.path(input_dir, "matrix.mtx.gz")
    )
    
    cat("  Creating symbolic links...\n")
    for (i in seq_along(original_files)) {
        # Remove existing symlink if it exists
        if (file.exists(standard_files[i])) {
            file.remove(standard_files[i])
        }
        # Create symlink
        file.symlink(basename(original_files[i]), standard_files[i])
    }
    
    # Try to read and process the matrix
    tryCatch({
        cat("  Reading 10X matrix...\n")
        filter_matrix <- Read10X(input_dir)
        
        # Check if Read10X returned a list (multimodal data)
        if (is.list(filter_matrix) && !is.data.frame(filter_matrix)) {
            cat("  Detected multimodal data with", length(filter_matrix), "modalities\n")
            # Extract names of modalities
            modality_names <- names(filter_matrix)
            cat("  Modalities found:", paste(modality_names, collapse=", "), "\n")
            
            # Process each modality separately
            for (modality in modality_names) {
                cat("  Processing modality:", modality, "\n")
                
                modal_matrix <- filter_matrix[[modality]]
                # Clean modality name for filename (replace spaces with underscores, make lowercase)
                clean_modality <- tolower(gsub(" ", "_", modality))
                output_h5 <- file.path(input_dir, paste0(sample, "_", clean_modality, "_filtered_feature_bc_matrix.h5"))
                
                cat("  Writing HDF5 file for", modality, "...\n")
                write10xCounts(
                    output_h5,
                    modal_matrix,
                    type = "HDF5",
                    genome = "mm10",
                    version = "3",
                    overwrite = TRUE,
                    gene.id = rownames(modal_matrix),
                    gene.symbol = rownames(modal_matrix)
                )
                
                cat("  ✓ Success! Output:", basename(output_h5), "\n")
            }
        } else {
            # Single modality data
            output_h5 <- file.path(input_dir, paste0(sample, "_filtered_feature_bc_matrix.h5"))
            
            cat("  Writing HDF5 file...\n")
            write10xCounts(
                output_h5,
                filter_matrix,
                type = "HDF5",
                genome = "mm10",
                version = "3",
                overwrite = TRUE,
                gene.id = rownames(filter_matrix),
                gene.symbol = rownames(filter_matrix)
            )
            
            cat("  ✓ Success! Output:", basename(output_h5), "\n")
        }
        
        # Clean up symbolic links
        cat("  Cleaning up symbolic links...\n")
        for (link_file in standard_files) {
            if (file.exists(link_file)) {
                file.remove(link_file)
            }
        }
        
    }, error = function(e) {
        cat("  ✗ Error processing sample:", conditionMessage(e), "\n")
        # Clean up symbolic links even on error
        for (link_file in standard_files) {
            if (file.exists(link_file)) {
                file.remove(link_file)
            }
        }
    })
    
    cat("\n")
}

cat("Processing complete!\n")
EOF