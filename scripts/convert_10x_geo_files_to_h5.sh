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

cat("Converting 10X matrices to HDF5 format for ArchR...\n")
cat("Base directory:", base_dir, "\n\n")

# Vector to store HDF5 file paths for ArchR import
h5_files <- c()

# Loop through each sample
for (sample in samples[1]) {
    current <- current + 1
    
    cat(sprintf("[%d/%d] Processing %s...\n", current, total_samples, sample))
    
    # Define input directory path
    input_dir <- file.path(base_dir, sample)
    
    # Check if directory exists
    if (!dir.exists(input_dir)) {
        cat("  ✗ Directory not found:", input_dir, "\n\n")
        next
    }
    
    # Original files with sample prefix
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
    
    # Standard 10X names (create symlinks)
    standard_files <- c(
        barcodes = file.path(input_dir, "barcodes.tsv.gz"),
        features = file.path(input_dir, "features.tsv.gz"),
        matrix = file.path(input_dir, "matrix.mtx.gz")
    )
    
    cat("  Creating symbolic links...\n")
    for (i in seq_along(original_files)) {
        if (file.exists(standard_files[i])) {
            file.remove(standard_files[i])
        }
        file.symlink(basename(original_files[i]), standard_files[i])
    }
    
    # Try to read and convert the matrix
    tryCatch({
        cat("  Reading 10X matrix...\n")
        data_10x <- Read10X(input_dir)
        print(paste0("the 10x matrix is of class: ", class(data_10x)))

        # Check if it's multimodal and extract Gene Expression
        if (is.list(data_10x) && !is.data.frame(data_10x)) {
            cat("  Detected multimodal data\n")
            
            if ("Gene Expression" %in% names(data_10x)) {
                cat("  Extracting 'Gene Expression' modality...\n")
                gene_expr_matrix <- data_10x[["Gene Expression"]]
                cat("  Extracting 'Peaks' modality...\n")
                peaks_matrix <- data_10x[["Peaks"]]
            } else {
                cat("  ✗ 'Gene Expression' modality not found\n")
                cat("  Available modalities:", paste(names(data_10x), collapse=", "), "\n\n")
                next
            }
        } else {
            gene_expr_matrix <- data_10x
        }

        # Combine gene_expr and peaks matrices
        combined_matrix <- rbind(gene_expr_matrix, peaks_matrix)
        combined_ids <- c(
                        rownames(gene_expr_matrix),
                        rownames(peaks_matrix)
                        )

        combined_symbols <- c(
                        rownames(gene_expr_matrix), # Or a vector of gene symbols if you have them
                        rownames(peaks_matrix)
                        )
        feature_types <- c(
                            rep("Gene Expression", nrow(gene_expr_matrix)),
                            rep("Peaks", nrow(peaks_matrix))
                            )
        
        # Write HDF5 with gene.type specified as "Gene Expression"
        cat("  Writing HDF5 file with proper feature_type annotation...\n")
        output_h5 <- file.path(input_dir, "filtered_feature_bc_matrix.h5")
        
        # write10xCounts(
        #     path = output_h5,
        #     x = gene_expr_matrix,
        #     type = "HDF5",
        #     genome = "hg38",
        #     version = "3",
        #     overwrite = TRUE,
        #     gene.type = rep("Gene Expression", nrow(gene_expr_matrix)),
        #     gene.id = rownames(gene_expr_matrix),
        #     gene.symbol = rownames(gene_expr_matrix)
        # )
        write10xCounts(
        path = output_h5,
        x = combined_matrix,
        type = "HDF5",
        genome = "hg38",
        version = "3",
        overwrite = TRUE,
        gene.type = feature_types,
        gene.id = combined_ids,
        gene.symbol = combined_symbols
        )
        
        cat("  ✓ Success! Output:", output_h5, "\n")
        h5_files <- c(h5_files, output_h5)
        
        # Clean up symbolic links
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

cat("HDF5 conversion complete!\n")
cat(sprintf("Successfully created %d HDF5 files\n\n", length(h5_files)))

# Now import into ArchR
if (length(h5_files) > 0) {
    cat("Starting ArchR import...\n\n")
    
    # Name the files by sample
    names(h5_files) <- basename(dirname(h5_files))
    
    # Import the feature matrices using ArchR
    ArrowFiles <- import10xFeatureMatrix(
        input = h5_files,
        names = names(h5_files),
        featureType = "Gene Expression"
    )
    
    cat("\n✓ ArchR import complete!\n")
    cat("Arrow files created:\n")
    for (af in ArrowFiles) {
        cat("  -", af, "\n")
    }
} else {
    cat("✗ No HDF5 files were created successfully. Cannot import to ArchR.\n")
}

cat("\nDone!\n")
EOF