#!/bin/bash
#SBATCH --job-name=run_01b_Gaiti_multiome_barcodeComparison
#SBATCH --output=run_01b_Gaiti_multiome_barcodeComparison_%j.out
#SBATCH --error=run_01b_Gaiti_multiome_barcodeComparison_%j.err
#SBATCH --partition=himem
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=12G
#SBATCH --time=01:00:00

# Load necessary modules (adjust as needed for your system)
module load R/4.4.1

# Run R script
Rscript - <<'EOF'
# load libraries
library(rhdf5)
library(utils)

# --- CONFIGURATION ---
files <- list(
  # Sample 1
  S1_Arrow = "/cluster/projects/wouterslab/ArchR103_4/GSM8820713_6467_A.arrow",
  S1_H5    = "/cluster/projects/wouterslab/ArchR103_4/data/Gaiti_multiome/GSM8820713_6467_A/6467_A_filtered_feature_bc_matrix.h5",
  S1_TSV   = "/cluster/projects/wouterslab/ArchR103_4/data/Gaiti_multiome/GSM8820713_6467_A/GSM8820713_6467_A_barcodes.tsv.gz",
  
  # Sample 2
  S2_Arrow = "/cluster/projects/wouterslab/ArchR103_4/GSM8820714_6467_B.arrow",
  S2_H5    = "/cluster/projects/wouterslab/ArchR103_4/data/Gaiti_multiome/GSM8820714_6467_B/6467_B_filtered_feature_bc_matrix.h5",
  S2_TSV   = "/cluster/projects/wouterslab/ArchR103_4/data/Gaiti_multiome/GSM8820714_6467_B/GSM8820714_6467_B_barcodes.tsv.gz"
)

output_dir <- "Gaiti_multiome_all_vs_all_comparison"
if (!dir.exists(output_dir)) dir.create(output_dir)

# --- STEP 1: LOAD AND NORMALIZE BARCODES ---

read_bc <- function(path, name) {
  if (!file.exists(path)) {
    warning(paste("File skipped (not found):", path))
    return(NULL)
  }
  
  message(paste("Loading:", name))
  
  # Determine type by extension
  if (grepl("\\.arrow$", path)) {
    # ArchR: Read Metadata, strip "Sample#" prefix
    raw <- h5read(path, "Metadata/CellNames")
    return(sub("^[^#]*#", "", raw))
    
  } else if (grepl("\\.h5$", path)) {
    # 10x H5: Read matrix/barcodes
    raw <- h5read(path, "matrix/barcodes")
    return(as.character(raw))
    
  } else if (grepl("tsv|txt", path)) {
    # 10x TSV: Read text file
    raw <- read.table(gzfile(path), header = FALSE, stringsAsFactors = FALSE)$V1
    return(raw)
  }
}

# Load all valid files into a list
bc_list <- list()
for (n in names(files)) {
  res <- read_bc(files[[n]], n)
  if (!is.null(res)) bc_list[[n]] <- res
}

# --- STEP 2: GENERATE ALL PAIRWISE COMPARISONS ---

# Generate all unique combinations of 2 files
pairs <- combn(names(bc_list), 2, simplify = FALSE)

message(paste("\nStarting", length(pairs), "pairwise comparisons...\n"))

summary_stats <- data.frame(
  Comparison = character(),
  Overlap = integer(),
  Unique_To_A = integer(),
  Unique_To_B = integer(),
  stringsAsFactors = FALSE
)

for (pair in pairs) {
  name_a <- pair[1]
  name_b <- pair[2]
  
  set_a <- bc_list[[name_a]]
  set_b <- bc_list[[name_b]]
  
  # Calculate sets
  overlap     <- intersect(set_a, set_b)
  unique_a    <- setdiff(set_a, set_b)
  unique_b    <- setdiff(set_b, set_a)
  
  # Log to console
  cat(sprintf("[%s vs %s] -> Overlap: %d | Unique %s: %d | Unique %s: %d\n", 
              name_a, name_b, length(overlap), name_a, length(unique_a), name_b, length(unique_b)))
  
  # Record stats
  summary_stats[nrow(summary_stats) + 1, ] <- list(
    paste(name_a, "vs", name_b, sep="_"),
    length(overlap),
    length(unique_a),
    length(unique_b)
  )
  
  # --- STEP 3: WRITE FILES ---
  # Construct a clean filename prefix
  file_prefix <- file.path(output_dir, paste0("COMPARE_", name_a, "_vs_", name_b))
  
  # Only write if data exists
  if(length(overlap) > 0)  writeLines(overlap,  paste0(file_prefix, "_OVERLAP.txt"))
  if(length(unique_a) > 0) writeLines(unique_a, paste0(file_prefix, "_UNIQUE_TO_", name_a, ".txt"))
  if(length(unique_b) > 0) writeLines(unique_b, paste0(file_prefix, "_UNIQUE_TO_", name_b, ".txt"))
}

# Save a master summary CSV
write.csv(summary_stats, file.path(output_dir, "summary_stats.csv"), row.names = FALSE)
message("\nDone! Summary CSV and individual text files saved.")

EOF

echo "barcode analysis complete"