####################################################################
# Write 10x HDF5 Files with Genomic Intervals
####################################################################
# This script provides functions to write 10x-compatible HDF5 files
# with genomic interval information, ensuring RangedSummarizedExperiment
# creation upon import.

library(DropletUtils)
library(GenomicRanges)
library(rhdf5)

#' Write 10x counts with genomic intervals included
#' 
#' This wrapper ensures that genomic interval information is included
#' in the HDF5 file, which is required for creating RangedSummarizedExperiment
#' 
#' @param path Output file path (ending in .h5)
#' @param counts Sparse count matrix (genes x cells)
#' @param gene.id Gene IDs (e.g., Ensembl IDs)
#' @param gene.symbol Gene symbols (e.g., gene names)
#' @param barcodes Cell barcodes
#' @param genome Genome version (e.g., "hg38", "mm10")
#' @param gene.ranges GRanges object with genomic coordinates for genes (REQUIRED for intervals)
#' @param version 10x file format version (default "3")
#' @param overwrite Whether to overwrite existing file
#' 
#' @return Path to created file (invisibly)
#' 
#' @examples
#' # Create gene ranges from annotation
#' library(EnsDb.Hsapiens.v86)
#' gene_ranges <- genes(EnsDb.Hsapiens.v86)
#' 
#' # Write with intervals
#' write10xCounts_with_intervals(
#'   path = "output.h5",
#'   counts = my_counts,
#'   gene.id = rownames(my_counts),
#'   gene.symbol = gene_symbols,
#'   barcodes = colnames(my_counts),
#'   genome = "hg38",
#'   gene.ranges = gene_ranges
#' )
#' @export
write10xCounts_with_intervals <- function(
  path,
  counts,
  gene.id,
  gene.symbol,
  barcodes,
  genome = "hg38",
  gene.ranges = NULL,
  version = "3",
  overwrite = TRUE
) {
  
  if (is.null(gene.ranges)) {
    stop("gene.ranges is required to include genomic intervals in the HDF5 file.\n",
         "Provide a GRanges object with genomic coordinates for your genes.\n",
         "Example: genes(EnsDb.Hsapiens.v86) for human, genes(EnsDb.Mmusculus.v79) for mouse")
  }
  
  if (!is(gene.ranges, "GRanges")) {
    stop("gene.ranges must be a GRanges object")
  }
  
  # Write the basic file first
  write10xCounts(
    path = path,
    x = counts,
    gene.id = gene.id,
    gene.symbol = gene.symbol,
    barcodes = barcodes,
    genome = genome,
    version = version,
    type = "HDF5",
    overwrite = overwrite
  )
  
  message("Basic HDF5 file created, now adding genomic intervals...")
  
  # Match genes to ranges
  # Try multiple matching strategies
  gene_names <- gene.symbol
  if (is.null(gene_names)) {
    gene_names <- gene.id
  }
  
  # Try to match by symbol first
  if ("symbol" %in% colnames(mcols(gene.ranges))) {
    match_idx <- match(gene_names, gene.ranges$symbol)
  } else if ("gene_name" %in% colnames(mcols(gene.ranges))) {
    match_idx <- match(gene_names, gene.ranges$gene_name)
  } else if ("SYMBOL" %in% colnames(mcols(gene.ranges))) {
    match_idx <- match(gene_names, gene.ranges$SYMBOL)
  } else {
    # Try matching by names of GRanges
    match_idx <- match(gene_names, names(gene.ranges))
  }
  
  # If still not matched, try by gene ID
  if (all(is.na(match_idx))) {
    if ("gene_id" %in% colnames(mcols(gene.ranges))) {
      match_idx <- match(gene.id, gene.ranges$gene_id)
    } else {
      match_idx <- match(gene.id, names(gene.ranges))
    }
  }
  
  num_matched <- sum(!is.na(match_idx))
  message(sprintf("Matched %d/%d genes to genomic ranges (%.1f%%)", 
                  num_matched, length(gene.id), 100*num_matched/length(gene.id)))
  
  if (num_matched == 0) {
    stop("Could not match any genes to the provided ranges. ",
         "Check that gene names/IDs match between your count matrix and gene.ranges")
  }
  
  # Create interval strings
  intervals <- rep("NA", length(gene.id))
  for (i in seq_along(gene.id)) {
    if (!is.na(match_idx[i])) {
      gr <- gene.ranges[match_idx[i]]
      # Format: chr:start-end
      intervals[i] <- paste0(seqnames(gr), ":", start(gr), "-", end(gr))
    }
  }
  
  # Check how many are still NA
  num_na <- sum(intervals == "NA")
  if (num_na > 0) {
    warning(sprintf("%d genes have no genomic interval and will be assigned to chrNA:1-1", num_na))
  }
  
  # Add intervals to the HDF5 file
  message("Writing interval data to HDF5...")
  
  # The interval field should already exist from write10xCounts, but it may be empty
  # We need to overwrite it
  tryCatch({
    # Try to delete existing interval if it exists
    h5delete(path, "/matrix/features/interval")
  }, error = function(e) {
    # Interval doesn't exist yet, that's fine
  })
  
  # Write the interval data
  h5write(intervals, path, "/matrix/features/interval")
  
  H5close()
  
  message("Successfully added genomic intervals to HDF5 file!")
  message("File should now create a RangedSummarizedExperiment when imported.")
  
  invisible(path)
}


#' Helper function to get gene ranges for common genomes
#' 
#' @param genome Genome name: "hg38", "hg19", "mm10", "mm39"
#' @return GRanges object with gene annotations
#' 
#' @examples
#' gene_ranges <- get_gene_ranges("hg38")
#' @export
get_gene_ranges <- function(genome = c("hg38", "hg19", "mm10", "mm39")) {
  
  genome <- match.arg(genome)
  
  message("Loading gene annotations for ", genome, "...")
  
  if (genome == "hg38") {
    if (!requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE)) {
      stop("Package 'EnsDb.Hsapiens.v86' is required for hg38.\n",
           "Install with: BiocManager::install('EnsDb.Hsapiens.v86')")
    }
    library(EnsDb.Hsapiens.v86)
    gene_ranges <- genes(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
    
  } else if (genome == "hg19") {
    if (!requireNamespace("EnsDb.Hsapiens.v75", quietly = TRUE)) {
      stop("Package 'EnsDb.Hsapiens.v75' is required for hg19.\n",
           "Install with: BiocManager::install('EnsDb.Hsapiens.v75')")
    }
    library(EnsDb.Hsapiens.v75)
    gene_ranges <- genes(EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75)
    
  } else if (genome == "mm10") {
    if (!requireNamespace("EnsDb.Mmusculus.v79", quietly = TRUE)) {
      stop("Package 'EnsDb.Mmusculus.v79' is required for mm10.\n",
           "Install with: BiocManager::install('EnsDb.Mmusculus.v79')")
    }
    library(EnsDb.Mmusculus.v79)
    gene_ranges <- genes(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79)
    
  } else if (genome == "mm39") {
    if (!requireNamespace("EnsDb.Mmusculus.v79", quietly = TRUE)) {
      stop("Package 'EnsDb.Mmusculus.v79' is required for mm39.\n",
           "Install with: BiocManager::install('EnsDb.Mmusculus.v79')")
    }
    library(EnsDb.Mmusculus.v79)
    gene_ranges <- genes(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79)
  }
  
  message("Loaded ", length(gene_ranges), " gene ranges")
  return(gene_ranges)
}


#' Complete example workflow
#' @examples
#' \dontrun{
#' # 1. Get gene ranges for your genome
#' gene_ranges <- get_gene_ranges("hg38")
#' 
#' # 2. Write the file with intervals
#' write10xCounts_with_intervals(
#'   path = "output_with_intervals.h5",
#'   counts = my_counts,
#'   gene.id = gene_ids,
#'   gene.symbol = gene_symbols,
#'   barcodes = cell_barcodes,
#'   genome = "hg38",
#'   gene.ranges = gene_ranges,
#'   overwrite = TRUE
#' )
#' 
#' # 3. Verify it worked
#' library(rhdf5)
#' intervals <- h5read("output_with_intervals.h5", "/matrix/features/interval")
#' cat("Number of genes with valid intervals:", sum(intervals != "NA"), "\n")
#' }
example_workflow <- function() {
  message("See function documentation for example usage")
}

# Print usage message when sourced
message("
========================================
write10x_with_intervals.R loaded
========================================

Main function: write10xCounts_with_intervals()

Quick start:
  1. Get gene ranges:
     gene_ranges <- get_gene_ranges('hg38')
  
  2. Write counts with intervals:
     write10xCounts_with_intervals(
       path = 'output.h5',
       counts = my_counts,
       gene.id = gene_ids,
       gene.symbol = gene_symbols,
       barcodes = cell_barcodes,
       genome = 'hg38',
       gene.ranges = gene_ranges
     )

Available genomes: hg38, hg19, mm10, mm39
========================================
")
