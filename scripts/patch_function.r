library(rhdf5)

#' Patch the package's internal function to fix array issue
#' 
#' This directly modifies the .import10xToSE function in the package namespace
#' 
#' @param package_name Name of the package (default "ArchR")
patch_import10x_function <- function(package_name = "ArchR") {
  
  # Try to get the namespace
  if (!package_name %in% loadedNamespaces()) {
    stop("Package '", package_name, "' is not loaded. Please load it first with library(", package_name, ")")
  }
  
  ns <- asNamespace(package_name)
  
  # Get the original function
  if (!exists(".import10xToSE", envir = ns, inherits = FALSE)) {
    stop("Function .import10xToSE not found in package namespace")
  }
  
  original_func <- get(".import10xToSE", envir = ns)
  
  # Create new version with fix
  new_func <- function(h5 = NULL, type10x = NULL, name = NULL, ranges = NULL) {
    
    # Shape
    shape <- rhdf5::h5read(h5, "/matrix/shape")
    
    # Read features10x
    features10x <- rhdf5::h5read(h5, "/matrix/features")
    features10x <- lapply(seq_along(features10x), function(x){
      if(length(features10x[[x]]) == shape[1]){
        if(object.size(features10x[[x]]) > object.size(S4Vectors::Rle(features10x[[x]]))){
          df <- S4Vectors::DataFrame(x = S4Vectors::Rle(as.vector(features10x[[x]])))
        }else{
          df <- S4Vectors::DataFrame(x = as.vector(features10x[[x]]))
        }
        colnames(df) <- names(as.vector(features10x))[x]
        df
      }else{
        NULL
      }
    })
    features10x <- Reduce("cbind", features10x[!unlist(lapply(features10x, is.null))])
    
    # Determine Idx
    if(!is.null(type10x)){
      idx <- which(paste0(features10x$feature_type) %in% type10x)
    }else{
      idx <- seq_len(nrow(features10x))
    }
    if(length(idx)==0){
      stop(
        paste0(
          h5,
          "\nMissing `type10x`! Feature Types in h5:\n",
          "\t", paste0(unique(features10x$feature_type),collapse="; ")
        )
      )
    }
    
    # Subset
    features10x <- features10x[idx, , drop=FALSE]
    
    # Interval processing (keeping original logic)
    if("interval" %in% colnames(features10x)){
      idxNA <- which(features10x$interval=="NA")
      
      if(length(idxNA) > 0){
        if(!is.null(ranges) & is(ranges, "GRanges")){
          idx1 <- paste0(GenomicRanges::seqnames(ranges)) %in% c(1:22, "X", "Y", "MT")
          if(length(idx1) > 0){
            ranges2 <- GenomicRanges::GRanges(
              seqnames = ifelse(idx1, paste0("chr", GenomicRanges::seqnames(ranges)), paste0(GenomicRanges::seqnames(ranges))),
              ranges = GenomicRanges::ranges(ranges)
            )
            S4Vectors::mcols(ranges2) <- S4Vectors::mcols(ranges)
          }
          
          features10xNA <- features10x[which(features10x$interval=="NA"),,drop=FALSE]
          namesNA <- features10xNA$name
          idxFix <- match(namesNA, S4Vectors::mcols(ranges2)[, grep("name", colnames(S4Vectors::mcols(ranges)), ignore.case=TRUE)])
          if(length(idxFix[!is.na(idxFix)]) > 0){
            message("Correcting missing intervals...")
            idx2 <- which(!is.na(idxFix))
            rangesFix <- ranges2[idxFix[idx2]]
            BiocGenerics::strand(rangesFix) <- "*"
            features10xNA$interval[idx2] <- paste0(rangesFix)
            features10x[which(features10x$interval=="NA"), ] <- features10xNA
          }
          
          features10xNA <- features10x[which(features10x$interval=="NA"),,drop=FALSE]
          if(nrow(features10xNA) > 0){
            features10xNA$interval <- paste0("chrNA:1-1")
            features10x[which(features10x$interval=="NA"), ] <- features10xNA
          }
          
          features10x$ranges <- GenomicRanges::GRanges(paste0(features10x$interval))
          features10x$interval <- NULL
        }else{
          features10xNA <- features10x[which(features10x$interval=="NA"),,drop=FALSE]
          if(nrow(features10xNA) > 0){
            features10xNA$interval <- paste0("chrNA:1-1")
            features10x[which(features10x$interval=="NA"), ] <- features10xNA
          }
          features10x$ranges <- GenomicRanges::GRanges(paste0(features10x$interval))
          features10x$interval <- NULL
        }
      }else{
        features10x$ranges <- GenomicRanges::GRanges(paste0(features10x$interval))
        features10x$interval <- NULL
      }
    }
    
    # Read Matrix - WITH FIX FOR ARRAY ISSUE
    mat <- Matrix::sparseMatrix(
      i = as.vector(rhdf5::h5read(h5, "/matrix/indices")),
      p = as.vector(rhdf5::h5read(h5, "/matrix/indptr")),
      x = as.numeric(as.vector(rhdf5::h5read(h5, "/matrix/data"))),
      dims = shape,
      index1 = FALSE
    )
    
    barcodes <- rhdf5::h5read(h5, "/matrix/barcodes")
    if(!is.null(name)){
      colnames(mat) <- paste0(name, "#", barcodes)
    }else{
      colnames(mat) <- barcodes
    }
    
    # Subset
    mat <- mat[idx, , drop = FALSE]
    gc()
    
    # Summarized Experiment
    if("ranges" %in% colnames(features10x)){
      mat <- SummarizedExperiment::SummarizedExperiment(
        assays = S4Vectors::SimpleList(data = mat),
        rowRanges = features10x$ranges
      )
      SummarizedExperiment::rowData(mat) <- features10x
      SummarizedExperiment::rowData(mat)$ranges <- NULL
    }else{
      mat <- SummarizedExperiment::SummarizedExperiment(
        assays = S4Vectors::SimpleList(data = mat),
        rowData = features10x
      )
    }
    
    rownames(mat) <- SummarizedExperiment::rowData(mat)$name
    
    # Sort (using package's internal function if available)
    if(exists(".sortRSE", envir = ns, inherits = FALSE)){
      sortFunc <- get(".sortRSE", envir = ns)
      mat <- sortFunc(mat)
    }
    
    return(mat)
  }
  
  # Assign the new function to the namespace
  unlockBinding(".import10xToSE", ns)
  assign(".import10xToSE", new_func, envir = ns)
  lockBinding(".import10xToSE", ns)
  
  message("Successfully patched .import10xToSE function in ", package_name, " package!")
  message("You can now use import10xFeatureMatrix() normally.")
}

# Usage:
# 1. Load your package first
# library(ArchR)  # or whatever package contains import10xFeatureMatrix
#
# 2. Apply the patch
# patch_import10x_function("ArchR")
#
# 3. Use the function normally
# rse <- import10xFeatureMatrix(
#   input = "filtered_feature_bc_matrix.h5",
#   names = "sample1"
# )