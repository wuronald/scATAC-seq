# motif bubble plot
# By: Ronald Wu
# last updated: Oct 15 2025

# load library
library(tidyverse)
library(here)

#' Motif Bubble Plot
#'
#' Creates a bubble plot visualizing transcription factor motif enrichment analysis results.
#' The size and color of bubbles represent the significance (FDR) of motif enrichment,
#' while the x-axis position represents the log2 fold enrichment.
#'
#' @param df A dataframe containing motif enrichment results. Must include columns:
#'   - `dataset`: Character column specifying the condition/dataset name
#'   - `TF`: Character column with transcription factor names
#'   - `log2(FE)`: Numeric column with log2 fold enrichment values
#'   - `q-value/FDR`: Numeric column with FDR-adjusted p-values
#'
#' @param cond A character string specifying which condition/dataset to plot.
#'   Must match a value in the `dataset` column of `df`.
#'
#' @param title A character string for the plot title. Default is "Motif Bubble Plot HOMER:".
#'   The condition name will be appended to this title.
#'
#' @return A ggplot object displaying the motif enrichment bubble plot.
#'   Significantly enriched motifs (FDR < 0.05) are highlighted in color.
#'   Bubble size is inversely proportional to FDR value (smaller p-values = larger bubbles).
#'
#' @details
#' The function performs the following steps:
#' \itemize{
#'   \item Filters the dataframe to the specified condition
#'   \item Counts and prints the number of significantly enriched motifs (FDR < 0.05)
#'   \item Calculates dynamic FDR breaks based on the data range
#'   \item Generates a ggplot with TFs ordered by log2 fold enrichment
#'   \item Color codes points based on statistical significance
#' }
#'
#' @examples
#' \dontrun{
#' MotifBubblePlot(motifs_df, cond = "knownmotif_N_BMP_edgeR")
#' MotifBubblePlot(motifs_df, cond = "knownmotif_N_VEH_edgeR", 
#'                 title = "ATAC-seq Motif Enrichment: ")
#' }
#'
#' @export
MotifBubblePlot <- function(df, cond, title ="Motif Bubble Plot HOMER:"){
  
  # Prefilter dataframe by condition
  df <- df %>% filter(dataset == cond)
  
  num_motifs <- df %>% pull(dataset) %>% length()
  print(paste("The number of motifs for: ", cond, "is", num_motifs))
  
  fdr_values <- df$`q-value/FDR`
  fdr_min <- min(fdr_values, na.rm = TRUE)
  fdr_max <- max(fdr_values, na.rm = TRUE)
  
  # Create transformed FDR column for sizing
  # Apply power transformation (square root) to increase visual separation
  # Add a small epsilon to avoid issues with zero values
  df <- df %>%
    mutate(fdr_transformed = (`q-value/FDR` + 1e-10) ^ 0.5)
  
  # Generate breaks at logarithmic intervals based on order of magnitude
  # This creates evenly-spaced breaks on a log scale, which is appropriate for
  # FDR values that often span multiple orders of magnitude
  if (fdr_min == 0) {
    # If minimum is 0, create breaks starting at 0 and then log-spaced for remaining values
    breaks <- c(0, 10^seq(ceiling(log10(fdr_max)) - 3, floor(log10(fdr_max))))
    breaks <- breaks[breaks <= fdr_max]
  } else {
    # Create log-spaced breaks across the full range of FDR values
    breaks <- 10^seq(floor(log10(fdr_min)), floor(log10(fdr_max)))
    breaks <- breaks[breaks >= fdr_min & breaks <= fdr_max]
  }
  
  # Fallback: if fewer than 3 breaks were generated, use linear spacing instead
  # This ensures the legend has sufficient granularity
  if (length(breaks) < 3) {
    breaks <- seq(fdr_min, fdr_max, length.out = 4)
  }
  
  # Create the bubble plot using transformed FDR values for sizing
  # Display original FDR values in the legend for interpretability
  g <- df %>% 
    ggplot(aes(x = reorder(TF, `log2(FE)`), y = `log2(FE)`, size = fdr_transformed, color = `q-value/FDR` < 0.05))+
    geom_point(stat = 'identity')+
    scale_size("FDR",
               trans = 'reverse',  # Reverse scale so smaller p-values = larger bubbles
               breaks = (breaks + 1e-10) ^ 0.5,  # Transform breaks to match data transformation
               labels = round(breaks, 4),  # Show original FDR values in legend
               range = c(1, 10))+  # Set bubble size range
    labs(title = paste(title, cond),
         x = "Transcription Factor Binding Motif",
         y = "log2 Fold Enrichment",
         color = "FDR < 0.05") +
    ylim(0, 2)+
    coord_flip()  # Flip coordinates for better TF name readability
  
  return(g)
}

