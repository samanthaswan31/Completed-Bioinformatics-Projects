#import libraries
library('dplyr')
library('tidyr')
library('ggplot2')
library('SummarizedExperiment')
library('DESeq2')
library('biomaRt')
library('fgsea')

# LOAD IN DATA
load_deseq <- function(filepath) { 
  # validate that the file exists
  if (!file.exists(filepath)) {
    stop("File does not exist.")}
  
  # get file extension
  ext <- tools::file_ext(filepath)
  
  # switch based on file extension
  deseq_data <- switch(ext,
                      csv = vroom::vroom(filepath, delim = ","),
                      stop("Invalid file: only .csv files are supported."))
  
  # rename the first column to "gene"
  colnames(deseq_data)[1] <- "gene"
  
  return(deseq_data)
  
} # closes load_deseq

#GENERATE HISTOGRAM FOR PVALUE AND PADJ COUNTS
plot_p_histogram <- function(deseq_data, pchoice) {
  
  ggplot(deseq_data, aes_string(x = pchoice)) +
    geom_histogram(color="black", fill="darkolivegreen") +
    labs(title = paste("Distribution of",pchoice), 
         x = pchoice, y = "Count")
  
} # closes plot_p_histogram

#GENERATE VOLCANO PLOT FOR PVALUE AND PADJ COUNTS
plot_p_volcano <- function(deseq_data, pchoice, pfilter) {
  
  # Ensure pchoice is a valid column name
  if (!pchoice %in% colnames(deseq_data)) {
    stop("Invalid column name provided in `pchoice`. Please use 'padj' or 'pvalue'.")
  }
  
  # Create a new column to indicate significance
  deseq_data$significant <- deseq_data$padj < pfilter
  
  # Dynamically use the column specified by pchoice
  plot <- ggplot(deseq_data, 
                 aes(x = log2FoldChange, 
                     y = -log10(.data[[pchoice]]), 
                     color = significant)) +
    geom_point() +  # Add points for each data point
    labs(x = "log2 fold-change", 
         y = paste0("-log10(", pchoice, ")"), 
         title = "Volcano Plot Colored By Relative Significance") +
    scale_color_manual(
      values = c("FALSE" = "black", "TRUE" = "red"),  # Set colors for significance
      name = paste0("P-adj <", pfilter)  # Update legend title dynamically
    ) +
    theme_light() +  # Use a minimal theme
    #scale_y_log10() +  # Apply logarithmic scaling to y-axis
    theme(plot.title = element_text(hjust = 0.5))
  
  return(plot)
}