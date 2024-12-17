#import libraries
library(janitor)
library(tidyr)
library(naniar)
library(ggplot2)
library(dplyr)
library(tidyverse)
library('RColorBrewer')
library(gplots)
library("igraph")
library(reshape2)


#LOAD NORMALIZED COUNTS

load_nadata <- function(filepath) { 
  # validate that the file exists
  if (!file.exists(filepath)) {
    stop("File does not exist.")}
  
  # get file extension
  ext <- tools::file_ext(filepath)
  
  # switch based on file extension
  na_data <- switch(ext,
                       csv = vroom::vroom(filepath, delim = ","),
                       stop("Invalid file: only .csv files are supported."))
  
  colnames(na_data)[1] <- "GeneID"
  
  
  return(na_data)
  
} # closes load_deseq

#CALCULATE CORRELATION MATRIX

nacorr_matrix <- function(na_data, genes_list) {
  # genes_list is passed directly from the reactive input
  
  # Filter for the selected genes
  normdata <- na_data %>% dplyr::filter(GeneID %in% genes_list)
  
  # Create an empty list to store expression data for each gene
  expression_data <- list()
  
  # Extract expression values for each gene
  for (gene in genes_list) {
    gene_expr <- normdata %>% dplyr::filter(GeneID == gene) %>% dplyr::select(-GeneID) %>% unlist() %>% as.numeric()
    expression_data[[gene]] <- gene_expr
  }
  
  # Combine all gene expression data into a data frame (columns are the samples, rows are genes)
  gene_expression_df <- as.data.frame(expression_data)
  
  # Calculate the correlation matrix between the genes across samples
  cor_matrix <- cor(gene_expression_df, use = "pairwise.complete.obs")
  
  # Return the correlation matrix
  return(cor_matrix)
}

#PLOT HEATMAP

nacorr_heatmap <- function(corr_matrix) {
    
    
    #create heatmap, two dynamic inputs
    pheatmap(corr_matrix, 
                        cluster_rows = FALSE,  # Do not cluster genes, focus on samples
                        cluster_cols = FALSE,   # Cluster samples
                        scale = "none",        # No scaling
                        show_rownames = TRUE,  # Show gene names
                        show_colnames = TRUE,  # Show sample names
                        legend = TRUE,
                        main = "Log Transformed Normalized Counts Post-Filter",
                        legend_title = "Log-Transformed Counts")
    
    
}

#PLOT NETWORK ANALYSIS
# Function to visualize correlation network with customized vertex color and labels
visualize_correlation_network <- function(correlation_matrix) {
  # Extract gene names from the correlation matrix (row names are the gene names)
  gene_names <- rownames(correlation_matrix)
  
  # Create the graph from the correlation matrix (undirected, weighted)
  graph <- graph.adjacency(correlation_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  # Generate a color palette (one color for each gene)
  unique_colors <- rainbow(length(gene_names))
  
  # Assign each vertex a unique color
  V(graph)$color <- unique_colors
  
  # Label vertices with gene names
  V(graph)$label <- gene_names
  
  # Plot the network
  network_vis <- plot(graph, 
                      vertex.label = V(graph)$label, 
                      vertex.label.cex = 0.8,  # Adjust label size
                      vertex.label.color = "black", 
                      edge.width = 1, 
                      main = "Correlation Network")
  
  return(network_vis)
}

# Function to calculate network metrics for each gene
calculate_na_summary <- function(correlation_matrix) {
  
  # Convert all correlation values to absolute values
  correlation_matrix <- abs(correlation_matrix)
  
  # Create the graph from the correlation matrix
  graph <- graph.adjacency(correlation_matrix, mode = "undirected", weighted = TRUE, diag = FALSE)
  
  # Calculate degree for each vertex (gene)
  degree_values <- degree(graph, mode = "all")
  
  # Calculate closeness centrality for each vertex (gene)
  closeness_values <- closeness(graph, mode = "all")
  
  # Calculate betweenness centrality for each vertex (gene)
  betweenness_values <- betweenness(graph, directed = FALSE)
  
  # Combine all metrics into a data frame
  na_summary <- data.frame(
    GeneID = names(degree_values),
    Degree = degree_values,
    Closeness_Centrality = closeness_values,
    Betweenness_Centrality = betweenness_values
  )
  
  na_summary <- na_summary[, -1]
  
  # Return the summary table
  return(na_summary)
}
