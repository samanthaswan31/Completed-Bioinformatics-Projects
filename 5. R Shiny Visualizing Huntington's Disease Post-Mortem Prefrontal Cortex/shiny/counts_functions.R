# load required libraries
library(tidyverse)
library(ggplot2) 
library(dplyr) 
library('RColorBrewer')
library(gplots)
library(pheatmap)

# LOAD IN DATA

load_counts <- function(filepath) { 
  # validate that the file exists
  if (!file.exists(filepath)) {
    stop("File does not exist.")}
  
  # get file extension
  ext <- tools::file_ext(filepath)
  
  # switch based on file extension
  normcount <- switch(ext,
                 csv = vroom::vroom(filepath, delim = ","),
                 stop("Invalid file: only .csv files are supported."))
  
  # rename the first column to "gene"
  colnames(normcount)[1] <- "gene"
  
  return(normcount)
  
} # closes load_counts

# REACTIVE FILTERING FUNCTIONS

# genes > X percentile variance
# genes > X number of sample that are non-zero

filter_counts <- function(normcount, zero_threshold, percentile_threshold) {
  
  # filter by row sum of non-zero > zero_threshold
  
  zero_filt_gene <- normcount %>%
    rowwise() %>%
    # sum across all genes comparing to var_threshold
    filter(sum(c_across(-gene) != 0) >= zero_threshold) %>%
    ungroup()
  
  # filter by percentile variance > percentile_threshold
  
  # calculate gene_var, store in 'variance' col
  gene_var <- zero_filt_gene %>%
    mutate(variance = apply(dplyr::select(., -gene), 1, var))
  
  # calculate variance threshold from percentile_threshold input
  var_threshold <- quantile(gene_var$variance, probs = percentile_threshold / 100)
  
  #filter
  filtered_genes <- gene_var %>%
    filter(variance >= var_threshold)
  
  return(filtered_genes)
  
} # closes zero_var_filter

# GENERATE SUMMARY TABLE 

summarize_counts <- function(normcount, filtered_genes) {
  
  # generate needed statistics
  total_genes <- nrow(normcount)
  passing_genes <- nrow(filtered_genes)
  not_passing_genes <- total_genes - passing_genes
  passing_percent <- round((passing_genes / total_genes) * 100, 2)
  not_passing_percent <- round((not_passing_genes / total_genes) * 100, 2)
  
  # combine into one tibble
  summary_tibble <- tibble(
    Total_Genes = total_genes,
    Passing_Filter = passing_genes,
    Passing_Percentage = passing_percent,
    Not_Passing_Filter = not_passing_genes,
    Not_Passing_Percentage = not_passing_percent)
  
  return(summary_tibble)
  
} # closes summarize_counts

# GENERATE PLOTTING DATA

process_normplot <- function(normcount, filtered_genes) {
  
  # calculate median, variance, zeros and add to table
  norm_plot_data <- normcount %>%
    rowwise() %>%
    mutate(
      median = median(c_across(-gene)),
      variance = var(c_across(-gene)),
      zeros = sum(c_across(-gene) == 0)) %>%
    ungroup()
  
  # generate new 'Significance' column: if gene$alldata in gene$filtered data = significant
  norm_plot_data <- norm_plot_data %>%
    mutate(Significance = if_else(paste(gene) %in% paste(filtered_genes$gene),"Significant", "Not Significant"))

  return(norm_plot_data)  
} # closes process_normplot

# GENERATE MEDIAN V VAR DIAGNOSTIC SCATTER PLOT

med_var_plot <- function(norm_plot_data) {
  
  # graph, title and color informatively
  diag_plot1 <- norm_plot_data %>%
    ggplot(aes(x = median, y = variance, color = Significance))+
    geom_point()+
    scale_y_log10()+ # clearer curve
    theme_minimal()+
    scale_color_manual(values = c("Significant" = "royalblue", "Not Significant" = "gray"))+
    xlab("Median Normalized Count")+
    ylab("Variance (log10)")+
    labs(title = "Median Normalized Count vs Variance")
  
  return(diag_plot1)
  
} # closes med_var_plot

# GENERATE MEDIAN V ZERO PLOT

med_zero_plot <- function(norm_plot_data) {
  
  diag_plot2 <- norm_plot_data %>%
    ggplot(aes(x = median, y = zeros, color = Significance))+
    geom_point()+
    theme_minimal()+
    scale_color_manual(values = c("Significant" = "firebrick2", "Not Significant" = "gray"))+
    xlab("Median Normalized Count (log10)")+
    ylab("Zero Count")+
    scale_x_log10()+
    labs(title = "Median Normalized Count vs Zero Count")
  
  return(diag_plot2)
  
} # closes med_zero_plot

# PERFORM AND PLOT PRINCIPAL COMPONENT ANALYSIS PLOTS (PC1 V PC2)

pca_norm <- function(normcount, sample_data, pc_choiceone, pc_choicetwo, pca_color_variable) {
  
  # format data by selecting only samples and transposing
  format_data <- normcount %>%
    as.data.frame() %>%
    dplyr::select(-gene) %>%
    t() # want to perform PCA by sample
  
  # perform PCA
  pca_data <- prcomp(
    format_data,
    center = TRUE, # center data by subtracting mean
    scale = FALSE) # do not want to scale to get axis to match
  
  # generate the PC scores by sample
  pca_scores <- as.data.frame(pca_data$x)
  # create row "sample" with rownames to match for meta
  pca_scores$sample_title <- rownames(pca_scores)
  
  # dynamically select only the PCs of interest for graphing
  pca_scores <- pca_scores %>%
    dplyr::select(all_of(paste0("PC", pc_choiceone)), 
                  all_of(paste0("PC", pc_choicetwo)), 
                  sample_title)
  
  # rename meta sample title for joining
  meta <- sample_data %>%
    as.data.frame()
  
  # combine meta and PCA data by sample
  combined_data <- full_join(
    x = pca_scores,
    y = meta, 
    by = "sample_title")
  
  # for labels calculate variance explained by each PC
  variance_explained <- (pca_data$sdev^2 / sum(pca_data$sdev^2)) * 100
  pc1_variance <- round(variance_explained[pc_choiceone]) # PC input one
  pc2_variance <- round(variance_explained[pc_choicetwo]) # PC input two
  
  # check if the color variable is categorical or continuous
  if (is.factor(combined_data[[pca_color_variable]]) || is.character(combined_data[[pca_color_variable]])) {
    # If the variable is categorical, ensure it's treated as a factor
    pca_plot <- combined_data %>%
      ggplot(aes(x = !!sym(paste0("PC", pc_choiceone)), 
                 y = !!sym(paste0("PC", pc_choicetwo)), 
                 fill = as.factor(!!sym(pca_color_variable)))) +
      geom_point() +
      scale_fill_discrete() + # Discrete color scale for categorical variables
      labs(title = "PCA Plot",
           x = paste("PC", pc_choiceone, ":", pc1_variance, "% variance"),
           y = paste("PC", pc_choicetwo, ":", pc2_variance, "% variance"))
  } else {
    # If the variable is continuous, use a gradient
    pca_plot <- combined_data %>%
      ggplot(aes(x = !!sym(paste0("PC", pc_choiceone)), 
                 y = !!sym(paste0("PC", pc_choicetwo)), 
                 fill = !!sym(pca_color_variable))) +
      geom_point() +
      scale_fill_gradient(low = "blue", high = "red") + # Continuous color scale
      labs(title = "PCA Plot",
           x = paste("PC", pc_choiceone, ":", pc1_variance, "% variance"),
           y = paste("PC", pc_choicetwo, ":", pc2_variance, "% variance"))
  }
  
  return(pca_plot)
} # closes pca_norm

# PLOT HEATMAP

plot_heatmap <- function(filtered_genes,colorblind, num_colors) {
  
  # generate heatmap_data by cleaning, transforming to matrix, log transforming
  heatmap_data <- filtered_genes %>% dplyr::select(-gene, -variance) %>% {log(. + .1)} # .1 gives inf values a value for graphing
  
  # create specified custom color palette with user input
  color_palette <- brewer.pal(num_colors, colorblind)
  
  #create heatmap, two dynamic inputs
  heatmap <- pheatmap(heatmap_data, 
                      cluster_rows = FALSE,  # Do not cluster genes, focus on samples
                      cluster_cols = TRUE,   # Cluster samples
                      scale = "none",        # No scaling
                      show_rownames = TRUE,  # Show gene names
                      show_colnames = TRUE,  # Show sample names
                      legend = TRUE,
                      main = "Log Transformed Normalized Counts Post-Filter",
                      legend_title = "Log-Transformed Counts",
                      color = color_palette)
  
  return(heatmap)
  
} # closes plot_heatmap

# PLOT CORRELATION PLOT

plot_corrmap <- function(filtered_genes, num_colors) {
  
  # generate corr_data by cleaning, transforming to matrix, log transforming
  corr_data <- filtered_genes %>% dplyr::select(-gene, -variance) %>% {log(. + .1)} # .1 gives inf values a value for graphing
  
  # calculate correlations between samples
  cor_matrix <- cor(corr_data)
  
  # create the heatmap with clustering enabled
  corr_heatmap <- pheatmap(cor_matrix,
           cluster_rows = TRUE, # want to cluster samples
           cluster_cols = TRUE,
           scale = "none", # already log transformed values
           show_rownames = TRUE,
           show_colnames = TRUE,
           legend = TRUE,
           main = "Log Transformed Normalized Counts Correlation by Grouped Samples",
           legend_title = "Relative Correlation",
           fontsize_col = 6,
           fontsize_row = 4)
  
  return(corr_heatmap)
  
} # closes plot_corrmap