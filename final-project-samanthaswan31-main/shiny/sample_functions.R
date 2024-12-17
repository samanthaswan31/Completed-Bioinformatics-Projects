# load required libraries
library(tidyverse)
library(ggplot2)
library(dplyr)

# LOAD IN DATA, DISPLAY "SAMPLE"

# load and clean sample data
sample_summary <- function(filepath) {
  
  # read and clean age at death data
  # validate that the file exists
    if (!file.exists(filepath)) {
      stop("File does not exist.")}
  
  # get file extension
  ext <- tools::file_ext(filepath)
  
  # switch based on file extension
  series_matrix <- switch(ext,
                      csv = vroom::vroom(filepath, delim = ","),
                      stop("Invalid file: only .tsv or .csv files are supported."))
  
  # read and clean age at death data
  age_at_death <- series_matrix[44, ]
  age_at_death <- age_at_death %>%
    pivot_longer(
      cols = everything(),
      names_to = "sample",
      values_to = "age_at_death") %>%
    mutate(age_at_death = as.numeric(stringr::str_extract(age_at_death, "\\d+$"))) %>%
    dplyr::filter(!row_number() == 1)
  
  # read and clean diagnosis data
  diagnosis <- series_matrix[42, ]
  diagnosis <- diagnosis %>%
    pivot_longer(
      cols = everything(),
      names_to = "sample",
      values_to = "diagnosis") %>%
    mutate(diagnosis = as.character(stringr::str_extract(diagnosis, "[^:]+$"))) %>%
    dplyr::filter(!row_number() == 1)
  
  # read and clean geo accession data
  geo_accession <- series_matrix[33, ]
  geo_accession <- geo_accession %>%
    pivot_longer(
      cols = everything(),
      names_to = "sample",
      values_to = "geo_accession") %>%
    dplyr::filter(!row_number() == 1)
  
  # read and clean sample title data
  sample_title <- series_matrix[32, ]
  sample_title <- sample_title %>%
    pivot_longer(
      cols = everything(),
      names_to = "sample",
      values_to = "sample_title") %>%
    dplyr::filter(!row_number() == 1)
  
  # read and clean pmi data
  pmi <- series_matrix[43, ]
  pmi <- pmi %>%
    pivot_longer(
      cols = everything(),
      names_to = "sample",
      values_to = "pmi") %>%
    mutate(pmi = as.numeric(stringr::str_extract(pmi, "\\d+$"))) %>%
    dplyr::filter(!row_number() == 1)
  
  # read and clean rin data
  rin <- series_matrix[45, ]
  rin <- rin %>%
    pivot_longer(
      cols = everything(),
      names_to = "sample",
      values_to = "rin") %>%
    mutate(rin = as.numeric(stringr::str_extract(rin, "\\d+(\\.\\d+)?$"))) %>%
    dplyr::filter(!row_number() == 1)
  
  # read and clean total reads data
  total_reads <- series_matrix[46, ]
  total_reads <- total_reads %>%
    pivot_longer(
      cols = everything(),
      names_to = "sample",
      values_to = "total_reads") %>%
    mutate(total_reads = as.numeric(stringr::str_extract(total_reads, "\\d+$"))) %>%
    dplyr::filter(!row_number() == 1)
  
  # read and clean age of onset data
  onset_age <- series_matrix[47, ]
  onset_age <- onset_age %>%
    pivot_longer(
      cols = everything(),
      names_to = "sample",
      values_to = "onset_age") %>%
    mutate(onset_age = as.numeric(stringr::str_extract(onset_age, "\\d+$"))) %>%
    dplyr::filter(!row_number() == 1)
  
  # read and clean disease duration data
  duration <- series_matrix[48, ]
  duration <- duration %>%
    pivot_longer(
      cols = everything(),
      names_to = "sample",
      values_to = "duration") %>%
    mutate(duration = as.numeric(stringr::str_extract(duration, "\\d+$"))) %>%
    dplyr::filter(!row_number() == 1)
  
  # read and clean cag repeat count data
  cag <- series_matrix[49, ]
  cag <- cag %>%
    pivot_longer(
      cols = everything(),
      names_to = "sample",
      values_to = "cag") %>%
    mutate(cag = as.numeric(stringr::str_extract(cag, "\\d+$"))) %>%
    dplyr::filter(!row_number() == 1)
  
  # read and clean vonsattel grade data
  vonsattel_grade <- series_matrix[50, ]
  vonsattel_grade <- vonsattel_grade %>%
    pivot_longer(
      cols = everything(),
      names_to = "sample",
      values_to = "vonsattel_grade") %>%
    mutate(vonsattel_grade = as.numeric(stringr::str_extract(vonsattel_grade, "\\d+$"))) %>%
    dplyr::filter(!row_number() == 1)
  
  # read and clean h-v striatal score data
  hv_striatal_score <- series_matrix[51, ]
  hv_striatal_score <- hv_striatal_score %>%
    pivot_longer(
      cols = everything(),
      names_to = "sample",
      values_to = "hv_striatal_score") %>%
    mutate(hv_striatal_score = as.numeric(stringr::str_extract(hv_striatal_score, "\\d+$"))) %>%
    dplyr::filter(!row_number() == 1)
  
  # read and clean h-v cortical score data
  hv_cortical_score <- series_matrix[52, ]
  hv_cortical_score <- hv_cortical_score %>%
    pivot_longer(
      cols = everything(),
      names_to = "sample",
      values_to = "hv_cortical_score") %>%
    mutate(hv_cortical_score = as.numeric(stringr::str_extract(hv_cortical_score, "\\d+$"))) %>%
    dplyr::filter(!row_number() == 1)
  
  # combine all cleaned data into a final tibble
  final_tibble <- age_at_death %>%
    full_join(diagnosis, by = "sample") %>%
    full_join(geo_accession, by = "sample") %>%
    full_join(sample_title, by = "sample") %>%
    full_join(pmi, by = "sample") %>%
    full_join(rin, by = "sample") %>%
    full_join(total_reads, by = "sample") %>%
    full_join(onset_age, by = "sample") %>%
    full_join(duration, by = "sample") %>%
    full_join(cag, by = "sample") %>%
    full_join(vonsattel_grade, by = "sample") %>%
    full_join(hv_striatal_score, by = "sample") %>%
    full_join(hv_cortical_score, by = "sample") %>%
    dplyr::select(-sample)
  return(final_tibble)
} #close load_sample_data


# GENERATE THE SAMPLE , DISPLAY IN "SUMMARY"

dataset_summary <- function(final_tibble) {
  # generate column_types to store class of each column
  column_types <- sapply(final_tibble, class)
  
  # convert the 'column_types' into a data frame (to bind with summary_tibble)
  column_types_df <- data.frame(column_name = names(column_types), column_type = column_types)
  
  # pivot summary_tibble to long format
  summary_tibble <- final_tibble %>%
    mutate_all(as.character) %>%
    pivot_longer(cols = everything(),
                 names_to = "column_name",
                 values_to = "value")
  
  # convert 'value' column to numeric
  summary_tibble <- summary_tibble %>%
    mutate(value = as.numeric(value))
  
  # calculate the mean/sd for each column
  condensed_summary <- summary_tibble %>%
    group_by(column_name) %>%
    summarise(mean_value = mean(value, na.rm = TRUE),
              standard_dev = sd(value, na.rm = TRUE))
  
  # format the 'mean_value' and 'standard_dev' columns to prevent scientific notation
  condensed_summary$mean_value <- format(condensed_summary$mean_value, scientific = FALSE)
  condensed_summary$standard_dev <- format(condensed_summary$standard_dev, scientific = FALSE)
  
  # bind the 'column_types' to the condensed_summary
  condensed_summary <- condensed_summary %>%
    left_join(column_types_df, by = "column_name")
  
  # manually fix the character column values
  condensed_summary[condensed_summary$column_name == "diagnosis", "mean_value"] <- "Neurologically normal/Huntington's Disease"
  condensed_summary[condensed_summary$column_name == "diagnosis", "standard_dev"] <- "Neurologically normal/Huntington's Disease"
  
  condensed_summary[condensed_summary$column_name == "geo_accession", "mean_value"] <- "GSM158____"
  condensed_summary[condensed_summary$column_name == "geo_accession", "standard_dev"] <- "GSM158____"
  
  condensed_summary[condensed_summary$column_name == "sample_title", "mean_value"] <- "C_00__/H_____"
  condensed_summary[condensed_summary$column_name == "sample_title", "standard_dev"] <- "C_00__/H_____"
  
  return(condensed_summary)
}

# PLOT THE SAMPLE DATA, DISPLAY IN "GRAPH"

sample_plot <- function(final_tibble, x_name, y_name, plot_type) {
  
  # Generate a dynamic title based on the plot type, x_name, and y_name
  title <- paste(plot_type, "of", x_name, "grouped by", y_name)
  
  # Check which type of plot the user wants
  if (plot_type == "Histogram") {
    # Create a histogram for the selected variable (x_name)
    plot <- ggplot(data = final_tibble) +
      geom_histogram(aes(x = !!sym(x_name), fill = !!sym(y_name)), bins = 30, position = "dodge") +
      xlab(x_name) +
      ylab("Frequency") +
      ggtitle(title) +  # Add the dynamic title
      theme_minimal()
    
  } else if (plot_type == "Density") {
    # Create a density plot for the selected variable (x_name)
    plot <- ggplot(data = final_tibble) +
      geom_density(aes(x = !!sym(x_name), fill = !!sym(y_name)), alpha = 0.5) +
      xlab(x_name) +
      ylab("Density") +
      ggtitle(title) +  # Add the dynamic title
      theme_minimal()
    
  } else if (plot_type == "Violin") {
    # Create a violin plot for x_name grouped by y_name
    plot <- ggplot(data = final_tibble) +
      geom_violin(aes(x = !!sym(y_name), y = !!sym(x_name), fill = !!sym(y_name))) +
      xlab(y_name) +
      ylab(x_name) +
      ggtitle(title) +  # Add the dynamic title
      theme_minimal()
  }
  
  return(plot)
}