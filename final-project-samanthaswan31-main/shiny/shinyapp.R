# load libraries
library(shiny)
library(bslib)
library(DT)
library(tidyverse)
library(ggplot2)
library(dplyr)
library('RColorBrewer')
library(gplots)
library(pheatmap)
library("igraph")
library(reshape2)

#increase upload limit to 30 MB (for norm counts matrix)
options(shiny.maxRequestSize = 30*1024^2)

# source scripts
source("sample_functions.R") # samples
source("counts_functions.R") #counts
source("deseq_functions.R")# differential expression
source("networkanalysis_functions.R") #network analysis

# user interface
ui <- page_navbar(
  # title navbar
  title = "Analysis of mRNA-Seq Expression in Brain Tissue in Huntington's Disease",
  
  # use navset_card_pill layout, set each tab
  navset_card_pill(
    
    # SAMPLES
    nav_panel("Samples", 
              sidebarLayout(
                sidebarPanel(
                  p("The Samples tab is designed to facilitate an intuitive and comprehensive exploration of the dataset by providing users with multiple avenues for interaction. Here's a deeper dive into the functionality of each section:"),
                  br(),
                  p("Sample Data - The Sample Data tab allows users to upload and inspect sample-level metadata. After selecting a file, the data will be displayed in a sortable and filterable table. This provides a clear overview of the sample attributes, such as age, diagnosis, and other clinical or experimental variables, helping users quickly understand the structure of their dataset."),
                  br(),
                  p("Dataset Summary - The Dataset Summary tab provides an aggregated view of the dataset, including key statistics such as the number of samples, variable distributions, and potential data inconsistencies or missing values. This tab offers users a high-level summary, assisting in the data cleaning and quality assessment process before diving deeper into analysis or visualization."),
                  br(),
                  p("Graph - The Graph tab serves as a dynamic data visualization tool, where users can explore relationships between different variables in the dataset. By selecting variables for the x and y axes, users can create various types of plots including Histograms, Density plots, and Violin plots. This interactive graphing feature empowers users to examine the distribution and relationships of variables, facilitating more insightful data interpretation. The customization options—such as grouping by categorical variables like diagnosis or sample_title—allow for tailored visualizations to highlight patterns in the data, aiding in hypothesis generation or trend identification."),
                  br(),
                  p("Together, these sections of the Samples tab provide a seamless workflow for data inspection, summarization, and visualization, giving users the tools they need to better understand their data, identify trends or anomalies, and communicate findings effectively. This tab ensures that users have a well-rounded understanding of their sample data, from raw metadata to insightful visual representation."),
                  br(),
                ), # closes general sidebarPanel for samples
                
                mainPanel( 
                  # add tabs for "Sample Data", "Dataset Summary", "Graph"
                  tabsetPanel(
                    
                    # sample data tab with sidebar for metadata input
                    tabPanel("Sample Data",
                             sidebarLayout(
                               sidebarPanel(
                                 # add metadata input
                                 fileInput("select_sample", "Load Sample Data:")
                               ), # closes 'Sample Data' sidebarPanel
                               mainPanel(
                                 DTOutput("sample_table") # sample data table
                               ) # closes 'Sample Data' mainPanel
                             ) # closes 'Sample Data' sidebarLayout
                    ), # closes 'Sample Data' tabPanel
                    
                    # dataset summary tab with its own sidebar
                    tabPanel("Dataset Summary",
                             sidebarLayout(
                               sidebarPanel(
                                 # placeholder for dataset summary controls (if needed)
                               ), # closes 'Dataset Summary' sidebarPanel
                               mainPanel(
                                 DTOutput("summary_tibble") # dataset summary table
                               ) # closes 'Dataset Summary' mainPanel
                             ) # closes 'Dataset Summary' sidebarLayout
                    ), # closes 'Dataset Summary' tabPanel
                    
                    # graph tab with its own sidebar for plot controls
                    tabPanel("Graph",
                             sidebarLayout(
                               sidebarPanel(
                                 # create radiobutton for x_name
                                 radioButtons("x_name", "Choose the Variable to Plot:",
                                              c("age_at_death", "pmi", "rin", "total_reads", "hv_cortical_score", "hv_striatal_score", "vonsattel_grade", "cag", "duration", "onset_age")),
                                 
                                 # create radiobutton for y_name
                                 radioButtons("y_name", "Choose the Variable to Group By:",
                                              c("sample_title", "geo_accession", "diagnosis")),
                                 
                                 # choose plot type
                                 selectInput("plot_type", "Choose Plot Type",
                                             choices = c("Histogram", "Density", "Violin"))
                               ), # closes 'Graph' sidebarPanel
                               mainPanel(
                                 plotOutput("sample_graph") # graph output
                               ) # closes 'Graph' mainPanel
                             ) # closes 'Graph' sidebarLayout
                    ) # closes 'Graph' tabPanel
                    
                  ) # closes 'samples' tabsetPanel
                ) # closes 'samples' mainPanel
              ) # closes 'samples' sidebarLayout
    ), # closes nav_panel("Samples")
    
    # COUNTS 
    
    nav_panel("Counts",
              sidebarLayout(
                sidebarPanel(
                  p("The Counts tab offers tools to explore and visualize normalized count data, providing users with a comprehensive view of the FPKM normalized dataset and allowing for deeper analysis. Here's an overview of the functionality of each section:"),
                  br(),
                  p("Normalized Counts Data - The Normalized Counts Data tab allows users to upload and view the FPKM normalized count data. After selecting a file, the data will be displayed in a table format, providing an organized view of the count data for each sample and gene. This enables users to evaluate the distribution of counts across their samples and identify any normalization issues or trends in gene expression."),
                  br(),
                  p("Diagnostic Scatter Plots - The Diagnostic Scatter Plots tab offers a way to visualize data quality and distribution. Users can adjust thresholds for minimum non-zero counts and percentile variance, allowing for interactive exploration of data integrity. This section includes scatter plots that display the relationship between median expression, variance, and the presence of zero counts, helping users assess the quality of their data and identify potential outliers."),
                  br(),
                  p("Together, these sections of the Counts tab facilitate an in-depth examination of the normalized count data, enabling users to identify trends, assess data quality, and ensure that the data is appropriately prepared for further analysis."),
                  br(),
                ), # closes general sidebarPanel for counts
                mainPanel(
                  # Add tabs for "Samples", "Diagnostic Scatter Plots", "Heatmap", "PCA Projection Plot"
                  tabsetPanel(
                    
                    # Counts tab with sidebar for the file input
                    tabPanel("Normalized Counts Data",
                             sidebarLayout(
                               sidebarPanel(
                                 # Add normalized counts input
                                 fileInput("select_normcount", "Load Normalized Counts Data:")
                               ), # closes 'Normalized Counts Data' sidebarPanel
                               mainPanel(
                                 fluidPage(
                                   DTOutput("norm_filtered_table"), # normalized count table
                                   DTOutput("norm_count_summary") # reactive summary text
                                 ) # closes 'Normalized Counts Data' mainPanel
                               ) # closes 'Normalized Counts Data' sidebarLayout
                             ) # closes 'Normalized Counts Data' tabPanel
                    ), # closes 'Normalized Counts Data' tabPanel
                    
                    # diagnostic scatterplots tab with its own sidebar
                    tabPanel("Diagnostic Scatter Plots",
                             sidebarLayout(
                               sidebarPanel(
                                 # zero threshold slider
                                 sliderInput("zero_threshold", "Select the minimum number of samples that are non-zero:",
                                             min = 0, max = 69, value = 35),
                                 
                                 # percentile threshold slider
                                 sliderInput("percentile_threshold", "Select the minimum percentile of variance:",
                                             min = 0, max = 100, value = 50)
                               ), # closes 'Diagnostic Scatter Plots' sidebarPanel
                               mainPanel(
                                 fluidPage(
                                   plotOutput("diagplot_mv"), # median versus variance plot
                                   plotOutput("diagplot_mz")  # median versus zero count plot
                                 ) # closes 'Diagnostic Scatter Plots' mainPanel
                               ) # closes 'Diagnostic Scatter Plots' sidebarLayout
                             ) # closes 'Diagnostic Scatter Plots' tabPanel
                    ), # closes 'Diagnostic Scatter Plots' tabPanel
                    
                    # heatmap tab with its own sidebar
                    tabPanel("Heatmaps",
                             sidebarLayout(
                               sidebarPanel(
                                 # colorblind input
                                 selectInput("colorblind", "Select Color Mode:",
                                              choices = list("Colorblind" = "Set1", "Non-colorblind" = "RdYlBu"),
                                              selected = "RdYlBu"), # default = non / continuous color palette
                                 
                                 # number of colors input
                                 sliderInput("num_colors", "Select number of colors to plot with:",
                                             min = 1, max = 9, value = 4)
                               ), # closes 'Heatmaps' sidebarPanel
                               mainPanel(
                                 fluidPage(
                                   plotOutput("heat_map"),
                                   plotOutput("correlation_heatmap")
                                 ) # closes 'Counts' fluidPage
                               ) # closes 'Heatmaps' mainPanel
                             ) # closes 'Heatmaps' sidebarLayout
                    ), # closes 'Heatmaps' tabPanel
                    
                    # PCA Projection Plot tab with its own sidebar
                    tabPanel("PCA Projection Plot",
                             sidebarLayout(
                               sidebarPanel(
                                 
                                 # Add slider inputs for PCA-specific controls here
                                 sliderInput("pc_choiceone", "Select a PC for the x-axis:",
                                             min = 1, max = 10, value = 1),  # default = PC1
                                 
                                 sliderInput("pc_choicetwo", "Select a PC for the y-axis:",
                                             min = 1, max = 10, value = 2),   # default = PC2
                                 
                                 radioButtons("pca_color_variable", "Choose the Variable to Color By:",
                                              c("age_at_death", "pmi", "rin", "total_reads", "hv_cortical_score", "hv_striatal_score", "vonsattel_grade", "cag", "duration", "onset_age","sample_title", "geo_accession", "diagnosis"))
                                 
                               ), # closes 'PCA Projection Plot' sidebarPanel
                               
                               mainPanel(
                                 plotOutput("pca_plot") # PCA projection plot
                               ) # closes 'PCA Projection Plot' mainPanel
                             ) # closes 'PCA Projection Plot' sidebarLayout
                    ) # closes 'PCA Projection Plot' tabPanel
                    
                  ) # closes 'counts' tabsetPanel
                ) # closes 'counts' mainPanel
              ) # closes 'counts' sidebarLayout
    ), # closes nav_panel("Counts")
    
    # DIFFERENTIAL EXPRESSION
    nav_panel("Differential Expression", 
              sidebarLayout(
                sidebarPanel(
                  p("The Differential Expression tab enables users to explore, analyze, and visualize the results of differential expression analysis. This section provides tools to view the results, create plots, and refine statistical thresholds. Here's a breakdown of the functionality in each section:"),
                  br(),
                  p("Differential Expression Data - The Differential Expression Data tab allows users to upload and view the results of differential expression analysis. After selecting a file, the DE results will be displayed in a table format, giving users a comprehensive view of gene expression changes between conditions. This allows for a closer inspection of p-values, adjusted p-values, log fold changes, and other relevant statistics for each gene, helping users identify significant expression differences."),
                  br(),
                  p("Differential Expression Graphs - The Differential Expression Graphs tab provides several visualization options to explore the results of the differential expression analysis. Users can choose between plotting raw p-values or adjusted p-values (padj) and filter the results based on a significance threshold. This section includes both histograms, which provide a visual summary of the distribution of p-values, and volcano plots, which allow users to visualize the relationship between statistical significance and effect size for each gene. These visualizations help in identifying patterns, trends, and outliers in the data."),
                  br(),
                  p("Together, these sections of the Differential Expression tab allow users to effectively examine and interpret the results of their differential expression analysis, empowering them to identify key genes, trends, and potential biological insights in their dataset."),
                  br(),
                ),#closes DE general sidebarPanel
                mainPanel(
                  # tabs for DE data, plots
                  tabsetPanel(
                    
                  # DE data tab with sidebar for file input
                  tabPanel("Differential Expression Data",
                           sidebarLayout(
                             sidebarPanel(
                               #add de data input
                               fileInput("select_deseqdata", "Load Differential Expression Counts Data:")
                             ),#closes DE data sidebarPanel
                             mainPanel(
                               fluidPage(
                                 DTOutput("deseq_table")
                               )#closes DE data fluidPage
                             )#closes DE data mainPanel
                           )#closes DE data sidebarLayout
                    ),#closes DE data tabPanel
                  
                  #DE plot tab 
                  tabPanel("Differential Expression Graphs",
                           sidebarLayout(
                             sidebarPanel(
                               radioButtons("p_histogram_choice", "Choose the P-statistic to Plot:",
                                            c("pvalue", "padj")),
                               sliderInput("pfilter_de", "Select the p-statistic Significance Level",
                                           min = 0, max = 1, value = .1)
                             ),#closes DE graph sidebarPanel
                             mainPanel(
                               fluidPage(
                                 plotOutput("p_histogram"),#padj and pvalue histogram count plots
                                 plotOutput("p_volcano") #padj and pvalue volcano count plots
                               )#closes DE graph fluidPage
                             )#closes DE graph mainPanel
                           )#closes DE graph sidebarLayout
                    
                  )#closes DE graph tabPanel
                )#closes DE tabsetPanel
                )#closes DE mainPanel
              )#closes DE sidebarLayout
              ), # closes nav_panel("Differential Expression")
    
    # NETWORK ANALYSIS
    nav_panel("Network Analysis", 
              sidebarLayout(
                sidebarPanel(
                  p("The Network Analysis tab provides users with tools to explore and visualize relationships between genes or other entities in the dataset. This section allows for data exploration, network visualization, and analysis of network metrics. Here's an overview of each section:"),
                                br(),
                                p("Network Analysis Data - The Network Analysis Data tab allows users to upload and view the normalized counts data that will be used for network analysis. After selecting a file, the data will be displayed in a table, allowing users to inspect the network data, including gene expression values or other relevant metrics. This provides the foundation for building the network and performing further analyses."),
                                br(),
                                p("Heatmap - The Heatmap tab enables users to visualize the relationships between genes of interest by creating a heatmap. Users can enter a list of gene names, and the system will generate a heatmap that highlights the expression patterns of those genes across the dataset. This visualization is especially useful for identifying clusters of genes with similar expression profiles, which could indicate functional or regulatory associations."),
                                br(),
                                p("Network Visualization - The Network Visualization tab allows users to explore the structure of the gene network by plotting the relationships between genes based on correlation. Users can set a minimum correlation threshold to determine which edges (connections) are drawn between genes. This interactive network visualization helps in identifying clusters of genes that are highly correlated and may share biological functions or pathways."),
                                br(),
                                p("Network Metrics - The Network Metrics tab provides quantitative measures of the network's structure and properties. Users can view various metrics that describe the network's connectivity, centrality, and other characteristics. These metrics are essential for understanding the organization and behavior of the network, as well as identifying key genes or hubs that play important roles within the network."),
                                br(),
                                p("Together, these sections of the Network Analysis tab allow users to explore and visualize gene relationships, generate heatmaps, and analyze network properties in a comprehensive manner. This tab provides the tools needed to understand complex biological networks and uncover potential insights into gene regulation, interactions, and pathways."),
                                br(),
                ),#closes NA sidebarPanel
                mainPanel(
                  #add tabs
                  tabsetPanel(
                    
                    tabPanel("Network Analysis Data",
                              sidebarLayout(
                                sidebarPanel(
                                  fileInput("select_nacount", "Load Normalized Counts Data:")
                                ),#closes NA data sidebarPanel
                                mainPanel(
                                  DTOutput("nacount_table")
                                )#closes NA data mainPanel
                              )#closes NA data sidebarLayout
                              ),#closes NA data tabPanel
                    
                    tabPanel("Heatmap",
                             sidebarLayout(
                               sidebarPanel(
                                 textAreaInput("select_genes", "Enter Genes of Interest", 
                                               value = "",  # Default empty value
                                               placeholder = "Enter one gene name per line", 
                                               rows = 4,  # Display space for two lines
                                               width = "100%"),
                                 actionButton("plot_button", "Plot Heatmap")
                               ),#closes heatmap sidebarPanel
                               mainPanel(
                                 fluidPage(
                                   plotOutput("na_heatmap")
                                 )#closes heatmap fluidpage
                               )#closes heatmap mainPanel
                             )#closes heatmap sidebarlayout
                    ),#closes heatmap tabpanel
                    
                    tabPanel("Network Visualization",
                             sidebarLayout(
                               sidebarPanel(
                                 sliderInput("min_corr", "Select Minimum Correlation for Drawing an Edge:",
                                             min = 0, max = 1, value = .1)
                               ),#closes sidebarPanel
                               mainPanel(
                                 fluidPage(
                                   plotOutput("na_visualization")
                                 )#closes vis fluidPage
                               )#closes vis mainPanel
                             )#closes vis sidebarLayout
                             ),#close NA tabPanel
                    
                    tabPanel("Network Metrics",
                             sidebarLayout(
                               sidebarPanel(
                                 
                               ),#closes NM sidebarPanel
                               mainPanel(
                                 fluidPage(
                                   DTOutput("na_metrics")
                                 )#closes NM
                               )#closes NM mainPanel
                             )#closes NM sidebarLayout
                             )#closes NM tabPanel

                    
                    
                  )#closes NA tabsetPanel
                )#closes NA mainPanel
              )#closes NA sidebarLayout
              ), # closes nav_panel("Network Analysis")
    
    # create nav_menu for supplementary materials/links
    navbarMenu("Supplementary Materials",
               nav_panel("Credits", "This project was made for BF 591 at Boston University for the Fall 2024 Final. Thank you to the BF 591 team for a great semester and the support to create this Shiny App."), # closes nav_panel("Credits and Links")
               "----", 
               "Description:", 
               nav_item(
                 a("NCBI Listing", href = "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810")
               ) # closes nav_item
    ) # closes navbarMenu
  ) # closes navset_card_pill
  
) # closes UI

# server logic
server <- function(input, output, session) {
  
  # SAMPLE
  # reactive to handle 'sample_data' loading
  sample_data <- reactive({
    # check if a file is selected
    req(input$select_sample) # ensure input exists
    file_path <- input$select_sample$datapath
    sample_summary(file_path) # call sample_summary from 'sample_functions.R'
  })
  # render 'sample_table'
  output$sample_table <- renderDT({
    sample_data() # display the tibble
  })
  # render 'sample_graph'
  output$sample_graph <- renderPlot({
    sample_plot(sample_data(), input$x_name, input$y_name, input$plot_type)
  })
  # render 'summary_tibble'
  output$summary_tibble <- renderDT({
    dataset_summary(sample_data()) # pass reactive sample_data() to dataset_summary from 'sample_functions.R'
  })
  
  # COUNTS
  # reactive to handle 'norm_data' data loading
  norm_data <- reactive({
    req(input$select_normcount)  # ensure input exists
    load_counts(input$select_normcount$datapath) # call load_counts from 'counts_functions.R'
  })
  
  # zero_threshold and percentile_threshold slider input, reactive
  norm_filtered <- reactive({
    req(norm_data()) # ensure input exists
    filter_counts(norm_data(), input$zero_threshold, input$percentile_threshold) # call filter_counts from 'counts_functions.R'
  })
  
  # display filtered data
  output$norm_filtered_table <- renderDT({
    req(norm_filtered)
    norm_filtered() # directly pass norm_filtered
  })
  
  # summary output, reactive
  output$norm_count_summary <- renderDT({
    req(norm_data())
    summarize_counts(norm_data(), norm_filtered()) # call summarize_counts from 'counts_functions.R'
  })
  
  # generate 'norm_plot_data' to be used for both diagnostic scatter plots
  # reactive to handle 'norm_data' data loading
  norm_plot_data <- reactive({
    req(norm_data())  # ensure input exists
    process_normplot(norm_data(), norm_filtered()) # call process_normplot from 'counts_functions.R'
  })
  
  # render median count v variance diagnostic scatter plot
  output$diagplot_mv <- renderPlot({
    med_var_plot(norm_plot_data()) # call med_var_plot from 'counts_functions.R'
  })  
  
  # render median count v zero count diagnostic scatter plot
  output$diagplot_mz <- renderPlot({
    med_zero_plot(norm_plot_data()) # call med_zero_plot from 'counts_functions.R'
  })
  
  # perform PCA and plot PC1 versus PC2 with dynamic color input
  output$pca_plot <- renderPlot({
    pca_norm(norm_data(), sample_data(), input$pc_choiceone, input$pc_choicetwo, input$pca_color_variable)
  })
  
  # plot heatmap of filtered counts
  output$heat_map <- renderPlot({
    plot_heatmap(norm_filtered(), input$colorblind, input$num_colors)
  })
  
  # plot heatmap of correlations between samples
  output$correlation_heatmap <- renderPlot({
    plot_corrmap(norm_filtered())
  })
  
  # DIFFERENTIAL EXPRESSION
  
  #load in data, make reactive
  deseq_data <- reactive({
    #check if file is selected
    req(input$select_deseqdata)
    file_path <- input$select_deseqdata$datapath
    load_deseq(file_path) #call from functions
  })
  
  #display plot of data
  output$deseq_table <- renderDT({
    deseq_data()
  })
  
  #plot pvalue and padj histograms
  output$p_histogram <- renderPlot({
    plot_p_histogram(deseq_data(), input$p_histogram_choice)
  })
  
  #plot pvalue and padj volcano plots
  output$p_volcano <- renderPlot({
    plot_p_volcano(deseq_data(), input$p_histogram_choice, input$pfilter_de)
  })

#NETWORK ANALYSIS
# Load in data, make reactive
na_data <- reactive({
  req(input$select_nacount)  # Ensure file is selected
  file_path <- input$select_nacount$datapath
  load_nadata(file_path)  # Your function to load the data
})
#display output table
output$nacount_table <- renderDT({
  na_data()
})

#calculate corr matrix
cor_matrix <- eventReactive(input$plot_button, {
  selected_genes <- strsplit(input$select_genes, "\n")[[1]] %>% trimws()
  
  # Ensure at least one gene is selected
  if (length(selected_genes) == 0) {
    showNotification("Please enter at least one gene name.", type = "error")
    return(NULL)
  }
  
  # Call nacorr_matrix with the selected genes
  nacorr_matrix(na_data(), selected_genes)  # Pass the selected genes to nacorr_matrix
})

# Render the heatmap plot when the plot button is clicked and the correlation matrix is ready
output$na_heatmap <- renderPlot({
  # Ensure we have a correlation matrix before rendering the plot
  req(cor_matrix())  # Wait for cor_matrix() to be available
  
  # Generate heatmap from the correlation matrix
  nacorr_heatmap(cor_matrix())
})
#output the visualization of network
output$na_visualization <- renderPlot({
  req(cor_matrix())  # Wait for cor_matrix() to be available
  req(input$min_corr)  # Ensure the slider value is available
  
  # Access the minimum correlation threshold from the slider
  min_corr_threshold <- input$min_corr
  
  # Filter the correlation matrix based on the selected minimum correlation value
  filtered_matrix <- cor_matrix()
  filtered_matrix[abs(filtered_matrix) < min_corr_threshold] <- NA  # Set correlations below threshold to NA
  
  # Pass the filtered matrix to the visualization function
  visualize_correlation_network(filtered_matrix)
})
#calculate network metrics
na_metric_summary <- reactive({
  req(cor_matrix())  # Wait for cor_matrix() to be available
  req(input$min_corr)
  
  calculate_na_summary(cor_matrix())
})
#output tibble
output$na_metrics <- renderDT({
  na_metric_summary()
})

}#closes server
shinyApp(ui, server) # closes shinyApp