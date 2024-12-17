# R Shiny Visualizing Huntington's Disease Post-Mortem Prefrontal Cortex

This Shiny app is designed to analyze the [GSE64810 RNA-Seq dataset](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE64810), which investigates transcriptional changes in postmortem human caudate nucleus samples associated with Huntington's Disease (HD). The dataset, available through the NCBI Gene Expression Omnibus (GEO), includes RNA-Seq data from 80 samples spanning:

**HD Progression Stages:**

Controls (n = 20),
Early-stage HD (Vonsattel Grade 0–1, n = 20),
Mid-stage HD (Grade 2, n = 20),
Late-stage HD (Grade 3–4, n = 20).
Sample Source: Postmortem caudate nucleus tissues.

Platform: Illumina HiSeq 2000 (RNA-Seq raw counts data).

**App Features**

The app allows users to interactively analyze and visualize this dataset with the following modules:

**Sample Metadata Exploration:**

View and summarize patient metadata, including disease stage, age, and sample quality metrics.

**Gene Expression Analysis:**

Filter genes by variance percentile and non-zero counts.
Visualize gene expression distributions across samples using histograms, boxplots, and scatter plots.
Perform Principal Component Analysis (PCA) to explore sample clustering.

**Differential Expression Analysis:**

Identify genes differentially expressed across HD stages and controls.
Examine results with sortable tables, interactive volcano plots, and MA plots.
Customize thresholds for log-fold change and statistical significance.

**Gene Network Visualization:**

Explore gene relationships with correlation heatmaps.
Generate interactive network graphs to visualize co-expression patterns.


This app provides an accessible and user-friendly platform for exploring RNA-Seq data, enabling users to uncover key transcriptional changes and insights into the molecular mechanisms driving Huntington's Disease progression.
