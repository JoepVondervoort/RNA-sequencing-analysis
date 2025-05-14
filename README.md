# RNA-seq Analysis Pipeline
**Overview**
This repository contains an R-based computational pipeline for analyzing RNA sequencing (RNA-seq) data with a focus on differential gene expression analysis using DESeq2. The pipeline is designed for comprehensive analysis of RNA-seq data, particularly for studies examining the effects of small activating RNAs (saRNAs) on gene expression.
The workflow includes quality control, differential expression analysis, visualization of results through various plots (MA plots, volcano plots, heatmaps, PCA), and functional enrichment analysis via Gene Ontology (GO) and KEGG pathways.

**Features**

Comprehensive RNA-seq Analysis: Complete workflow from count data to biological interpretation
Differential Expression Analysis: Robust identification of differentially expressed genes using DESeq2
Target Gene Analysis: Special focus on SCN1A gene expression changes
Off-Target Effect Analysis: Identification and visualization of potential off-target effects from saRNAs
Advanced Visualization: Publication-quality plots for MA, volcano, PCA, heatmaps, and expression profiles
Functional Enrichment: GO term and KEGG pathway analysis with direction-specific enrichment (up/down-regulated)
Structured Output: Organized directory structure with date-stamped results for reproducibility
Detailed Reporting: Comprehensive analysis reports and statistics

**Prerequisites**
R Packages Required:
Core Analysis
R# CRAN packages
required_packages <- c("BiocManager", "readr", "dplyr", "ggplot2", "pheatmap", 
                      "stringr", "RColorBrewer", "tibble", "ggrepel", "tidyr")

# Bioconductor packages
bioc_packages <- c("DESeq2", "org.Hs.eg.db", "AnnotationDbi", "clusterProfiler", 
                  "enrichplot", "DOSE", "pathview")
                  
**Installation**

Clone this repository to your local machine:

bashgit clone https://github.com/JoepVondervoort/rna-seq-analysis.git

cd rna-seq-analysis

Installing R Dependencies

R# Install required CRAN packages

for (pkg in c("BiocManager", "readr", "dplyr", "ggplot2", "pheatmap", 
             "stringr", "RColorBrewer", "tibble", "ggrepel", "tidyr")) {
             
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
    
}

# Install Bioconductor packages

BiocManager::install(c("DESeq2", "org.Hs.eg.db", "AnnotationDbi", 
                     "clusterProfiler", "enrichplot", "DOSE", "pathview"))
**Usage**
Input Data Format
The pipeline expects two primary input files:

A feature counts file with gene expression data
A metadata file describing experimental conditions

**Basic Execution**

R# Set your input and output paths
input_base_path <- "path/to/input/data"
output_base_path <- "path/to/output/directory"

**# Run the main analysis**

source("scripts/rna_seq_analysis.R")
Expected Output
The pipeline generates a structured output directory with the following components:

RNA_seq_analysis_YYYYMMDD/

├── results/                       # Raw DESeq2 results in CSV format

├── figures/

│   ├── MA_plots/                  # MA plots for all comparisons

│   ├── volcano_plots/             # Volcano plots for all comparisons

│   ├── heatmaps/                  # Heatmaps and sample distance plots

│   ├── PCA/                       # PCA plots

│   ├── gene_ontology/             # GO and KEGG pathway analysis results

│   └── SCN1A/                     # SCN1A-specific expression analysis

└── debug_logs/                    # Detailed logs and debugging information



Example
A minimal working example:
R# Example with included test data
input_base_path <- "example/data"
output_base_path <- "example/results"
source("scripts/rna_seq_analysis.R")
Customization
The pipeline can be customized by modifying the following parameters:

Significance Thresholds: Adjust p-value and fold change thresholds
GO Analysis Parameters: Modify GO term filtering criteria
Visualization Settings: Customize plot aesthetics and dimensions
Off-Target Gene Lists: Update potential off-target genes for each saRNA
