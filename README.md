# RNA-seq Analysis Pipeline with DESeq2

## Overview
This comprehensive R pipeline performs differential gene expression analysis on RNA-seq data using DESeq2. It was specifically designed to analyze the effects of small activating RNAs (saRNAs) on gene expression, with a particular focus on SCN1A gene regulation in the context of epilepsy research. The pipeline can integrate data from two separate experiments by combining their feature count tables.

## Key Features

### Core Analysis
- **Multi-Experiment Integration**: Combines feature count tables from two independent experiments into a unified analysis
- **DESeq2 Differential Expression Analysis**: Performs pairwise comparisons between treatment conditions (saRNA7, saRNA10, saRNA12, saRNA14) and scrambled control (SCR)
- **Customizable Analysis Parameters**: Centralized parameter configuration for p-value thresholds, fold change cutoffs, and visualization settings
- **Comprehensive Quality Control**: Including PCA plots, sample distance heatmaps, and diagnostic outputs

### Visualization Suite
- **MA Plots**: With automatic labeling of significant genes and off-target effects
- **Volcano Plots**: Individual and combined multi-panel views with customizable scaling
- **Expression Heatmaps**: Top differentially expressed genes with hierarchical clustering
- **PCA Analysis**: Multiple principal component views with scree plots
- **Summary Barplots**: Overview of up/down-regulated genes across conditions

### Specialized Analyses
- **SCN1A-Focused Analysis**: Dedicated visualizations for the target gene including expression plots, fold change comparisons, and diagnostic outputs
- **lncRNA Expression Analysis**: Identifies and analyzes long non-coding RNAs with scatter plots and expression rankings
- **Venn Diagrams**: Compares differentially expressed genes across conditions with detailed gene list exports
- **Gene Ontology Enrichment**: GO term analysis (BP, MF, CC) and KEGG pathway enrichment with multiple visualization options
- **Off-Target Gene Tracking**: Monitors and highlights genes with similar sequences to the target

### Output Organization
- Automatically creates timestamped output directories
- Generates comprehensive analysis reports
- Exports all results in multiple formats (PDF, PNG, CSV)
- Maintains detailed debug logs and parameter records

## Requirements

### R Packages
- **Core**: DESeq2, tidyverse (dplyr, ggplot2, tidyr)
- **Bioconductor**: org.Hs.eg.db, AnnotationDbi, clusterProfiler, enrichplot, DOSE, pathview
- **Visualization**: pheatmap, RColorBrewer, ggrepel, VennDiagram, gridExtra
- **Additional**: stringr, tibble

### Input Data
The pipeline expects:
1. **Two feature count tables** from separate experiments (e.g., biological replicates or technical replicates)
2. **Preprocessed combined count matrix** (`combined_featurecounts_filtered.csv`) that merges the two experiments
3. Sample metadata (can be automatically generated from sample names)

#### Input File Structure:
```
input_directory/
├── Featurecounts/
│   ├── experiment1_featureCounts.raw
│   ├── experiment2_featureCounts.raw
│   ├── combined_featurecounts_filtered.csv  # Preprocessed merged data
│   └── metadata.csv (optional)
```

## Usage

1. **Prepare Input Data**: 
   - Ensure your two feature count tables are properly formatted
   - Merge them into a combined count matrix (preprocessing step)
   - Place files in the expected directory structure

2. **Set Parameters**: Modify the `analysis_params` list to adjust:
   - Significance thresholds (padj, log2FC)
   - Visualization parameters
   - Gene lists for off-target monitoring

3. **Configure Paths**:
   ```r
   input_base_path <- "path/to/your/data"
   output_base_path <- "path/for/results"
   ```

4. **Run Analysis**:
   ```r
   source("RNA_seq_analysis_pipeline.R")
   main()
   ```

## Output Structure
```
RNA_seq_analysis_YYYYMMDD_HHMMSS/
├── results/                 # DESeq2 results tables
├── figures/
│   ├── MA_plots/           # MA plots with gene labels
│   ├── volcano_plots/      # Individual and combined volcano plots
│   ├── heatmaps/          # Expression heatmaps
│   ├── PCA/               # PCA and sample clustering
│   ├── gene_ontology/     # GO enrichment results
│   ├── SCN1A/             # Target gene analysis
│   ├── lncRNA/            # lncRNA expression analysis
│   └── venn_diagrams/     # DEG overlaps between conditions
├── debug_logs/            # Analysis logs and diagnostics
└── RNA_seq_analysis_report.txt  # Comprehensive summary report
```

## Features in Detail

### Multi-Experiment Handling
- Seamlessly combines data from two independent experiments
- Handles batch effects through DESeq2's normalization
- Maintains experiment tracking throughout analysis

### Automated Gene Annotation
- Maps Ensembl IDs to gene symbols and Entrez IDs
- Handles missing annotations gracefully

### Statistical Robustness
- Multiple testing correction (Benjamini-Hochberg)
- Outlier detection and diagnostic reporting
- Size factor normalization checks

### Customization Options
- Color schemes for all visualizations
- Adjustable plot dimensions and resolution
- Flexible gene labeling strategies

### Error Handling
- Comprehensive try-catch blocks for each analysis module
- Continues analysis even if individual components fail
- Detailed error logging for troubleshooting

## Citation
If you use this pipeline in your research, please cite:
- DESeq2: Love MI, Huber W, Anders S (2014). "Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2." Genome Biology, 15, 550.
- Additional tools as appropriate for your analysis

## Notes
- Designed for human RNA-seq data (uses org.Hs.eg.db)
- Optimized for experiments with multiple treatment conditions vs. control
- Includes specific features for saRNA and gene activation studies
- Requires preprocessing step to merge feature count tables from two experiments

