# RNA-seq analysis with DESeq2 in R

# Install and load required libraries
required_packages <- c("BiocManager", "readr", "dplyr", "ggplot2", "pheatmap", 
                      "stringr", "RColorBrewer", "tibble", "ggrepel", "tidyr")
for (pkg in required_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg, repos = "https://cloud.r-project.org")
}

# Install Bioconductor packages
bioc_packages <- c("DESeq2", "org.Hs.eg.db", "AnnotationDbi", "clusterProfiler", 
                  "enrichplot", "DOSE", "pathview")
for (pkg in bioc_packages) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg)
}

if (!requireNamespace("VennDiagram", quietly = TRUE)) {
  install.packages("VennDiagram", repos = "https://cloud.r-project.org")
}

# Load all required libraries
library(DESeq2)
library(readr)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(stringr)
library(RColorBrewer)
library(tibble)
library(ggrepel)
library(tidyr)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(clusterProfiler)
library(enrichplot)
library(DOSE)
library(VennDiagram)

# Add this check after loading libraries:
check_go_packages <- function() {
  required_go <- c("clusterProfiler", "enrichplot", "DOSE", "pathview", "org.Hs.eg.db")
  missing <- c()
  
  for (pkg in required_go) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      missing <- c(missing, pkg)
    }
  }
  
  if (length(missing) > 0) {
    cat("Missing GO analysis packages:", paste(missing, collapse = ", "), "\n")
    cat("Installing missing packages...\n")
    BiocManager::install(missing)
  } else {
    cat("All GO analysis packages are installed.\n")
  }
}

# Run the check
check_go_packages()

# ============================================================================
# ANALYSIS PARAMETERS - MODIFY THESE TO CHANGE THRESHOLDS ACROSS ALL ANALYSES
# ============================================================================

# Create a configuration list for all analysis parameters
analysis_params <- list(
  # Significance thresholds
  padj_threshold = 0.01,      # Adjusted p-value threshold (use 0.05 for less stringent, 0.01 for more stringent)
  log2fc_threshold = 0.5,     # Log2 fold change threshold (use 1.0 for 2-fold change)
  
  # Gene Ontology parameters
  go_pvalue_cutoff = 0.05,
  go_qvalue_cutoff = 0.2,
  
  # Heatmap parameters
  heatmap_log2fc_threshold = 1.0,  # More stringent threshold for heatmaps
  heatmap_top_genes = 50,          # Number of top genes to show in heatmaps
  
  # MA plot parameters
  ma_plot_ylim = c(-3, 3),         # Y-axis limits for MA plots
  
  # Volcano plot parameters
  volcano_xlim = c(-4.5, 4.5),     # X-axis limits for volcano plots
  
  # Off-target genes for each condition
  off_target_genes = list(
    "saRNA7" = c("HDAC8", "RASSF8-AS1", "NET1", "KHDRBS2", "ARHGAP10", "TRIM2", "BTBD8", "ASAP2"),
    "saRNA10" = c("ICA1"),
    "saRNA12" = c("NUDT3", "PLPPR5"),
    "saRNA14" = c("MAP3K21", "DNAH2", "THADA", "DOCK3", "LINC03025")
  ),
  
  # Color palette
  saRNA_colors = c("SCR" = "#404040", 
                   "saRNA7" = "#E69F00", 
                   "saRNA10" = "#56B4E9", 
                   "saRNA12" = "#009E73", 
                   "saRNA14" = "#CC79A7")
)

# Print current parameters for documentation
cat("=== ANALYSIS PARAMETERS ===\n")
cat("Adjusted p-value threshold:", analysis_params$padj_threshold, "\n")
cat("Log2 fold change threshold:", analysis_params$log2fc_threshold, "\n")
cat("Heatmap log2FC threshold:", analysis_params$heatmap_log2fc_threshold, "\n")
cat("GO p-value cutoff:", analysis_params$go_pvalue_cutoff, "\n")
cat("GO q-value cutoff:", analysis_params$go_qvalue_cutoff, "\n")
cat("===========================\n\n")

# ============================================================================
# ANALYSIS FUNCTIONS
# ============================================================================

# Updated create_output_dirs function with Venn diagram directory
create_output_dirs <- function(base_path) {
  # Add current date and time to folder name
  current_datetime <- format(Sys.time(), "%Y%m%d_%H%M%S")
  output_folder <- paste0("RNA_seq_analysis_", current_datetime)
  
  # Make sure base_path doesn't end with a slash
  base_path <- gsub("[/\\]+$", "", base_path)
  
  # Define directory structure
  dirs <- c(
    "results",
    "figures/MA_plots",
    "figures/volcano_plots",
    "figures/heatmaps",
    "figures/PCA",
    "figures/gene_ontology",
    "figures/SCN1A",
    "figures/lncRNA",
    "figures/venn_diagrams",  # Add this line
    "debug_logs"
  )
  
  # Print the full path we're creating (for debugging)
  full_path <- file.path(base_path, output_folder)
  cat("Creating directory:", full_path, "\n")
  dir.create(full_path, recursive = TRUE, showWarnings = TRUE)
  
  # Create directories
  for (dir in dirs) {
    dir_path <- file.path(base_path, output_folder, dir)
    cat("Creating:", dir_path, "\n")
    created <- dir.create(dir_path, recursive = TRUE, showWarnings = TRUE)
    cat("  Success:", created, "\n")
  }
  
  # Return paths structure
  return(list(
    base = file.path(base_path, output_folder),
    results = file.path(base_path, output_folder, "results"),
    ma_plots = file.path(base_path, output_folder, "figures", "MA_plots"),
    volcano_plots = file.path(base_path, output_folder, "figures", "volcano_plots"),
    heatmaps = file.path(base_path, output_folder, "figures", "heatmaps"),
    pca = file.path(base_path, output_folder, "figures", "PCA"),
    go = file.path(base_path, output_folder, "figures", "gene_ontology"),
    scn1a = file.path(base_path, output_folder, "figures", "SCN1A"),
    lncrna = file.path(base_path, output_folder, "figures", "lncRNA"),
    venn = file.path(base_path, output_folder, "figures", "venn_diagrams"),  # Add this line
    debug = file.path(base_path, output_folder, "debug_logs")
  ))
}

# Function to create metadata with corrected condition labels
create_metadata <- function(filtered_combined_counts) {
  # Extract sample names
  sample_names <- colnames(filtered_combined_counts)[-1]  # Exclude Geneid column
  
  # Define conditions from sample names - focus on saRNA numbers
  conditions <- sapply(sample_names, function(name) {
    if (grepl("7", name)) {
      "saRNA7"
    } else if (grepl("10", name)) {
      "saRNA10"
    } else if (grepl("12", name)) {
      "saRNA12"
    } else if (grepl("14", name)) {
      "saRNA14"
    } else {
      "SCR"
    }
  })
  
  # Create metadata DataFrame
  metadata <- data.frame(
    SampleName = sample_names,
    Condition = conditions,
    stringsAsFactors = FALSE
  )
  
  # Print debugging info
  cat("Created metadata with", nrow(metadata), "samples\n")
  print(table(metadata$Condition))
  
  return(metadata)
}

# Function to map Ensembl IDs to gene symbols
map_ensembl_to_symbol <- function(ensembl_ids) {
  # Extract just the ENSEMBL ID without version numbers if present
  clean_ids <- gsub("\\..*$", "", ensembl_ids)
  
  # Map to gene symbols
  cat("Mapping", length(ensembl_ids), "Ensembl IDs to gene symbols...\n")
  gene_symbols <- suppressWarnings(
    AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = clean_ids,
      column = "SYMBOL",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
  )
  
  # Check mapping success
  mapped_count <- sum(!is.na(gene_symbols))
  cat("Successfully mapped", mapped_count, "Ensembl IDs (", 
      round(mapped_count/length(ensembl_ids)*100, 1), "%)\n")
  
  # Replace NAs with original IDs
  gene_symbols[is.na(gene_symbols)] <- ensembl_ids[is.na(gene_symbols)]
  
  return(gene_symbols)
}

# Function to map Ensembl IDs to Entrez IDs for gene ontology analysis
map_ensembl_to_entrez <- function(ensembl_ids) {
  # Extract just the ENSEMBL ID without version numbers if present
  clean_ids <- gsub("\\..*$", "", ensembl_ids)
  
  # Map to Entrez IDs
  cat("Mapping", length(ensembl_ids), "Ensembl IDs to Entrez IDs...\n")
  entrez_ids <- suppressWarnings(
    AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = clean_ids,
      column = "ENTREZID",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
  )
  
  # Check mapping success
  mapped_count <- sum(!is.na(entrez_ids))
  cat("Successfully mapped", mapped_count, "Ensembl IDs to Entrez IDs (", 
      round(mapped_count/length(ensembl_ids)*100, 1), "%)\n")
  
  return(entrez_ids)
}

# Function to run DESeq2 analysis - Updated to use analysis_params
run_deseq2_analysis <- function(counts_df, metadata_df, debug_path, params) {
  # Convert data to matrix and DataFrame
  cat("Preparing count matrix...\n")
  countData <- as.matrix(counts_df[,-1])
  rownames(countData) <- counts_df[[1]]
  
  # Write count data summary to debug file
  debug_file <- file.path(debug_path, "count_data_summary.txt")
  cat("Count data summary\n", file = debug_file)
  cat("Date:", Sys.Date(), "\n", file = debug_file, append = TRUE)
  cat("Total genes:", nrow(countData), "\n", file = debug_file, append = TRUE)
  cat("Total samples:", ncol(countData), "\n", file = debug_file, append = TRUE)
  cat("Analysis parameters:\n", file = debug_file, append = TRUE)
  cat("  padj threshold:", params$padj_threshold, "\n", file = debug_file, append = TRUE)
  cat("  log2FC threshold:", params$log2fc_threshold, "\n", file = debug_file, append = TRUE)
  cat("Count statistics:\n", file = debug_file, append = TRUE)
  
  count_stats <- data.frame(
    Min = apply(countData, 2, min),
    Mean = round(apply(countData, 2, mean), 2),
    Median = apply(countData, 2, median),
    Max = apply(countData, 2, max),
    Zeros = apply(countData == 0, 2, sum),
    ZerosPercent = round(apply(countData == 0, 2, sum) / nrow(countData) * 100, 2)
  )
  
  write.table(count_stats, file = debug_file, append = TRUE, sep = "\t", quote = FALSE)
  
  # Convert Condition to factor and ensure SCR is the reference level
  metadata_df$Condition <- factor(metadata_df$Condition, 
                                 levels = c("SCR", "saRNA7", "saRNA10", "saRNA12", "saRNA14"))
  
  # Use base R instead of tibble's column_to_rownames
  colData <- metadata_df
  rownames(colData) <- colData$SampleName
  colData$SampleName <- NULL
  
  # Verify sample matching
  cat("Checking sample matching between count data and metadata...\n")
  count_samples <- colnames(countData)
  meta_samples <- rownames(colData)
  
  matching_samples <- intersect(count_samples, meta_samples)
  cat("Matching samples:", length(matching_samples), "out of", length(count_samples), "in count data\n")
  
  if (length(matching_samples) < length(count_samples)) {
    missing_samples <- setdiff(count_samples, meta_samples)
    cat("WARNING: Missing", length(missing_samples), "samples in metadata:", 
        paste(missing_samples, collapse=", "), "\n")
  }
  
  # Update count data to include only matching samples
  if (length(matching_samples) < length(count_samples)) {
    countData <- countData[, matching_samples]
    cat("Count data updated to include only matching samples.\n")
  }
  
  # Create DESeqDataSet
  cat("Creating DESeqDataSet with", nrow(countData), "genes and", ncol(countData), "samples\n")
  dds <- DESeqDataSetFromMatrix(
    countData = countData, 
    colData = colData, 
    design = ~Condition
  )
  
  # Run DESeq2 analysis
  cat("Running DESeq2 analysis...\n")
  dds <- DESeq(dds)
  
  # Get size factors for normalization check
  size_factors <- sizeFactors(dds)
  cat("Size factors range:", min(size_factors), "to", max(size_factors), "\n")
  write.table(data.frame(Sample = names(size_factors), SizeFactor = size_factors),
             file = file.path(debug_path, "size_factors.txt"), sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Perform pairwise comparisons against SCR
  results_list <- list()
  conditions <- levels(colData$Condition)
  conditions <- conditions[conditions != "SCR"]
  
  cat("Performing comparisons for conditions:", paste(conditions, collapse = ", "), "\n")
  
  # Map Ensembl IDs to gene symbols and Entrez IDs
  cat("Mapping gene identifiers...\n")
  gene_symbols <- map_ensembl_to_symbol(rownames(dds))
  entrez_ids <- map_ensembl_to_entrez(rownames(dds))
  
  for (cond in conditions) {
    cat("  Processing comparison:", cond, "vs SCR\n")
    res <- results(dds, contrast = c("Condition", cond, "SCR"))
    resOrdered <- res[order(res$pvalue), ]
    
    # Convert to data frame and add gene symbols and Entrez IDs
    res_df <- as.data.frame(resOrdered)
    res_df$gene_symbol <- gene_symbols[match(rownames(res_df), rownames(dds))]
    res_df$entrez <- entrez_ids[match(rownames(res_df), rownames(dds))]
    
    # Mark SCN1A gene(s)
    res_df$is_SCN1A <- res_df$gene_symbol == "SCN1A"
    
    results_list[[cond]] <- res_df
    
    # Summary of results using centralized thresholds
    sig_genes <- sum(results_list[[cond]]$padj < params$padj_threshold, na.rm = TRUE)
    up_genes <- sum(results_list[[cond]]$padj < params$padj_threshold & 
                   results_list[[cond]]$log2FoldChange > params$log2fc_threshold, na.rm = TRUE)
    down_genes <- sum(results_list[[cond]]$padj < params$padj_threshold & 
                     results_list[[cond]]$log2FoldChange < -params$log2fc_threshold, na.rm = TRUE)
    
    cat("    Found", sig_genes, "significantly differentially expressed genes (padj <", 
        params$padj_threshold, ")\n")
    cat("    Up-regulated genes (log2FC >", params$log2fc_threshold, "):", up_genes, "\n")
    cat("    Down-regulated genes (log2FC <", -params$log2fc_threshold, "):", down_genes, "\n")
    
    # Check if SCN1A is in results
    scn1a_rows <- which(res_df$is_SCN1A)
    if (length(scn1a_rows) > 0) {
      cat("    SCN1A found in results with log2FoldChange:", 
          round(res_df$log2FoldChange[scn1a_rows], 3), 
          "and padj:", format(res_df$padj[scn1a_rows], scientific = TRUE, digits = 3), "\n")
    } else {
      cat("    SCN1A not found in results\n")
    }
    
    # Write detailed statistics to file
    stats_file <- file.path(debug_path, paste0(cond, "_vs_SCR_stats.txt"))
    cat(paste0("Statistics for ", cond, " vs SCR\n"), file = stats_file)
    cat("Date:", Sys.Date(), "\n", file = stats_file, append = TRUE)
    cat("Analysis parameters:\n", file = stats_file, append = TRUE)
    cat("  padj threshold:", params$padj_threshold, "\n", file = stats_file, append = TRUE)
    cat("  log2FC threshold:", params$log2fc_threshold, "\n", file = stats_file, append = TRUE)
    cat("Total genes analyzed:", nrow(res_df), "\n", file = stats_file, append = TRUE)
    cat("Genes with valid p-value:", sum(!is.na(res_df$pvalue)), "\n", file = stats_file, append = TRUE)
    cat("Genes with valid padj:", sum(!is.na(res_df$padj)), "\n", file = stats_file, append = TRUE)
    cat("Significant genes (padj <", params$padj_threshold, "):", sig_genes, "\n", file = stats_file, append = TRUE)
    cat("Up-regulated genes (padj <", params$padj_threshold, ", log2FC >", params$log2fc_threshold, "):", up_genes, "\n", file = stats_file, append = TRUE)
    cat("Down-regulated genes (padj <", params$padj_threshold, ", log2FC <", -params$log2fc_threshold, "):", down_genes, "\n", file = stats_file, append = TRUE)
    cat("SCN1A genes found:", length(scn1a_rows), "\n", file = stats_file, append = TRUE)
    
    if (length(scn1a_rows) > 0) {
      cat("SCN1A details:\n", file = stats_file, append = TRUE)
      for (i in scn1a_rows) {
        cat(paste0("  Ensembl ID: ", rownames(res_df)[i], "\n"), file = stats_file, append = TRUE)
        cat(paste0("  log2FC: ", res_df$log2FoldChange[i], "\n"), file = stats_file, append = TRUE)
        cat(paste0("  p-value: ", res_df$pvalue[i], "\n"), file = stats_file, append = TRUE)
        cat(paste0("  padj: ", res_df$padj[i], "\n"), file = stats_file, append = TRUE)
      }
    }
  }
  
  return(list(dds = dds, results_list = results_list, gene_symbols = gene_symbols, entrez_ids = entrez_ids))
}

# Function to generate PCA plots
generate_pca_plots <- function(dds, pca_output_path, color_palette) {
  cat("Generating PCA plots...\n")
  
  # Get variance stabilized data
  vsd <- vst(dds, blind = TRUE)
  
  # Calculate PCA
  pcaData <- plotPCA(vsd, intgroup = "Condition", returnData = TRUE)
  percentVar <- round(100 * attr(pcaData, "percentVar"))
  
  # Create enhanced PCA plot
  p1 <- ggplot(pcaData, aes(x = PC1, y = PC2, color = Condition)) +
    geom_point(size = 4, alpha = 0.8) +
    scale_color_manual(values = color_palette) +
    xlab(paste0("PC1: ", percentVar[1], "% variance")) +
    ylab(paste0("PC2: ", percentVar[2], "% variance")) +
    labs(
      title = "Principal Component Analysis",
      subtitle = "Samples colored by condition",
      caption = paste0("Analysis date: ", Sys.Date())
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      panel.border = element_rect(color = "grey80", fill = NA),
      axis.title = element_text(face = "bold")
    ) +
    coord_fixed()
  
  # Save plot
  ggsave(file.path(pca_output_path, "PCA_plot_PC1_PC2.pdf"), 
         plot = p1, width = 10, height = 8)
  ggsave(file.path(pca_output_path, "PCA_plot_PC1_PC2.png"), 
         plot = p1, width = 10, height = 8, dpi = 300)
  
  # Create PCA plot with sample labels
  p2 <- p1 + 
    geom_text_repel(
      aes(label = rownames(pcaData)),
      size = 3,
      box.padding = 0.5,
      point.padding = 0.3,
      max.overlaps = 30
    )
  
  ggsave(file.path(pca_output_path, "PCA_plot_PC1_PC2_labeled.pdf"), 
         plot = p2, width = 12, height = 9)
  ggsave(file.path(pca_output_path, "PCA_plot_PC1_PC2_labeled.png"), 
         plot = p2, width = 12, height = 9, dpi = 300)
  
  # Create PCA plot for PC1 vs PC3
  if (ncol(assay(vsd)) > 3) {
    # Calculate full PCA manually to get more PCs
    pca <- prcomp(t(assay(vsd)))
    percentVar_all <- pca$sdev^2 / sum(pca$sdev^2) * 100
    
    # Create data frame for plotting
    pcaData_all <- data.frame(
      PC1 = pca$x[,1],
      PC2 = pca$x[,2],
      PC3 = pca$x[,3],
      Condition = vsd$Condition,
      row.names = colnames(vsd)
    )
    
    # PC1 vs PC3 plot
    p3 <- ggplot(pcaData_all, aes(x = PC1, y = PC3, color = Condition)) +
      geom_point(size = 4, alpha = 0.8) +
      scale_color_manual(values = color_palette) +
      xlab(paste0("PC1: ", round(percentVar_all[1], 1), "% variance")) +
      ylab(paste0("PC3: ", round(percentVar_all[3], 1), "% variance")) +
      labs(
        title = "Principal Component Analysis (PC1 vs PC3)",
        subtitle = "Alternative view of sample relationships"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        panel.border = element_rect(color = "grey80", fill = NA),
        axis.title = element_text(face = "bold")
      ) +
      coord_fixed()
    
    ggsave(file.path(pca_output_path, "PCA_plot_PC1_PC3.pdf"), 
           plot = p3, width = 10, height = 8)
    ggsave(file.path(pca_output_path, "PCA_plot_PC1_PC3.png"), 
           plot = p3, width = 10, height = 8, dpi = 300)
    
    # Create scree plot
    variance_df <- data.frame(
      PC = paste0("PC", 1:min(20, length(percentVar_all))),
      Variance = percentVar_all[1:min(20, length(percentVar_all))],
      Cumulative = cumsum(percentVar_all[1:min(20, length(percentVar_all))])
    )
    
    variance_df$PC <- factor(variance_df$PC, levels = variance_df$PC)
    
    p_scree <- ggplot(variance_df, aes(x = PC, y = Variance)) +
      geom_bar(stat = "identity", fill = "steelblue", alpha = 0.8) +
      geom_line(aes(y = Cumulative, group = 1), color = "red", size = 1.2) +
      geom_point(aes(y = Cumulative), color = "red", size = 3) +
      scale_y_continuous(
        name = "Variance Explained (%)",
        sec.axis = sec_axis(~., name = "Cumulative Variance (%)")
      ) +
      labs(
        title = "PCA Scree Plot",
        subtitle = "Variance explained by each principal component",
        x = "Principal Component"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title = element_text(face = "bold"),
        panel.grid.major = element_line(color = "grey90"),
        panel.border = element_rect(color = "grey80", fill = NA)
      )
    
    ggsave(file.path(pca_output_path, "PCA_scree_plot.pdf"), 
           plot = p_scree, width = 10, height = 8)
    ggsave(file.path(pca_output_path, "PCA_scree_plot.png"), 
           plot = p_scree, width = 10, height = 8, dpi = 300)
  }
  
  # Save PCA data
  write.csv(pcaData, file.path(pca_output_path, "PCA_data.csv"), row.names = TRUE)
  
  cat("Created PCA plots\n")
}

# Function to generate SCN1A-specific expression plots
generate_scn1a_plots <- function(dds, results_list, scn1a_output_path, color_palette) {
  cat("Generating SCN1A expression plots...\n")
  
  # Get normalized counts
  normalized_counts <- counts(dds, normalized = TRUE)
  
  # Find SCN1A gene(s)
  scn1a_indices <- which(results_list[[1]]$gene_symbol == "SCN1A")
  
  if (length(scn1a_indices) == 0) {
    cat("  SCN1A not found in the dataset, skipping SCN1A plots\n")
    return()
  }
  
  # Get SCN1A expression data
  scn1a_ensembl <- rownames(results_list[[1]])[scn1a_indices]
  
  for (i in seq_along(scn1a_ensembl)) {
    ensembl_id <- scn1a_ensembl[i]
    
    # Create expression data frame
    expr_data <- data.frame(
      Sample = colnames(normalized_counts),
      Expression = normalized_counts[ensembl_id, ],
      Condition = colData(dds)$Condition,
      stringsAsFactors = FALSE
    )
    
    # Order conditions
    expr_data$Condition <- factor(expr_data$Condition, 
                                 levels = c("SCR", "saRNA7", "saRNA10", "saRNA12", "saRNA14"))
    
    # Create boxplot
    p1 <- ggplot(expr_data, aes(x = Condition, y = Expression, fill = Condition)) +
      geom_boxplot(alpha = 0.8, outlier.shape = NA) +
      geom_jitter(width = 0.2, size = 3, alpha = 0.6) +
      scale_fill_manual(values = color_palette) +
      labs(
        title = paste0("SCN1A Expression (", ensembl_id, ")"),
        subtitle = "Normalized counts by condition",
        x = "",
        y = "Normalized Expression",
        caption = paste0("Analysis date: ", Sys.Date())
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        legend.position = "none",
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.border = element_rect(color = "grey80", fill = NA)
      )
    
    # Save boxplot
    ggsave(file.path(scn1a_output_path, paste0("SCN1A_boxplot_", i, ".pdf")), 
           plot = p1, width = 10, height = 8)
    ggsave(file.path(scn1a_output_path, paste0("SCN1A_boxplot_", i, ".png")), 
           plot = p1, width = 10, height = 8, dpi = 300)
    
    # Create bar plot with mean and error bars
    expr_summary <- expr_data %>%
      group_by(Condition) %>%
      summarise(
        Mean = mean(Expression),
        SD = sd(Expression),
        SE = sd(Expression) / sqrt(n()),
        n = n()
      )
    
    p2 <- ggplot(expr_summary, aes(x = Condition, y = Mean, fill = Condition)) +
      geom_bar(stat = "identity", alpha = 0.8) +
      geom_errorbar(aes(ymin = Mean - SE, ymax = Mean + SE), 
                    width = 0.2, size = 0.8) +
      geom_point(data = expr_data, aes(x = Condition, y = Expression),
                position = position_jitter(width = 0.1), 
                size = 2, alpha = 0.6) +
      scale_fill_manual(values = color_palette) +
      labs(
        title = paste0("SCN1A Mean Expression (", ensembl_id, ")"),
        subtitle = "Mean ± SEM with individual data points",
        x = "",
        y = "Normalized Expression",
        caption = paste0("n = ", paste(expr_summary$n, collapse = ", "))
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        plot.caption = element_text(hjust = 1, size = 10),
        legend.position = "none",
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.border = element_rect(color = "grey80", fill = NA)
      )
    
    ggsave(file.path(scn1a_output_path, paste0("SCN1A_barplot_", i, ".pdf")), 
           plot = p2, width = 10, height = 8)
    ggsave(file.path(scn1a_output_path, paste0("SCN1A_barplot_", i, ".png")), 
           plot = p2, width = 10, height = 8, dpi = 300)
    
    # Create a heatmap showing SCN1A expression across all samples
    # Get VST normalized data for better visualization
    vsd <- vst(dds, blind = TRUE)
    scn1a_vst <- assay(vsd)[ensembl_id, , drop = FALSE]
    
    # Create annotation for heatmap
    annotation_col <- data.frame(
      Condition = colData(dds)$Condition,
      row.names = colnames(scn1a_vst)
    )
    
    # Set colors
    ann_colors <- list(Condition = color_palette)
    
    # Create heatmap
    pdf(file.path(scn1a_output_path, paste0("SCN1A_heatmap_", i, ".pdf")), 
        width = 12, height = 4)
    pheatmap(scn1a_vst,
            cluster_rows = FALSE,
            cluster_cols = TRUE,
            annotation_col = annotation_col,
            annotation_colors = ann_colors,
            main = paste0("SCN1A Expression Heatmap (", ensembl_id, ")"),
            show_colnames = TRUE,
            show_rownames = TRUE,
            fontsize = 10,
            fontsize_col = 8,
            color = colorRampPalette(c("navy", "white", "red"))(50),
            border_color = NA)
    dev.off()
    
    # Save expression data
    write.csv(expr_data, 
             file.path(scn1a_output_path, paste0("SCN1A_expression_data_", i, ".csv")),
             row.names = FALSE)
    
    # Create fold change plot
    fc_data <- data.frame(
      Condition = c("saRNA7", "saRNA10", "saRNA12", "saRNA14"),
      log2FC = numeric(4),
      padj = numeric(4),
      stringsAsFactors = FALSE
    )
    
    for (j in 1:nrow(fc_data)) {
      cond <- fc_data$Condition[j]
      fc_data$log2FC[j] <- results_list[[cond]]$log2FoldChange[scn1a_indices[i]]
      fc_data$padj[j] <- results_list[[cond]]$padj[scn1a_indices[i]]
    }
    
    fc_data$Significant <- fc_data$padj < 0.05
    
    p3 <- ggplot(fc_data, aes(x = Condition, y = log2FC, fill = Significant)) +
      geom_bar(stat = "identity", alpha = 0.8) +
      geom_hline(yintercept = 0, linetype = "solid", color = "black") +
      geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey50") +
      scale_fill_manual(values = c("FALSE" = "grey60", "TRUE" = "red")) +
      labs(
        title = paste0("SCN1A Log2 Fold Changes (", ensembl_id, ")"),
        subtitle = "Compared to SCR control",
        x = "",
        y = "Log2 Fold Change",
        caption = "Red bars: padj < 0.05 | Dashed lines: ±0.5 log2FC"
      ) +
      theme_minimal(base_size = 14) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
        plot.subtitle = element_text(hjust = 0.5, size = 12),
        plot.caption = element_text(hjust = 1, size = 10),
        legend.position = "none",
        axis.title.y = element_text(face = "bold"),
        axis.text.x = element_text(size = 12, face = "bold"),
        panel.grid.major.x = element_blank(),
        panel.border = element_rect(color = "grey80", fill = NA)
      ) +
      # Add padj values as text
      geom_text(aes(label = paste0("p=", format(padj, scientific = TRUE, digits = 2))),
               vjust = ifelse(fc_data$log2FC > 0, -0.5, 1.5),
               size = 3.5)
    
    ggsave(file.path(scn1a_output_path, paste0("SCN1A_fold_changes_", i, ".pdf")), 
           plot = p3, width = 10, height = 8)
    ggsave(file.path(scn1a_output_path, paste0("SCN1A_fold_changes_", i, ".png")), 
           plot = p3, width = 10, height = 8, dpi = 300)
  }
  
  cat("Created SCN1A expression plots\n")
}

# Additional check for gridExtra package for combined volcano plots
if (!requireNamespace("gridExtra", quietly = TRUE)) {
  install.packages("gridExtra", repos = "https://cloud.r-project.org")
}
library(gridExtra)

# Function to generate volcano plots with off-target gene labels and combined plot
generate_volcano_plots <- function(results_list, volcano_output_path, params) {
  cat("Generating volcano plots with improved styling and off-target gene labeling...\n")
  cat("  Using padj threshold:", params$padj_threshold, "\n")
  cat("  Using log2FC threshold:", params$log2fc_threshold, "\n")
  
  # Store individual plots for combined figure
  plot_list <- list()
  plot_list_uniform <- list()  # For uniform scale plots
  
  # Find the maximum y-value across all conditions for uniform scaling
  max_y_all <- 0
  for (cond in names(results_list)) {
    res_df <- results_list[[cond]]
    current_max <- max(-log10(res_df$padj), na.rm = TRUE)
    if (current_max > max_y_all) {
      max_y_all <- current_max
    }
  }
  
  # Create nice breaks for uniform scale
  uniform_y_breaks <- pretty(c(0, max_y_all), n = 8)
  
  for (cond in names(results_list)) {
    cat(paste0("  Creating volcano plot for ", cond, " vs SCR\n"))
    res_df <- results_list[[cond]]
    
    # Get total counts for reporting
    total_genes <- nrow(res_df)
    genes_with_padj <- sum(!is.na(res_df$padj))
    scn1a_count <- sum(res_df$is_SCN1A, na.rm = TRUE)
    
    # Count significant genes with the centralized thresholds
    sig_genes <- sum(!is.na(res_df$padj) & res_df$padj < params$padj_threshold, na.rm = TRUE)
    up_genes <- sum(!is.na(res_df$padj) & res_df$padj < params$padj_threshold & 
                   res_df$log2FoldChange > params$log2fc_threshold, na.rm = TRUE)
    down_genes <- sum(!is.na(res_df$padj) & res_df$padj < params$padj_threshold & 
                     res_df$log2FoldChange < -params$log2fc_threshold, na.rm = TRUE)
    
    cat("    Total genes: ", total_genes, "\n")
    cat("    Genes with padj value: ", genes_with_padj, "\n")
    cat("    Significant genes (padj < ", params$padj_threshold, "): ", sig_genes, "\n")
    cat("    Up-regulated genes: ", up_genes, "\n")
    cat("    Down-regulated genes: ", down_genes, "\n")
    cat("    SCN1A genes found: ", scn1a_count, "\n")
    
    if (scn1a_count > 0) {
      scn1a_rows <- which(res_df$is_SCN1A)
      for (i in scn1a_rows) {
        cat("    SCN1A details: log2FC =", res_df$log2FoldChange[i], ", padj =", res_df$padj[i], "\n")
      }
    }
    
    # Add categories for significance using centralized thresholds
    res_df$Expression <- "Unchanged"
    res_df$Expression[!is.na(res_df$padj) & res_df$padj < params$padj_threshold & 
                     res_df$log2FoldChange < -params$log2fc_threshold] <- "Down-regulated"
    res_df$Expression[!is.na(res_df$padj) & res_df$padj < params$padj_threshold & 
                     res_df$log2FoldChange > params$log2fc_threshold] <- "Up-regulated"
    
    # Add a special category for SCN1A
    if (scn1a_count > 0) {
      # Mark SCN1A regardless of significance for visibility
      res_df$Expression[res_df$is_SCN1A] <- "SCN1A"
    }
    
    # Add flag for off-target genes based on the current condition
    res_df$is_OffTarget <- FALSE
    if (cond %in% names(params$off_target_genes)) {
      # Mark the specific off-target genes for this condition
      res_df$is_OffTarget <- res_df$gene_symbol %in% params$off_target_genes[[cond]]
      # Log similar sequence genes found
      off_target_found <- res_df$gene_symbol[res_df$is_OffTarget]
      if (length(off_target_found) > 0) {
        cat("    Similar sequence genes found:", paste(off_target_found, collapse=", "), "\n")
        # For genes that are both similar sequence and would be classified by expression, prioritize similar sequence
        res_df$Expression[res_df$is_OffTarget] <- "Off-Target"
      } else {
        cat("    No similar sequence genes found in the dataset\n")
      }
    }
    
    # Convert to factor with defined levels
    res_df$Expression <- factor(res_df$Expression, 
                      levels = c("Down-regulated", "Up-regulated", "SCN1A", "Off-Target", "Unchanged"))
    
    # Calculate max y value for this specific plot
    max_y_current <- max(-log10(res_df$padj), na.rm = TRUE)
    
    # Create nice breaks for individual plot
    if (max_y_current > 100) {
      # For very large values, use breaks of 50
      individual_y_breaks <- seq(0, ceiling(max_y_current/50)*50, by = 50)
    } else if (max_y_current > 50) {
      # For medium-large values, use breaks of 20
      individual_y_breaks <- seq(0, ceiling(max_y_current/20)*20, by = 20)
    } else if (max_y_current > 20) {
      # For medium values, use breaks of 10
      individual_y_breaks <- seq(0, ceiling(max_y_current/10)*10, by = 10)
    } else {
      # For smaller values, use breaks of 5
      individual_y_breaks <- seq(0, ceiling(max_y_current/5)*5, by = 5)
    }
    
    # Function to create the base plot
    create_base_plot <- function(y_breaks, y_limit) {
      p <- ggplot(res_df, aes(x = log2FoldChange, y = -log10(padj), color = Expression)) +
        # Add points
        geom_point(size = 1.5, alpha = 0.8) +
        # Set colors - using a nice purple for similar sequence genes
        scale_color_manual(values = c("Down-regulated" = "#3366CC", 
                                    "Up-regulated" = "#CC0000", 
                                    "SCN1A" = "green3",
                                    "Off-Target" = "#9966CC",  # Purple color for off-target genes
                                    "Unchanged" = "grey70"),
                         labels = c("Down-regulated", "Up-regulated", 
                                   "SCN1A", "Similar Sequence Gene", "Unchanged")) +
        # Add reference lines with the defined thresholds
        geom_vline(xintercept = c(-params$log2fc_threshold, params$log2fc_threshold), 
                  linetype = "dotted", color = "darkgray", size = 0.5) +
        geom_hline(yintercept = -log10(params$padj_threshold), linetype = "dotted", color = "darkgray", size = 0.5) +
        # Set limits with custom breaks
        scale_y_continuous(limits = c(0, y_limit),
                           breaks = y_breaks,
                           expand = expansion(mult = c(0, 0.02))) +
        scale_x_continuous(limits = params$volcano_xlim, 
                           breaks = c(-4, -2, -params$log2fc_threshold, 0, 
                                     params$log2fc_threshold, 2, 4),
                           labels = function(x) {
                             ifelse(abs(x) == params$log2fc_threshold,
                                    paste0(ifelse(x < 0, "-", ""), params$log2fc_threshold),
                                    as.character(x))
                           }) +
        # Add detailed labels
        labs(
          title = paste0(cond, " → SCR"),
          subtitle = paste0("Total: ", total_genes, " | Up: ", up_genes, " | Down: ", down_genes),
          x = "log2(FC)",
          y = "-log10(padj)"
        ) +
        # Set theme
        theme_minimal(base_size = 10) +
        theme(
          panel.grid.major = element_line(color = "grey92"),
          panel.grid.minor = element_line(color = "grey96", linetype = "dotted"),
          panel.border = element_rect(color = "black", fill = NA, size = 0.5),
          plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
          plot.subtitle = element_text(hjust = 0.5, size = 9),
          legend.position = "none",  # Remove legend for individual plots in combined view
          axis.title = element_text(face = "bold", size = 10),
          axis.text = element_text(size = 8),
          plot.margin = margin(5, 5, 5, 5)
        )
      
      return(p)
    }
    
    # Create plot with individual scaling
    p <- create_base_plot(individual_y_breaks, max(individual_y_breaks))
    
    # Create plot with uniform scaling
    p_uniform <- create_base_plot(uniform_y_breaks, max(uniform_y_breaks))
    
    # Select top genes for labeling (fewer for combined plot)
    top_down <- res_df %>%
      filter(Expression == "Down-regulated") %>%
      arrange(padj) %>%
      head(3)  # Reduced from 5 to 3 for cleaner combined plot
    
    top_up <- res_df %>%
      filter(Expression == "Up-regulated") %>%
      arrange(padj) %>%
      head(3)  # Reduced from 5 to 3
    
    # Combine and label them
    top_genes <- rbind(top_down, top_up)
    
    # Function to add labels to plot
    add_labels <- function(plot) {
      if (nrow(top_genes) > 0) {
        plot <- plot + geom_text_repel(
          data = top_genes,
          aes(label = gene_symbol),
          size = 2.5,  # Smaller text for combined plot
          fontface = "bold",
          box.padding = 0.3,
          point.padding = 0.2,
          force = 8,
          max.overlaps = 20,
          min.segment.length = 0.1,
          segment.color = "grey50",
          segment.size = 0.2,
          color = ifelse(top_genes$Expression == "Up-regulated", "#CC0000", "#3366CC")
        )
      }
      
      # Always make sure SCN1A is shown if present
      if (scn1a_count > 0) {
        scn1a_data <- res_df[res_df$is_SCN1A, ]
        
        plot <- plot + geom_point(
          data = scn1a_data,
          size = 2.5,
          color = "green4",
          shape = 19
        ) +
        geom_text_repel(
          data = scn1a_data,
          aes(label = "SCN1A"),
          size = 3,
          fontface = "bold",
          box.padding = 0.5,
          force = 10,
          color = "green4",
          segment.color = "green4",
          segment.size = 0.3
        )
      }
      
      # Add labels for off-target genes
      if (cond %in% names(params$off_target_genes)) {
        off_target_data <- res_df[res_df$is_OffTarget, ]
        if (nrow(off_target_data) > 0) {
          plot <- plot + geom_point(
            data = off_target_data,
            size = 2.5,
            color = "#9966CC",
            shape = 19
          ) +
          geom_text_repel(
            data = off_target_data,
            aes(label = gene_symbol),
            size = 2.5,
            fontface = "bold",
            box.padding = 0.4,
            force = 8,
            color = "#9966CC",
            segment.color = "#9966CC",
            segment.size = 0.3,
            max.overlaps = 20
          )
        }
      }
      
      return(plot)
    }
    
    # Add labels to both plots
    p <- add_labels(p)
    p_uniform <- add_labels(p_uniform)
    
    # Store plots for combined figures
    plot_list[[cond]] <- p
    plot_list_uniform[[cond]] <- p_uniform
    
    # Save individual plot with legend
    p_individual <- p + 
      theme(legend.position = "right",
            plot.title = element_text(size = 14),
            plot.subtitle = element_text(size = 10),
            plot.caption = element_text(hjust = 1, size = 8, color = "grey50"),
            plot.margin = margin(20, 20, 20, 20)) +
      labs(caption = paste0("Analysis date: ", Sys.Date(), " | Purple points: genes with similar sequence to target"))
    
    ggsave(file.path(volcano_output_path, paste0("Volcano_plot_", cond, "_vs_SCR.pdf")), 
           plot = p_individual, width = 10, height = 8)
    ggsave(file.path(volcano_output_path, paste0("Volcano_plot_", cond, "_vs_SCR.png")), 
           plot = p_individual, width = 10, height = 8, dpi = 300)
  }
  
  # Create combined multi-panel volcano plot with individual scaling
  cat("\n  Creating combined multi-panel volcano plot...\n")
  
  library(gridExtra)
  library(grid)
  
  # Create a legend plot
  legend_data <- data.frame(
    x = 1:5,
    y = 1:5,
    Expression = factor(c("Down-regulated", "Up-regulated", "SCN1A", "Off-Target", "Unchanged"),
                       levels = c("Down-regulated", "Up-regulated", "SCN1A", "Off-Target", "Unchanged"))
  )
  
  legend_plot <- ggplot(legend_data, aes(x = x, y = y, color = Expression)) +
    geom_point(size = 3) +
    scale_color_manual(values = c("Down-regulated" = "#3366CC", 
                                "Up-regulated" = "#CC0000", 
                                "SCN1A" = "green3",
                                "Off-Target" = "#9966CC",
                                "Unchanged" = "grey70"),
                     labels = c("Down-regulated", "Up-regulated", 
                               "SCN1A", "Similar Sequence Gene", "Unchanged"),
                     name = "Gene Category") +
    theme_void() +
    theme(legend.position = "bottom",
          legend.title = element_text(face = "bold", size = 12),
          legend.text = element_text(size = 10),
          legend.key.size = unit(0.8, "cm"))
  
  # Extract just the legend
  g_legend <- function(p) {
    tmp <- ggplot_gtable(ggplot_build(p))
    leg <- which(sapply(tmp$grobs, function(x) x$name) == "guide-box")
    legend <- tmp$grobs[[leg]]
    return(legend)
  }
  
  legend <- g_legend(legend_plot)
  
  # Arrange plots in a 2x2 grid with shared legend (individual scaling)
  combined_plot <- arrangeGrob(
    plot_list[["saRNA7"]],
    plot_list[["saRNA10"]],
    plot_list[["saRNA12"]],
    plot_list[["saRNA14"]],
    ncol = 2,
    nrow = 2,
    top = textGrob("Volcano Plots: All Conditions vs SCR (Individual Scales)", 
                   gp = gpar(fontsize = 16, fontface = "bold")),
    bottom = textGrob(paste0("Thresholds: padj < ", params$padj_threshold, 
                            ", |log2FC| > ", params$log2fc_threshold,
                            " | Analysis date: ", Sys.Date()),
                     gp = gpar(fontsize = 10, col = "grey50"))
  )
  
  # Combine with legend
  final_plot <- arrangeGrob(combined_plot, legend, 
                           ncol = 1, 
                           heights = c(10, 1))
  
  # Save combined plot with individual scaling
  ggsave(file.path(volcano_output_path, "Volcano_plots_combined_individual_scales.pdf"), 
         plot = final_plot, width = 14, height = 14, dpi = 300)
  ggsave(file.path(volcano_output_path, "Volcano_plots_combined_individual_scales.png"), 
         plot = final_plot, width = 14, height = 14, dpi = 600)
  
  # Create combined plot with uniform scaling
  combined_plot_uniform <- arrangeGrob(
    plot_list_uniform[["saRNA7"]],
    plot_list_uniform[["saRNA10"]],
    plot_list_uniform[["saRNA12"]],
    plot_list_uniform[["saRNA14"]],
    ncol = 2,
    nrow = 2,
    top = textGrob("Volcano Plots: All Conditions vs SCR (Uniform Scale)", 
                   gp = gpar(fontsize = 16, fontface = "bold")),
    bottom = textGrob(paste0("Thresholds: padj < ", params$padj_threshold, 
                            ", |log2FC| > ", params$log2fc_threshold,
                            " | Y-axis max: ", max(uniform_y_breaks),
                            " | Analysis date: ", Sys.Date()),
                     gp = gpar(fontsize = 10, col = "grey50"))
  )
  
  # Combine with legend
  final_plot_uniform <- arrangeGrob(combined_plot_uniform, legend, 
                                   ncol = 1, 
                                   heights = c(10, 1))
  
  # Save combined plot with uniform scaling
  ggsave(file.path(volcano_output_path, "Volcano_plots_combined_uniform_scale.pdf"), 
         plot = final_plot_uniform, width = 14, height = 14, dpi = 300)
  ggsave(file.path(volcano_output_path, "Volcano_plots_combined_uniform_scale.png"), 
         plot = final_plot_uniform, width = 14, height = 14, dpi = 600)
  
  cat("Created volcano plots with improved y-axis breaks\n")
  cat("Generated combined plots with both individual and uniform scales\n")
}

# Function to generate enhanced summary barplot - Updated to use params
generate_summary_barplot <- function(results_list, base_output_path, params) {
  cat("Generating enhanced summary barplot...\n")
  cat("  Using padj threshold:", params$padj_threshold, "\n")
  cat("  Using log2FC threshold:", params$log2fc_threshold, "\n")
  
  # Collect significant gene counts using centralized thresholds
  sig_counts <- data.frame(
    Condition = names(results_list),
    Down_regulated = sapply(results_list, function(res) 
      sum(!is.na(res$padj) & res$padj < params$padj_threshold & 
          res$log2FoldChange < -params$log2fc_threshold, na.rm = TRUE)),
    Up_regulated = sapply(results_list, function(res) 
      sum(!is.na(res$padj) & res$padj < params$padj_threshold & 
          res$log2FoldChange > params$log2fc_threshold, na.rm = TRUE))
  )
  
  # Calculate totals
  sig_counts$Total <- sig_counts$Up_regulated + sig_counts$Down_regulated
  
  # Save count data
  write.csv(sig_counts, file.path(base_output_path, "figures", "DEG_counts.csv"), row.names = FALSE)
  
  # Reshape for plotting
  sig_counts_long <- tidyr::pivot_longer(sig_counts, 
                                       cols = c("Down_regulated", "Up_regulated"),
                                       names_to = "Regulation",
                                       values_to = "Count")
  
  # Ensure conditions are in the correct order
  sig_counts_long$Condition <- factor(sig_counts_long$Condition, 
                                     levels = c("SCR", "saRNA7", "saRNA10", "saRNA12", "saRNA14"))
  
  # Set factor levels to ensure down-regulated comes first
  sig_counts_long$Regulation <- factor(sig_counts_long$Regulation, 
                                      levels = c("Down_regulated", "Up_regulated"))
  
  # Set custom colors
  reg_colors <- c("Down_regulated" = "#0072B2", "Up_regulated" = "#D55E00") 
  
  # Create enhanced barplot
  p <- ggplot() +
    # Add bars with improved styling
    geom_bar(
      data = sig_counts_long, 
      aes(x = Condition, y = Count, fill = Regulation),
      stat = "identity", 
      position = position_dodge(width = 0.75), 
      width = 0.7,
      color = "black",
      size = 0.3
    ) +
    # Add text labels on bars
    geom_text(
      data = sig_counts_long,
      aes(x = Condition, y = Count, label = Count, group = Regulation),
      position = position_dodge(width = 0.75),
      vjust = -0.5,
      color = "black",
      size = 4,
      fontface = "bold"
    ) +
    # Use custom colors
    scale_fill_manual(
      values = reg_colors,
      labels = c("Down-regulated", "Up-regulated"),
      name = "Gene Regulation"
    ) +
    # Improved labels
    labs(
      title = "Differentially Expressed Genes by Condition",
      subtitle = paste0("padj < ", params$padj_threshold, ", |log2FC| > ", params$log2fc_threshold),
      x = "",
      y = "Number of Genes",
      caption = paste0("Analysis date: ", Sys.Date())
    ) +
    # Enhanced theme
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20, margin = margin(b = 10)),
      plot.subtitle = element_text(hjust = 0.5, size = 14, color = "gray30", margin = margin(b = 20)),
      plot.caption = element_text(hjust = 1, size = 10, color = "gray50", margin = margin(t = 15)),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(size = 12),
      legend.background = element_rect(fill = "white", color = "gray80"),
      legend.key.size = unit(1, "cm"),
      axis.title.y = element_text(face = "bold", size = 16, margin = margin(r = 10)),
      axis.text.x = element_text(size = 14, face = "bold", color = "black", margin = margin(t = 5)),
      axis.text.y = element_text(size = 12),
      panel.grid.major.y = element_line(color = "gray90"),
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "gray80", fill = NA, size = 1),
      panel.background = element_rect(fill = "gray98"),
      plot.background = element_rect(fill = "white", color = NA),
      plot.margin = margin(20, 20, 50, 20) # Increased bottom margin for total labels
    ) +
    # Customized y axis starting at 0
    scale_y_continuous(
      expand = expansion(mult = c(0, 0.1)),
      limits = c(0, NA),
      breaks = scales::pretty_breaks(n = 8)
    )
  
  # Add total counts separately as an annotation below the x-axis
  p <- p +
    annotate(
      "text",
      x = sig_counts$Condition,
      y = rep(0, nrow(sig_counts)),
      label = paste("Total:", sig_counts$Total),
      vjust = 7.5,
      hjust = 0.5,
      size = 3.8,
      fontface = "bold"
    )
  
  # Save as PDF with higher resolution
  ggsave(file.path(base_output_path, "figures", "DEG_summary_barplot.pdf"), 
         plot = p, width = 12, height = 9, dpi = 300)
  
  # Save as PNG with even higher resolution for presentations
  ggsave(file.path(base_output_path, "figures", "DEG_summary_barplot.png"), 
         plot = p, width = 12, height = 9, dpi = 600)
  
  # Also save a version with white background for publications
  p_pub <- p + theme(
    panel.background = element_rect(fill = "white"),
    plot.background = element_rect(fill = "white"),
    panel.grid.major.y = element_line(color = "gray95"),
    panel.border = element_rect(color = "black", size = 0.5, fill = NA)
  )
  ggsave(file.path(base_output_path, "figures", "DEG_summary_barplot_publication.pdf"), 
         plot = p_pub, width = 12, height = 9, dpi = 300)
  
  cat("Created enhanced summary barplot\n")
}

# Function to run Gene Ontology enrichment analysis - Updated to use params
run_go_analysis <- function(results_list, go_output_path, params) {
  cat("Running Gene Ontology enrichment analysis with direction separation...\n")
  cat("  Using padj threshold:", params$padj_threshold, "\n")
  cat("  Using log2FC threshold:", params$log2fc_threshold, "\n")
  cat("  GO p-value cutoff:", params$go_pvalue_cutoff, "\n")
  cat("  GO q-value cutoff:", params$go_qvalue_cutoff, "\n")
  
  # Add debug info
  cat("  Number of conditions to analyze:", length(results_list), "\n")
  cat("  Conditions:", paste(names(results_list), collapse = ", "), "\n")
  
  for (cond in names(results_list)) {
    cat("  Processing GO analysis for:", cond, "vs SCR\n")
    res_df <- results_list[[cond]]
    
    # Separate up-regulated and down-regulated genes using centralized thresholds
    up_genes <- res_df[!is.na(res_df$padj) & res_df$padj < params$padj_threshold & 
                      res_df$log2FoldChange > params$log2fc_threshold, ]
    down_genes <- res_df[!is.na(res_df$padj) & res_df$padj < params$padj_threshold & 
                        res_df$log2FoldChange < -params$log2fc_threshold, ]
    
    # Create dirs for up and down regulated genes
    up_dir <- file.path(go_output_path, paste0(cond, "_upregulated"))
    down_dir <- file.path(go_output_path, paste0(cond, "_downregulated"))
    dir.create(up_dir, recursive = TRUE, showWarnings = FALSE)
    dir.create(down_dir, recursive = TRUE, showWarnings = FALSE)
    
    # Process up-regulated genes
    if (nrow(up_genes) >= 10) {
      cat("    Processing", nrow(up_genes), "up-regulated genes\n")
      
      # Get Entrez IDs of up-regulated genes (removing NAs)
      up_entrez <- up_genes$entrez[!is.na(up_genes$entrez)]
      
      if (length(up_entrez) >= 10) {
        # Save the gene list used for GO analysis
        gene_list_file <- file.path(up_dir, "upregulated_genes_for_GO_analysis.txt")
        gene_data <- data.frame(
          EnsemblID = rownames(up_genes),
          GeneSymbol = up_genes$gene_symbol,
          EntrezID = up_genes$entrez,
          log2FoldChange = up_genes$log2FoldChange,
          padj = up_genes$padj
        )
        write.table(gene_data, gene_list_file, sep="\t", quote=FALSE, row.names=FALSE)
        
        # Run GO enrichment for up-regulated genes
        process_go_categories(up_entrez, res_df, cond, "upregulated", up_dir, params)
      } else {
        cat("    Not enough mapped Entrez IDs for up-regulated genes in", cond, "(need at least 10)\n")
      }
    } else {
      cat("    Not enough up-regulated genes for", cond, "(need at least 10)\n")
    }
    
    # Process down-regulated genes
    if (nrow(down_genes) >= 10) {
      cat("    Processing", nrow(down_genes), "down-regulated genes\n")
      
      # Get Entrez IDs of down-regulated genes (removing NAs)
      down_entrez <- down_genes$entrez[!is.na(down_genes$entrez)]
      
      if (length(down_entrez) >= 10) {
        # Save the gene list used for GO analysis
        gene_list_file <- file.path(down_dir, "downregulated_genes_for_GO_analysis.txt")
        gene_data <- data.frame(
          EnsemblID = rownames(down_genes),
          GeneSymbol = down_genes$gene_symbol,
          EntrezID = down_genes$entrez,
          log2FoldChange = down_genes$log2FoldChange,
          padj = down_genes$padj
        )
        write.table(gene_data, gene_list_file, sep="\t", quote=FALSE, row.names=FALSE)
        
        # Run GO enrichment for down-regulated genes
        process_go_categories(down_entrez, res_df, cond, "downregulated", down_dir, params)
      } else {
        cat("    Not enough mapped Entrez IDs for down-regulated genes in", cond, "(need at least 10)\n")
      }
    } else {
      cat("    Not enough down-regulated genes for", cond, "(need at least 10)\n")
    }
    
    # Original combined analysis for backward compatibility
    sig_genes <- res_df[!is.na(res_df$padj) & res_df$padj < params$padj_threshold & 
                       abs(res_df$log2FoldChange) > params$log2fc_threshold, ]
    
    # Skip if not enough significant genes
    if (nrow(sig_genes) < 10) {
      cat("    Not enough significant genes for", cond, "(need at least 10), skipping combined GO analysis\n")
      next
    }
    
    # Get Entrez IDs of significant genes (removing NAs)
    sig_entrez <- sig_genes$entrez[!is.na(sig_genes$entrez)]
    
    if (length(sig_entrez) < 10) {
      cat("    Not enough mapped Entrez IDs for", cond, "(need at least 10), skipping combined GO analysis\n")
      next
    }
    
    cat("    Running combined GO analysis with", length(sig_entrez), "genes\n")
    
    # Save the gene list used for GO analysis
    gene_list_file <- file.path(go_output_path, paste0(cond, "_genes_for_GO_analysis.txt"))
    gene_data <- data.frame(
      EnsemblID = rownames(sig_genes),
      GeneSymbol = sig_genes$gene_symbol,
      EntrezID = sig_genes$entrez,
      log2FoldChange = sig_genes$log2FoldChange,
      padj = sig_genes$padj
    )
    write.table(gene_data, gene_list_file, sep="\t", quote=FALSE, row.names=FALSE)
    
    # Run GO for combined genes (up and down)
    process_go_categories(sig_entrez, res_df, cond, "combined", go_output_path, params)
  }
  
  cat("Gene Ontology analysis completed with direction separation\n")
}

# Helper function to process GO categories - Updated to use params
process_go_categories <- function(entrez_ids, res_df, cond, direction, output_path, params) {
  # Run GO enrichment analyses for BP, MF, and CC categories
  go_categories <- c("BP", "MF", "CC")
  go_results <- list()
  
  for (cat in go_categories) {
    tryCatch({
      go_res <- enrichGO(
        gene = entrez_ids,
        universe = res_df$entrez[!is.na(res_df$entrez)],
        OrgDb = org.Hs.eg.db,
        ont = cat,
        pAdjustMethod = "BH",
        pvalueCutoff = params$go_pvalue_cutoff,
        qvalueCutoff = params$go_qvalue_cutoff,
        readable = TRUE
      )
      
      # Simplify GO terms to reduce redundancy
      go_res_simp <- simplify(go_res, cutoff = 0.7, by = "p.adjust", select_fun = min)
      
      go_results[[cat]] <- go_res_simp
      
      cat("    Found", nrow(go_res_simp), "enriched", cat, "terms for", direction, "\n")
      
      if (nrow(go_res_simp) > 0) {
        # Save the results as CSV
        write.csv(
          as.data.frame(go_res_simp),
          file.path(output_path, paste0(cond, "_", direction, "_", cat, "_GO_enrichment.csv")),
          row.names = FALSE
        )
        
        # Generate plots
        # Bar plot
        p1 <- barplot(go_res_simp, showCategory = 15, 
                      title = paste0(cond, " - ", direction, " ", cat, " GO Terms"))
        ggsave(
          file.path(output_path, paste0(cond, "_", direction, "_", cat, "_GO_barplot.pdf")),
          plot = p1, width = 10, height = 8
        )
        
        # Dot plot
        p2 <- dotplot(go_res_simp, showCategory = 15, 
                      title = paste0(cond, " - ", direction, " ", cat, " GO Terms"))
        ggsave(
          file.path(output_path, paste0(cond, "_", direction, "_", cat, "_GO_dotplot.pdf")),
          plot = p2, width = 10, height = 8
        )
        
        # Enrichment map
        if (nrow(go_res_simp) >= 5) {
          p3 <- emapplot(go_res_simp, showCategory = min(20, nrow(go_res_simp)))
          ggsave(
            file.path(output_path, paste0(cond, "_", direction, "_", cat, "_GO_enrichment_map.pdf")),
            plot = p3, width = 12, height = 10
          )
        }
        
        # Gene-Concept Network
        p4 <- cnetplot(go_res_simp, showCategory = min(10, nrow(go_res_simp)))
        ggsave(
          file.path(output_path, paste0(cond, "_", direction, "_", cat, "_GO_cnetplot.pdf")),
          plot = p4, width = 12, height = 10
        )
        
        cat("    Generated GO analysis plots for", cond, direction, cat, "\n")
      } else {
        cat("    No significant GO terms found for", cond, direction, cat, "\n")
      }
    }, error = function(e) {
      cat("    Error in GO analysis for", cond, direction, cat, ":", e$message, "\n")
    })
  }
  
  # Run KEGG pathway analysis if we have enough genes
  if (length(entrez_ids) >= 10) {
    tryCatch({
      cat("    Running KEGG pathway analysis for", direction, "genes...\n")
      kegg_res <- enrichKEGG(
        gene = entrez_ids,
        universe = res_df$entrez[!is.na(res_df$entrez)],
        organism = "hsa",
        pvalueCutoff = params$go_pvalue_cutoff,
        pAdjustMethod = "BH"
      )
      
      if (!is.null(kegg_res) && nrow(kegg_res) > 0) {
        cat("    Found", nrow(kegg_res), "enriched KEGG pathways for", direction, "\n")
        
        # Save the results as CSV
        write.csv(
          as.data.frame(kegg_res),
          file.path(output_path, paste0(cond, "_", direction, "_KEGG_pathway_enrichment.csv")),
          row.names = FALSE
        )
        
        # Generate plots
        # Bar plot
        p1 <- barplot(kegg_res, showCategory = 15, 
                      title = paste0(cond, " - ", direction, " KEGG Pathways"))
        ggsave(
          file.path(output_path, paste0(cond, "_", direction, "_KEGG_barplot.pdf")),
          plot = p1, width = 10, height = 8
        )
        
        # Dot plot
        p2 <- dotplot(kegg_res, showCategory = 15, 
                      title = paste0(cond, " - ", direction, " KEGG Pathways"))
        ggsave(
          file.path(output_path, paste0(cond, "_", direction, "_KEGG_dotplot.pdf")),
          plot = p2, width = 10, height = 8
        )
        
        cat("    Generated KEGG pathway analysis plots for", cond, direction, "\n")
      } else {
        cat("    No significant KEGG pathways found for", cond, direction, "\n")
      }
    }, error = function(e) {
      cat("    Error in KEGG analysis for", cond, direction, ":", e$message, "\n")
    })
  }
}

# Function to generate enhanced MA plots - Updated to use params
generate_ma_plots <- function(dds, results_list, ma_output_path, params) {
  cat("Generating enhanced MA plots with off-target gene labeling...\n")
  cat("  Using padj threshold:", params$padj_threshold, "\n")
  cat("  Y-axis limits:", params$ma_plot_ylim, "\n")
  
  for (cond in names(results_list)) {
    cat("  Creating MA plot for", cond, "vs SCR\n")
    res_obj <- results(dds, contrast = c("Condition", cond, "SCR"))
    
    # Convert to data frame for ggplot
    res_df <- as.data.frame(res_obj)
    
    # Add gene symbols
    res_df$gene_symbol <- results_list[[cond]]$gene_symbol[match(rownames(res_df), rownames(results_list[[cond]]))]
    
    # Add significance categories for coloring using centralized threshold
    res_df$Significance <- "Not Significant"
    res_df$Significance[!is.na(res_df$padj) & res_df$padj < params$padj_threshold] <- "Significant"
    
    # Mark SCN1A
    res_df$is_SCN1A <- results_list[[cond]]$is_SCN1A[match(rownames(res_df), rownames(results_list[[cond]]))]
    
    # Add flag for off-target genes based on the current condition
    res_df$is_OffTarget <- FALSE
    if (cond %in% names(params$off_target_genes)) {
      # Mark the specific off-target genes for this condition
      res_df$is_OffTarget <- res_df$gene_symbol %in% params$off_target_genes[[cond]]
      # Log similar sequence genes found
      off_target_found <- res_df$gene_symbol[res_df$is_OffTarget]
      if (length(off_target_found) > 0) {
        cat("    Similar sequence genes found:", paste(off_target_found, collapse=", "), "\n")
      } else {
        cat("    No similar sequence genes found in the dataset\n")
      }
    }
    
    # Debugging info
    cat("    Total genes:", nrow(res_df), "\n")
    cat("    Significant genes (padj <", params$padj_threshold, "):", 
        sum(res_df$Significance == "Significant", na.rm = TRUE), "\n")
    cat("    SCN1A present:", sum(res_df$is_SCN1A, na.rm = TRUE), "\n")
    
    # Select top genes for labeling
    top_genes <- res_df %>%
      filter(Significance == "Significant") %>%
      arrange(padj) %>%
      head(10)
    
    # Create ggplot MA plot
    p <- ggplot(res_df, aes(x = log10(baseMean), y = log2FoldChange)) +
      # Add points with color by significance
      geom_point(aes(color = Significance), size = 1.2, alpha = 0.7) +
      # Add horizontal line at y=0
      geom_hline(yintercept = 0, linetype = "solid", color = "darkgray") +
      # Add horizontal lines at log2FC thresholds
      geom_hline(yintercept = c(-params$log2fc_threshold, params$log2fc_threshold), 
                linetype = "dotted", color = c("blue", "red")) +
      # Set colors similar to volcano plot
      scale_color_manual(values = c("Significant" = "#3366CC", "Not Significant" = "grey70")) +
      # Add labels
      labs(
        title = paste0("MA Plot: ", cond, " vs SCR"),
        subtitle = paste0("Total genes: ", nrow(res_df), " | Significant (padj < ", 
                         params$padj_threshold, "): ", 
                         sum(res_df$Significance == "Significant", na.rm = TRUE)),
        x = "log10(Mean Expression)",
        y = "log2(Fold Change)",
        caption = paste0("Purple points: genes with similar sequence to target | ",
                        "Dotted lines: log2FC = ±", params$log2fc_threshold)
      ) +
      # Match theme of volcano plots
      theme_light(base_size = 12) +
      theme(
        panel.grid.major = element_line(color = "white"),
        panel.grid.minor = element_line(color = "white"),
        panel.border = element_rect(color = "grey80", fill = NA),
        panel.background = element_rect(fill = "grey95"),
        plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        plot.caption = element_text(hjust = 1, size = 8, color = "grey50"),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        legend.background = element_rect(fill = "white", color = "grey80"),
        axis.title = element_text(face = "bold", size = 12)
      ) +
      # Set limits using centralized parameters
      ylim(params$ma_plot_ylim)
    
    # Add labels for top genes
    if (nrow(top_genes) > 0) {
      p <- p + geom_text_repel(
        data = top_genes,
        aes(label = gene_symbol),
        size = 3.5,
        fontface = "bold",
        box.padding = 0.5,
        point.padding = 0.3,
        force = 12,
        max.overlaps = 30,
        min.segment.length = 0.2,
        segment.color = "grey50",
        segment.size = 0.3,
        color = "#3366CC"
      )
    }
    
    # Highlight SCN1A if present
    scn1a_data <- res_df[res_df$is_SCN1A, ]
    if (nrow(scn1a_data) > 0) {
      cat("    SCN1A log2FoldChange:", scn1a_data$log2FoldChange, "\n")
      cat("    SCN1A baseMean:", scn1a_data$baseMean, "\n")
      cat("    SCN1A padj:", scn1a_data$padj, "\n")
      
      p <- p + geom_point(
        data = scn1a_data,
        size = 3,
        color = "green4",
        shape = 19
      ) +
      geom_text_repel(
        data = scn1a_data,
        aes(label = "SCN1A"),
        size = 4,
        fontface = "bold",
        box.padding = 1,
        force = 15,
        color = "green4",
        segment.color = "green4",
        segment.size = 0.4
      )
    }
    
    # Highlight off-target genes if present
    if (cond %in% names(params$off_target_genes)) {
      off_target_data <- res_df[res_df$is_OffTarget, ]
      if (nrow(off_target_data) > 0) {
        cat("    Including similar sequence genes in MA plot\n")
        
        p <- p + geom_point(
          data = off_target_data,
          size = 3,
          color = "#9966CC",  # Purple color for similar sequence genes
          shape = 19
        ) +
        geom_text_repel(
          data = off_target_data,
          aes(label = gene_symbol),
          size = 3.5,
          fontface = "bold",
          box.padding = 0.8,
          force = 12,
          color = "#9966CC",
          segment.color = "#9966CC",
          segment.size = 0.4,
          max.overlaps = 30
        )
      }
    }
    
    # Save the plot
    ggsave(file.path(ma_output_path, paste0("MA_plot_", cond, "_vs_SCR.pdf")), 
           plot = p, width = 10, height = 8)
    ggsave(file.path(ma_output_path, paste0("MA_plot_", cond, "_vs_SCR.png")), 
           plot = p, width = 10, height = 8, dpi = 300)
  }
  cat("Created enhanced MA plots with off-target gene labeling\n")
}

# Function to generate heatmaps - Updated to use params
generate_heatmaps <- function(dds, results_list, heatmaps_output_path, params) {
  cat("Generating heatmaps...\n")
  cat("  Using padj threshold:", params$padj_threshold, "\n")
  cat("  Using log2FC threshold for heatmaps:", params$heatmap_log2fc_threshold, "\n")
  cat("  Top genes to show:", params$heatmap_top_genes, "\n")
  
  # Sample distance heatmap
  cat("  Creating sample distance heatmap...\n")
  vsd <- vst(dds, blind = TRUE)
  sample_dists <- dist(t(assay(vsd)))
  sample_dist_matrix <- as.matrix(sample_dists)
  rownames(sample_dist_matrix) <- colnames(vsd)
  colnames(sample_dist_matrix) <- colnames(vsd)
  
  # Create annotation data frame for the heatmap
  annotation_col <- as.data.frame(colData(dds)[, "Condition", drop = FALSE])
  
  # Set colors for the annotation
  ann_colors <- list(Condition = params$saRNA_colors)
  
  pdf(file.path(heatmaps_output_path, "Sample_distance_heatmap.pdf"), width = 10, height = 9)
  pheatmap(sample_dist_matrix, 
           clustering_distance_rows = sample_dists,
           clustering_distance_cols = sample_dists, 
           main = "Sample Distance Heatmap",
           annotation_col = annotation_col,
           annotation_colors = ann_colors,
           fontsize = 10,
           fontsize_row = 8,
           fontsize_col = 8,
           color = colorRampPalette(rev(brewer.pal(9, "Blues")))(255))
  dev.off()
  
  # Generate DEG heatmaps for each comparison
  for (cond in names(results_list)) {
    cat("  Creating DEG heatmap for", cond, "vs SCR...\n")
    res_df <- results_list[[cond]]
    
    # Select differentially expressed genes using centralized thresholds
    sig_genes <- rownames(res_df[!is.na(res_df$padj) & 
                                res_df$padj < params$padj_threshold & 
                                abs(res_df$log2FoldChange) > params$heatmap_log2fc_threshold, ])
    sig_count <- length(sig_genes)
    
    cat("    Found", sig_count, "significant DEGs (padj <", params$padj_threshold, 
        ", |log2FC| >", params$heatmap_log2fc_threshold, ")\n")
    
    # If there are significant genes, create the heatmap
    if (sig_count > 0) {
      # Limit to top genes if more than specified
      if (sig_count > params$heatmap_top_genes) {
        sig_genes <- sig_genes[order(res_df[sig_genes, "padj"])][1:params$heatmap_top_genes]
        cat("    Using top", params$heatmap_top_genes, "genes for heatmap\n")
      }
      
      # Extract normalized counts for these genes
      top_gene_counts <- assay(vsd)[sig_genes, ]
      
      # Scale the rows for better visualization
      top_gene_counts_scaled <- t(scale(t(top_gene_counts)))
      
      # Create annotation data frame
      sample_info <- as.data.frame(colData(dds)[, "Condition", drop = FALSE])
      
      # Get gene symbols for rownames
      gene_labels <- res_df$gene_symbol[match(rownames(top_gene_counts_scaled), rownames(res_df))]
      rownames(top_gene_counts_scaled) <- gene_labels
      
      # Check for SCN1A in the selected genes
      scn1a_present <- "SCN1A" %in% gene_labels
      if (scn1a_present) {
        cat("    SCN1A is among the top DEGs\n")
      }
      
      # Set colors for the annotation
      ann_colors <- list(Condition = params$saRNA_colors)
      
      # Create heatmap
      pdf(file.path(heatmaps_output_path, paste0("Top_DEGs_heatmap_", cond, "_vs_SCR.pdf")), 
          width = 12, height = 14)
      pheatmap(top_gene_counts_scaled,
              annotation_col = sample_info,
              annotation_colors = ann_colors,
              main = paste0("Top ", length(sig_genes), " Differentially Expressed Genes: ", 
                           cond, " vs SCR\n",
                           "(padj < ", params$padj_threshold, ", |log2FC| > ", 
                           params$heatmap_log2fc_threshold, ")"),
              cluster_rows = TRUE,
              cluster_cols = TRUE,
              show_rownames = TRUE,
              fontsize = 12,
              fontsize_row = 9,
              fontsize_col = 10,
              color = colorRampPalette(rev(brewer.pal(9, "RdBu")))(255),
              border_color = NA,
              clustering_method = "ward.D2")
      dev.off()
      
      # Save the gene list used in the heatmap
      gene_list_file <- file.path(heatmaps_output_path, paste0("Top_DEGs_", cond, "_vs_SCR.txt"))
      gene_data <- data.frame(
        EnsemblID = sig_genes,
        GeneSymbol = gene_labels,
        log2FoldChange = res_df[sig_genes, "log2FoldChange"],
        pvalue = res_df[sig_genes, "pvalue"],
        padj = res_df[sig_genes, "padj"]
      )
      write.table(gene_data, gene_list_file, sep="\t", quote=FALSE, row.names=FALSE)
      cat("    Gene list saved to:", gene_list_file, "\n")
    } else {
      cat("    No significant genes found for", cond, ", skipping heatmap\n")
    }
  }
  cat("Created heatmaps\n")
}

# Enhanced function to generate Venn diagrams comparing DEGs across conditions
generate_venn_diagrams <- function(results_list, venn_output_path, params) {
  cat("Generating Venn diagrams for DEG comparisons...\n")
  cat("  Using padj threshold:", params$padj_threshold, "\n")
  cat("  Using log2FC threshold:", params$log2fc_threshold, "\n")
  
  # Create output directory
  dir.create(venn_output_path, recursive = TRUE, showWarnings = FALSE)
  
  # Extract DEG lists for each condition - now returns gene symbols instead of IDs
  get_degs <- function(res_df, direction = "all", return_symbols = TRUE) {
    if (direction == "up") {
      degs_idx <- which(!is.na(res_df$padj) & 
                       res_df$padj < params$padj_threshold & 
                       res_df$log2FoldChange > params$log2fc_threshold)
    } else if (direction == "down") {
      degs_idx <- which(!is.na(res_df$padj) & 
                       res_df$padj < params$padj_threshold & 
                       res_df$log2FoldChange < -params$log2fc_threshold)
    } else {  # all
      degs_idx <- which(!is.na(res_df$padj) & 
                       res_df$padj < params$padj_threshold & 
                       abs(res_df$log2FoldChange) > params$log2fc_threshold)
    }
    
    if (return_symbols) {
      # Return gene symbols for better readability
      return(res_df$gene_symbol[degs_idx])
    } else {
      # Return Ensembl IDs
      return(rownames(res_df)[degs_idx])
    }
  }
  
  # Compare saRNA10, saRNA12, and saRNA14 (conditions that upregulate SCN1A)
  conditions_to_compare <- c("saRNA10", "saRNA12", "saRNA14")
  
  # Create color palette with better colors
  venn_colors <- c("#FF6B6B", "#4ECDC4", "#45B7D1")  # Coral, Turquoise, Sky Blue
  
  # 1. All DEGs comparison
  deg_lists_all <- list()
  deg_lists_all_ids <- list()  # Keep IDs for data processing
  for (cond in conditions_to_compare) {
    deg_lists_all[[cond]] <- get_degs(results_list[[cond]], "all", return_symbols = TRUE)
    deg_lists_all_ids[[cond]] <- get_degs(results_list[[cond]], "all", return_symbols = FALSE)
  }
  
  # Create enhanced Venn diagram for all DEGs
  venn.plot <- venn.diagram(
    x = deg_lists_all,
    category.names = conditions_to_compare,
    filename = NULL,
    output = TRUE,
    
    # Colors with transparency
    fill = venn_colors,
    alpha = 0.5,
    
    # Numbers styling
    cex = 2,
    fontface = "bold",
    fontfamily = "sans",
    
    # Category names styling
    cat.cex = 2,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-20, 20, 180),
    cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans",
    cat.col = venn_colors,
    
    # Border styling
    lwd = 2,
    lty = 'solid',
    
    # Title
    main = paste0("All DEGs (padj < ", params$padj_threshold, ", |log2FC| > ", params$log2fc_threshold, ")"),
    main.cex = 2.5,
    main.fontface = "bold",
    main.fontfamily = "sans",
    
    # Layout
    margin = 0.1,
    
    # Additional styling
    ext.text = TRUE,
    ext.line.lwd = 2,
    ext.dist = -0.15,
    ext.length = 0.85,
    ext.pos = 30
  )
  
  # Save the plot
  pdf(file.path(venn_output_path, "Venn_all_DEGs_saRNA10_12_14.pdf"), width = 12, height = 12)
  grid.draw(venn.plot)
  dev.off()
  
  # Also save as PNG with high resolution
  png(file.path(venn_output_path, "Venn_all_DEGs_saRNA10_12_14.png"), 
      width = 3600, height = 3600, res = 300)
  grid.draw(venn.plot)
  dev.off()
  
  # 2. Up-regulated genes only
  deg_lists_up <- list()
  deg_lists_up_ids <- list()
  for (cond in conditions_to_compare) {
    deg_lists_up[[cond]] <- get_degs(results_list[[cond]], "up", return_symbols = TRUE)
    deg_lists_up_ids[[cond]] <- get_degs(results_list[[cond]], "up", return_symbols = FALSE)
  }
  
  venn.plot.up <- venn.diagram(
    x = deg_lists_up,
    category.names = conditions_to_compare,
    filename = NULL,
    output = TRUE,
    
    # Use warmer colors for up-regulated
    fill = c("#FF6B6B", "#FFD93D", "#6BCB77"),  # Coral, Yellow, Green
    alpha = 0.5,
    
    # Numbers styling
    cex = 2,
    fontface = "bold",
    fontfamily = "sans",
    
    # Category names styling
    cat.cex = 2,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-20, 20, 180),
    cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans",
    cat.col = c("#FF6B6B", "#FFD93D", "#6BCB77"),
    
    # Border styling
    lwd = 2,
    lty = 'solid',
    
    # Title
    main = paste0("Up-regulated DEGs (padj < ", params$padj_threshold, ", log2FC > ", params$log2fc_threshold, ")"),
    main.cex = 2.5,
    main.fontface = "bold",
    
    margin = 0.1
  )
  
  pdf(file.path(venn_output_path, "Venn_upregulated_DEGs_saRNA10_12_14.pdf"), width = 12, height = 12)
  grid.draw(venn.plot.up)
  dev.off()
  
  png(file.path(venn_output_path, "Venn_upregulated_DEGs_saRNA10_12_14.png"), 
      width = 3600, height = 3600, res = 300)
  grid.draw(venn.plot.up)
  dev.off()
  
  # 3. Down-regulated genes only
  deg_lists_down <- list()
  deg_lists_down_ids <- list()
  for (cond in conditions_to_compare) {
    deg_lists_down[[cond]] <- get_degs(results_list[[cond]], "down", return_symbols = TRUE)
    deg_lists_down_ids[[cond]] <- get_degs(results_list[[cond]], "down", return_symbols = FALSE)
  }
  
  venn.plot.down <- venn.diagram(
    x = deg_lists_down,
    category.names = conditions_to_compare,
    filename = NULL,
    output = TRUE,
    
    # Use cooler colors for down-regulated
    fill = c("#4B89DC", "#9B59B6", "#5D9CEC"),  # Blues and purple
    alpha = 0.5,
    
    # Numbers styling
    cex = 2,
    fontface = "bold",
    fontfamily = "sans",
    
    # Category names styling
    cat.cex = 2,
    cat.fontface = "bold",
    cat.default.pos = "outer",
    cat.pos = c(-20, 20, 180),
    cat.dist = c(0.055, 0.055, 0.085),
    cat.fontfamily = "sans",
    cat.col = c("#4B89DC", "#9B59B6", "#5D9CEC"),
    
    # Border styling
    lwd = 2,
    lty = 'solid',
    
    # Title
    main = paste0("Down-regulated DEGs (padj < ", params$padj_threshold, ", log2FC < -", params$log2fc_threshold, ")"),
    main.cex = 2.5,
    main.fontface = "bold",
    
    margin = 0.1
  )
  
  pdf(file.path(venn_output_path, "Venn_downregulated_DEGs_saRNA10_12_14.pdf"), width = 12, height = 12)
  grid.draw(venn.plot.down)
  dev.off()
  
  png(file.path(venn_output_path, "Venn_downregulated_DEGs_saRNA10_12_14.png"), 
      width = 3600, height = 3600, res = 300)
  grid.draw(venn.plot.down)
  dev.off()
  
  # Calculate overlaps and save gene lists with enhanced formatting
  cat("\n  Venn diagram statistics:\n")
  
  # Enhanced save function that creates well-formatted Excel-like CSV files
  save_gene_list_enhanced <- function(gene_symbols, gene_ids, filename, results_list, description) {
    if (length(gene_symbols) == 0) {
      cat("    No genes to save for:", description, "\n")
      return()
    }
    
    # Create data frame with gene information
    gene_data <- data.frame(
      GeneSymbol = gene_symbols,
      stringsAsFactors = FALSE
    )
    
    # Add Ensembl IDs if available
    if (!is.null(gene_ids) && length(gene_ids) == length(gene_symbols)) {
      gene_data$EnsemblID <- gene_ids
    }
    
    # Add expression data from each condition
    for (cond in conditions_to_compare) {
      res_df <- results_list[[cond]]
      
      # Match by gene symbol
      match_idx <- match(gene_symbols, res_df$gene_symbol)
      
      gene_data[[paste0(cond, "_log2FC")]] <- round(res_df$log2FoldChange[match_idx], 3)
      gene_data[[paste0(cond, "_padj")]] <- format(res_df$padj[match_idx], scientific = TRUE, digits = 3)
      
      # Add significance indicator
      sig_indicator <- rep("", length(gene_symbols))
      sig_indicator[!is.na(res_df$padj[match_idx]) & 
                   res_df$padj[match_idx] < params$padj_threshold & 
                   abs(res_df$log2FoldChange[match_idx]) > params$log2fc_threshold] <- "*"
      gene_data[[paste0(cond, "_sig")]] <- sig_indicator
    }
    
    # Calculate mean log2FC across conditions
    log2fc_cols <- grep("_log2FC$", colnames(gene_data))
    gene_data$Mean_log2FC <- round(rowMeans(gene_data[, log2fc_cols], na.rm = TRUE), 3)
    
    # Sort by mean log2FC (descending)
    gene_data <- gene_data[order(abs(gene_data$Mean_log2FC), decreasing = TRUE), ]
    
    # Add rank
    gene_data <- cbind(Rank = 1:nrow(gene_data), gene_data)
    
    # Write to file with description header
    output_file <- file.path(venn_output_path, filename)
    
    # Write header information
    cat("# ", description, "\n", file = output_file)
    cat("# Analysis date: ", format(Sys.Date(), "%Y-%m-%d"), "\n", file = output_file, append = TRUE)
    cat("# Total genes: ", nrow(gene_data), "\n", file = output_file, append = TRUE)
    cat("# Significance threshold: padj < ", params$padj_threshold, ", |log2FC| > ", params$log2fc_threshold, "\n", file = output_file, append = TRUE)
    cat("# * indicates significant in that condition\n", file = output_file, append = TRUE)
    cat("#\n", file = output_file, append = TRUE)
    
    # Write the data
    write.table(gene_data, output_file, sep = "\t", quote = FALSE, row.names = FALSE, append = TRUE)
    
    cat("    Saved", nrow(gene_data), "genes to:", basename(filename), "\n")
    
    # Print top 5 genes
    if (nrow(gene_data) > 0) {
      cat("      Top genes:", paste(head(gene_data$GeneSymbol, 5), collapse = ", "), "\n")
    }
  }
  
  # Process all DEGs overlaps
  all_overlap <- Reduce(intersect, deg_lists_all)
  all_overlap_ids <- Reduce(intersect, deg_lists_all_ids)
  cat("    Common DEGs in all three conditions:", length(all_overlap), "\n")
  
  # Check if SCN1A is in the common list
  if ("SCN1A" %in% all_overlap) {
    cat("    *** SCN1A is among the common DEGs! ***\n")
  }
  
  # Save common genes
  save_gene_list_enhanced(all_overlap, all_overlap_ids, 
                         "Common_DEGs_all_three_conditions.txt", 
                         results_list,
                         "Genes differentially expressed in all three conditions (saRNA10, saRNA12, saRNA14)")
  
  # Pairwise overlaps
  overlap_10_12 <- intersect(deg_lists_all$saRNA10, deg_lists_all$saRNA12)
  overlap_10_12_ids <- intersect(deg_lists_all_ids$saRNA10, deg_lists_all_ids$saRNA12)
  
  overlap_10_14 <- intersect(deg_lists_all$saRNA10, deg_lists_all$saRNA14)
  overlap_10_14_ids <- intersect(deg_lists_all_ids$saRNA10, deg_lists_all_ids$saRNA14)
  
  overlap_12_14 <- intersect(deg_lists_all$saRNA12, deg_lists_all$saRNA14)
  overlap_12_14_ids <- intersect(deg_lists_all_ids$saRNA12, deg_lists_all_ids$saRNA14)
  
  cat("    saRNA10 ∩ saRNA12:", length(overlap_10_12), "\n")
  cat("    saRNA10 ∩ saRNA14:", length(overlap_10_14), "\n")
  cat("    saRNA12 ∩ saRNA14:", length(overlap_12_14), "\n")
  
  # Save pairwise overlaps (excluding genes in all three)
  overlap_10_12_only <- setdiff(overlap_10_12, all_overlap)
  overlap_10_12_only_ids <- setdiff(overlap_10_12_ids, all_overlap_ids)
  save_gene_list_enhanced(overlap_10_12_only, overlap_10_12_only_ids,
                         "Common_DEGs_saRNA10_saRNA12_only.txt", 
                         results_list,
                         "Genes differentially expressed in saRNA10 and saRNA12 (but not saRNA14)")
  
  overlap_10_14_only <- setdiff(overlap_10_14, all_overlap)
  overlap_10_14_only_ids <- setdiff(overlap_10_14_ids, all_overlap_ids)
  save_gene_list_enhanced(overlap_10_14_only, overlap_10_14_only_ids,
                         "Common_DEGs_saRNA10_saRNA14_only.txt", 
                         results_list,
                         "Genes differentially expressed in saRNA10 and saRNA14 (but not saRNA12)")
  
  overlap_12_14_only <- setdiff(overlap_12_14, all_overlap)
  overlap_12_14_only_ids <- setdiff(overlap_12_14_ids, all_overlap_ids)
  save_gene_list_enhanced(overlap_12_14_only, overlap_12_14_only_ids,
                         "Common_DEGs_saRNA12_saRNA14_only.txt", 
                         results_list,
                         "Genes differentially expressed in saRNA12 and saRNA14 (but not saRNA10)")
  
  # Save unique genes for each condition
  unique_10 <- setdiff(deg_lists_all$saRNA10, union(deg_lists_all$saRNA12, deg_lists_all$saRNA14))
  unique_10_ids <- setdiff(deg_lists_all_ids$saRNA10, union(deg_lists_all_ids$saRNA12, deg_lists_all_ids$saRNA14))
  
  unique_12 <- setdiff(deg_lists_all$saRNA12, union(deg_lists_all$saRNA10, deg_lists_all$saRNA14))
  unique_12_ids <- setdiff(deg_lists_all_ids$saRNA12, union(deg_lists_all_ids$saRNA10, deg_lists_all_ids$saRNA14))
  
  unique_14 <- setdiff(deg_lists_all$saRNA14, union(deg_lists_all$saRNA10, deg_lists_all$saRNA12))
  unique_14_ids <- setdiff(deg_lists_all_ids$saRNA14, union(deg_lists_all_ids$saRNA10, deg_lists_all_ids$saRNA12))
  
  cat("    Unique to saRNA10:", length(unique_10), "\n")
  cat("    Unique to saRNA12:", length(unique_12), "\n")
  cat("    Unique to saRNA14:", length(unique_14), "\n")
  
  save_gene_list_enhanced(unique_10, unique_10_ids,
                         "Unique_DEGs_saRNA10_only.txt", 
                         results_list,
                         "Genes differentially expressed only in saRNA10")
  
  save_gene_list_enhanced(unique_12, unique_12_ids,
                         "Unique_DEGs_saRNA12_only.txt", 
                         results_list,
                         "Genes differentially expressed only in saRNA12")
  
  save_gene_list_enhanced(unique_14, unique_14_ids,
                         "Unique_DEGs_saRNA14_only.txt", 
                         results_list,
                         "Genes differentially expressed only in saRNA14")
  
  # Create enhanced summary plot using ggplot2
  library(ggplot2)
  library(tidyr)
  
  # Prepare data for enhanced bar plot
  venn_summary <- data.frame(
    Category = c("All three conditions", 
                 "saRNA10 & saRNA12 only", 
                 "saRNA10 & saRNA14 only", 
                 "saRNA12 & saRNA14 only",
                 "saRNA10 only", 
                 "saRNA12 only", 
                 "saRNA14 only"),
    Count = c(length(all_overlap),
              length(overlap_10_12_only),
              length(overlap_10_14_only),
              length(overlap_12_14_only),
              length(unique_10), 
              length(unique_12), 
              length(unique_14)),
    Type = c("Common", "Pairwise", "Pairwise", "Pairwise", 
             "Unique", "Unique", "Unique")
  )
  
  # Set factor levels for ordering
  venn_summary$Category <- factor(venn_summary$Category, 
                                 levels = rev(venn_summary$Category))
  
  # Create enhanced bar plot with color coding
  p_summary <- ggplot(venn_summary, aes(x = Category, y = Count, fill = Type)) +
    geom_bar(stat = "identity", width = 0.7) +
    geom_text(aes(label = Count), hjust = -0.2, size = 5, fontface = "bold") +
    coord_flip() +
    scale_fill_manual(values = c("Common" = "#2ECC71", 
                                "Pairwise" = "#3498DB", 
                                "Unique" = "#E74C3C")) +
    labs(
      title = "DEG Overlap Summary",
      subtitle = paste0("Comparing saRNA10, saRNA12, and saRNA14\n",
                       "(padj < ", params$padj_threshold, ", |log2FC| > ", params$log2fc_threshold, ")"),
      x = "",
      y = "Number of Genes",
      fill = "Overlap Type"
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 20),
      plot.subtitle = element_text(hjust = 0.5, size = 14, color = "gray30"),
      axis.text.y = element_text(size = 12, face = "bold"),
      axis.text.x = element_text(size = 12),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      plot.margin = margin(10, 20, 10, 10)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.15)),
                      breaks = scales::pretty_breaks(n = 8))
  
  ggsave(file.path(venn_output_path, "Venn_summary_barplot.pdf"), 
         plot = p_summary, width = 12, height = 10)
  ggsave(file.path(venn_output_path, "Venn_summary_barplot.png"), 
         plot = p_summary, width = 12, height = 10, dpi = 300)
  
  # Create a comprehensive summary table
  summary_table <- data.frame(
    Comparison = c("Total DEGs", "Up-regulated", "Down-regulated"),
    saRNA10 = c(length(deg_lists_all$saRNA10), 
                length(deg_lists_up$saRNA10), 
                length(deg_lists_down$saRNA10)),
    saRNA12 = c(length(deg_lists_all$saRNA12), 
                length(deg_lists_up$saRNA12), 
                length(deg_lists_down$saRNA12)),
    saRNA14 = c(length(deg_lists_all$saRNA14), 
                length(deg_lists_up$saRNA14), 
                length(deg_lists_down$saRNA14)),
    All_Three = c(length(all_overlap), 
                  length(Reduce(intersect, deg_lists_up)), 
                  length(Reduce(intersect, deg_lists_down)))
  )
  
  write.table(summary_table, 
             file.path(venn_output_path, "Venn_summary_statistics.txt"),
             sep = "\t", quote = FALSE, row.names = FALSE)
  
  # Save separate lists for up and down regulated overlaps
  cat("\n  Saving directional overlap lists...\n")
  
  # Up-regulated overlaps
  up_all_overlap <- Reduce(intersect, deg_lists_up)
  if (length(up_all_overlap) > 0) {
    up_all_overlap_ids <- rownames(results_list$saRNA10)[match(up_all_overlap, results_list$saRNA10$gene_symbol)]
    save_gene_list_enhanced(up_all_overlap, up_all_overlap_ids,
                           "Common_upregulated_all_three_conditions.txt", 
                           results_list,
                           "Genes up-regulated in all three conditions")
  }
  
  # Down-regulated overlaps
  down_all_overlap <- Reduce(intersect, deg_lists_down)
  if (length(down_all_overlap) > 0) {
    down_all_overlap_ids <- rownames(results_list$saRNA10)[match(down_all_overlap, results_list$saRNA10$gene_symbol)]
    save_gene_list_enhanced(down_all_overlap, down_all_overlap_ids,
                           "Common_downregulated_all_three_conditions.txt", 
                           results_list,
                           "Genes down-regulated in all three conditions")
  }
  
  cat("\nCreated enhanced Venn diagrams with comprehensive gene lists\n")
  cat("All gene lists include gene symbols for easy identification\n")
}

# Function to generate lncRNA expression scatter plots - Updated to use params
generate_lncrna_scatter_plots <- function(dds, results_list, lncrna_output_path, params) {
  cat("Generating long non-coding RNA (lncRNA) expression scatter plots...\n")
  cat("  Using padj threshold:", params$padj_threshold, "\n")
  
  # Create output directory for lncRNA plots
  dir.create(lncrna_output_path, recursive = TRUE, showWarnings = FALSE)
  
  # Get normalized expression data
  vsd <- vst(dds, blind = TRUE)
  normalized_counts <- assay(vsd)
  
  # Function to identify lncRNAs based on gene symbols
  identify_lncrnas <- function(gene_symbols) {
    # Common lncRNA patterns
    lncrna_patterns <- c(
      "-AS[0-9]*$",      # Antisense RNAs (e.g., SCN1A-AS1)
      "^LINC[0-9]+",     # Long intergenic non-coding RNAs
      "^XIST$",          # X-inactive specific transcript
      "^TSIX$",          # TSIX transcript
      "^HOTAIR$",        # HOX transcript antisense RNA
      "^MALAT1$",        # Metastasis associated lung adenocarcinoma transcript 1
      "^NEAT1$",         # Nuclear enriched abundant transcript 1
      "^H19$",           # H19 imprinted maternally expressed transcript
      "^MEG[0-9]+$",     # Maternally expressed genes
      "^MIR.*HG$",       # MicroRNA host genes
      "^SNHG[0-9]+$",    # Small nucleolar RNA host genes
      "^DANCR$",         # Differentiation antagonizing non-protein coding RNA
      "-IT[0-9]*$",      # Intronic transcript
      "^AC[0-9]+\\.[0-9]+", # Many lncRNAs have these identifiers
      "^AL[0-9]+\\.[0-9]+",
      "^AP[0-9]+\\.[0-9]+",
      "^BX[0-9]+\\.[0-9]+",
      "^CTA-",           # Cancer-testis antigen family (many are lncRNAs)
      "^CTB-",
      "^CTC-",
      "^CTD-",
      "^RP[0-9]+-",      # Many lncRNAs start with RP
      "^ENSG.*\\..*"     # Ensembl IDs without mapped symbols (often lncRNAs)
    )
    
    # Check if gene symbol matches any lncRNA pattern
    is_lncrna <- rep(FALSE, length(gene_symbols))
    for (pattern in lncrna_patterns) {
      is_lncrna <- is_lncrna | grepl(pattern, gene_symbols, ignore.case = FALSE)
    }
    
    return(is_lncrna)
  }
  
  # For each condition comparison, create scatter plot
  for (cond in names(results_list)) {
    cat(paste0("  Creating lncRNA scatter plot for ", cond, " vs SCR\n"))
    
    res_df <- results_list[[cond]]
    
    # Identify lncRNAs
    res_df$is_lncrna <- identify_lncrnas(res_df$gene_symbol)
    
    # Filter for lncRNAs only
    lncrna_df <- res_df[res_df$is_lncrna, ]
    
    cat("    Found", nrow(lncrna_df), "lncRNAs\n")
    
    if (nrow(lncrna_df) < 5) {
      cat("    Too few lncRNAs found, skipping plot for", cond, "\n")
      next
    }
    
    # Get mean expression for each lncRNA across samples
    # First get sample indices for SCR and current condition
    scr_samples <- which(colData(dds)$Condition == "SCR")
    cond_samples <- which(colData(dds)$Condition == cond)
    
    # Calculate mean expression for SCR and condition
    lncrna_expression <- data.frame(
      gene_id = rownames(lncrna_df),
      gene_symbol = lncrna_df$gene_symbol,
      SCR_mean = rowMeans(normalized_counts[rownames(lncrna_df), scr_samples, drop = FALSE]),
      Condition_mean = rowMeans(normalized_counts[rownames(lncrna_df), cond_samples, drop = FALSE]),
      log2FoldChange = lncrna_df$log2FoldChange,
      padj = lncrna_df$padj,
      stringsAsFactors = FALSE
    )
    
    # Add significance categories using centralized threshold
    lncrna_expression$Significance <- "Not Significant"
    lncrna_expression$Significance[!is.na(lncrna_expression$padj) & 
                                  lncrna_expression$padj < params$padj_threshold] <- "Significant"
    
    # Check if SCN1A-AS1 is present
    lncrna_expression$is_SCN1A_AS1 <- grepl("SCN1A-AS", lncrna_expression$gene_symbol)
    scn1a_as1_present <- any(lncrna_expression$is_SCN1A_AS1)
    
    if (scn1a_as1_present) {
      cat("    SCN1A-AS1 found in lncRNA dataset\n")
      scn1a_as1_data <- lncrna_expression[lncrna_expression$is_SCN1A_AS1, ]
      cat("    SCN1A-AS1 details: log2FC =", scn1a_as1_data$log2FoldChange[1], 
          ", padj =", scn1a_as1_data$padj[1], "\n")
    }
    
    # Calculate overall expression (mean of SCR and condition means) for ranking
    lncrna_expression$overall_expression <- (lncrna_expression$SCR_mean + 
                                           lncrna_expression$Condition_mean) / 2
    
    # Select top expressed lncRNAs
    top_expressed <- lncrna_expression %>%
      arrange(desc(overall_expression)) %>%
      head(10)
    
    # Select most differentially expressed lncRNAs
    top_de <- lncrna_expression %>%
      filter(Significance == "Significant") %>%
      arrange(padj) %>%
      head(5)
    
    # Combine for labeling (remove duplicates)
    genes_to_label <- unique(c(top_expressed$gene_id, top_de$gene_id))
    
    # Always include SCN1A-AS1 if present
    if (scn1a_as1_present) {
      genes_to_label <- unique(c(genes_to_label, 
                               lncrna_expression$gene_id[lncrna_expression$is_SCN1A_AS1]))
    }
    
    # Create scatter plot
    p <- ggplot(lncrna_expression, aes(x = SCR_mean, y = Condition_mean)) +
      # Add diagonal reference line (y = x)
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray50") +
      # Add 2-fold change lines
      geom_abline(intercept = 0, slope = 2, linetype = "dotted", color = "red", alpha = 0.5) +
      geom_abline(intercept = 0, slope = 0.5, linetype = "dotted", color = "blue", alpha = 0.5) +
      # Add points
      geom_point(aes(color = Significance), size = 2, alpha = 0.7) +
      # Color scale
      scale_color_manual(values = c("Significant" = "#E41A1C", "Not Significant" = "grey60")) +
      # Labels
      labs(
        title = paste0("lncRNA Expression: ", cond, " vs SCR"),
        subtitle = paste0("Total lncRNAs: ", nrow(lncrna_expression), 
                         " | Significant (padj < ", params$padj_threshold, "): ", 
                         sum(lncrna_expression$Significance == "Significant")),
        x = "Mean Expression in SCR (VST normalized)",
        y = paste0("Mean Expression in ", cond, " (VST normalized)"),
        caption = "Dotted lines: 2-fold change boundaries | Dashed line: y = x"
      ) +
      # Theme
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 12, color = "gray30"),
        plot.caption = element_text(hjust = 1, size = 9, color = "gray50"),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "gray80", fill = NA),
        axis.title = element_text(face = "bold", size = 12)
      ) +
      # Use log scale for better visualization
      scale_x_continuous(trans = "log10", 
                        labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      scale_y_continuous(trans = "log10", 
                        labels = scales::trans_format("log10", scales::math_format(10^.x))) +
      annotation_logticks()
    
    # Add labels for selected genes
    label_data <- lncrna_expression[lncrna_expression$gene_id %in% genes_to_label, ]
    
    if (nrow(label_data) > 0) {
      p <- p + geom_text_repel(
        data = label_data,
        aes(label = gene_symbol),
        size = 3.5,
        box.padding = 0.5,
        point.padding = 0.3,
        force = 10,
        max.overlaps = 30,
        min.segment.length = 0.2,
        segment.color = "grey50",
        segment.size = 0.3,
        color = "black",
        fontface = ifelse(label_data$Significance == "Significant", "bold", "plain")
      )
    }
    
    # Highlight SCN1A-AS1 if present
    if (scn1a_as1_present) {
      scn1a_as1_plot_data <- lncrna_expression[lncrna_expression$is_SCN1A_AS1, ]
      
      p <- p + 
        geom_point(
          data = scn1a_as1_plot_data,
          size = 4,
          color = "green4",
          shape = 19
        ) +
        geom_text_repel(
          data = scn1a_as1_plot_data,
          aes(label = gene_symbol),
          size = 4,
          fontface = "bold",
          box.padding = 1,
          force = 15,
          color = "green4",
          segment.color = "green4",
          segment.size = 0.5
        )
    }
    
    # Save plots
    ggsave(file.path(lncrna_output_path, paste0("lncRNA_scatter_", cond, "_vs_SCR.pdf")), 
           plot = p, width = 10, height = 9, dpi = 300)
    ggsave(file.path(lncrna_output_path, paste0("lncRNA_scatter_", cond, "_vs_SCR.png")), 
           plot = p, width = 10, height = 9, dpi = 600)
    
    # Create a bar plot showing top expressed lncRNAs
    top_20_lncrnas <- lncrna_expression %>%
      arrange(desc(overall_expression)) %>%
      head(20) %>%
      mutate(gene_symbol = factor(gene_symbol, levels = rev(gene_symbol)))
    
    p_bar <- ggplot(top_20_lncrnas, aes(x = gene_symbol, y = overall_expression)) +
      geom_bar(stat = "identity", aes(fill = Significance), width = 0.7) +
      geom_hline(yintercept = median(lncrna_expression$overall_expression), 
                 linetype = "dashed", color = "gray40", alpha = 0.7) +
      coord_flip() +
      scale_fill_manual(values = c("Significant" = "#E41A1C", "Not Significant" = "grey70")) +
      labs(
        title = paste0("Top 20 Expressed lncRNAs: ", cond, " vs SCR"),
        subtitle = paste0("Dashed line: median lncRNA expression | Significance: padj < ", 
                         params$padj_threshold),
        x = "",
        y = "Mean Expression (VST normalized)",
        fill = "DE Status"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
        plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray30"),
        legend.position = "right",
        axis.text.y = element_text(face = ifelse(top_20_lncrnas$is_SCN1A_AS1, "bold", "plain"),
                                  color = ifelse(top_20_lncrnas$is_SCN1A_AS1, "green4", "black")),
        panel.grid.major.y = element_blank(),
        panel.border = element_rect(color = "gray80", fill = NA)
      )
    
    # Save bar plot
    ggsave(file.path(lncrna_output_path, paste0("lncRNA_top_expressed_", cond, "_vs_SCR.pdf")), 
           plot = p_bar, width = 10, height = 8, dpi = 300)
    ggsave(file.path(lncrna_output_path, paste0("lncRNA_top_expressed_", cond, "_vs_SCR.png")), 
           plot = p_bar, width = 10, height = 8, dpi = 600)
    
    # Save lncRNA data
    write.csv(lncrna_expression, 
             file.path(lncrna_output_path, paste0("lncRNA_expression_data_", cond, "_vs_SCR.csv")),
             row.names = FALSE)
    
    # Log summary statistics
    cat("    Top 5 expressed lncRNAs:\n")
    for (i in 1:min(5, nrow(top_expressed))) {
      cat("      ", top_expressed$gene_symbol[i], ": mean expression =", 
          round(top_expressed$overall_expression[i], 2), "\n")
    }
    
    if (sum(lncrna_expression$Significance == "Significant") > 0) {
      cat("    Top differentially expressed lncRNAs:\n")
      sig_lncrnas <- lncrna_expression %>%
        filter(Significance == "Significant") %>%
        arrange(padj) %>%
        head(5)
      
      for (i in 1:nrow(sig_lncrnas)) {
        cat("      ", sig_lncrnas$gene_symbol[i], ": log2FC =", 
            round(sig_lncrnas$log2FoldChange[i], 2), 
            ", padj =", format(sig_lncrnas$padj[i], scientific = TRUE, digits = 3), "\n")
      }
    }
  }
  
  cat("Created lncRNA expression scatter plots\n")
}

# Function to create a detailed analysis report - Updated to use params
create_analysis_report <- function(dds, results_list, paths, params) {
  cat("Creating analysis report...\n")
  
  # Create report file
  report_file <- file.path(paths$base, "RNA_seq_analysis_report.txt")
  
  # Write report header
  cat("RNA-seq Analysis Report\n", file = report_file)
  cat("=====================\n", file = report_file, append = TRUE)
  cat("Date: ", Sys.Date(), "\n\n", file = report_file, append = TRUE)
  
  # Analysis parameters
  cat("Analysis Parameters\n", file = report_file, append = TRUE)
  cat("------------------\n", file = report_file, append = TRUE)
  cat("Adjusted p-value threshold: ", params$padj_threshold, "\n", file = report_file, append = TRUE)
  cat("Log2 fold change threshold: ", params$log2fc_threshold, "\n", file = report_file, append = TRUE)
  cat("Heatmap log2FC threshold: ", params$heatmap_log2fc_threshold, "\n", file = report_file, append = TRUE)
  cat("GO p-value cutoff: ", params$go_pvalue_cutoff, "\n", file = report_file, append = TRUE)
  cat("GO q-value cutoff: ", params$go_qvalue_cutoff, "\n\n", file = report_file, append = TRUE)
  
  # Sample information
  cat("Sample Information\n", file = report_file, append = TRUE)
  cat("-----------------\n", file = report_file, append = TRUE)
  cat("Total samples: ", ncol(dds), "\n", file = report_file, append = TRUE)
  
  # Conditions summary
  condition_counts <- table(dds$Condition)
  cat("Samples by condition:\n", file = report_file, append = TRUE)
  for (cond in names(condition_counts)) {
    cat("  ", cond, ": ", condition_counts[cond], "\n", file = report_file, append = TRUE)
  }
  cat("\n", file = report_file, append = TRUE)
  
  # Analysis summary
  cat("Analysis Summary\n", file = report_file, append = TRUE)
  cat("---------------\n", file = report_file, append = TRUE)
  cat("Total genes analyzed: ", nrow(dds), "\n", file = report_file, append = TRUE)
  
  # Results by condition
  cat("\nDifferential Expression Results\n", file = report_file, append = TRUE)
  cat("-----------------------------\n", file = report_file, append = TRUE)
  
  for (cond in names(results_list)) {
    res_df <- results_list[[cond]]
    
    # Count significant genes using centralized thresholds
    sig_genes <- sum(!is.na(res_df$padj) & res_df$padj < params$padj_threshold, na.rm = TRUE)
    up_genes <- sum(!is.na(res_df$padj) & res_df$padj < params$padj_threshold & 
                   res_df$log2FoldChange > params$log2fc_threshold, na.rm = TRUE)
    down_genes <- sum(!is.na(res_df$padj) & res_df$padj < params$padj_threshold & 
                     res_df$log2FoldChange < -params$log2fc_threshold, na.rm = TRUE)
    
    cat("\n", cond, " vs SCR:\n", file = report_file, append = TRUE)
    cat("  Significant genes (padj < ", params$padj_threshold, "): ", sig_genes, "\n", 
        file = report_file, append = TRUE)
    cat("  Up-regulated genes (padj < ", params$padj_threshold, ", log2FC > ", 
        params$log2fc_threshold, "): ", up_genes, "\n", file = report_file, append = TRUE)
    cat("  Down-regulated genes (padj < ", params$padj_threshold, ", log2FC < -", 
        params$log2fc_threshold, "): ", down_genes, "\n", file = report_file, append = TRUE)
    
    # SCN1A status
    scn1a_rows <- which(res_df$is_SCN1A)
    if (length(scn1a_rows) > 0) {
      cat("  SCN1A status:\n", file = report_file, append = TRUE)
      for (i in scn1a_rows) {
        cat("    Ensembl ID: ", rownames(res_df)[i], "\n", file = report_file, append = TRUE)
        cat("    log2FC: ", res_df$log2FoldChange[i], "\n", file = report_file, append = TRUE)
        cat("    padj: ", res_df$padj[i], "\n", file = report_file, append = TRUE)
        status <- "Not significant"
        if (!is.na(res_df$padj[i]) && res_df$padj[i] < params$padj_threshold) {
          if (res_df$log2FoldChange[i] > params$log2fc_threshold) {
            status <- paste0("Significantly up-regulated (threshold: log2FC > ", 
                           params$log2fc_threshold, ")")
          } else if (res_df$log2FoldChange[i] < -params$log2fc_threshold) {
            status <- paste0("Significantly down-regulated (threshold: log2FC < -", 
                           params$log2fc_threshold, ")")
          } else {
            status <- paste0("Significant but |log2FC| < ", params$log2fc_threshold)
          }
        }
        cat("    Status: ", status, "\n", file = report_file, append = TRUE)
      }
    } else {
      cat("  SCN1A not found in results\n", file = report_file, append = TRUE)
    }
    
    # Top genes
    top_up <- res_df %>%
      filter(!is.na(padj) & padj < params$padj_threshold & 
             log2FoldChange > params$log2fc_threshold) %>%
      arrange(padj) %>%
      head(5)
    
    top_down <- res_df %>%
      filter(!is.na(padj) & padj < params$padj_threshold & 
             log2FoldChange < -params$log2fc_threshold) %>%
      arrange(padj) %>%
      head(5)
    
    cat("  Top 5 up-regulated genes:\n", file = report_file, append = TRUE)
    if (nrow(top_up) > 0) {
      for (i in 1:nrow(top_up)) {
        cat("    ", top_up$gene_symbol[i], ": log2FC = ", round(top_up$log2FoldChange[i], 2), 
            ", padj = ", format(top_up$padj[i], scientific = TRUE, digits = 3), "\n", 
            file = report_file, append = TRUE)
      }
    } else {
      cat("    None found\n", file = report_file, append = TRUE)
    }
    
    cat("  Top 5 down-regulated genes:\n", file = report_file, append = TRUE)
    if (nrow(top_down) > 0) {
      for (i in 1:nrow(top_down)) {
        cat("    ", top_down$gene_symbol[i], ": log2FC = ", round(top_down$log2FoldChange[i], 2), 
            ", padj = ", format(top_down$padj[i], scientific = TRUE, digits = 3), "\n", 
            file = report_file, append = TRUE)
      }
    } else {
      cat("    None found\n", file = report_file, append = TRUE)
    }
  }
  
  # Output file locations
  cat("\nOutput File Locations\n", file = report_file, append = TRUE)
  cat("-------------------\n", file = report_file, append = TRUE)
  cat("Results directory: ", paths$results, "\n", file = report_file, append = TRUE)
  cat("MA plots directory: ", paths$ma_plots, "\n", file = report_file, append = TRUE)
  cat("Volcano plots directory: ", paths$volcano_plots, "\n", file = report_file, append = TRUE)
  cat("PCA plots directory: ", paths$pca, "\n", file = report_file, append = TRUE)
  cat("Heatmaps directory: ", paths$heatmaps, "\n", file = report_file, append = TRUE)
  cat("Gene Ontology directory: ", paths$go, "\n", file = report_file, append = TRUE)
  cat("SCN1A analysis directory: ", paths$scn1a, "\n", file = report_file, append = TRUE)
  cat("lncRNA analysis directory: ", paths$lncrna, "\n", file = report_file, append = TRUE)
  cat("Venn diagrams directory: ", paths$venn, "\n", file = report_file, append = TRUE)  # Add this line
  cat("Debug logs directory: ", paths$debug, "\n", file = report_file, append = TRUE)
  
  cat("Analysis report created at", report_file, "\n")
  return(report_file)
}

# Function to analyze the 22 common genes in detail
# Function to analyze the 22 common genes in detail
analyze_common_degs <- function(results_list, venn_path, base_path, params) {
  cat("\n=== ANALYZING 22 COMMON DEGs ===\n")
  
  # Load the common genes list
  common_genes_file <- file.path(venn_path, "Common_DEGs_all_three_conditions.txt")
  
  # Check if file exists
  if (!file.exists(common_genes_file)) {
    cat("ERROR: Common genes file not found at:", common_genes_file, "\n")
    cat("Looking for files in venn directory:\n")
    print(list.files(venn_path))
    return(NULL)
  }
  
  # Read the file, skipping comment lines
  common_genes_data <- read.table(common_genes_file, 
                                 sep = "\t", 
                                 header = TRUE, 
                                 comment.char = "#",
                                 stringsAsFactors = FALSE)
  
  cat("Found", nrow(common_genes_data), "common genes\n")
  print(head(common_genes_data))
  
  # 1. Create a detailed summary table
  gene_details <- data.frame(
    Gene = common_genes_data$GeneSymbol,
    saRNA10_FC = common_genes_data$saRNA10_log2FC,
    saRNA12_FC = common_genes_data$saRNA12_log2FC,
    saRNA14_FC = common_genes_data$saRNA14_log2FC,
    Mean_FC = common_genes_data$Mean_log2FC,
    Direction = ifelse(common_genes_data$Mean_log2FC > 0, "Up", "Down"),
    stringsAsFactors = FALSE
  )
  
  # Sort by mean fold change
  gene_details <- gene_details[order(abs(gene_details$Mean_FC), decreasing = TRUE), ]
  
  cat("\nTop common genes by fold change:\n")
  print(head(gene_details, 10))
  
  # 2. Create visualization of common genes
  library(ggplot2)
  library(tidyr)
  
  # Reshape for plotting
  gene_details_long <- gene_details %>%
    pivot_longer(cols = c(saRNA10_FC, saRNA12_FC, saRNA14_FC),
                names_to = "Condition",
                values_to = "log2FC") %>%
    mutate(Condition = gsub("_FC", "", Condition))
  
  # Create heatmap-style plot
  p1 <- ggplot(gene_details_long, aes(x = Condition, y = Gene, fill = log2FC)) +
    geom_tile() +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", 
                        midpoint = 0, limits = c(-2, 2),
                        oob = scales::squish) +
    theme_minimal() +
    theme(axis.text.y = element_text(size = 10),
          axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Expression Changes in 22 Common DEGs",
         subtitle = "Genes differentially expressed in all three conditions",
         x = "", y = "", fill = "log2FC")
  
  ggsave(file.path(output_path, "Common_22_genes_heatmap.pdf"),
         plot = p1, width = 8, height = 10)
  
  # 3. Try GO analysis with relaxed parameters
  cat("\nAttempting GO analysis with relaxed parameters...\n")
  
  # Get Entrez IDs for common genes
  ensembl_ids <- common_genes_data$EnsemblID
  clean_ids <- gsub("\\..*$", "", ensembl_ids)
  
  entrez_ids <- suppressWarnings(
    AnnotationDbi::mapIds(
      org.Hs.eg.db,
      keys = clean_ids,
      column = "ENTREZID",
      keytype = "ENSEMBL",
      multiVals = "first"
    )
  )
  
  entrez_ids <- entrez_ids[!is.na(entrez_ids)]
  cat("Mapped to", length(entrez_ids), "Entrez IDs\n")
  
  if (length(entrez_ids) >= 5) {
    # Try GO with very relaxed parameters
    tryCatch({
      go_bp <- enrichGO(
        gene = entrez_ids,
        OrgDb = org.Hs.eg.db,
        ont = "BP",
        pAdjustMethod = "BH",
        pvalueCutoff = 0.1,  # Relaxed
        qvalueCutoff = 0.3,  # Relaxed
        readable = TRUE
      )
      
      if (nrow(go_bp) > 0) {
        cat("Found", nrow(go_bp), "enriched BP terms\n")
        write.csv(as.data.frame(go_bp), 
                 file.path(output_path, "Common_22_genes_GO_BP.csv"),
                 row.names = FALSE)
      } else {
        cat("No enriched GO terms found (expected with only 22 genes)\n")
      }
    }, error = function(e) {
      cat("GO analysis error:", e$message, "\n")
    })
  }
  
  # 4. Manual pathway annotation
  cat("\nPerforming manual annotation of key genes...\n")
  
  key_genes <- c("SCN1A", "CACNA1E", "VGF", "ABCA12", "CLDN11")
  annotations <- data.frame(
    Gene = key_genes,
    Function = c(
      "Voltage-gated sodium channel - epilepsy gene",
      "Voltage-gated calcium channel - neuronal signaling",
      "VGF nerve growth factor - neuropeptide precursor",
      "ATP-binding cassette transporter - lipid transport",
      "Claudin-11 - tight junction protein"
    ),
    Category = c("Ion channel", "Ion channel", "Neuropeptide", 
                "Transporter", "Cell junction")
  )
  
  cat("\nKey genes in common set:\n")
  print(annotations)
  
  # 5. String network analysis suggestion
  cat("\n\nSUGGESTIONS FOR FURTHER ANALYSIS:\n")
  cat("1. STRING network analysis at: https://string-db.org/\n")
  cat("   - Input the 22 gene symbols to see protein interactions\n")
  cat("2. Manual curation of gene functions\n")
  cat("3. Compare with published gene sets\n")
  
  # Save gene list for external tools
  write.table(unique(common_genes_data$GeneSymbol),
             file.path(output_path, "Common_22_genes_symbols_only.txt"),
             quote = FALSE, row.names = FALSE, col.names = FALSE)
  
  return(gene_details)
}

# Alternative: Analyze the larger overlaps
analyze_expanded_gene_sets <- function(output_path) {
  cat("\n=== ANALYZING LARGER GENE SETS ===\n")
  
  # Analyze genes common to saRNA12 & saRNA14 (373 genes)
  # This larger set is more amenable to GO analysis
  
  cat("\nFor more robust GO analysis, consider analyzing:\n")
  cat("1. saRNA12 ∩ saRNA14: 373 genes (better statistical power)\n")
  cat("2. All upregulated genes across conditions\n")
  cat("3. All downregulated genes across conditions\n")
}

# Function to save results
save_results <- function(results_list, output_path) {
  cat("Saving DESeq2 results to CSV files...\n")
  
  for (cond in names(results_list)) {
    output_file <- file.path(output_path, paste0("DESeq2_results_", cond, "_vs_SCR.csv"))
    write_csv(results_list[[cond]], output_file)
    cat("Results saved to:", output_file, "\n")
  }
}

# Main function - Updated to use centralized parameters
main <- function() {
  # Define input and output paths separately
  input_base_path <- "D:/20250415 RNA seq analysis"
  output_base_path <- "D:/"
  cat("Input path:", input_base_path, "\n")
  cat("Output path:", output_base_path, "\n")
  
  # Create organized output directories with date in path
  paths <- create_output_dirs(output_base_path)
  cat("Created output directories in:", paths$base, "\n")
  
  # Save analysis parameters to a file for reference
  params_file <- file.path(paths$base, "analysis_parameters.txt")
  cat("Analysis Parameters\n", file = params_file)
  cat("==================\n", file = params_file, append = TRUE)
  cat("Date:", Sys.Date(), "\n", file = params_file, append = TRUE)
  cat("padj_threshold:", analysis_params$padj_threshold, "\n", file = params_file, append = TRUE)
  cat("log2fc_threshold:", analysis_params$log2fc_threshold, "\n", file = params_file, append = TRUE)
  cat("heatmap_log2fc_threshold:", analysis_params$heatmap_log2fc_threshold, "\n", file = params_file, append = TRUE)
  cat("heatmap_top_genes:", analysis_params$heatmap_top_genes, "\n", file = params_file, append = TRUE)
  cat("go_pvalue_cutoff:", analysis_params$go_pvalue_cutoff, "\n", file = params_file, append = TRUE)
  cat("go_qvalue_cutoff:", analysis_params$go_qvalue_cutoff, "\n", file = params_file, append = TRUE)
  cat("ma_plot_ylim:", paste(analysis_params$ma_plot_ylim, collapse = ", "), "\n", file = params_file, append = TRUE)
  cat("volcano_xlim:", paste(analysis_params$volcano_xlim, collapse = ", "), "\n", file = params_file, append = TRUE)
  
  # Start a log file
  log_file <- file.path(paths$debug, "analysis_log.txt")
  cat("RNA-seq analysis log\n", file = log_file)
  cat("Date:", Sys.Date(), "\n", file = log_file, append = TRUE)
  cat("R version:", R.version.string, "\n", file = log_file, append = TRUE)
  cat("Analysis parameters:\n", file = log_file, append = TRUE)
  cat("  padj threshold:", analysis_params$padj_threshold, "\n", file = log_file, append = TRUE)
  cat("  log2FC threshold:", analysis_params$log2fc_threshold, "\n", file = log_file, append = TRUE)
  
  cat("Starting RNA-seq analysis with DESeq2\n", file = log_file, append = TRUE)
  
  # Check if the preprocessed files already exist - use input_base_path
  counts_path <- file.path(input_base_path, "Featurecounts", "combined_featurecounts_filtered.csv")
  metadata_path <- file.path(input_base_path, "Featurecounts", "metadata.csv")
  
  if (file.exists(counts_path) && file.exists(metadata_path)) {
    cat("Using existing preprocessed files:\n")
    cat("  Counts file:", counts_path, "\n")
    cat("  Metadata file:", metadata_path, "\n")
    
    cat("Loading preprocessed files", file = log_file, append = TRUE)
    cat("  Counts file:", counts_path, "\n", file = log_file, append = TRUE)
    cat("  Metadata file:", metadata_path, "\n", file = log_file, append = TRUE)
    
    # Load preprocessed data
    counts_df <- read_csv(counts_path)
    
    # Log data dimensions
    cat("  Loaded count data:", nrow(counts_df), "genes,", ncol(counts_df) - 1, "samples\n", 
        file = log_file, append = TRUE)
    
    # Create corrected metadata from count data
    metadata_df <- create_metadata(counts_df)
    cat("Created metadata with corrected condition labels focusing on saRNAs 7, 10, 12, 14, and SCR\n")
    
    # Log metadata
    cat("  Created metadata with", nrow(metadata_df), "samples\n", file = log_file, append = TRUE)
    cat("  Conditions:\n", file = log_file, append = TRUE)
    print(table(metadata_df$Condition), file = log_file, append = TRUE)
  } else {
    # Try to find the raw input files to process
    feature_counts_dir <- file.path(input_base_path, "Featurecounts")
    raw_files <- list.files(feature_counts_dir, pattern = ".*_featureCounts\\.raw$", full.names = TRUE)
    
    if (length(raw_files) > 0) {
      cat("Found", length(raw_files), "raw feature count files.\n")
      cat("Raw files found:", paste(basename(raw_files), collapse=", "), "\n")
      cat("You'll need to process these files first to create the combined count matrix.\n")
    } else {
      error_msg <- paste("Required preprocessed files not found at", counts_path, 
                         paste0("and no raw feature count files found in ", feature_counts_dir))
      cat(error_msg, "\n")
      cat("ERROR:", error_msg, "\n", file = log_file, append = TRUE)
      stop(error_msg)
    }
  }

  # Run DESeq2 analysis with parameters
  cat("Running DESeq2 analysis...\n", file = log_file, append = TRUE)
  analysis_results <- run_deseq2_analysis(counts_df, metadata_df, paths$debug, analysis_params)
  dds <- analysis_results$dds
  results_list <- analysis_results$results_list
  gene_symbols <- analysis_results$gene_symbols
  entrez_ids <- analysis_results$entrez_ids
  
  # Save results
  cat("Saving results...\n", file = log_file, append = TRUE)
  save_results(results_list, paths$results)
  
  # Generate all plots and perform analysis steps, with error handling
  tryCatch({
    cat("Generating MA plots with similar sequence gene labeling...\n", file = log_file, append = TRUE)
    generate_ma_plots(dds, results_list, paths$ma_plots, analysis_params)
  }, error = function(e) {
    warning_msg <- paste("Error creating MA plots:", e$message)
    warning(warning_msg)
    cat("WARNING:", warning_msg, "\n", file = log_file, append = TRUE)
  })
  
  tryCatch({
    cat("Generating PCA plots...\n", file = log_file, append = TRUE)
    generate_pca_plots(dds, paths$pca, analysis_params$saRNA_colors)
  }, error = function(e) {
    warning_msg <- paste("Error creating PCA plots:", e$message)
    warning(warning_msg)
    cat("WARNING:", warning_msg, "\n", file = log_file, append = TRUE)
  })
  
  tryCatch({
    cat("Generating volcano plots with similar sequence gene labeling...\n", file = log_file, append = TRUE)
    generate_volcano_plots(results_list, paths$volcano_plots, analysis_params)
  }, error = function(e) {
    warning_msg <- paste("Error creating volcano plots:", e$message)
    warning(warning_msg)
    cat("WARNING:", warning_msg, "\n", file = log_file, append = TRUE)
  })
  
  tryCatch({
    cat("Generating heatmaps...\n", file = log_file, append = TRUE)
    generate_heatmaps(dds, results_list, paths$heatmaps, analysis_params)
  }, error = function(e) {
    warning_msg <- paste("Error creating heatmaps:", e$message)
    warning(warning_msg)
    cat("WARNING:", warning_msg, "\n", file = log_file, append = TRUE)
  })
  
  tryCatch({
    cat("Generating SCN1A expression plots...\n", file = log_file, append = TRUE)
    generate_scn1a_plots(dds, results_list, paths$scn1a, analysis_params$saRNA_colors)
  }, error = function(e) {
    warning_msg <- paste("Error creating SCN1A expression plots:", e$message)
    warning(warning_msg)
    cat("WARNING:", warning_msg, "\n", file = log_file, append = TRUE)
  })
  
  tryCatch({
    cat("Generating summary barplot...\n", file = log_file, append = TRUE)
    generate_summary_barplot(results_list, paths$base, analysis_params)
  }, error = function(e) {
    warning_msg <- paste("Error creating summary barplot:", e$message)
    warning(warning_msg)
    cat("WARNING:", warning_msg, "\n", file = log_file, append = TRUE)
  })
  
  tryCatch({
    cat("Generating lncRNA expression scatter plots...\n", file = log_file, append = TRUE)
    generate_lncrna_scatter_plots(dds, results_list, paths$lncrna, analysis_params)
  }, error = function(e) {
    warning_msg <- paste("Error creating lncRNA scatter plots:", e$message)
    warning(warning_msg)
    cat("WARNING:", warning_msg, "\n", file = log_file, append = TRUE)
  })


  # Generate Venn diagrams
  tryCatch({
    cat("Generating Venn diagrams for DEG comparisons...\n", file = log_file, append = TRUE)
    generate_venn_diagrams(results_list, paths$venn, analysis_params)
  }, error = function(e) {
    warning_msg <- paste("Error creating Venn diagrams:", e$message)
    warning(warning_msg)
    cat("WARNING:", warning_msg, "\n", file = log_file, append = TRUE)
  })
  
  # Run Gene Ontology analysis
  tryCatch({
    cat("Running Gene Ontology analysis...\n", file = log_file, append = TRUE)
    run_go_analysis(results_list, paths$go, analysis_params)
  }, error = function(e) {
    warning_msg <- paste("Error running Gene Ontology analysis:", e$message)
    warning(warning_msg)
    cat("WARNING:", warning_msg, "\n", file = log_file, append = TRUE)
  })
  # Further Venn diagram analysis 
  tryCatch({
  cat("Analyzing common DEGs from Venn diagram...\n", file = log_file, append = TRUE)
  analyze_common_degs(results_list, paths$venn, paths$base, analysis_params)
}, error = function(e) {
  warning_msg <- paste("Error analyzing common DEGs:", e$message)
  warning(warning_msg)
  cat("WARNING:", warning_msg, "\n", file = log_file, append = TRUE)
})
  
  # Create analysis report
  tryCatch({
    cat("Creating analysis report...\n", file = log_file, append = TRUE)
    report_file <- create_analysis_report(dds, results_list, paths, analysis_params)
    cat("Analysis report created at:", report_file, "\n", file = log_file, append = TRUE)
  }, error = function(e) {
    warning_msg <- paste("Error creating analysis report:", e$message)
    warning(warning_msg)
    cat("WARNING:", warning_msg, "\n", file = log_file, append = TRUE)
  })
  
  cat("Analysis completed successfully!\n")
  cat("Analysis completed successfully!\n", file = log_file, append = TRUE)
  cat("End time:", Sys.time(), "\n", file = log_file, append = TRUE)
}

# Add this diagnostic function to check what's happening with SCN1A
diagnose_scn1a_expression <- function(dds, results_list, output_path) {
  cat("\n=== DIAGNOSING SCN1A EXPRESSION ISSUE ===\n")
  
  # Find SCN1A
  scn1a_idx <- which(results_list[[1]]$gene_symbol == "SCN1A")
  if (length(scn1a_idx) == 0) {
    cat("SCN1A not found!\n")
    return()
  }
  
  scn1a_ensembl <- rownames(results_list[[1]])[scn1a_idx]
  cat("SCN1A Ensembl ID:", scn1a_ensembl, "\n\n")
  
  # Get raw counts
  raw_counts <- counts(dds, normalized = FALSE)[scn1a_ensembl, ]
  norm_counts <- counts(dds, normalized = TRUE)[scn1a_ensembl, ]
  
  # Create detailed data frame
  sample_data <- data.frame(
    Sample = colnames(dds),
    Condition = colData(dds)$Condition,
    Raw_Count = raw_counts,
    Normalized_Count = norm_counts,
    Size_Factor = sizeFactors(dds),
    stringsAsFactors = FALSE
  )
  
  # Calculate means by condition
  condition_summary <- sample_data %>%
    group_by(Condition) %>%
    summarise(
      n = n(),
      Raw_Mean = mean(Raw_Count),
      Raw_SD = sd(Raw_Count),
      Norm_Mean = mean(Normalized_Count),
      Norm_SD = sd(Normalized_Count),
      Mean_Size_Factor = mean(Size_Factor)
    )
  
  cat("CONDITION SUMMARY:\n")
  print(condition_summary)
  cat("\n")
  
  # Calculate fold changes manually
  scr_raw_mean <- condition_summary$Raw_Mean[condition_summary$Condition == "SCR"]
  scr_norm_mean <- condition_summary$Norm_Mean[condition_summary$Condition == "SCR"]
  
  fc_comparison <- data.frame(
    Condition = c("saRNA7", "saRNA10", "saRNA12", "saRNA14"),
    Raw_FC = numeric(4),
    Norm_FC = numeric(4),
    log2FC_Raw = numeric(4),
    log2FC_Norm = numeric(4),
    log2FC_DESeq2 = numeric(4),
    padj_DESeq2 = numeric(4)
  )
  
  for (i in 1:nrow(fc_comparison)) {
    cond <- fc_comparison$Condition[i]
    cond_raw_mean <- condition_summary$Raw_Mean[condition_summary$Condition == cond]
    cond_norm_mean <- condition_summary$Norm_Mean[condition_summary$Condition == cond]
    
    fc_comparison$Raw_FC[i] <- cond_raw_mean / scr_raw_mean
    fc_comparison$Norm_FC[i] <- cond_norm_mean / scr_norm_mean
    fc_comparison$log2FC_Raw[i] <- log2(fc_comparison$Raw_FC[i])
    fc_comparison$log2FC_Norm[i] <- log2(fc_comparison$Norm_FC[i])
    fc_comparison$log2FC_DESeq2[i] <- results_list[[cond]]$log2FoldChange[scn1a_idx]
    fc_comparison$padj_DESeq2[i] <- results_list[[cond]]$padj[scn1a_idx]
  }
  
  cat("\nFOLD CHANGE COMPARISON:\n")
  print(fc_comparison)
  cat("\n")
  
  # Check for outliers
  cat("CHECKING FOR OUTLIERS:\n")
  for (cond in unique(sample_data$Condition)) {
    cond_data <- sample_data[sample_data$Condition == cond, ]
    if (nrow(cond_data) > 2) {
      # Calculate z-scores
      z_scores <- abs((cond_data$Normalized_Count - mean(cond_data$Normalized_Count)) / sd(cond_data$Normalized_Count))
      outliers <- which(z_scores > 2.5)
      if (length(outliers) > 0) {
        cat(cond, "- Potential outliers (z-score > 2.5):\n")
        print(cond_data[outliers, ])
      }
    }
  }
  
  # Create diagnostic plots
  library(ggplot2)
  library(gridExtra)
  
  # Plot 1: Raw counts by sample
  p1 <- ggplot(sample_data, aes(x = Sample, y = Raw_Count, fill = Condition)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8)) +
    labs(title = "SCN1A Raw Counts by Sample", y = "Raw Count") +
    facet_wrap(~Condition, scales = "free_x")
  
  # Plot 2: Size factors by condition
  p2 <- ggplot(sample_data, aes(x = Condition, y = Size_Factor)) +
    geom_boxplot(aes(fill = Condition)) +
    geom_jitter(width = 0.2, alpha = 0.5) +
    theme_minimal() +
    labs(title = "Size Factors by Condition", y = "Size Factor") +
    theme(legend.position = "none")
  
  # Plot 3: Normalized vs Raw counts
  p3 <- ggplot(sample_data, aes(x = Raw_Count, y = Normalized_Count, color = Condition)) +
    geom_point(size = 3, alpha = 0.7) +
    geom_smooth(method = "lm", se = FALSE) +
    theme_minimal() +
    labs(title = "Raw vs Normalized SCN1A Counts", 
         x = "Raw Count", 
         y = "Normalized Count")
  
  # Plot 4: Fold change comparison
  fc_plot_data <- fc_comparison %>%
    select(Condition, log2FC_Raw, log2FC_Norm, log2FC_DESeq2) %>%
    pivot_longer(cols = starts_with("log2FC"), 
                names_to = "Method", 
                values_to = "log2FC")
  
  fc_plot_data$Method <- factor(fc_plot_data$Method, 
                               levels = c("log2FC_Raw", "log2FC_Norm", "log2FC_DESeq2"),
                               labels = c("Raw Counts", "Normalized", "DESeq2"))
  
  p4 <- ggplot(fc_plot_data, aes(x = Condition, y = log2FC, fill = Method)) +
    geom_bar(stat = "identity", position = "dodge") +
    geom_hline(yintercept = 0, linetype = "solid") +
    geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", alpha = 0.5) +
    theme_minimal() +
    labs(title = "SCN1A log2 Fold Changes: Comparison of Methods",
         y = "log2(Fold Change)") +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  # Save diagnostic plots
  diagnostic_plot <- arrangeGrob(p1, p2, p3, p4, ncol = 2, nrow = 2)
  ggsave(file.path(output_path, "SCN1A_diagnostic_plots.pdf"), 
         plot = diagnostic_plot, width = 16, height = 12)
  
  # Save detailed sample data
  write.csv(sample_data, 
           file.path(output_path, "SCN1A_sample_level_data.csv"),
           row.names = FALSE)
  
  # Save comparison data
  write.csv(fc_comparison,
           file.path(output_path, "SCN1A_fold_change_comparison.csv"),
           row.names = FALSE)
  
  cat("\nDiagnostic files saved to:", output_path, "\n")
  
  # Print warnings
  if (any(condition_summary$n < 3)) {
    cat("\nWARNING: Some conditions have very few replicates:\n")
    print(condition_summary[condition_summary$n < 3, c("Condition", "n")])
  }
  
  # Check if DESeq2 results are very different from normalized means
  large_discrepancy <- abs(fc_comparison$log2FC_Norm - fc_comparison$log2FC_DESeq2) > 0.5
  if (any(large_discrepancy)) {
    cat("\nWARNING: Large discrepancies between normalized means and DESeq2 results:\n")
    print(fc_comparison[large_discrepancy, ])
  }
}

# Run the script
main()