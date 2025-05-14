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

# Function to create output directories with date in path
create_output_dirs <- function(base_path) {
  # Add current date to folder name
  current_date <- format(Sys.Date(), "%Y%m%d")
  output_folder <- paste0("RNA_seq_analysis_", current_date)
  
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

# Function to run DESeq2 analysis
run_deseq2_analysis <- function(counts_df, metadata_df, debug_path) {
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
    
    # Summary of results
    sig_genes <- sum(results_list[[cond]]$padj < 0.05, na.rm = TRUE)
    up_genes <- sum(results_list[[cond]]$padj < 0.05 & results_list[[cond]]$log2FoldChange > 0.5, na.rm = TRUE)
    down_genes <- sum(results_list[[cond]]$padj < 0.05 & results_list[[cond]]$log2FoldChange < -0.5, na.rm = TRUE)
    cat("    Found", sig_genes, "significantly differentially expressed genes (padj < 0.05)\n")
    cat("    Up-regulated genes:", up_genes, "\n")
    cat("    Down-regulated genes:", down_genes, "\n")
    
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
    cat("Total genes analyzed:", nrow(res_df), "\n", file = stats_file, append = TRUE)
    cat("Genes with valid p-value:", sum(!is.na(res_df$pvalue)), "\n", file = stats_file, append = TRUE)
    cat("Genes with valid padj:", sum(!is.na(res_df$padj)), "\n", file = stats_file, append = TRUE)
    cat("Significant genes (padj < 0.05):", sig_genes, "\n", file = stats_file, append = TRUE)
    cat("Up-regulated genes (padj < 0.05, log2FC > 1):", up_genes, "\n", file = stats_file, append = TRUE)
    cat("Down-regulated genes (padj < 0.05, log2FC < -1):", down_genes, "\n", file = stats_file, append = TRUE)
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

# Function to save results
save_results <- function(results_list, output_path) {
  cat("Saving DESeq2 results to CSV files...\n")
  
  for (cond in names(results_list)) {
    output_file <- file.path(output_path, paste0("DESeq2_results_", cond, "_vs_SCR.csv"))
    write_csv(results_list[[cond]], output_file)
    cat("Results saved to:", output_file, "\n")
  }
}

# Function to generate enhanced MA plots with off-target gene labeling
generate_ma_plots <- function(dds, results_list, ma_output_path) {
  cat("Generating enhanced MA plots with off-target gene labeling...\n")
  
  # Define potential off-target genes for each condition (based on sequence similarity)
  off_target_genes <- list(
    "saRNA7" = c("HDAC8", "RASSF8-AS1", "NET1", "KHDRBS2", "ARHGAP10", "TRIM2", "BTBD8", "ASAP2"),
    "saRNA10" = c("ICA1"),
    "saRNA12" = c("NUDT3", "PLPPR5"),
    "saRNA14" = c("MAP3K21", "DNAH2", "THADA", "DOCK3", "LINC03025")
  )
  
  for (cond in names(results_list)) {
    cat("  Creating MA plot for", cond, "vs SCR\n")
    res_obj <- results(dds, contrast = c("Condition", cond, "SCR"))
    
    # Convert to data frame for ggplot
    res_df <- as.data.frame(res_obj)
    
    # Add gene symbols
    res_df$gene_symbol <- results_list[[cond]]$gene_symbol[match(rownames(res_df), rownames(results_list[[cond]]))]
    
    # Add significance categories for coloring
    res_df$Significance <- "Not Significant"
    res_df$Significance[!is.na(res_df$padj) & res_df$padj < 0.05] <- "Significant"
    
    # Mark SCN1A
    res_df$is_SCN1A <- results_list[[cond]]$is_SCN1A[match(rownames(res_df), rownames(results_list[[cond]]))]
    
    # Add flag for off-target genes based on the current condition
    res_df$is_OffTarget <- FALSE
    if (cond %in% names(off_target_genes)) {
      # Mark the specific off-target genes for this condition
      res_df$is_OffTarget <- res_df$gene_symbol %in% off_target_genes[[cond]]
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
    cat("    Significant genes (padj < 0.05):", sum(res_df$Significance == "Significant", na.rm = TRUE), "\n")
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
      # Add horizontal lines at log2FC = +/- 1
      geom_hline(yintercept = c(-1, 1), linetype = "dotted", color = c("blue", "red")) +
      # Set colors similar to volcano plot
      scale_color_manual(values = c("Significant" = "#3366CC", "Not Significant" = "grey70")) +
      # Add labels
      labs(
        title = paste0("MA Plot: ", cond, " vs SCR"),
        subtitle = paste0("Total genes: ", nrow(res_df), " | Significant: ", 
                         sum(res_df$Significance == "Significant", na.rm = TRUE)),
        x = "log10(Mean Expression)",
        y = "log2(Fold Change)",
        caption = "Purple points: genes with similar sequence to target"
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
      # Set limits for better visibility
      ylim(-3, 3)
    
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
    if (cond %in% names(off_target_genes)) {
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
    
    # Also save the traditional DESeq2 MA plot for comparison
    pdf(file.path(ma_output_path, paste0("DESeq2_MA_plot_", cond, "_vs_SCR.pdf")), width = 10, height = 8)
    plotMA(res_obj, main = paste0("DESeq2 MA Plot: ", cond, " vs SCR"), ylim = c(-3, 3), alpha = 0.05)
    has_scn1a <- FALSE
    has_off_target <- FALSE
    
    # Add SCN1A
    if (nrow(scn1a_data) > 0) {
      for (i in 1:nrow(scn1a_data)) {
        points(scn1a_data$baseMean[i], scn1a_data$log2FoldChange[i], 
               col = "green4", pch = 19, cex = 1.5)
        text(scn1a_data$baseMean[i], scn1a_data$log2FoldChange[i], 
             labels = "SCN1A", pos = 3, col = "green4", font = 2)
      }
      has_scn1a <- TRUE
    }
    
    # Add off-target genes to the traditional plot if present
    if (cond %in% names(off_target_genes)) {
      off_target_data <- res_df[res_df$is_OffTarget, ]
      if (nrow(off_target_data) > 0) {
        for (i in 1:nrow(off_target_data)) {
          points(off_target_data$baseMean[i], off_target_data$log2FoldChange[i], 
                 col = "#9966CC", pch = 19, cex = 1.5)
          text(off_target_data$baseMean[i], off_target_data$log2FoldChange[i], 
               labels = off_target_data$gene_symbol[i], pos = 3, col = "#9966CC", font = 2)
        }
        has_off_target <- TRUE
      }
    }
    
    # Add a custom legend
    if (has_scn1a || has_off_target) {
      legend_labels <- c()
      legend_colors <- c()
      legend_pch <- c()
      
      if (has_scn1a) {
        legend_labels <- c(legend_labels, "SCN1A")
        legend_colors <- c(legend_colors, "green4")
        legend_pch <- c(legend_pch, 19)
      }
      
      if (has_off_target) {
        legend_labels <- c(legend_labels, "Similar Sequence Gene")
        legend_colors <- c(legend_colors, "#9966CC")
        legend_pch <- c(legend_pch, 19)
      }
      
      legend("topright", legend = legend_labels, 
             col = legend_colors, pch = legend_pch, 
             cex = 0.8, pt.cex = 1.2, bty = "n")
    }
    
    dev.off()
  }
  cat("Created enhanced MA plots with off-target gene labeling\n")
}

# Function to generate volcano plots with off-target gene labels
generate_volcano_plots <- function(results_list, volcano_output_path) {
  cat("Generating volcano plots with improved styling and off-target gene labeling...\n")
  
  # Define potential off-target genes for each condition (based on sequence similarity)
  off_target_genes <- list(
    "saRNA7" = c("HDAC8", "RASSF8-AS1", "NET1", "KHDRBS2", "ARHGAP10", "TRIM2", "BTBD8", "ASAP2"),
    "saRNA10" = c("ICA1"),
    "saRNA12" = c("NUDT3", "PLPPR5"),
    "saRNA14" = c("MAP3K21", "DNAH2", "THADA", "DOCK3", "LINC03025")
  )
  
  for (cond in names(results_list)) {
    cat(paste0("  Creating volcano plot for ", cond, " vs SCR\n"))
    res_df <- results_list[[cond]]
    
    # Verification and debugging info
    total_genes <- nrow(res_df)
    genes_with_padj <- sum(!is.na(res_df$padj))
    sig_genes <- sum(!is.na(res_df$padj) & res_df$padj < 0.05, na.rm = TRUE)
    up_genes <- sum(!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange > 0.5, na.rm = TRUE)
    down_genes <- sum(!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange < -0.5, na.rm = TRUE)
    scn1a_count <- sum(res_df$is_SCN1A, na.rm = TRUE)
    
    cat("    Total genes: ", total_genes, "\n")
    cat("    Genes with padj value: ", genes_with_padj, "\n")
    cat("    Significant genes (padj < 0.05): ", sig_genes, "\n")
    cat("    Up-regulated genes: ", up_genes, "\n")
    cat("    Down-regulated genes: ", down_genes, "\n")
    cat("    SCN1A genes found: ", scn1a_count, "\n")
    
    if (scn1a_count > 0) {
      scn1a_rows <- which(res_df$is_SCN1A)
      for (i in scn1a_rows) {
        cat("    SCN1A details: log2FC =", res_df$log2FoldChange[i], ", padj =", res_df$padj[i], "\n")
      }
    }
    
    # Add categories for significance
    res_df$Expression <- "Unchanged"
    res_df$Expression[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange < -0.5] <- "Down-regulated"
    res_df$Expression[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange > 0.5] <- "Up-regulated"
    
    # Add a special category for SCN1A
    if (scn1a_count > 0) {
      # Mark SCN1A regardless of significance for visibility
      res_df$Expression[res_df$is_SCN1A] <- "SCN1A"
    }
    
    # Add flag for off-target genes based on the current condition
    res_df$is_OffTarget <- FALSE
    if (cond %in% names(off_target_genes)) {
      # Mark the specific off-target genes for this condition
      res_df$is_OffTarget <- res_df$gene_symbol %in% off_target_genes[[cond]]
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
    
    # Create the volcano plot with verification in title
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
      # Add reference lines
      geom_vline(xintercept = c(-0.5, 0.5), linetype = "dotted", color = c("blue", "red")) +
      geom_hline(yintercept = -log10(0.05), linetype = "dotted", color = "blue") +
      # Set limits
      ylim(0, max(-log10(res_df$padj), na.rm = TRUE) * 1.05) +
      xlim(-4.5, 4.5) +
      # Add detailed labels
      labs(
        title = paste0(gsub("saRNA", "sample", cond), " -> scrambledA"),
        subtitle = paste0("Total: ", total_genes, " | Up: ", up_genes, " | Down: ", down_genes),
        x = "log2(FC)",
        y = "-log10(padj)",
        caption = paste0("Analysis date: ", Sys.Date(), " | Purple points: genes with similar sequence to target")
      ) +
      # Set theme
      theme_light(base_size = 12) +
      theme(
        panel.grid.major = element_line(color = "white"),
        panel.grid.minor = element_line(color = "white"),
        panel.border = element_rect(color = "grey80", fill = NA),
        panel.background = element_rect(fill = "grey95"),
        plot.title = element_text(hjust = 0.5, size = 14),
        plot.subtitle = element_text(hjust = 0.5, size = 10),
        plot.caption = element_text(hjust = 1, size = 8, color = "grey50"),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        legend.background = element_rect(fill = "white", color = "grey80"),
        axis.title = element_text(face = "bold", size = 12),
        axis.text = element_text(size = 10)
      )
    
    # Select top 5 genes for each category to label
    top_down <- res_df %>%
      filter(Expression == "Down-regulated") %>%
      arrange(padj) %>%
      head(5)
    
    top_up <- res_df %>%
      filter(Expression == "Up-regulated") %>%
      arrange(padj) %>%
      head(5)
    
    # Combine and label them
    top_genes <- rbind(top_down, top_up)
    
    if (nrow(top_genes) > 0) {
      cat("    Top 5 down-regulated genes: ", paste(top_down$gene_symbol, collapse=", "), "\n")
      cat("    Top 5 up-regulated genes: ", paste(top_up$gene_symbol, collapse=", "), "\n")
      
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
        color = ifelse(top_genes$Expression == "Up-regulated", "#CC0000", "#3366CC")
      )
    }
    
    # Always make sure SCN1A is shown if present
    if (scn1a_count > 0) {
      scn1a_data <- res_df[res_df$is_SCN1A, ]
      cat("    Including SCN1A in volcano plot\n")
      
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
    
    # Add labels for off-target genes
    if (cond %in% names(off_target_genes)) {
      off_target_data <- res_df[res_df$is_OffTarget, ]
      if (nrow(off_target_data) > 0) {
        cat("    Including similar sequence genes in volcano plot\n")
        
        p <- p + geom_point(
          data = off_target_data,
          size = 3,
          color = "#9966CC",
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
    ggsave(file.path(volcano_output_path, paste0("Volcano_plot_", cond, "_vs_SCR.pdf")), 
           plot = p, width = 10, height = 8)
    ggsave(file.path(volcano_output_path, paste0("Volcano_plot_", cond, "_vs_SCR.png")), 
           plot = p, width = 10, height = 8, dpi = 300)
    
    # Save a debug file with key statistics
    debug_file <- file.path(volcano_output_path, paste0("debug_", cond, "_vs_SCR.txt"))
    cat(paste0("Debug information for ", cond, " vs SCR\n"), file = debug_file)
    cat(paste0("Analysis date: ", Sys.Date(), "\n"), file = debug_file, append = TRUE)
    cat(paste0("Total genes: ", total_genes, "\n"), file = debug_file, append = TRUE)
    cat(paste0("Genes with padj value: ", genes_with_padj, "\n"), file = debug_file, append = TRUE)
    cat(paste0("Significant genes (padj < 0.05): ", sig_genes, "\n"), file = debug_file, append = TRUE)
    cat(paste0("Up-regulated genes: ", up_genes, "\n"), file = debug_file, append = TRUE)
    cat(paste0("Down-regulated genes: ", down_genes, "\n"), file = debug_file, append = TRUE)
    cat(paste0("SCN1A genes found: ", scn1a_count, "\n"), file = debug_file, append = TRUE)
    if (scn1a_count > 0) {
      for (i in scn1a_rows) {
        cat(paste0("SCN1A details: log2FC = ", res_df$log2FoldChange[i], ", padj = ", res_df$padj[i], "\n"), 
            file = debug_file, append = TRUE)
      }
    }
    cat(paste0("Top 5 down-regulated genes: ", paste(top_down$gene_symbol, collapse=", "), "\n"), 
        file = debug_file, append = TRUE)
    cat(paste0("Top 5 up-regulated genes: ", paste(top_up$gene_symbol, collapse=", "), "\n"), 
        file = debug_file, append = TRUE)
    
    # Add off-target gene info to debug file
    if (cond %in% names(off_target_genes)) {
      cat(paste0("Genes with similar sequence to target: ", 
                paste(off_target_genes[[cond]], collapse=", "), "\n"), 
          file = debug_file, append = TRUE)
      
      off_target_data <- res_df[res_df$is_OffTarget, ]
      if (nrow(off_target_data) > 0) {
        cat(paste0("Similar sequence genes found in dataset: ", nrow(off_target_data), "\n"), 
            file = debug_file, append = TRUE)
        
        for (i in 1:nrow(off_target_data)) {
          cat(paste0("  ", off_target_data$gene_symbol[i], ": log2FC = ", 
                    round(off_target_data$log2FoldChange[i], 2), 
                    ", padj = ", format(off_target_data$padj[i], scientific = TRUE, digits = 3), "\n"), 
              file = debug_file, append = TRUE)
        }
      } else {
        cat("No off-target genes found in dataset\n", file = debug_file, append = TRUE)
      }
    }
  }
  cat("Created volcano plots with off-target gene labeling\n")
}

# Function to generate heatmaps
generate_heatmaps <- function(dds, results_list, heatmaps_output_path, saRNA_colors) {
  cat("Generating heatmaps...\n")
  
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
  ann_colors <- list(Condition = saRNA_colors)
  
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
    
    # Select top 50 differentially expressed genes sorted by adjusted p-value
    sig_genes <- rownames(res_df[!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 1, ])
    sig_count <- length(sig_genes)
    
    cat("    Found", sig_count, "significant DEGs\n")
    
    # If there are significant genes, create the heatmap
    if (sig_count > 0) {
      # Limit to top 50 if more than 50
      if (sig_count > 50) {
        sig_genes <- sig_genes[order(res_df[sig_genes, "padj"])][1:50]
        cat("    Using top 50 genes for heatmap\n")
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
      ann_colors <- list(Condition = saRNA_colors)
      
      # Create heatmap
      pdf(file.path(heatmaps_output_path, paste0("Top_DEGs_heatmap_", cond, "_vs_SCR.pdf")), 
          width = 12, height = 14)
      pheatmap(top_gene_counts_scaled,
              annotation_col = sample_info,
              annotation_colors = ann_colors,
              main = paste0("Top ", length(sig_genes), " Differentially Expressed Genes: ", cond, " vs SCR"),
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

# Function to generate enhanced SCN1A expression plots
generate_scn1a_plots <- function(dds, results_list, scn1a_output_path, saRNA_colors) {
  cat("Generating enhanced SCN1A expression plot...\n")
  
  # Find SCN1A
  scn1a_ids <- unique(unlist(lapply(results_list, function(res) rownames(res)[res$is_SCN1A])))
  
  if (length(scn1a_ids) > 0) {
    cat("  Found", length(scn1a_ids), "SCN1A transcripts:", paste(scn1a_ids, collapse=", "), "\n")
    vsd <- vst(dds, blind = TRUE)
    
    # Extract SCN1A normalized counts
    scn1a_data <- data.frame()
    
    for (id in scn1a_ids) {
      if (id %in% rownames(vsd)) {
        temp_data <- data.frame(
          Sample = colnames(vsd),
          Condition = dds$Condition,
          Expression = assay(vsd)[id, ],
          Gene = paste0("SCN1A (", id, ")")
        )
        scn1a_data <- rbind(scn1a_data, temp_data)
      }
    }
    
    if (nrow(scn1a_data) > 0) {
      # Save raw data to file
      write.csv(scn1a_data, file.path(scn1a_output_path, "SCN1A_expression_data.csv"), row.names = FALSE)
      
      # Ensure conditions are in the correct order
      scn1a_data$Condition <- factor(scn1a_data$Condition, 
                                    levels = c("SCR", "saRNA7", "saRNA10", "saRNA12", "saRNA14"))
      
      # Calculate statistics for error bars
      scn1a_stats <- scn1a_data %>%
        group_by(Condition, Gene) %>%
        summarize(
          Mean = mean(Expression),
          Median = median(Expression),
          SD = sd(Expression),
          N = n(),
          SE = SD / sqrt(N),
          CI95 = qt(0.975, df = N-1) * SE,
          .groups = "drop"
        )
      
      # Create enhanced boxplot with custom styling
      p <- ggplot() +
        # Add boxplot layer with improved styling
        geom_boxplot(
          data = scn1a_data, 
          aes(x = Condition, y = Expression, fill = Condition),
          width = 0.6, 
          alpha = 0.8,
          outlier.shape = NA
        ) +
        # Add individual points with jitter
        geom_jitter(
          data = scn1a_data,
          aes(x = Condition, y = Expression, color = Condition),
          width = 0.15, 
          height = 0,
          size = 3,
          alpha = 0.7,
          shape = 16
        ) +
        # Add error bars based on 95% CI
        geom_errorbar(
          data = scn1a_stats,
          aes(x = Condition, ymin = Mean - CI95, ymax = Mean + CI95, group = Condition),
          width = 0.2,
          color = "black",
          size = 0.5
        ) +
        # Add mean points
        geom_point(
          data = scn1a_stats,
          aes(x = Condition, y = Mean),
          size = 3,
          shape = 23,
          fill = "white",
          color = "black",
          stroke = 1.5
        ) +
        # Set colors
        scale_fill_manual(values = saRNA_colors) +
        scale_color_manual(values = saRNA_colors) +
        # Improved labels
        labs(
          title = "SCN1A Expression Analysis",
          subtitle = paste0("Based on ", length(scn1a_ids), " transcript(s)"),
          x = "",
          y = "Normalized Expression (VST)",
          caption = paste0("Analysis date: ", Sys.Date())
        ) +
        # Enhanced theme
        theme_minimal(base_size = 14) +
        theme(
          plot.title = element_text(hjust = 0.5, face = "bold", size = 20, margin = margin(b = 10)),
          plot.subtitle = element_text(hjust = 0.5, size = 14, color = "gray30", margin = margin(b = 20)),
          plot.caption = element_text(hjust = 1, size = 10, color = "gray50", margin = margin(t = 15)),
          legend.position = "none",
          axis.title.y = element_text(face = "bold", size = 16, margin = margin(r = 10)),
          axis.text.x = element_text(size = 14, face = "bold", color = "black"),
          axis.text.y = element_text(size = 12),
          panel.grid.major.y = element_line(color = "gray90"),
          panel.grid.major.x = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_rect(color = "gray80", fill = NA, size = 1),
          panel.background = element_rect(fill = "gray98"),
          plot.background = element_rect(fill = "white", color = NA),
          plot.margin = margin(20, 20, 20, 20)
        ) +
        # Customized y axis
        scale_y_continuous(
          expand = expansion(mult = c(0.05, 0.1)),
          breaks = scales::pretty_breaks(n = 8)
        )
      
      # If multiple SCN1A transcripts, facet the plot
      if (length(unique(scn1a_data$Gene)) > 1) {
        p <- p + facet_wrap(~Gene, scales = "free_y")
      }
      
      # Add a subtle highlight to SCR condition for reference
      p <- p + geom_rect(
        data = data.frame(xmin = 0.5, xmax = 1.5, ymin = -Inf, ymax = Inf),
        aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
        fill = "gray90",
        alpha = 0.3,
        inherit.aes = FALSE
      )
      
      # Save as PDF with higher resolution
      ggsave(file.path(scn1a_output_path, "SCN1A_expression.pdf"), 
             plot = p, width = 12, height = 9, dpi = 300)
      
      # Save as PNG with even higher resolution for presentations
      ggsave(file.path(scn1a_output_path, "SCN1A_expression.png"), 
             plot = p, width = 12, height = 9, dpi = 600)
             
      # Also save a version with white background for publications
      p_pub <- p + theme(
        panel.background = element_rect(fill = "white"),
        plot.background = element_rect(fill = "white"),
        panel.grid.major.y = element_line(color = "gray95"),
        panel.border = element_rect(color = "black", size = 0.5, fill = NA)
      )
      ggsave(file.path(scn1a_output_path, "SCN1A_expression_publication.pdf"), 
             plot = p_pub, width = 12, height = 9, dpi = 300)
             
      cat("  Created enhanced SCN1A expression plots\n")
      
      # Create a statistical summary
      write.csv(scn1a_stats, file.path(scn1a_output_path, "SCN1A_expression_stats.csv"), row.names = FALSE)
      cat("  SCN1A expression statistics saved\n")
    }
  } else {
    cat("  SCN1A not found in the dataset, skipping SCN1A expression plot\n")
  }
}

# Function to generate enhanced summary barplot with improved layout
generate_summary_barplot <- function(results_list, base_output_path) {
  cat("Generating enhanced summary barplot...\n")
  
  # Collect significant gene counts
  sig_counts <- data.frame(
    Condition = names(results_list),
    Down_regulated = sapply(results_list, function(res) sum(!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange < -0.5, na.rm = TRUE)),
    Up_regulated = sapply(results_list, function(res) sum(!is.na(res$padj) & res$padj < 0.05 & res$log2FoldChange > 0.5, na.rm = TRUE))
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
      subtitle = paste0("padj < 0.05, |log2FC| > 0.5"),
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
  # First get the max count for proper positioning
  max_count <- max(sig_counts_long$Count)
  
  # Create a data frame for annotations
  annotations <- data.frame(
    x = levels(sig_counts_long$Condition),
    y = rep(0, length(levels(sig_counts_long$Condition))),
    label = paste("Total:", sig_counts$Total[match(levels(sig_counts_long$Condition), sig_counts$Condition)])
  )
  annotations$label[is.na(annotations$label)] <- ""
  
  # Add the annotation text below the plot
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

# Function to run Gene Ontology enrichment analysis with separate up/down regulation
run_go_analysis <- function(results_list, go_output_path) {
  cat("Running Gene Ontology enrichment analysis with direction separation...\n")
  
  # Set up parameters for enrichment analysis
  pvalue_cutoff <- 0.05
  qvalue_cutoff <- 0.2
  
  for (cond in names(results_list)) {
    cat("  Processing GO analysis for:", cond, "vs SCR\n")
    res_df <- results_list[[cond]]
    
    # Separate up-regulated and down-regulated genes
    up_genes <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange > 0.5, ]
    down_genes <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange < -0.5, ]
    
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
        process_go_categories(up_entrez, res_df, cond, "upregulated", up_dir, pvalue_cutoff, qvalue_cutoff)
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
        process_go_categories(down_entrez, res_df, cond, "downregulated", down_dir, pvalue_cutoff, qvalue_cutoff)
      } else {
        cat("    Not enough mapped Entrez IDs for down-regulated genes in", cond, "(need at least 10)\n")
      }
    } else {
      cat("    Not enough down-regulated genes for", cond, "(need at least 10)\n")
    }
    
    # Original combined analysis for backward compatibility
    sig_genes <- res_df[!is.na(res_df$padj) & res_df$padj < 0.05 & abs(res_df$log2FoldChange) > 0.5, ]
    
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
    process_go_categories(sig_entrez, res_df, cond, "combined", go_output_path, pvalue_cutoff, qvalue_cutoff)
  }
  
  cat("Gene Ontology analysis completed with direction separation\n")
}

# Helper function to process GO categories
process_go_categories <- function(entrez_ids, res_df, cond, direction, output_path, pvalue_cutoff, qvalue_cutoff) {
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
        pvalueCutoff = pvalue_cutoff,
        qvalueCutoff = qvalue_cutoff,
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
        pvalueCutoff = pvalue_cutoff,
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

# Function to create a detailed analysis report
create_analysis_report <- function(dds, results_list, paths) {
  cat("Creating analysis report...\n")
  
  # Create report file
  report_file <- file.path(paths$base, "RNA_seq_analysis_report.txt")
  
  # Write report header
  cat("RNA-seq Analysis Report\n", file = report_file)
  cat("=====================\n", file = report_file, append = TRUE)
  cat("Date: ", Sys.Date(), "\n\n", file = report_file, append = TRUE)
  
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
    
    # Count significant genes
    sig_genes <- sum(!is.na(res_df$padj) & res_df$padj < 0.05, na.rm = TRUE)
    up_genes <- sum(!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange > 1, na.rm = TRUE)
    down_genes <- sum(!is.na(res_df$padj) & res_df$padj < 0.05 & res_df$log2FoldChange < -1, na.rm = TRUE)
    
    cat("\n", cond, " vs SCR:\n", file = report_file, append = TRUE)
    cat("  Significant genes (padj < 0.05): ", sig_genes, "\n", file = report_file, append = TRUE)
    cat("  Up-regulated genes (padj < 0.05, log2FC > 1): ", up_genes, "\n", file = report_file, append = TRUE)
    cat("  Down-regulated genes (padj < 0.05, log2FC < -1): ", down_genes, "\n", file = report_file, append = TRUE)
    
    # SCN1A status
    scn1a_rows <- which(res_df$is_SCN1A)
    if (length(scn1a_rows) > 0) {
      cat("  SCN1A status:\n", file = report_file, append = TRUE)
      for (i in scn1a_rows) {
        cat("    Ensembl ID: ", rownames(res_df)[i], "\n", file = report_file, append = TRUE)
        cat("    log2FC: ", res_df$log2FoldChange[i], "\n", file = report_file, append = TRUE)
        cat("    padj: ", res_df$padj[i], "\n", file = report_file, append = TRUE)
        status <- "Not significant"
        if (!is.na(res_df$padj[i]) && res_df$padj[i] < 0.05) {
          if (res_df$log2FoldChange[i] > 1) {
            status <- "Significantly up-regulated"
          } else if (res_df$log2FoldChange[i] < -1) {
            status <- "Significantly down-regulated"
          } else {
            status <- "Significant but |log2FC| < 1"
          }
        }
        cat("    Status: ", status, "\n", file = report_file, append = TRUE)
      }
    } else {
      cat("  SCN1A not found in results\n", file = report_file, append = TRUE)
    }
    
    # Top genes
    top_up <- res_df %>%
      filter(!is.na(padj) & padj < 0.05 & log2FoldChange > 0.5) %>%
      arrange(padj) %>%
      head(5)
    
    top_down <- res_df %>%
      filter(!is.na(padj) & padj < 0.05 & log2FoldChange < -0.5) %>%
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
  cat("Debug logs directory: ", paths$debug, "\n", file = report_file, append = TRUE)
  
  cat("Analysis report created at", report_file, "\n")
  return(report_file)
}

# Main function
main <- function() {
  # Define input and output paths separately
  input_base_path <- "D:/20250415 RNA seq analysis"
  output_base_path <- "D:/"
  cat("Input path:", input_base_path, "\n")
  cat("Output path:", output_base_path, "\n")
  
  # Create organized output directories with date in path
  paths <- create_output_dirs(output_base_path)
  cat("Created output directories in:", paths$base, "\n")
  
  # Define a professional color palette
  saRNA_colors <- c("SCR" = "#404040", 
                   "saRNA7" = "#E69F00", 
                   "saRNA10" = "#56B4E9", 
                   "saRNA12" = "#009E73", 
                   "saRNA14" = "#CC79A7")
  
  # Start a log file
  log_file <- file.path(paths$debug, "analysis_log.txt")
  cat("RNA-seq analysis log\n", file = log_file)
  cat("Date:", Sys.Date(), "\n", file = log_file, append = TRUE)
  cat("R version:", R.version.string, "\n", file = log_file, append = TRUE)
  
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

    # Run DESeq2 analysis
  cat("Running DESeq2 analysis...\n", file = log_file, append = TRUE)
  analysis_results <- run_deseq2_analysis(counts_df, metadata_df, paths$debug)
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
    generate_ma_plots(dds, results_list, paths$ma_plots)
  }, error = function(e) {
    warning_msg <- paste("Error creating MA plots:", e$message)
    warning(warning_msg)
    cat("WARNING:", warning_msg, "\n", file = log_file, append = TRUE)
  })
  
  tryCatch({
    cat("Generating PCA plots...\n", file = log_file, append = TRUE)
    generate_pca_plots(dds, paths$pca, saRNA_colors)
  }, error = function(e) {
    warning_msg <- paste("Error creating PCA plots:", e$message)
    warning(warning_msg)
    cat("WARNING:", warning_msg, "\n", file = log_file, append = TRUE)
  })
  
  tryCatch({
    cat("Generating volcano plots with similar sequence gene labeling...\n", file = log_file, append = TRUE)
    generate_volcano_plots(results_list, paths$volcano_plots)
  }, error = function(e) {
    warning_msg <- paste("Error creating volcano plots:", e$message)
    warning(warning_msg)
    cat("WARNING:", warning_msg, "\n", file = log_file, append = TRUE)
  })
  
  tryCatch({
    cat("Generating heatmaps...\n", file = log_file, append = TRUE)
    generate_heatmaps(dds, results_list, paths$heatmaps, saRNA_colors)
  }, error = function(e) {
    warning_msg <- paste("Error creating heatmaps:", e$message)
    warning(warning_msg)
    cat("WARNING:", warning_msg, "\n", file = log_file, append = TRUE)
  })
  
  tryCatch({
    cat("Generating SCN1A expression plots...\n", file = log_file, append = TRUE)
    generate_scn1a_plots(dds, results_list, paths$scn1a, saRNA_colors)
  }, error = function(e) {
    warning_msg <- paste("Error creating SCN1A expression plots:", e$message)
    warning(warning_msg)
    cat("WARNING:", warning_msg, "\n", file = log_file, append = TRUE)
  })
  
  tryCatch({
    cat("Generating summary barplot...\n", file = log_file, append = TRUE)
    generate_summary_barplot(results_list, paths$base)
  }, error = function(e) {
    warning_msg <- paste("Error creating summary barplot:", e$message)
    warning(warning_msg)
    cat("WARNING:", warning_msg, "\n", file = log_file, append = TRUE)
  })
  
  # Run Gene Ontology analysis
  tryCatch({
    cat("Running Gene Ontology analysis...\n", file = log_file, append = TRUE)
    run_go_analysis(results_list, paths$go)
  }, error = function(e) {
    warning_msg <- paste("Error running Gene Ontology analysis:", e$message)
    warning(warning_msg)
    cat("WARNING:", warning_msg, "\n", file = log_file, append = TRUE)
  })
  
  # Create analysis report
  tryCatch({
    cat("Creating analysis report...\n", file = log_file, append = TRUE)
    report_file <- create_analysis_report(dds, results_list, paths)
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

# Run the script
main()