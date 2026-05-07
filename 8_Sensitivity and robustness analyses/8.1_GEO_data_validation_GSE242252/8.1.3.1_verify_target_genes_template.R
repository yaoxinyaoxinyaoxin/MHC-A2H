#!/usr/bin/env Rscript
# ==============================================================================
# [Script]: 8.1.3.1_verify_target_genes_template.R
# [Method]: Sensitivity and Robustness Analysis (GEO Validation)
# [Step]: 8.1.3_
# 
# [Function]:
# Parameterized script for GEO data processing and gene validation.
# 
# [Input]:
#   --input_path   : Input directory or file path
#   --out_dir      : Output directory path
# 
# [Output]:
# Processed matrices or validation results.
# ==============================================================================

rm(list = ls())
gc()

# 1. Load Dependencies
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(data.table)
  library(tibble)
  if(requireNamespace("limma", quietly = TRUE)) library(limma)
  if(requireNamespace("ggplot2", quietly = TRUE)) library(ggplot2)
  if(requireNamespace("ggpubr", quietly = TRUE)) library(ggpubr)
})

# 2. Parse Command Line Arguments
option_list <- list(
  make_option(c("--input_path"), type="character", default=NULL, help="Input file or directory path"),
  make_option(c("--gene_list"), type="character", default=NULL, help="Path to target genes list (if applicable)"),
  make_option(c("--out_dir"), type="character", default="./GEO_Results", help="Output directory path")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_path)) {
  print_help(opt_parser)
  stop("Missing required input arguments.")
}

INPUT_PATH <- opt$input_path
GENE_LIST <- opt$gene_list
OUTPUT_DIR <- opt$out_dir

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

base_dir <- "./"
matrix_dir <- "./"
output_dir <- base_dir

# Set working directory

# Create logs directory
log_dir <- file.path(base_dir, "logs")
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)

# Setup logging
log_file <- file.path(log_dir, paste0("GSE242252_Gene_Verification_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
sink(log_file, split = TRUE)

cat(paste0("[", Sys.time(), "] Script started.\n"))

# Load libraries
cat(paste0("[", Sys.time(), "] Loading libraries...\n"))
suppressPackageStartupMessages({
  library(tidyverse)
  library(EnsDb.Hsapiens.v75)
  library(ggpubr)
  library(readr)
  library(openxlsx)
  library(data.table)
})

# ------------------------------------------------------------------------------
# 2. Read Gene List & Map IDs
# 2. ID
# ------------------------------------------------------------------------------
cat(paste0("[", Sys.time(), "] Reading gene list...\n"))

gene_list_file <- "Gene_List_Single_Outcome_Herpes_Zoster.txt"
if (!file.exists(gene_list_file)) {
  stop(paste("Gene list file not found:", gene_list_file))
}

# Read gene symbols
genes_df <- read_tsv(gene_list_file, col_names = FALSE, show_col_types = FALSE)
gene_symbols <- genes_df[[1]]
cat(paste0("Number of genes in list: ", length(gene_symbols), "\n"))

# Map to Ensembl IDs
# Ensembl ID
cat(paste0("[", Sys.time(), "] Mapping symbols to Ensembl IDs using EnsDb.Hsapiens.v75...\n"))

# Use mapIds for conversion
# mapIds
# Note: Using v75 as requested. Check if valid keys exist.
tryCatch({
  gene_ids <- mapIds(EnsDb.Hsapiens.v75, 
                     keys = gene_symbols, 
                     column = "GENEID", 
                     keytype = "SYMBOL", 
                     multiVals = "first")
  
  # Create a mapping dataframe
  gene_map <- data.frame(
    Symbol = gene_symbols,
    Ensembl_ID = gene_ids,
    stringsAsFactors = FALSE
  )
  
  # Filter out unmapped genes
  valid_map <- gene_map %>% dplyr::filter(!is.na(Ensembl_ID))
  
  cat(paste0("Successfully mapped ", nrow(valid_map), " out of ", length(gene_symbols), " genes.\n"))
  
  if (nrow(valid_map) == 0) {
    stop("No valid Ensembl IDs found for the provided gene symbols.")
  }
  
  # Save mapping info
  write.xlsx(gene_map, file.path(output_dir, "Gene_ID_Mapping_Log.xlsx"))
  
}, error = function(e) {
  cat("Error in ID mapping: ", e$message, "\n")
  stop(e)
})

target_ids <- valid_map$Ensembl_ID

# ------------------------------------------------------------------------------
# 3. Extract Expression Data
# 3. 
# ------------------------------------------------------------------------------
cat(paste0("[", Sys.time(), "] Extracting expression data from matrix files...\n"))

# Define file paths for each group
files <- list(
  Control = file.path(matrix_dir, "GSE242252_Expression_Control.txt"),
  Acute_HZ = file.path(matrix_dir, "GSE242252_Expression_Acute_HZ.txt"),
  Resolved_1Y = file.path(matrix_dir, "GSE242252_Expression_HZ_Resolved_1Y.txt")
)

combined_data <- data.frame()

for (group in names(files)) {
  fpath <- files[[group]]
  cat(paste0("Processing ", group, " from ", basename(fpath), "...\n"))
  
  if (!file.exists(fpath)) {
    warning(paste("File not found:", fpath))
    next
  }
  
  # Read header to get sample names
  # Using data.table::fread for efficiency, looking for target IDs
  # fread
  
  # Read only rows matching target IDs would be ideal, but grep might be slow if file is huge.
  # Given the file size is likely manageable, read and filter.
  # Using a system command to grep headers + matching lines is faster for huge files.
  # , grep
  
  # Construct grep pattern
  # grep
  # We need to match GeneID column. Pattern: "^ENSG..."
  # Use a temporary file to store filtered content
  
  temp_file <- file.path(output_dir, paste0("temp_", group, ".txt"))
  
  # Get header
  header_cmd <- paste0("head -n 1 \"", fpath, "\" > \"", temp_file, "\"")
  system(header_cmd)
  
  # Grep gene IDs
  # Create a pattern file for grep -F
  pattern_file <- file.path(output_dir, "grep_patterns.txt")
  writeLines(target_ids, pattern_file)
  
  # Grep matching lines and append to temp file
  # grep -F -f pattern_file file >> temp_file
  grep_cmd <- paste0("grep -F -f \"", pattern_file, "\" \"", fpath, "\" >> \"", temp_file, "\"")
  system(grep_cmd)
  
  # Read the filtered data
  dt <- fread(temp_file)
  
  # Clean up temp files
  unlink(temp_file)
  unlink(pattern_file)
  
  if (nrow(dt) > 0) {
    # Reshape to long format
    # Columns are: GeneID, Sample1, Sample2...
    dt_long <- dt %>%
      pivot_longer(cols = -GeneID, names_to = "Sample", values_to = "Expression") %>%
      mutate(Group = group)
    
    combined_data <- bind_rows(combined_data, dt_long)
  } else {
    cat(paste0("Warning: No matching genes found in ", group, "\n"))
  }
}

if (nrow(combined_data) == 0) {
  stop("No expression data extracted.")
}

# Join with gene symbols
final_data <- combined_data %>%
  left_join(valid_map, by = c("GeneID" = "Ensembl_ID")) %>%
  dplyr::select(Symbol, GeneID, Sample, Group, Expression)

# Ensure Group order
final_data$Group <- factor(final_data$Group, levels = c("Control", "Acute_HZ", "Resolved_1Y"))

cat(paste0("Total data points: ", nrow(final_data), "\n"))

# ------------------------------------------------------------------------------
# 4. Statistical Analysis
# 4. 
# ------------------------------------------------------------------------------
cat(paste0("[", Sys.time(), "] Performing statistical analysis...\n"))

# Define comparisons
my_comparisons <- list(
  c("Control", "Acute_HZ"),
  c("Control", "Resolved_1Y"),
  c("Acute_HZ", "Resolved_1Y")
)

# Function to perform stats for each gene
perform_stats <- function(df) {
  # Kruskal-Wallis test (global)
  kw_res <- kruskal.test(Expression ~ Group, data = df)
  
  # Pairwise Wilcoxon test
  pairwise_res <- pairwise.wilcox.test(df$Expression, df$Group, p.adjust.method = "BH")
  
  # Extract p-values
  p_acute_vs_ctrl <- pairwise_res$p.value["Acute_HZ", "Control"]
  p_resolved_vs_ctrl <- pairwise_res$p.value["Resolved_1Y", "Control"]
  p_resolved_vs_acute <- pairwise_res$p.value["Resolved_1Y", "Acute_HZ"]
  
  # Calculate Means
  means <- df %>%
    group_by(Group) %>%
    summarise(Mean_Expr = mean(Expression), .groups = "drop") %>%
    pivot_wider(names_from = Group, values_from = Mean_Expr, names_prefix = "Mean_")
  
  tibble(
    KW_pvalue = kw_res$p.value,
    P_Acute_vs_Control = p_acute_vs_ctrl,
    P_Resolved_vs_Control = p_resolved_vs_ctrl,
    P_Resolved_vs_Acute = p_resolved_vs_acute
  ) %>% bind_cols(means)
}

# Run stats per gene
stats_results <- final_data %>%
  group_by(Symbol, GeneID) %>%
  do(perform_stats(.)) %>%
  ungroup()

# Save stats
stats_file <- file.path(output_dir, "Gene_Verification_Statistics.xlsx")
write.xlsx(stats_results, stats_file)
cat(paste0("Statistics saved to: ", stats_file, "\n"))

# ------------------------------------------------------------------------------
# 5. Visualization
# 5. 
# ------------------------------------------------------------------------------
cat(paste0("[", Sys.time(), "] Generating boxplots...\n"))

# Plotting function
plot_gene <- function(gene_sym, data_subset) {
  p <- ggboxplot(data_subset, x = "Group", y = "Expression",
                 color = "Group", palette = "jco",
                 add = "jitter",
                 title = gene_sym,
                 ylab = "Expression Level", xlab = "Group") +
    stat_compare_means(comparisons = my_comparisons, method = "wilcox.test", label = "p.signif") +
    stat_compare_means(label.y.npc = "bottom") + # Add global p-value
    theme_minimal() +
    theme(legend.position = "none",
          plot.title = element_text(hjust = 0.5, face = "bold"))
  return(p)
}

# Generate plots for all genes
pdf_file <- file.path(output_dir, "Gene_Expression_Boxplots.pdf")
pdf(pdf_file, width = 8, height = 6)

genes <- unique(final_data$Symbol)
for (g in genes) {
  sub_df <- final_data %>% dplyr::filter(Symbol == g)
  print(plot_gene(g, sub_df))
}

dev.off()
cat(paste0("Plots saved to: ", pdf_file, "\n"))


# ------------------------------------------------------------------------------
# 6. Finalizing
# 6. 
# ------------------------------------------------------------------------------
end_time <- Sys.time()
cat(paste0("[", end_time, "] Analysis complete.\n"))
cat(paste0("Results directory: ", output_dir, "\n"))

sink()

# Create Readme for this analysis
# Readme
readme_content <- paste0(
  "# Gene Expression Verification Analysis / \n\n",
  "## Date / \n", Sys.Date, "\n\n",
  "## Description / \n",
  "This analysis verifies the expression of specific HZ-related genes across three cohorts: Control, Acute HZ, and Resolved 1Y.\n",
  "HZ（、、1）. \n\n",
  "## Input / \n",
  "- Gene List: Gene_List_Single_Outcome_Herpes_Zoster.txt\n",
  "- Expression Data: Split_Matrix_20260118_190833/\n\n",
  "## Output / \n",
  "- Gene_ID_Mapping_Log.xlsx: Mapping between symbols and Ensembl IDs.\n",
  "- Gene_Verification_Statistics.xlsx: Statistical test results (Wilcoxon/Kruskal).\n",
  "- Gene_Expression_Boxplots.pdf: Boxplots for each gene.\n",
  "- logs/: Execution logs.\n\n",
  "## Methods / \n",
  "- ID Mapping: EnsDb.Hsapiens.v75\n",
  "- Statistics: Wilcoxon rank-sum test (pairwise), Kruskal-Wallis test (global).\n"
)

writeLines(readme_content, file.path(output_dir, "README_Verification.md"))
