#!/usr/bin/env Rscript
# ==============================================================================
# [Script]: 8.1.2.1_split_expression_matrix_template.R
# [Method]: Sensitivity and Robustness Analysis (GEO Validation)
# [Step]: 8.1.2_(,,1)
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

options(stringsAsFactors = FALSE)

# Load libraries with checks
if (!require("data.table")) install.packages("data.table")
if (!require("dplyr")) install.packages("dplyr")
if (!require("stringr")) install.packages("stringr")
if (!require("readr")) install.packages("readr")

library(data.table)
library(dplyr)
library(stringr)
library(readr)

# Start timer
start_time <- Sys.time()
cat("Script started at:", as.character(start_time), "\n")

# -----------------------------------------------------------------------------
# 2. Define Paths
# -----------------------------------------------------------------------------
project_dir <- "./"
input_matrix_path <- file.path(project_dir, "2_/GSE242252_Expression_Matrix_20260118_184611/GSE242252_Gene_Expression_Matrix.txt")
input_metadata_path <- file.path(project_dir, "2_/GSE242252_Expression_Matrix_20260118_184611/GSE242252_Sample_Metadata_README.txt")
output_base_dir <- file.path(project_dir, "3_(,,1)")

# Create timestamped output directory
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- file.path(output_base_dir, paste0("Split_Matrix_", timestamp))

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Create log folder
log_dir <- file.path(output_dir, "logs")
if (!dir.exists(log_dir)) {
  dir.create(log_dir)
}

cat("Output directory created at:", output_dir, "\n")

# -----------------------------------------------------------------------------
# 3. Read Data
# -----------------------------------------------------------------------------
cat("Reading expression matrix...\n")
# Using fread for efficiency
# fread
expr_df <- fread(input_matrix_path, data.table = FALSE)

cat("Reading metadata...\n")
meta_df <- fread(input_metadata_path, data.table = FALSE)

# Check data structure (Head)
cat("Expression Matrix Head:\n")
print(head(expr_df[, 1:5]))
cat("Metadata Head:\n")
print(head(meta_df))

# -----------------------------------------------------------------------------
# 4. Data Cleaning (Gene IDs)
# -----------------------------------------------------------------------------
cat("Cleaning Gene IDs...\n")

# Remove version number from GeneID (assumed to be the first column)
# GeneID
# Check if first column is GeneID
gene_col_name <- colnames(expr_df)[1]
if (gene_col_name != "GeneID") {
  warning("First column is not named 'GeneID', assuming first column contains gene IDs.")
}

# Clean IDs: remove dot and digits after it
# ID: 
expr_df[[1]] <- gsub("\\.[0-9]+$", "", expr_df[[1]])

# Handle duplicates: aggregate by mean if duplicates exist after trimming
# : , 
if (any(duplicated(expr_df[[1]]))) {
  cat("Duplicate Gene IDs found after cleaning. Aggregating by mean...\n")
  # Melt, group by GeneID, summarize, then reshape back
  # , GeneID, , 
  # Note: This can be memory intensive. An alternative is distinct if values are identical, or just simple aggregation.
  # Using dplyr for aggregation
  expr_df <- expr_df %>%
    group_by(!!sym(colnames(expr_df)[1])) %>%
    summarise(across(where(is.numeric), mean)) %>%
    ungroup()
  cat("Aggregation complete.\n")
}

# Set row names
rownames(expr_df) <- expr_df[[1]]
expr_df <- expr_df[, -1] # Remove the GeneID column

# -----------------------------------------------------------------------------
# 5. Process Metadata & Define Groups
# -----------------------------------------------------------------------------
cat("Processing groups...\n")

# Ensure GSM_ID column exists in metadata
if (!"GSM_ID" %in% colnames(meta_df)) {
  stop("Metadata does not contain 'GSM_ID' column.")
}

# Ensure disease state column exists
if (!"disease state" %in% colnames(meta_df)) {
  stop("Metadata does not contain 'disease state' column.")
}

# Filter groups
# Group 1: Healthy Control
group_control_samples <- meta_df %>%
  filter(`disease state` == "Healthy control") %>%
  pull(GSM_ID)

# Group 2: Acute Phase (Herpes Zoster)
group_acute_samples <- meta_df %>%
  filter(`disease state` == "Herpes Zoster") %>%
  pull(GSM_ID)

# Group 3: 1 Year Post-Infection (Herpes Zoster resolved (1Y))
group_1y_samples <- meta_df %>%
  filter(`disease state` == "Herpes Zoster resolved (1Y)") %>%
  pull(GSM_ID)

cat("Samples found:\n")
cat("Control:", length(group_control_samples), "\n")
cat("Acute:", length(group_acute_samples), "\n")
cat("1 Year:", length(group_1y_samples), "\n")

# -----------------------------------------------------------------------------
# 6. Split Matrix
# -----------------------------------------------------------------------------
cat("Splitting matrices...\n")

# Helper function to subset
subset_matrix <- function(full_matrix, sample_ids, group_name) {
  # Intersect samples to ensure they exist in matrix
  valid_samples <- intersect(sample_ids, colnames(full_matrix))
  
  if (length(valid_samples) == 0) {
    warning(paste("No valid samples found for group:", group_name))
    return(NULL)
  }
  
  if (length(valid_samples) < length(sample_ids)) {
    missing <- setdiff(sample_ids, valid_samples)
    cat(paste("Warning: Missing samples in matrix for", group_name, ":", paste(missing, collapse=", "), "\n"))
  }
  
  sub_mat <- full_matrix[, valid_samples, drop = FALSE]
  # Add GeneID back as a column for output
  sub_mat <- data.frame(GeneID = rownames(sub_mat), sub_mat, check.names = FALSE)
  return(sub_mat)
}

mat_control <- subset_matrix(expr_df, group_control_samples, "Healthy Control")
mat_acute <- subset_matrix(expr_df, group_acute_samples, "Acute Phase")
mat_1y <- subset_matrix(expr_df, group_1y_samples, "1 Year Post")

# -----------------------------------------------------------------------------
# 7. Save Outputs
# -----------------------------------------------------------------------------
cat("Saving files...\n")

write_output <- function(df, filename) {
  if (!is.null(df)) {
    path <- file.path(output_dir, filename)
    fwrite(df, path, sep = "\t", quote = FALSE, row.names = FALSE)
    cat("Saved:", path, "\n")
  }
}

write_output(mat_control, "GSE242252_Expression_Control.txt")
write_output(mat_acute, "GSE242252_Expression_Acute_HZ.txt")
write_output(mat_1y, "GSE242252_Expression_HZ_Resolved_1Y.txt")

# -----------------------------------------------------------------------------
# 8. Generate Logs
# -----------------------------------------------------------------------------
end_time <- Sys.time()
run_time <- end_time - start_time

log_content <- c(
  paste("Project:", "8_"),
  paste("Date:", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  paste("Script:", "split_expression_matrix.R"),
  paste("Input Matrix:", input_matrix_path),
  paste("Input Metadata:", input_metadata_path),
  "--------------------------------------------------",
  "Execution Summary:",
  paste("Total Genes processed:", nrow(expr_df)),
  paste("Control Samples:", ncol(mat_control) - 1),
  paste("Acute HZ Samples:", ncol(mat_acute) - 1),
  paste("Resolved 1Y Samples:", ncol(mat_1y) - 1),
  paste("Total Run Time:", round(run_time, 2), units(run_time)),
  "--------------------------------------------------",
  "Outputs:",
  file.path(output_dir, "GSE242252_Expression_Control.txt"),
  file.path(output_dir, "GSE242252_Expression_Acute_HZ.txt"),
  file.path(output_dir, "GSE242252_Expression_HZ_Resolved_1Y.txt")
)

# Save Log
log_file <- file.path(log_dir, paste0("GSE242252_Analysis_Log_", timestamp, ".txt"))
writeLines(log_content, log_file)

# Copy this script to output directory
script_path <- commandArgs(trailingOnly = FALSE)[4] # Try to get current script path
# Since we are running interactively or via wrapper, we might not get the exact path easily this way.
# Instead, we will write the current content to a file in the output dir if we are generating it.
# In this environment, I will save the file explicitly.

cat("Process completed successfully.\n")
