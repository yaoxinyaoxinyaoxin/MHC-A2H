#!/usr/bin/env Rscript
# ==============================================================================
# [Script]: 8.1.1.1_process_GEO_data_template.R
# [Method]: Sensitivity and Robustness Analysis (GEO Validation)
# [Step]: 8.1.1_
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

project_root <- "./"
# Input Directories
raw_data_dir <- file.path(project_root, "5_GEO/GSE242252/1_/GSE242252_RAW")
series_matrix_path <- file.path(project_root, "5_GEO/GSE242252/1_/GSE242252_series_matrix.txt")
# Output Base Directory
output_base_dir <- file.path(project_root, "5_GEO/GSE242252/2_")

# Create Output Directory with Timestamp 
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- file.path(output_base_dir, paste0("GSE242252_Expression_Matrix_", timestamp))

if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
  message(paste("Created output directory:", output_dir))
}

# Initialize Log File 
log_file <- file.path(output_dir, "process_log.txt")
sink(log_file, split = TRUE) # Output to both console and file

cat("=================================================================\n")
cat("Script started at:", as.character(Sys.time()), "\n")
cat("Output Directory:", output_dir, "\n")
cat("=================================================================\n\n")

# =========================================================================
# 2. Process RAW Data 
# =========================================================================

cat("Step 1: Processing RAW data files...\n")

# List all .txt.gz files ( .txt.gz )
raw_files <- list.files(raw_data_dir, pattern = "\\.txt\\.gz$", full.names = TRUE)
cat("Found", length(raw_files), "RAW files.\n")

if (length(raw_files) == 0) {
  stop("No .txt.gz files found in RAW directory!")
}

# Function to read a single file 
# Returns a data.table with GeneID and Count
read_sample <- function(file_path) {
  # Extract GSM ID from filename (e.g., GSM7757087_CO_1-1Y_S20.txt.gz -> GSM7757087)
  filename <- basename(file_path)
  gsm_id <- str_extract(filename, "^GSM[0-9]+")
  
  # Read file using fread for speed ( fread )
  # Assuming 2 columns: GeneID, Count. No header.
  dt <- fread(file_path, header = FALSE, col.names = c("GeneID", "Count"))
  
  # Set the column name of Count to GSM ID ( Count  GSM ID)
  setnames(dt, "Count", gsm_id)
  
  return(dt)
}

# Read all files 
# We read the first file to establish the gene order, then cbind others if order matches
# This is much faster than merging 80 files
cat("Reading files and merging...\n")

# Read first file
first_file <- raw_files[1]
base_dt <- read_sample(first_file)
gene_ids <- base_dt$GeneID
n_genes <- length(gene_ids)
cat("Number of genes detected:", n_genes, "\n")

# Initialize matrix with the first sample 
expression_matrix <- base_dt

# Loop through the rest of the files 
for (i in 2:length(raw_files)) {
  current_file <- raw_files[i]
  temp_dt <- read_sample(current_file)
  
  # Validation: Check if GeneIDs match exactly (: GeneID )
  if (!identical(temp_dt$GeneID, gene_ids)) {
    warning(paste("GeneID mismatch in file:", basename(current_file), ". Attempting to merge safely."))
    # If mismatch, fall back to safe merge (slower but safe)
    expression_matrix <- merge(expression_matrix, temp_dt, by = "GeneID", all = TRUE)
  } else {
    # If match, just cbind the new column (, )
    # This assumes data.table structure. We remove GeneID from temp to avoid duplication
    new_col_name <- colnames(temp_dt)[2] # The GSM ID
    expression_matrix[, (new_col_name) := temp_dt[[2]]]
  }
  
  if (i %% 10 == 0) cat("Processed", i, "of", length(raw_files), "files...\n")
}

cat("Data merging completed.\n")
cat("Matrix dimensions:", nrow(expression_matrix), "genes x", ncol(expression_matrix)-1, "samples.\n")

# Save Expression Matrix 
output_matrix_file <- file.path(output_dir, "GSE242252_Gene_Expression_Matrix.txt")
fwrite(expression_matrix, output_matrix_file, sep = "\t")
cat("Saved expression matrix to:", output_matrix_file, "\n\n")

# =========================================================================
# 3. Process Metadata and Generate README ( README)
# =========================================================================

cat("Step 2: Processing metadata from series_matrix.txt...\n")

# Read series_matrix.txt ( series_matrix.txt)
# We use readLines to handle the variable header structure ( readLines )
lines <- readLines(series_matrix_path)

# Helper function to extract line by prefix (: )
get_line_content <- function(lines, prefix) {
  line <- lines[grep(paste0("^", prefix), lines)]
  if (length(line) == 0) return(NULL)
  # Remove prefix and quotes, split by tab (, )
  content <- str_split(line, "\t")[[1]]
  content <- content[-1] # Remove the label itself
  content <- str_remove_all(content, '"') # Remove quotes
  return(content)
}

# Extract relevant fields 
gsm_ids <- get_line_content(lines, "!Sample_geo_accession")
sample_titles <- get_line_content(lines, "!Sample_title")
sample_source <- get_line_content(lines, "!Sample_source_name_ch1")

# Extract characteristics (multiple lines) (, )
# Typically format is "!Sample_characteristics_ch1"
char_indices <- grep("^!Sample_characteristics_ch1", lines)
char_data_list <- list()

for (idx in char_indices) {
  line_content <- str_split(lines[idx], "\t")[[1]][-1]
  line_content <- str_remove_all(line_content, '"')
  # Usually format is "key: value"
  # We'll just keep them as columns Char_1, Char_2 etc. or try to parse key
  char_data_list[[paste0("Char_", idx)]] <- line_content
}

# Create Metadata Data Frame 
metadata_df <- data.frame(
  GSM_ID = gsm_ids,
  Sample_Title = sample_titles,
  Source = sample_source,
  stringsAsFactors = FALSE
)

# Add characteristics columns 
# Try to extract keys from the first sample to name columns
for (i in seq_along(char_data_list)) {
  vals <- char_data_list[[i]]
  # Parse key from first element e.g. "tissue: Whole blood" -> "tissue"
  first_val <- vals[1]
  if (grepl(":", first_val)) {
    col_name <- str_trim(str_split(first_val, ":")[[1]][1])
    # Remove key from values
    clean_vals <- str_trim(str_remove(vals, paste0("^", col_name, ":")))
    metadata_df[[col_name]] <- clean_vals
  } else {
    metadata_df[[paste0("Characteristic_", i)]] <- vals
  }
}

# Save Metadata/README (/README)
readme_file <- file.path(output_dir, "GSE242252_Sample_Metadata_README.txt")
write_tsv(metadata_df, readme_file)
cat("Saved metadata README to:", readme_file, "\n")

# Create a more human-readable README summary ( README )
summary_file <- file.path(output_dir, "README_Analysis_Summary.txt")
sink(summary_file)
cat("Project: GSE242252 Analysis\n")
cat("Date:", as.character(Sys.time()), "\n")
cat("------------------------------------------------\n")
cat("1. Data Processing:\n")
cat("   - Raw data merged from", length(raw_files), "files.\n")
cat("   - Total Genes:", nrow(expression_matrix), "\n")
cat("   - Total Samples:", ncol(expression_matrix) - 1, "\n")
cat("   - Output Matrix:", basename(output_matrix_file), "\n\n")
cat("2. Sample Naming Convention:\n")
cat("   - GSM_ID: GEO Sample Accession (Unique Identifier)\n")
cat("   - Sample_Title: Original sample name provided by submitter\n")
cat("   - Mapping Details are available in:", basename(readme_file), "\n\n")
cat("3. Sample Groups Overview:\n")
print(table(metadata_df$disease_state)) # Assuming 'disease state' or similar column exists
cat("------------------------------------------------\n")
sink()

# Restore log sink 
sink() 
sink(log_file, append = TRUE)

cat("\nStep 2 Completed.\n")

# =========================================================================
# 4. Finalization 
# =========================================================================

# Copy script to output directory 
current_script_path <- file.path(output_base_dir, "01_process_GSE242252_data.R")
if (file.exists(current_script_path)) {
  file.copy(current_script_path, file.path(output_dir, "01_process_GSE242252_data.R"))
  cat("Script copied to output directory.\n")
} else {
  warning("Could not find script file to copy: ", current_script_path)
}

cat("=================================================================\n")
cat("Analysis Finished Successfully.\n")
cat("Results saved in:", output_dir, "\n")
time_taken <- Sys.time()
cat("End Time:", as.character(time_taken), "\n")
cat("=================================================================\n")

# Close log 
sink()
