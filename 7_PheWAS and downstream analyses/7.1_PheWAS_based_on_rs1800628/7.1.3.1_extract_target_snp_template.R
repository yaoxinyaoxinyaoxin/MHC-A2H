#!/usr/bin/env Rscript
# ==============================================================================
# [Script]: 7.1.3.1_extract_target_snp_template.R
# [Method]: PheWAS Data Extraction
# [Step]: Extract Target SNP (e.g., rs1800628)
# 
# [Function]:
# Extract target SNP data from corresponding GWAS/eQTL databases for PheWAS.
# 
# [Input]:
#   --input_dir    : Directory or file containing the source summary statistics
#   --target_snp   : The rsID of the target SNP (default: rs1800628)
#   --out_dir      : Output directory path
# 
# [Output]:
# Extracted SNP statistics in CSV format.
# ==============================================================================

rm(list = ls())
gc()

# 1. Load Dependencies
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(parallel)
  if(requireNamespace("ieugwasr", quietly = TRUE)) library(ieugwasr)
})

# 2. Parse Command Line Arguments
option_list <- list(
  make_option(c("--input_dir"), type="character", default=NULL, help="Input directory/file with summary stats"),
  make_option(c("--target_snp"), type="character", default="rs1800628", help="Target SNP rsID"),
  make_option(c("--out_dir"), type="character", default="./PheWAS_Results", help="Output directory path")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_dir)) {
  print_help(opt_parser)
  stop("Missing required input directory/file.")
}

INPUT_DIR <- opt$input_dir
TARGET_SNP <- opt$target_snp
OUTPUT_DIR <- opt$out_dir

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

LOG_FILE <- file.path(OUTPUT_DIR, "extraction.log")
log_message <- function(msg) {
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  formatted_msg <- paste(timestamp, msg)
  cat(formatted_msg, "\n")
  write(formatted_msg, file = LOG_FILE, append = TRUE)
}

log_message(paste("Starting extraction for", TARGET_SNP))

start_time <- Sys.time()

# Unified working environment variables
input_dir <- "./"
output_base_dir <- INPUT_DIR
target_rsid <- "TARGET_SNP"
project_name <- "OneK1K_sc_eQTL_TARGET_SNP"

# Timestamp
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- file.path(output_base_dir, paste0(project_name, "_", timestamp))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Log directory
log_dir <- file.path(output_dir, "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

# readme
# Readme directory
readme_dir <- file.path(output_dir, "readme")
dir.create(readme_dir, recursive = TRUE, showWarnings = FALSE)

# Script directory
script_dir <- file.path(output_dir, "scripts")
dir.create(script_dir, recursive = TRUE, showWarnings = FALSE)

# Log file path
log_file <- file.path(log_dir, paste0(project_name, "_", timestamp, "_run_log.txt"))
sink(log_file, split = TRUE)

cat("========================================\n")
cat("Starting extraction of", target_rsid, "from OneK1K data\n")
cat("Start time:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Input directory:", input_dir, "\n")
cat("Output directory:", output_dir, "\n")
cat("========================================\n")

# Get file list
file_list <- list.files(input_dir, pattern = "\\.tsv\\.gz$", full.names = TRUE)
total_files <- length(file_list)
cat("Found", total_files, "files to process.\n")

total_extracted_records <- 0
processed_files <- 0
failed_files <- 0

# Process each file serially
for (i in seq_along(file_list)) {
  current_file <- file_list[i]
  file_name <- basename(current_file)
  
  cat(sprintf("[%d/%d] Processing file: %s\n", i, total_files, file_name))
  
  tryCatch({
    # Read full data
    dat <- fread(current_file)
    
    # RSID
    # Extract target RSID
    extracted_dat <- subset(dat, RSID == target_rsid)
    
    if (nrow(extracted_dat) > 0) {
      # , P_VALUE, P_VALUE
      # Deduplicate, keeping the minimum P_VALUE, and sort by P_VALUE
      extracted_dat <- extracted_dat %>%
        group_by(CELL_ID, CELL_TYPE, RSID, GENE, GENE_ID) %>%
        slice_min(order_by = P_VALUE, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        arrange(P_VALUE) # P  / Sort by P-value (ascending)
      
      # Generate output file name
      base_name <- tools::file_path_sans_ext(file_name) #  .gz
      base_name <- tools::file_path_sans_ext(base_name) #  .tsv
      out_file <- file.path(output_dir, paste0(base_name, ".csv"))
      
      # CSV
      # Save as CSV
      fwrite(extracted_dat, out_file)
      
      records_count <- nrow(extracted_dat)
      total_extracted_records <- total_extracted_records + records_count
      cat("  -> Extracted", records_count, "records and saved to", basename(out_file), "\n")
    } else {
      cat("  -> No data found for RSID:", target_rsid, "\n")
    }
    
    processed_files <- processed_files + 1
    
  }, error = function(e) {
    cat("  -> ERROR processing file:", file_name, "\n")
    cat("  -> Error message:", e$message, "\n")
    failed_files <<- failed_files + 1
  }, finally = {
    #  (, )
    # Clean memory (crucial as files are large)
    if (exists("dat")) rm(dat)
    if (exists("extracted_dat")) rm(extracted_dat)
    gc #  / Force garbage collection
  })
}

# End time
end_time <- Sys.time()
run_time <- end_time - start_time

cat("========================================\n")
cat("All files processed.\n")
cat("Total files:", total_files, "\n")
cat("Successfully processed:", processed_files, "\n")
cat("Failed:", failed_files, "\n")
cat("Total extracted records:", total_extracted_records, "\n")
cat("End time:", format(end_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("Run time:", run_time, "\n")
cat("========================================\n")

# readme
# Generate readme file
readme_content_cn <- paste0(
  "# : ", project_name, "\n",
  "# : ", format(Sys.time, "%Y-%m-%d"), "\n",
  "# :  OneK1K  ", target_rsid, "  eQTL ,  P . \n",
  "# : ", input_dir, "\n",
  "# : ", output_dir, "\n",
  "# :\n",
  "  - : ", total_files, "\n",
  "  - : ", processed_files, "\n",
  "  - : ", failed_files, "\n",
  "  - : ", total_extracted_records, "\n",
  "# : ", round(as.numeric(run_time, units = "mins"), 2), " \n"
)

readme_content_en <- paste0(
  "# Project Name: ", project_name, "\n",
  "# Date: ", format(Sys.time(), "%Y-%m-%d"), "\n",
  "# Description: Serially extract eQTL results for ", target_rsid, " from OneK1K dataset, deduplicate and sort by P-value ascending.\n",
  "# Input directory: ", input_dir, "\n",
  "# Output directory: ", output_dir, "\n",
  "# Statistics:\n",
  "  - Total files processed: ", total_files, "\n",
  "  - Successfully processed: ", processed_files, "\n",
  "  - Failed processing: ", failed_files, "\n",
  "  - Total extracted records: ", total_extracted_records, "\n",
  "# Run time: ", round(as.numeric(run_time, units = "mins"), 2), " mins\n"
)

readme_file <- file.path(readme_dir, paste0(project_name, "_", timestamp, "_readme.md"))
writeLines(c("###  (Chinese Version) ###", readme_content_cn, "\n", "###  (English Version) ###", readme_content_en), readme_file)

# Restore standard output
sink()

# scripts
# Copy the currently running script to the scripts directory in the output
current_script_path <- "./"
if(file.exists(current_script_path)) {
  file.copy(current_script_path, file.path(script_dir, "extract_OneK1K_TARGET_SNP.R"))
}
