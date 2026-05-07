#!/usr/bin/env Rscript
# ==============================================================================
# [Script]: 7.6.1.4_Local_LD_Clumping_template.R
# [Method]: cis-MR IV Selection Pipeline Step
# [Step]: 7.6.1.4_Local_LD_Clumping.R
# 
# [Function]:
# Parameterized script for cis-MR instrumental variable selection.
# 
# [Input]:
#   --input_path   : Input directory or file path
#   --out_dir      : Output directory path
#   --bfile_local  : Local LD reference panel (hg19)
#   --plink_bin    : Path to plink binary
# 
# [Output]:
# Processed results for the next pipeline step.
# ==============================================================================

rm(list = ls())
gc()

# 1. Load Dependencies
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(data.table)
  library(readr)
  library(stringr)
  library(fs)
  library(parallel)
  # Include specific libs that might be used
  if(requireNamespace("EnsDb.Hsapiens.v75", quietly = TRUE)) library(EnsDb.Hsapiens.v75)
  if(requireNamespace("TxDb.Hsapiens.UCSC.hg19.knownGene", quietly = TRUE)) library(TxDb.Hsapiens.UCSC.hg19.knownGene)
  if(requireNamespace("org.Hs.eg.db", quietly = TRUE)) library(org.Hs.eg.db)
  if(requireNamespace("AnnotationDbi", quietly = TRUE)) library(AnnotationDbi)
  if(requireNamespace("TwoSampleMR", quietly = TRUE)) library(TwoSampleMR)
  if(requireNamespace("ieugwasr", quietly = TRUE)) library(ieugwasr)
})

# 2. Parse Command Line Arguments
option_list <- list(
  make_option(c("--input_path"), type="character", default=NULL, help="Input file or directory path"),
  make_option(c("--gwas_dir"), type="character", default=NULL, help="GWAS summary stats directory"),
  make_option(c("--out_dir"), type="character", default="./Output", help="Output directory path"),
  make_option(c("--bfile_local"), type="character", default=NULL, help="Local LD reference panel (hg19) for clumping"),
  make_option(c("--plink_bin"), type="character", default=NULL, help="Path to plink binary for clumping"),
  make_option(c("--failed_list"), type="character", default=NULL, help="List of failed files for fallback")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Assign common parameters
INPUT_PATH <- opt$input_path
GWAS_DIR <- opt$gwas_dir
OUTPUT_DIR <- opt$out_dir
BFILE_LOCAL <- opt$bfile_local
PLINK_BIN <- opt$plink_bin
FAILED_LIST <- opt$failed_list

dir_create(OUTPUT_DIR)

LOG_FILE <- file.path(OUTPUT_DIR, "analysis.log")
log_message <- function(msg) {
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  formatted_msg <- paste(timestamp, msg)
  cat(formatted_msg, "\n")
  write(formatted_msg, file = LOG_FILE, append = TRUE)
}

log_message(paste("Starting", "7.6.1.4_Local_LD_Clumping.R"))
\n", start_time))

# 4.  (Define input and output paths)
source_dir <- INPUT_PATH
base_out_dir <- OUTPUT_DIR

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- file.path(base_out_dir, paste0("Analysis_Result_", timestamp))
csv_out_dir <- file.path(output_dir, "csv_results")
readme_dir <- file.path(output_dir, "readme")
log_dir <- file.path(output_dir, "logs")

#  (Create directory structure)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(csv_out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(readme_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)

#  (Configure log output)
log_file <- file.path(log_dir, "run_log.txt")
sink(log_file, split = TRUE)

cat(sprintf(" (Output directories initialized at): %s\n", output_dir))

# 5. LD (Set LD clumping parameters)
clump_kb <- 10000
clump_r2 <- 0.3
clump_p <- 1
ld_ref_panel <- "./"
plink_bin <- PLINK_BIN

# 6.  (Get all source files to process)
files <- list.files(source_dir, pattern = "\\.csv$", full.names = TRUE)
total_files <- length(files)
cat(sprintf(" %d  (Found %d source data files to process).\n", total_files, total_files))

#  (Initialize record lists)
summary_list <- list()
skipped_ld_files <- c()

# 7.  (Iterate and process each file)
for (i in seq_along(files)) {
  file_path <- files[i]
  file_name <- basename(file_path)
  
  cat(sprintf("[%d/%d]  (Processing file): %s\n", i, total_files, file_name))
  
  #  (Read data)
  #  fread  (Use fread for efficient reading)
  IV1 <- tryCatch({
    as.data.frame(fread(file_path))
  }, error = function(e) {
    cat(sprintf(" (Failed to read file) %s: %s\n", file_name, conditionMessage(e)))
    return(NULL)
  })
  
  if (is.null(IV1) || nrow(IV1) == 0) {
    cat(sprintf(" %s ,  (File is empty or failed to read, skipping)\n", file_name))
    next
  }
  
  orig_nrow <- nrow(IV1)
  final_nrow <- orig_nrow
  status_msg <- ""
  
  # : 1SNP vs SNP (Classification: 1 SNP vs Multiple SNPs)
  if (orig_nrow > 1) {
    # LD (LD clumping retry mechanism)
    ld_retry_count <- 0
    ld_max_retries <- 3
    ld_success <- FALSE
    
    while (!ld_success && ld_retry_count < ld_max_retries) {
      tryCatch({
        #  (Prepare input data for clumping)
        IV1_clump <- data.frame(
          rsid = IV1$SNP,
          pval = IV1$pval.exposure
        )
        
        # LD (Perform local LD clumping)
        IV1_clumped <- ieugwasr::ld_clump_local(
          dat = IV1_clump,
          clump_kb = clump_kb,
          clump_r2 = clump_r2,
          clump_p = clump_p,
          bfile = ld_ref_panel,
          plink_bin = plink_bin
        )
        
        #  (Filter original data based on clumping results)
        IV1 <- IV1[IV1$SNP %in% IV1_clumped$rsid, ]
        final_nrow <- nrow(IV1)
        
        cat(sprintf("  -> 3-LD %d SNP (Step 3 - Retained %d SNPs after LD clumping)\n", final_nrow, final_nrow))
        ld_success <- TRUE
        status_msg <- "LD Clumped"
        
      }, error = function(e) {
        ld_retry_count <<- ld_retry_count + 1
        if (ld_retry_count < ld_max_retries) {
          cat(sprintf("  -> LD%d,  (LD clump retry %d failed, error): %s\n", ld_retry_count, ld_retry_count, conditionMessage(e)))
          Sys.sleep(2)  # 2 (Wait 2 seconds before retrying)
        } else {
          cat(sprintf("  -> LD%d, LD,  (LD clump failed after %d retries, skipping LD step, error): %s\n", ld_max_retries, ld_max_retries, conditionMessage(e)))
          # LD (Record skipped LD files)
          skipped_ld_files <<- c(skipped_ld_files, file_name)
          cat(sprintf("  -> 3-LD,  %d SNP (Step 3 - Skipped LD clumping due to retry failure, retained %d SNPs)\n", orig_nrow, orig_nrow))
        }
      })
    }
    
    if (!ld_success) {
      status_msg <- "LD Failed (Skipped)"
    }
    
  } else {
    cat(sprintf("  -> 3-LD（SNP≤1）,  %d SNP (Step 3 - Skipped LD clumping since SNP count <= 1, retained %d SNPs)\n", orig_nrow, orig_nrow))
    status_msg <- "Skipped (<=1 SNP)"
  }
  
  #  (Save processed results)
  out_file <- file.path(csv_out_dir, file_name)
  fwrite(IV1, out_file, row.names = FALSE)
  
  #  (Record to summary table)
  summary_list[[length(summary_list) + 1]] <- data.frame(
    File_Name = file_name,
    Original_SNPs = orig_nrow,
    Final_SNPs = final_nrow,
    Status = status_msg,
    stringsAsFactors = FALSE
  )
}

# 8.  (Save summary results and record files)
if (length(summary_list) > 0) {
  summary_df <- bind_rows(summary_list)
  summary_file <- file.path(readme_dir, "LD_Clumping_Summary.csv")
  fwrite(summary_df, summary_file, row.names = FALSE)
  cat(sprintf("\n (Overall summary saved to): %s\n", summary_file))
}

#  (Save failed records)
skipped_file_path <- file.path(readme_dir, "Failed_LD_Files_Record.txt")
if (length(skipped_ld_files) > 0) {
  writeLines(skipped_ld_files, skipped_file_path)
  cat(sprintf(" %d LD,  (A total of %d files failed LD clumping, records saved to): %s\n", length(skipped_ld_files), length(skipped_ld_files), skipped_file_path))
} else {
  cat("SNPLDSNP. (All multi-SNP files successfully clumped or none existed.)\n")
  writeLines("No files failed LD clumping.", skipped_file_path)
}

# 9. README (Generate README file)
readme_content <- sprintf(
" (Project Name): LD (Local LD Clumping for Instrumental Variables)
 (Date): %s
 (Output Path): %s

##  (Task Description)
, LD（PLINK）. 
Local LD clumping of preliminarily filtered instrumental variables using PLINK.
 (Parameters): clump_kb = %d, clump_r2 = %f, clump_p = %d

##  (Workflow)
1. CSV (Read all CSV files from source directory).
2. SNP,  (If only 1 SNP, copy and save directly).
3. >1SNP, LD (If >1 SNPs, perform local LD clumping).
4. 3. ,  (Includes 3-retry mechanism. If retries fail, retain original file and record as failed).

##  (Processing Statistics)
-  (Total processed files): %d
- LD (Number of files failed LD clumping): %d
", 
  Sys.time(), output_dir, clump_kb, clump_r2, clump_p, 
  total_files, length(skipped_ld_files)
)

writeLines(readme_content, file.path(readme_dir, "README_LD_Clumping.md"))

# 10.  (Record end time and total duration)
end_time <- Sys.time()
cat(sprintf("\n (Script end time): %s\n", end_time))
cat(sprintf(" (Total duration): %s\n", round(difftime(end_time, start_time, units = "mins"), 2), " (mins)\n"))

#  (Close log output)
sink()
