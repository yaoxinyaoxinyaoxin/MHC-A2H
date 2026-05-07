#!/usr/bin/env Rscript
# ==============================================================================
# [Script]: 7.1.1.1_extract_target_snp_template.R
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

dir_gwas <- "./"
file_info <- "./"
out_base_dir <- INPUT_DIR
target_rsid <- "TARGET_SNP"

#  / Create output directory with timestamp
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path(out_base_dir, paste0("Extraction_", timestamp))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

#  results  / Create results subfolder for individual data
results_dir <- file.path(out_dir, "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

#  / Configure log file
log_dir <- file.path(out_dir, "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
log_file <- file.path(log_dir, paste0("UKB_protein_", target_rsid, "_extraction_", timestamp, ".log"))

write_log <- function(msg) {
  cat(paste0("[", Sys.time(), "] ", msg, "\n"), file = log_file, append = TRUE)
  cat(paste0("[", Sys.time(), "] ", msg, "\n"))
}

start_time <- Sys.time()
write_log("... / Starting analysis task...")

#  / Read phenotype information
write_log("... / Reading phenotype information file...")
info_df <- fread(file_info)

#  / Test mode configuration
test_mode <- FALSE
test_num <- 10

#  GWAS  / Get GWAS file list
gwas_files <- list.files(dir_gwas, pattern = "\\.tsv\\.gz$", full.names = TRUE)

if (test_mode) {
  write_log(paste0(",  ", test_num, " .  / Test mode enabled, only processing the first ", test_num, " files."))
  gwas_files <- head(gwas_files, test_num)
} else {
  write_log(paste0(" ", length(gwas_files), "  GWAS .  / Found ", length(gwas_files), " GWAS files."))
}

#  ( fread  subset) / Define extraction function (using fread and subset)
extract_rsid <- function(file_path) {
  gcst_id <- str_extract(basename(file_path), "GCST[0-9]+")
  res <- data.table()
  
  tryCatch({
    #  fread  / Read data using fread
    # , 
    dat <- fread(file_path, showProgress = FALSE)
    
    #  rs_id  / Ensure rs_id column exists
    if ("rs_id" %in% colnames(dat)) {
      #  subset  / Extract data using subset
      res <- subset(dat, rs_id == target_rsid)
      
      if (nrow(res) > 0) {
        #  / Merge phenotype info
        pheno_info <- info_df[STUDY == gcst_id, ]
        if (nrow(pheno_info) > 0) {
          res$GCST <- gcst_id
          res$Phenotype <- pheno_info$`DISEASE/TRAIT`[1]
        } else {
          res$GCST <- gcst_id
          res$Phenotype <- NA
        }
        
        #  results  / Save individual result to results subfolder
        indiv_out_file <- file.path(results_dir, paste0(gcst_id, "_", target_rsid, ".csv"))
        fwrite(res, indiv_out_file)
      }
    } else {
       write_log(paste0(": ", gcst_id, "  rs_id .  / Warning: rs_id column not found in file ", gcst_id))
    }
  }, error = function(e) {
    write_log(paste0(": ", gcst_id, " . / Error reading or processing file ", gcst_id, ": ", e$message))
  })
  
  return(res)
}

#  / Multi-threading processing with batching mechanism
num_cores <- 10 #  / Control number of threads
write_log(paste0(" ", num_cores, " ... / Using ", num_cores, " threads for parallel processing..."))

batch_size <- 50 #  / Set batch size
if (test_mode) batch_size <- 5 #  / Reduce batch size in test mode

total_files <- length(gwas_files)
num_batches <- ceiling(total_files / batch_size)
all_res_list <- list()

for (i in seq_len(num_batches)) {
  start_idx <- (i - 1) * batch_size + 1
  end_idx <- min(i * batch_size, total_files)
  batch_files <- gwas_files[start_idx:end_idx]
  
  write_log(paste0(" ", i, "/", num_batches, " ( ", start_idx, "-", end_idx, ")... / Processing batch ", i, "/", num_batches, "..."))
  
  #  mclapply  ( mc.preschedule = TRUE  FD leak)
  batch_res <- mclapply(batch_files, extract_rsid, mc.cores = num_cores, mc.preschedule = TRUE)
  
  #  data.table / Filter out empty data.tables
  batch_res <- Filter(function(x) nrow(x) > 0, batch_res)
  
  if (length(batch_res) > 0) {
    all_res_list <- append(all_res_list, batch_res)
  }
  
  #  / Proactive memory cleanup
  rm(batch_res, batch_files)
  gc()
}

#  / Summarize all results and write to file
write_log("... / Merging all extraction results...")
out_file <- file.path(out_dir, paste0("UKB_protein_", target_rsid, "_extraction_summary_", timestamp, ".csv"))

if (length(all_res_list) > 0) {
    final_res <- rbindlist(all_res_list, fill = TRUE)
    fwrite(final_res, out_file)
    total_records <- nrow(final_res)
    write_log(paste0(" ", total_records, " . : ", out_file, " / Total extracted ", total_records, " records. Results summarized to: ", out_file))
} else {
    total_records <- 0
    write_log(".  / No matching records found.")
}

#   / Save script to output directory (optional backup)
script_out_dir <- file.path(out_dir, "scripts")
dir.create(script_out_dir, recursive = TRUE, showWarnings = FALSE)
# , 
current_script <- tryCatch({
  normalizePath(sys.frame(1)$ofile)
}, error = function(e) {
  #  Rscript ,  commandArgs 
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    normalizePath(sub("--file=", "", file_arg[1]))
  } else {
    NA
  }
})
if (!is.na(current_script) && file.exists(current_script)) {
  file.copy(current_script, file.path(script_out_dir, "extract_TARGET_SNP.R"), overwrite = TRUE)
}

#  readme  / Generate readme and statistics
readme_dir <- file.path(out_dir, "readme")
dir.create(readme_dir, recursive = TRUE, showWarnings = FALSE)
readme_file <- file.path(readme_dir, paste0("UKB_protein_", target_rsid, "_", timestamp, "_readme.txt"))

end_time <- Sys.time()
run_time <- round(difftime(end_time, start_time, units = "mins"), 2)

readme_content <- c(
  "Project: UKB Plasma Protein GWAS extraction for TARGET_SNP",
  paste0("Date: ", Sys.time()),
  "Author: Trae",
  "Description: Extracted data for TARGET_SNP from 3049 UKB plasma protein GWAS summary statistics using fread and subset, with batching and memory cleanup.",
  paste0("Total GWAS files processed: ", length(gwas_files)),
  paste0("Total extracted records: ", total_records),
  paste0("Run time: ", run_time, " mins"),
  "Output files:",
  paste0("  - Individual results directory: ", results_dir),
  if(total_records > 0) paste0("  - Summary file: ", out_file) else "  - Summary file: None",
  "",
  "--------------------------------------------------",
  ": 3049 UKB  GWAS TARGET_SNP ",
  paste0(": ", Sys.time),
  ": Trae",
  ": fread  subset 3049UKBGWASTARGET_SNP, . . ",
  paste0("GWAS: ", length(gwas_files)),
  paste0(": ", total_records),
  paste0(": ", run_time, " "),
  ": " ,
  paste0("  - : ", results_dir),
  if(total_records > 0) paste0("  - : ", out_file) else "  - : "
)

writeLines(readme_content, readme_file)
write_log(".  / Analysis task completed.")
