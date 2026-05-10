#!/usr/bin/env Rscript
# ==============================================================================
# [Script]: 7.1.2.1_extract_target_snp_template.R
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

})

# 1. 
# 1. Unified workspace environment and directory setup
target_snp <- "TARGET_SNP"
work_dir <- "./"
gwas_dir <- INPUT_DIR
meta_file <- file.path(work_dir, "834_UKB_.csv")


# Create timestamped output directory
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
project_name <- "834UKB_TARGET_SNP"
out_dir <- file.path(work_dir, paste0("Analysis_", project_name, "_", timestamp))
ind_res_dir <- file.path(out_dir, "results")
log_dir <- file.path(out_dir, "logs")
tmp_dir <- file.path(out_dir, "temp")

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(ind_res_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(log_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(tmp_dir, showWarnings = FALSE, recursive = TRUE)

# Set log file
log_file <- file.path(log_dir, paste0(project_name, "_", Sys.Date(), "_", timestamp, "_extract_log.txt"))
log_msg <- function(msg) {
  time_str <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(paste0("[", time_str, "] ", msg, "\n"), file = log_file, append = TRUE)
  cat(paste0("[", time_str, "] ", msg, "\n"))
}

log_msg(" / Analysis started.")
start_time <- Sys.time()

# 2. 
# 2. Read phenotype metadata
log_msg(" / Reading metadata file.")
meta <- fread(meta_file)
#  / Check column names
if (!"STUDY ACCESSION" %in% names(meta)) {
  stop(" (STUDY ACCESSION) / Metadata is missing required column (STUDY ACCESSION).")
}
if (!"DISEASE/TRAIT" %in% names(meta)) {
  stop(" (DISEASE/TRAIT) / Metadata is missing required column (DISEASE/TRAIT).")
}

# 3. 
# 3. Get source data file list
log_msg("GWAS / Scanning GWAS source data directory.")
gwas_files <- list.files(gwas_dir, pattern = "\\.tsv\\.gz$", full.names = TRUE)
    log_msg(sprintf(" %d  / Found %d files in total.", length(gwas_files), length(gwas_files)))

if(length(gwas_files) == 0) {
  stop(" .tsv.gz , .  / No .tsv.gz files found, please check the path.")
}

# 4.  
# 4. Core processing function (with error handling)
process_file <- function(file_path) {
  #  worker  1
  # Set inner parallel threads to 1 for each worker
  setDTthreads(1)
  
  tryCatch({
    # GCST
    # Extract GCST number
    gcst_id <- str_extract(basename(file_path), "GCST[0-9]+")
    
    # awk, SNP, 
    # Temporary unzip and filter using awk to extract target SNP, reducing memory and time
    temp_unzip <- file.path(tmp_dir, paste0(gcst_id, "_", Sys.getpid(), "_", sample(1e6, 1), ".tsv"))
    
    # awk ,  rsid 
    # awk script, matching the target SNP from the rsid column
    awk_script <- sprintf('NR==1 {for(i=1;i<=NF;i++) {if($i=="rsid") col=i}; print $0} NR>1 && col!="" && $col=="%s" {print $0}', target_snp)
    cmd <- paste("gzcat", shQuote(file_path), "| awk -F'\\t'", shQuote(awk_script), ">", shQuote(temp_unzip))
    system(cmd)
    
    if (!file.exists(temp_unzip) || file.info(temp_unzip)$size == 0) {
      if (file.exists(temp_unzip)) file.remove(temp_unzip)
      return(NULL)
    }
    
    # Read the filtered data
    dat <- fread(temp_unzip)
    file.remove(temp_unzip) #  / Delete immediately after reading
    
    # SNP
    # If data contains target SNP
    if (nrow(dat) > 0) {
      # , GCST
      # Find phenotype information, primary key is GCST ID
      pheno <- meta[["DISEASE/TRAIT"]][meta[["STUDY ACCESSION"]] == gcst_id][1]
      
      # Add new columns
      dat$GCST <- gcst_id
      dat$Phenotype <- pheno
      
      # Save individual result
      out_name <- file.path(ind_res_dir, paste0(gcst_id, "_", target_snp, ".csv"))
      fwrite(dat, out_name)
      
      #  / Clean memory
      res <- copy(dat)
      rm(dat)
      gc()
      
      return(res)
    }
    
    return(NULL)
  }, error = function(e) {
    # Error logging
    err_msg <- paste("Error processing file", file_path, ":", e$message)
    cat(err_msg, file = file.path(log_dir, "error_log.txt"), append = TRUE, sep = "\n")
    return(NULL)
  })
}

# 5.  (, )
# 5. Batch parallel processing (Sequential between batches, parallel within batches)
# : 10worker10 / Outer parallel: 10 workers processing 10 files
outer_cores <- 10
batch_size <- 10  # 10 / Process 10 files per batch
num_batches <- ceiling(length(gwas_files) / batch_size)

log_msg(sprintf(",  %d ,  %d  / Starting batch processing: %d batches, %d files per batch.", 
                num_batches, batch_size, num_batches, batch_size))
log_msg(sprintf(": 10worker, worker 1 / Parallel setting: 10 workers, 1 core per worker."))

all_results_list <- list()

for (i in 1:num_batches) {
  batch_start <- (i - 1) * batch_size + 1
  batch_end <- min(i * batch_size, length(gwas_files))
  batch_files <- gwas_files[batch_start:batch_end]
  
  log_msg(sprintf(" %d/%d  ( %d  %d) / Processing batch %d/%d (files %d to %d)...", 
                  i, num_batches, batch_start, batch_end, i, num_batches, batch_start, batch_end))
  
  #  mclapply ,  mc.preschedule = TRUE  FD 
  # Use mclapply for parallel processing with mc.preschedule = TRUE to prevent FD leaks
  batch_res <- mclapply(batch_files, process_file, mc.cores = outer_cores, mc.preschedule = TRUE)
  
  #  NULL 
  # Remove NULL results
  batch_res <- Filter(Negate(is.null), batch_res)
  
  if (length(batch_res) > 0) {
    all_results_list <- append(all_results_list, batch_res)
  }
  
  # Clear memory
  rm(batch_res, batch_files)
  gc()
}

# 6. 
# 6. Summarize all results
log_msg(" / Summarizing all data.")
if (length(all_results_list) > 0) {
  final_dat <- rbindlist(all_results_list, use.names = TRUE, fill = TRUE)
  
  # ,  GCST  Phenotype 
  # Adjust column order to put GCST and Phenotype first
  cols <- c("GCST", "Phenotype", setdiff(names(final_dat), c("GCST", "Phenotype")))
  setcolorder(final_dat, cols)
  
  final_out_file <- file.path(out_dir, paste0(project_name, "_", Sys.Date(), "_", timestamp, "_", target_snp, "_summary.csv"))
  fwrite(final_dat, final_out_file)
  log_msg(sprintf(": %s / Summary result saved to: %s", final_out_file, final_out_file))
} else {
  log_msg("SNP / Target SNP data not found in any files.")
}

# 7. Generate statistics
log_msg(" / Generating statistics.")
stats_file <- file.path(log_dir, paste0(project_name, "_", Sys.Date(), "_", timestamp, "_stats.txt"))

end_time <- Sys.time()
run_time <- difftime(end_time, start_time, units = "mins")

stats_text <- c(
  "==================================================",
  " / Data Extraction Summary",
  "==================================================",
  sprintf(" / Project Name: %s", project_name),
  sprintf(" / Analysis Date: %s", Sys.Date()),
  sprintf(" / Analysis Time: %s", format(Sys.time(), "%Y-%m-%d %H:%M:%S")),
  sprintf("SNP / Target SNP: %s", target_snp),
  sprintf(" / Total files processed: %d", length(gwas_files)),
  sprintf("SNP / Files successfully extracted: %d", length(all_results_list)),
  sprintf(" / Total Time: %.2f mins", run_time),
  sprintf(" / Output Directory: %s", out_dir),
  "==================================================",
  " / Directory Structure:",
  "  - results/: SNP / Contains individual SNP extraction results for each phenotype.",
  "  - logs/:  / Execution log files.",
  "  - temp/:  / Temporary files (cleared after analysis).",
  "=================================================="
)

writeLines(stats_text, stats_file)

# Delete temporary folder
log_msg(" / Cleaning up temporary folder.")
unlink(tmp_dir, recursive = TRUE)

# , 
# Save a copy of this script to the output directory for traceability
file.copy("extract_TARGET_SNP.R", file.path(out_dir, "extract_TARGET_SNP.R"))

log_msg(sprintf(" / Execution Result:  %d  / Successfully extracted %d files.", length(all_results_list), length(all_results_list)))
log_msg(sprintf("！: %.2f  / Analysis completed! Total time: %.2f mins.", run_time, run_time))
