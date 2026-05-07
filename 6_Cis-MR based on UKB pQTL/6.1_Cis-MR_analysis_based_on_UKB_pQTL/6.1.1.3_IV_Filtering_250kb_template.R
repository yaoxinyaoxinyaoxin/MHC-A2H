#!/usr/bin/env Rscript
# ==============================================================================
# [Script]: 6.1.1.3_IV_Filtering_250kb_template.R
# [Method]: cis-MR IV Selection Pipeline Step
# [Step]: 6.1.1.3_IV_Filtering_250kb.R
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

log_message(paste("Starting", "6.1.1.3_IV_Filtering_250kb.R"))
\n")

num_cores <- detectCores() - 2
if(num_cores < 1) num_cores <- 1
cat(sprintf(",  %d ...\n", num_cores))

process_file <- function(file_path) {
  tryCatch({
    # 1.  (Read data)
    dat <- fread(file_path, data.table = FALSE)
    
    if(nrow(dat) == 0) return(NULL)
    
    #  (Column mapping)
    col_mapping <- c(
      "rs_id" = "SNP",
      "beta" = "beta.exposure",
      "standard_error" = "se.exposure",
      "effect_allele_frequency" = "eaf.exposure",
      "p_value" = "pval.exposure",
      "n" = "samplesize.exposure",
      "chromosome" = "chr.exposure",
      "base_pair_location" = "pos.exposure",
      "effect_allele" = "effect_allele.exposure",
      "other_allele" = "other_allele.exposure"
    )
    
    for (old_col in names(col_mapping)) {
      if (old_col %in% colnames(dat)) {
        colnames(dat)[colnames(dat) == old_col] <- col_mapping[old_col]
      }
    }
    
    #  rsid  (Remove rows with missing rsid)
    dat <- dat[!is.na(dat$SNP) & dat$SNP != "", ]
    if(nrow(dat) == 0) return(NULL)
    
    # 2. : P < 5e-8, MAF > 0.01 (Filter: P < 5e-8, MAF > 0.01)
    dat$MAF <- pmin(dat$eaf.exposure, 1 - dat$eaf.exposure)
    dat <- dat[dat$pval.exposure < p_threshold & dat$MAF > maf_threshold, ]
    
    if(nrow(dat) == 0) return(NULL)
    
    # 3.  F  F > 10 (Calculate F-statistic and filter F > 10)
    beta <- dat$beta.exposure
    se <- dat$se.exposure
    eaf <- dat$eaf.exposure
    N <- dat$samplesize.exposure
    k <- 1
    
    # R² (Calculate R² using improved formula)
    numerator <- 2 * beta^2 * eaf * (1 - eaf)
    denominator <- numerator + 2 * se^2 * N * eaf * (1 - eaf)
    R2 <- numerator / denominator
    
    # F: F = [R² × (N - k - 1)] / [(1 - R²) × k] (Calculate F-statistic using standard formula)
    F_stat <- (R2 * (N - k - 1)) / ((1 - R2) * k)
    
    dat$R2 <- R2
    dat$F_stat <- F_stat
    
    dat <- dat[dat$F_stat > f_threshold, ]
    
    if(nrow(dat) == 0) return(NULL)
    
    # 4.  (Save results)
    # LD, 
    out_file <- file.path(csv_dir, basename(file_path))
    fwrite(dat, out_file, row.names = FALSE)
    
    #  (Clean memory)
    rm(dat)
    return("SUCCESS")
  }, error = function(e) {
    cat(sprintf("Error processing file %s: %s\n", basename(file_path), e$message), 
        file = file.path(log_dir, "error.log"), append = TRUE)
    return(NULL)
  })
}

#  (Execute parallel processing)
# (chunking),  (Use chunking to optimize memory, regular cleanup)
chunk_size <- 500 # 500, 
file_chunks <- split(files, ceiling(seq_along(files) / chunk_size))

processed_files <- c()

for(i in seq_along(file_chunks)) {
  cat(sprintf("Processing chunk %d of %d...\n", i, length(file_chunks)))
  
  # mclapply (Use mclapply for parallel processing)
  res <- mclapply(file_chunks[[i]], process_file, mc.cores = num_cores)
  processed_files <- c(processed_files, unlist(res))
  
  #  (Regular memory cleanup)
  rm(res)
  gc()
}

num_success <- sum(processed_files == "SUCCESS", na.rm = TRUE)

cat(". \n")
cat(" - IV: ", num_success, "\n")

# -------------------------------------------------------------------------
#  README   (Generate README and Statistics)
# -------------------------------------------------------------------------
readme_filename <- sprintf("%s_%s_IV_Filtering_Readme.txt", project_name, timestamp)
stats_filename <- sprintf("%s_%s_IV_Filtering_Stats.csv", project_name, timestamp)

# README 
readme_content <- paste0(
  "============================================================\n",
  " (IV Filtering Analysis Report - No LD Clumping)\n",
  "============================================================\n\n",
  "1.  (Execution Time): ", Sys.time, "\n",
  "2.  (Script Function): SNP, (P、MAF、F), LD. \n",
  "3.  (Input Directory): ", source_dir, "\n",
  "4.  (Output Directory): ", out_dir, "\n\n",
  "5.  (Filtering Criteria):\n",
  "   - P-value < ", p_threshold, "\n",
  "   - MAF > ", maf_threshold, "\n",
  "   - F-statistic > ", f_threshold, "\n\n",
  "6.  (Result Statistics):\n",
  "   -  (Total Files Processed): ", length(files), "\n",
  "   - IV (Files with Valid IVs): ", num_success, "\n",
  "   - : LD, P、MAFFSNPs. \n\n",
  "============================================================\n"
)

#  README 
writeLines(readme_content, file.path(readme_dir, readme_filename))

stats_df <- data.frame(
  Metric = c("Total_Files_Processed", "Files_With_Valid_IVs", "P_Value_Threshold", "MAF_Threshold", "F_Statistic_Threshold"),
  Value = c(length(files), num_success, p_threshold, maf_threshold, f_threshold)
)
fwrite(stats_df, file.path(readme_dir, stats_filename))

cat(", : ", out_dir, "\n")
cat("Readme: ", readme_dir, "\n")
