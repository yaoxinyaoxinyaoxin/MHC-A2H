#!/usr/bin/env Rscript
# ==============================================================================
# [Script]: 6.1.1.2_extract_snps_by_gene_region_template.R
# [Method]: cis-MR IV Selection Pipeline Step
# [Step]: 6.1.1.2_extract_snps_by_gene_region.R
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

log_message(paste("Starting", "6.1.1.2_extract_snps_by_gene_region.R"))

log_message(". (Script finished.)")
