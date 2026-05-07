#!/usr/bin/env Rscript
# ==============================================================================
# [Script]: 7.6.1.5_Fallback_LD_Clumping_template.R
# [Method]: cis-MR IV Selection Pipeline Step
# [Step]: 7.6.1.5_Fallback_LD_Clumping.R
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

log_message(paste("Starting", "7.6.1.5_Fallback_LD_Clumping.R"))
\n", start_time))

# 4.  (Define paths and parameters)
failed_list_path <- "./"
source_dir <- INPUT_PATH
target_dir <- "./"
ld_ref_panel <- "./"
plink_bin <- PLINK_BIN
r2_threshold <- 0.3

# 5.  (Read failed files list)
if (!file.exists(failed_list_path)) {
  stop(sprintf(" (Cannot find failed record file): %s", failed_list_path))
}

failed_files <- readLines(failed_list_path)
failed_files <- failed_files[failed_files != ""] #  (Remove empty lines)

cat(sprintf(" %d  (Found %d files for fallback processing).\n", length(failed_files), length(failed_files)))

# 6.  (Iterate and process each file)
for (i in seq_along(failed_files)) {
  file_name <- failed_files[i]
  source_file <- file.path(source_dir, file_name)
  target_file <- file.path(target_dir, file_name)
  
  cat(sprintf("\n[%d/%d]  (Processing file): %s\n", i, length(failed_files), file_name))
  
  if (!file.exists(source_file)) {
    cat(sprintf("  -> :  (Warning: Source file does not exist): %s\n", source_file))
    next
  }
  
  #  (Read source data)
  IV <- as.data.frame(fread(source_file))
  orig_nrow <- nrow(IV)
  
  if (orig_nrow <= 1) {
    cat(sprintf("  -> SNP<=1,  (SNP count <= 1, no processing needed)\n"))
    next
  }
  
  # P (Sort by P-value ascending)
  IV <- IV[order(IV$pval.exposure), ]
  
  # SNP (Extract SNP list and save to temp file)
  tmp_snps <- tempfile(fileext = ".txt")
  writeLines(IV$SNP, tmp_snps)
  tmp_out <- tempfile()
  
  # PLINKR^2 (Compute R^2 using PLINK)
  # : 
  # --r2 R
  # --ld-window-kb 10000 10MB
  # --ld-window 99999 SNP
  # --ld-window-r2 0.3 R^2 > 0.3
  plink_cmd <- sprintf("'%s' --bfile '%s' --extract '%s' --r2 --ld-window-kb 10000 --ld-window 99999 --ld-window-r2 %f --out '%s'",
                       plink_bin, ld_ref_panel, tmp_snps, r2_threshold, tmp_out)
  
  # PLINK (Run PLINK command and ignore standard output)
  system(plink_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  ld_file <- paste0(tmp_out, ".ld")
  
  if (file.exists(ld_file) && file.info(ld_file)$size > 0) {
    # LD (Read LD results)
    # PLINK .ld : CHR_A, BP_A, SNP_A, CHR_B, BP_B, SNP_B, R2
    ld_data <- fread(ld_file)
    
    # SNP (Initialize lists for SNPs to keep and process)
    snps_to_process <- IV$SNP
    keep_snps <- character(0)
    
    #  (Iterative filtering)
    while (length(snps_to_process) > 0) {
      # PSNP (Get the current SNP with the smallest P-value)
      top_snp <- snps_to_process[1]
      keep_snps <- c(keep_snps, top_snp)
      
      #  top_snp R^2 > 0.3 SNP (Find all SNPs with R^2 > 0.3 with top_snp)
      # PLINK --ld-window-r2 0.3, 
      ld_pairs <- ld_data[(SNP_A == top_snp & SNP_B %in% snps_to_process) | 
                          (SNP_B == top_snp & SNP_A %in% snps_to_process)]
      
      snps_in_ld <- unique(c(ld_pairs$SNP_A, ld_pairs$SNP_B))
      
      #  top_snp LDSNP 
      # (Remove top_snp and SNPs in LD from the processing list)
      snps_to_process <- setdiff(snps_to_process, c(top_snp, snps_in_ld))
    }
    
    # SNP (Filter data based on the kept SNP list)
    IV_filtered <- IV[IV$SNP %in% keep_snps, ]
    final_nrow <- nrow(IV_filtered)
    cat(sprintf("  -> LD,  %d  %d SNP (Fallback LD clumping done, retained %d out of %d SNPs)\n", 
                orig_nrow, final_nrow, final_nrow, orig_nrow))
    
  } else {
    #  .ld ,  SNP , 
    #  SNP  R^2  <= 0.3. 
    # ,  SNP. 
    cat(sprintf("  -> PLINK R^2 > %f LDSNP,  %d SNP (No LD pairs with R^2 > %f found or SNPs not in reference panel, keeping all %d SNPs)\n", 
                r2_threshold, orig_nrow, r2_threshold, orig_nrow))
    IV_filtered <- IV
  }
  
  #  (Overwrite save to target directory)
  fwrite(IV_filtered, target_file, row.names = FALSE)
  
  #  (Clean up temporary files)
  unlink(c(tmp_snps, paste0(tmp_out, "*")))
}

# 7.  (Record end time and total duration)
end_time <- Sys.time()
cat(sprintf("\n (Fallback script end time): %s\n", end_time))
cat(sprintf(" (Total duration): %s\n", round(difftime(end_time, start_time, units = "mins"), 2), " (mins)\n"))
