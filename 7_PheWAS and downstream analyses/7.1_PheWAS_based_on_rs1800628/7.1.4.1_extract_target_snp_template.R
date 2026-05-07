#!/usr/bin/env Rscript
# ==============================================================================
# [Script]: 7.1.4.1_extract_target_snp_template.R
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

snp_id <- "TARGET_SNP"

gwas_list_file <- "./"


base_output_dir <- OUTPUT_DIR

# （NULL）
resume_folder <- NULL
# resume_folder <- "/Users/.../TARGET_SNP_PheWAS_20260323_051211"

n_cores <- parallel::detectCores() - 1

# APIGWAS
batch_size <- 10


############################################################
# 2 
############################################################

if(is.null(resume_folder)){
  
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  output_dir <- file.path(
    base_output_dir,
    paste0("TARGET_SNP_PheWAS_",timestamp)
  )
  
  dir.create(output_dir, recursive=TRUE)
  
}else{
  
  output_dir <- resume_folder
  
}

print(paste("Output:",output_dir))


############################################################
# 3 GWAS
############################################################

gwas_table <- fread(gwas_list_file)

gwas_ids <- gwas_table$`OpenGWAS ID`

cat("Total GWAS:",length(gwas_ids),"\n")


############################################################
# 4 
############################################################

if(!is.null(resume_folder)){
  
  existing_files <- list.files(
    output_dir,
    pattern="\\.csv$"
  )
  
  existing_ids <- str_extract(
    existing_files,
    "(?<=TARGET_SNP_).*(?=\\.csv)"
  )
  
  gwas_ids <- setdiff(gwas_ids, existing_ids)
  
  cat("Remaining:",length(gwas_ids),"\n")
  
}


############################################################
# 5 GWAS
############################################################

split_batches <- split(
  gwas_ids,
  ceiling(seq_along(gwas_ids)/batch_size)
)

cat("Total batches:",length(split_batches),"\n")


############################################################
# 6 API
############################################################
fetch_batch <- function(gwas_batch){
  
  max_retry <- 5
  retry <- 1
  
  while(retry <= max_retry){
    
    result <- try({
      
      dat <- associations(
        variants = snp_id,
        id = gwas_batch,
        proxies = 0
      )
      
      if(nrow(dat)==0) return(NULL)
      
      # TARGET_SNP
      dat <- dat[dat$rsid == snp_id, ]
      
      if(nrow(dat)==0) return(NULL)
      
      split_dat <- split(dat, dat$id)
      
      for(i in names(split_dat)){
        
        fwrite(
          split_dat[[i]],
          file=file.path(
            output_dir,
            paste0(snp_id,"_",i,".csv")
          )
        )
        
      }
      
      return(TRUE)
      
    }, silent=TRUE)
    
    if(!inherits(result,"try-error")){
      return(TRUE)
    }
    
    Sys.sleep(2^retry)
    retry <- retry + 1
  }
 
  
  # 5
  
  log_file <- file.path(output_dir,"failed_batches.txt")
  
  write(
    paste(Sys.time(), paste(gwas_batch,collapse=",")),
    file=log_file,
    append=TRUE
  )
  
}







############################################################
# 7 
############################################################

plan(multisession, workers = n_cores)

future_lapply(
  
  split_batches,
  fetch_batch
  
)


############################################################
# 8 
############################################################

cat("Finished\n")
