#!/usr/bin/env Rscript
# ==============================================================================
# [Script]: 2.1.1.4_sc_eQTL_preprocessing_for_SuSiE_template.R
# [Method]: SuSiE-based Colocalization 
# [Step]: sc-eQTL Data Preprocessing (eQTL)
# 
# [Function]:
# Prepares sc-eQTL summary statistics as inputs for susieR/coloc. 
#       Filters SNPs by UKB LD reference panel, aligns allele directions, 
#       calculates MAF, and filters SNPs with MAF >= 0.05. 
#       Preserves cell-type specific directory structure and checks for duplicate SNPs.
#        UKB LD  SNP ,  MAF < 0.05 
#        SNP. SNP. 
# 
# [Data Availability / ]:
# Requires gene_info.csv (with 'cell' and 'gene_id' columns) and raw sc-eQTL data files.
# ==============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)
Sys.setenv(TZ = "UTC")

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(future.apply)
  library(parallel)
})

option_list <- list(
  make_option(c("--gene_info_csv"), type="character", default=NULL,
              help="Path to CSV containing 'cell' and 'gene_id' columns", metavar="character"),
  make_option(c("--eqtl_src_dir"), type="character", default=NULL,
              help="Root directory containing cell-type subdirectories with raw sc-eQTL CSVs", metavar="character"),
  make_option(c("--ld_snp_info"), type="character", default=NULL,
              help="Path to UKB LD reference panel SNP info", metavar="character"),
  make_option(c("--out_dir"), type="character", default="./sc_eqtl_preprocessed",
              help="Output directory path", metavar="character"),
  make_option(c("--cores"), type="integer", default=1,
              help="Number of cores for parallel processing [default: %default]", metavar="integer"),
  make_option(c("--test_mode"), type="logical", default=FALSE,
              help="Enable test mode to run only the first 20 genes [default: %default]", metavar="logical")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (any(sapply(list(opt$gene_info_csv, opt$eqtl_src_dir, opt$ld_snp_info), is.null))) {
  print_help(opt_parser)
  stop("Missing required arguments: --gene_info_csv, --eqtl_src_dir, --ld_snp_info")
}

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_root <- file.path(opt$out_dir, paste0("sc_eQTL_preprocessing_", timestamp))

# Subdirectories
dir.create(file.path(out_root, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "stats"), showWarnings = FALSE)
dir.create(file.path(out_root, "logs"), showWarnings = FALSE)

# Initialize logging
log_file <- file.path(out_root, "logs", "run_log.txt")
sink(log_file, append = TRUE, split = TRUE)
cat("========================================\n")
cat("sc-eQTL Preprocessing for SuSiE/coloc\n")
cat("========================================\n")
cat("Start Time:", as.character(Sys.time()), "\n")
cat("Cores:", opt$cores, "\n")
cat("Test Mode:", opt$test_mode, "\n")

# Read LD panel
cat("Reading LD reference panel...\n")
ld_dt <- fread(opt$ld_snp_info, sep = "\t")
setnames(ld_dt, old = c("SNP_ID","Allele1","Allele2"), new = c("rsid","ld_a1","ld_a2"), skip_absent=TRUE)

# Read gene info
cat("Reading gene info CSV...\n")
gene_info <- fread(opt$gene_info_csv)
if (!all(c("cell", "gene_id") %in% colnames(gene_info))) {
  stop("Gene info file must contain 'cell' and 'gene_id' columns.")
}

to_process <- gene_info %>% dplyr::select(cell, gene_id) %>% dplyr::distinct()

if (opt$test_mode) {
  to_process <- head(to_process, 20)
  cat("Test mode enabled: only processing first 20 genes.\n")
}

cat("Total tasks:", nrow(to_process), "\n")

plan(multisession, workers = min(opt$cores, parallel::detectCores()))

process_gene <- function(ens_id, cell_type) {
  cell_dir <- file.path(out_root, "results", cell_type)
  if (!dir.exists(cell_dir)) dir.create(cell_dir, recursive = TRUE, showWarnings = FALSE)
  
  src_file <- file.path(opt$eqtl_src_dir, cell_type, paste0(ens_id, ".csv"))
  out_file <- file.path(cell_dir, paste0(ens_id, ".csv"))
  
  if (file.exists(out_file) || !file.exists(src_file)) return(FALSE)
  
  dt <- tryCatch({ fread(src_file) }, error = function(e) NULL)
  if (is.null(dt) || nrow(dt) == 0) return(FALSE)
  
  # Standardize cols
  cols <- colnames(dt)
  rsid <- intersect(cols, c("SNP","rsid","variant"))[1]
  chr <- intersect(cols, c("chr","chromosome","CHR"))[1]
  pos <- intersect(cols, c("pos","position","POS","BP"))[1]
  ea <- intersect(cols, c("effect_allele","A1","ALT","EA"))[1]
  oa <- intersect(cols, c("other_allele","A2","REF","NEA"))[1]
  beta <- intersect(cols, c("beta","effect_size","BETA"))[1]
  se <- intersect(cols, c("se","standard_error","SE"))[1]
  pval <- intersect(cols, c("pval","p_value","P","P_VALUE"))[1]
  eaf <- intersect(cols, c("eaf","effect_allele_frequency","EAF","freq"))[1]
  n_col <- intersect(cols, c("N","samplesize","n"))[1]
  
  if (any(is.na(c(rsid, chr, pos, ea, oa, beta, se, pval, n_col)))) return(FALSE)
  
  dt[[ea]] <- toupper(dt[[ea]])
  dt[[oa]] <- toupper(dt[[oa]])
  
  merged <- dt %>% inner_join(ld_dt, by = c(setNames("rsid", rsid)))
  
  status <- ifelse(merged[[ea]] == merged$ld_a1 & merged[[oa]] == merged$ld_a2, "consistent",
              ifelse(merged[[ea]] == merged$ld_a2 & merged[[oa]] == merged$ld_a1, "flipped", "inconsistent"))
              
  beta_val <- as.numeric(merged[[beta]])
  eaf_val <- if (!is.na(eaf)) as.numeric(merged[[eaf]]) else rep(NA, nrow(merged))
  
  out <- data.table(
    snp = merged[[rsid]], chr = merged[[chr]], pos = merged[[pos]],
    A1 = merged$ld_a1, A2 = merged$ld_a2,
    beta = ifelse(status == "flipped", -beta_val, beta_val),
    se = as.numeric(merged[[se]]), pval = as.numeric(merged[[pval]]),
    eaf = ifelse(status == "flipped" & !is.na(eaf_val), 1 - eaf_val, eaf_val),
    N = as.numeric(merged[[n_col]])
  )
  
  out <- out[status %in% c("consistent", "flipped") & !is.na(snp)]
  
  if ("eaf" %in% colnames(out) && !all(is.na(out$eaf))) {
    out$maf <- ifelse(out$eaf > 0.5, 1 - out$eaf, out$eaf)
    out <- out[out$maf >= 0.05, ]
    out$maf <- NULL
  }
  
  # Record duplicates but do not remove
  dup_snps <- out$snp[duplicated(out$snp)]
  if (length(dup_snps) > 0) {
    dup_info <- data.frame(
      cell = cell_type, gene = ens_id, n_duplicates = length(dup_snps),
      stringsAsFactors = FALSE
    )
    dup_file <- file.path(out_root, "stats", paste0("dup_", cell_type, "_", ens_id, ".csv"))
    fwrite(dup_info, dup_file)
  }
  
  if (nrow(out) > 0) fwrite(out, out_file)
  return(TRUE)
}

future_lapply(seq_len(nrow(to_process)), function(i) {
  process_gene(to_process$gene_id[i], to_process$cell[i])
}, future.seed = TRUE)

cat("sc-eQTL processing completed.\n")
cat("End Time:", as.character(Sys.time()), "\n")
sink()
