#!/usr/bin/env Rscript
# ==============================================================================
# [Script]: 2.1.1.2_eQTL_preprocessing_for_SuSiE_coloc_template.R
# [Method]: SuSiE-based Colocalization 
# [Step]: Bulk eQTL Data Preprocessing (Bulk eQTL)
# 
# [Function]:
# Prepares bulk eQTL summary statistics as inputs for susieR/coloc. 
#       Filters SNPs by UKB LD reference panel, strictly aligns allele directions, 
#       calculates MAF, and filters SNPs with MAF < 0.05. Uses parallel processing 
#       for multiple gene subsets.
#        UKB LD  SNP ,  MAF < 0.05 
#        SNP. . 
# 
# [Data Availability / ]:
# Requires gene_info.csv and raw eQTL data files.
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
              help="Path to CSV containing ENS ID, gene_name, Sample size", metavar="character"),
  make_option(c("--eqtl_src_dir"), type="character", default=NULL,
              help="Root directory containing raw eQTL CSVs", metavar="character"),
  make_option(c("--ld_snp_info"), type="character", default=NULL,
              help="Path to UKB LD reference panel SNP info", metavar="character"),
  make_option(c("--out_dir"), type="character", default="./eqtl_preprocessed",
              help="Output directory path", metavar="character"),
  make_option(c("--cores"), type="integer", default=1,
              help="Number of cores for parallel processing [default: %default]", metavar="integer")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (any(sapply(list(opt$gene_info_csv, opt$eqtl_src_dir, opt$ld_snp_info), is.null))) {
  print_help(opt_parser)
  stop("Missing required arguments.")
}

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_root <- file.path(opt$out_dir, paste0("eQTL_preprocessing_", timestamp))
dir.create(file.path(out_root, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "logs"), showWarnings = FALSE)

log_file <- file.path(out_root, "logs", "run_log.txt")
sink(log_file, append = TRUE, split = TRUE)
cat("Starting bulk eQTL Preprocessing...\n")

ld_dt <- fread(opt$ld_snp_info, sep = "\t")
setnames(ld_dt, old = c("SNP_ID","Allele1","Allele2"), new = c("rsid","ld_a1","ld_a2"), skip_absent=TRUE)

gene_info <- fread(opt$gene_info_csv)
# Detect ENS ID and N columns dynamically
ens_col <- intersect(colnames(gene_info), c("ens_gene_id","ENSG","ensg","ensembl_gene_id"))[1]
n_col <- intersect(colnames(gene_info), c("Sample size","Sample_size","N","n"))[1]

to_process <- gene_info

plan(multisession, workers = min(opt$cores, parallel::detectCores()))

process_gene <- function(ens_id, sample_size) {
  out_file <- file.path(out_root, "results", paste0(ens_id, ".csv"))
  if (file.exists(out_file)) return(TRUE)
  
  files <- list.files(opt$eqtl_src_dir, pattern = paste0("^", ens_id, ".*\\.csv$|^eqtl-a-", ens_id, ".*\\.csv$"), recursive = TRUE, full.names = TRUE)
  if (length(files) == 0) return(FALSE)
  
  dt <- fread(files[1])
  if (nrow(dt) == 0) return(FALSE)
  
  # Standardize cols
  cols <- colnames(dt)
  rsid <- intersect(cols, c("SNP","rsid","variant"))[1]
  chr <- intersect(cols, c("chr","chromosome"))[1]
  pos <- intersect(cols, c("pos","position"))[1]
  ea <- intersect(cols, c("effect_allele","A1","ALT"))[1]
  oa <- intersect(cols, c("other_allele","A2","REF"))[1]
  beta <- intersect(cols, c("beta","effect_size"))[1]
  se <- intersect(cols, c("se","standard_error"))[1]
  pval <- intersect(cols, c("pval","p_value","P"))[1]
  eaf <- intersect(cols, c("eaf","effect_allele_frequency"))[1]
  
  if (any(is.na(c(rsid, chr, pos, ea, oa, beta, se, pval)))) return(FALSE)
  
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
    N = sample_size
  )
  
  out <- out[status %in% c("consistent", "flipped") & !is.na(snp)]
  if ("eaf" %in% colnames(out) && !all(is.na(out$eaf))) {
    out$maf <- ifelse(out$eaf > 0.5, 1 - out$eaf, out$eaf)
    out <- out[out$maf >= 0.05, ]
    out$maf <- NULL
  }
  
  if (nrow(out) > 0) fwrite(out, out_file)
  return(TRUE)
}

future_lapply(seq_len(nrow(to_process)), function(i) {
  process_gene(to_process[[ens_col]][i], to_process[[n_col]][i])
}, future.seed = TRUE)

cat("eQTL processing completed.\n")
sink()
