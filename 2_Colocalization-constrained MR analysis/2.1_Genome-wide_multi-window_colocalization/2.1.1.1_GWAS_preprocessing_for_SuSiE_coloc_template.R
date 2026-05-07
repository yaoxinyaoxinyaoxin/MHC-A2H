#!/usr/bin/env Rscript
# ==============================================================================
# [Script]: 2.1.1.1_GWAS_preprocessing_for_SuSiE_coloc_template.R
# [Method]: SuSiE-based Colocalization 
# [Step]: GWAS Data Preprocessing (GWAS)
# 
# [Function]:
# Prepares GWAS summary statistics (e.g., HZ, Aging) for susieR/coloc input. 
#       Maps columns to standard names, strictly aligns alleles to a specified LD 
#       reference panel, calculates MAF, and filters out SNPs with MAF < 0.05.
#       LD, MAFMAF < 0.05SNP. 
# 
# [Data Availability / ]:
# Ensure the source GWAS and UKB LD reference panel are available.
# ==============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)
Sys.setenv(TZ = "UTC")

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(readr)
})

option_list <- list(
  make_option(c("--gwas_file"), type="character", default=NULL,
              help="Path to GWAS source data (.tsv.gz or .txt)", metavar="character"),
  make_option(c("--ld_snp_info"), type="character", default=NULL,
              help="Path to UKB LD reference panel SNP info", metavar="character"),
  make_option(c("--out_dir"), type="character", default="./gwas_preprocessed",
              help="Output directory path [default: %default]", metavar="character"),
  make_option(c("--phenotype"), type="character", default="Phenotype",
              help="Name of the phenotype (e.g., HZ, Aging) [default: %default]", metavar="character"),
  make_option(c("--force_N"), type="integer", default=NULL,
              help="Inject this sample size (N) if the column is missing in GWAS data", metavar="integer")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$gwas_file) || is.null(opt$ld_snp_info)) {
  print_help(opt_parser)
  stop("Missing required arguments: --gwas_file or --ld_snp_info.")
}

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_root <- file.path(opt$out_dir, paste0(opt$phenotype, "_GWAS_preprocessing_", timestamp))

dir.create(file.path(out_root, "results"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_root, "stats"), showWarnings = FALSE)
dir.create(file.path(out_root, "logs"), showWarnings = FALSE)

log_file <- file.path(out_root, "logs", paste0("run_log_", timestamp, ".txt"))
sink(log_file, append = TRUE, split = TRUE)

cat("========================================\n")
cat("GWAS Preprocessing for SuSiE/coloc\n")
cat("========================================\n")
cat("Phenotype:", opt$phenotype, "\n")
cat("GWAS File:", opt$gwas_file, "\n")
cat("LD Panel:", opt$ld_snp_info, "\n")

safe_upper <- function(x) { toupper(gsub("\n|\r", "", as.character(x))) }

detect_gwas_cols <- function(dt) {
  cols <- colnames(dt)
  map <- list()
  map$rsid <- intersect(cols, c("rsid","RSID","MarkerName","variant_id","ID","SNP"))[1]
  map$chr  <- intersect(cols, c("chromosome","CHR","Chromosome","chr"))[1]
  map$pos  <- intersect(cols, c("base_pair_location","BP","position","Position"))[1]
  map$ea   <- intersect(cols, c("effect_allele","EA","A1","Allele1","ALT"))[1]
  map$oa   <- intersect(cols, c("other_allele","NEA","A2","Allele2","REF"))[1]
  map$beta <- intersect(cols, c("beta","BETA","effect_size","ES"))[1]
  map$se   <- intersect(cols, c("standard_error","SE","se"))[1]
  map$pval <- intersect(cols, c("p_value","pval","p","P","Pvalue"))[1]
  map$eaf  <- intersect(cols, c("effect_allele_frequency","eaf","EAF","AF","MAF"))[1]
  map$n    <- intersect(cols, c("n","N","samplesize","sample_size"))[1]
  map
}

read_ld_panel <- function(ld_path) {
  ld <- fread(ld_path, sep = "\t")
  setnames(ld, old = c("SNP_ID","Chromosome","Position_BP","Allele1","Allele2","Allele_Direction"), 
               new = c("rsid","chr","pos","ld_a1","ld_a2","ld_dir"), skip_absent=TRUE)
  ld$rsid <- as.character(ld$rsid)
  ld$ld_a1 <- safe_upper(ld$ld_a1)
  ld$ld_a2 <- safe_upper(ld$ld_a2)
  ld
}

cat("Reading LD reference panel...\n")
ld_dt <- read_ld_panel(opt$ld_snp_info)

cat("Reading GWAS source data...\n")
gwas_dt <- tryCatch({ fread(opt$gwas_file) }, error = function(e) NULL)
if (is.null(gwas_dt) || nrow(gwas_dt) == 0) stop("Failed to read GWAS data.")

source_total <- nrow(gwas_dt)
map <- detect_gwas_cols(gwas_dt)

if (is.na(map$n) && !is.null(opt$force_N)) {
  cat("Missing Sample Size (N) column. Injecting N =", opt$force_N, "\n")
  gwas_dt$N <- opt$force_N
  map$n <- "N"
}

required <- c("rsid","chr","pos","ea","oa","beta","se","pval")
if (!all(!is.na(map[required]))) stop("Missing required columns in GWAS data.")

# Align Alleles
cat("Aligning alleles to LD panel...\n")
gwas_dt[[map$ea]] <- safe_upper(gwas_dt[[map$ea]])
gwas_dt[[map$oa]] <- safe_upper(gwas_dt[[map$oa]])
gwas_dt <- gwas_dt[gwas_dt[[map$ea]] %in% c("A","C","G","T") & gwas_dt[[map$oa]] %in% c("A","C","G","T"), ]

merged <- gwas_dt %>% inner_join(ld_dt, by = c(setNames("rsid", map$rsid)))

status <- ifelse(is.na(merged$ld_a1) | is.na(merged$ld_a2), "missing_ld",
           ifelse(is.na(merged[[map$ea]]) | is.na(merged[[map$oa]]), "missing_data",
             ifelse(merged[[map$ea]] == merged$ld_a1 & merged[[map$oa]] == merged$ld_a2, "consistent",
               ifelse(merged[[map$ea]] == merged$ld_a2 & merged[[map$oa]] == merged$ld_a1, "flipped", "inconsistent"))))

beta <- as.numeric(merged[[map$beta]])
eaf  <- suppressWarnings(as.numeric(merged[[map$eaf]]))

out <- data.table(
  snp = merged[[map$rsid]],
  chr = suppressWarnings(as.integer(merged[[map$chr]])),
  pos = suppressWarnings(as.integer(merged[[map$pos]])),
  A1  = merged$ld_a1,
  A2  = merged$ld_a2,
  beta = ifelse(status == "flipped", -beta, beta),
  se   = suppressWarnings(as.numeric(merged[[map$se]])),
  pval = suppressWarnings(as.numeric(merged[[map$pval]])),
  eaf  = ifelse(status == "flipped" & !is.na(eaf), 1 - eaf, eaf),
  N    = suppressWarnings(as.numeric(merged[[map$n]]))
)

out <- out[status %in% c("consistent","flipped")]
out <- out[!is.na(snp) & !is.na(chr) & !is.na(pos) & !is.na(beta) & !is.na(se) & !is.na(pval)]

# MAF Filtering
if (!is.null(out$eaf) && !all(is.na(out$eaf))) {
  out$maf <- ifelse(out$eaf > 0.5, 1 - out$eaf, out$eaf)
  out <- out[out$maf >= 0.05, ]
  out$maf <- NULL
}

out_file <- file.path(out_root, "results", paste0(opt$phenotype, "_GWAS_preprocessed.csv"))
fwrite(out, out_file)
cat("Saved preprocessed results to:", out_file, "\n")

cat("Process completed successfully.\n")
sink()
