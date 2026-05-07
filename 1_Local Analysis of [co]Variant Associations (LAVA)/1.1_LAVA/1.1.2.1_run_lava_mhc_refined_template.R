#!/usr/bin/env Rscript
# ==============================================================================
# [Script]: 1.1.2.1_run_lava_mhc_refined_template.R
# [Method]: LAVA (Local Analysis of [co]Variant Associations)
# [Step]: LAVA Analysis (LAVA)
# 
# [Function]:
# A generalized template to perform LAVA for evaluating local genetic 
#       correlations across multiple phenotypes. Designed to be coordinate-agnostic
#       by strictly utilizing rsID-driven SNP sets from a reference panel.
#       rsID, GWAS. 
# 
# [Data Availability / ]:
# LAVA UKB reference panels can be downloaded from the official LAVA GitHub.
# ==============================================================================

# 1. Environment & Dependencies / 
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(LAVA)
})

# 2. Command Line Arguments / 
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(c("--input_info"), type="character", default=NULL,
              help="Path to input.info.txt generated from preprocessing step", metavar="character"),
  make_option(c("--ref_prefix"), type="character", default=NULL,
              help="Prefix of the LAVA reference panel (e.g., lava-ukb-v1.1)", metavar="character"),
  make_option(c("--ld_blocks"), type="character", default=NULL,
              help="Path to the LD blocks definition file (e.g., LAVA_s2500_m25_f1_w200.blocks)", metavar="character"),
  make_option(c("--out_dir"), type="character", default="./lava_output",
              help="Output directory path [default: %default]", metavar="character"),
  make_option(c("--chr"), type="integer", default=6,
              help="Chromosome to analyze (e.g., 6 for MHC) [default: %default]", metavar="integer"),
  make_option(c("--start_bp"), type="integer", default=25000000,
              help="Start coordinate in bp (default 25Mb for MHC) [default: %default]", metavar="integer"),
  make_option(c("--end_bp"), type="integer", default=34000000,
              help="End coordinate in bp (default 34Mb for MHC) [default: %default]", metavar="integer"),
  make_option(c("--maf_filter"), type="numeric", default=0.05,
              help="Minimum MAF threshold for SNPs in the reference [default: %default]", metavar="numeric")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_info) || is.null(opt$ref_prefix) || is.null(opt$ld_blocks)) {
  print_help(opt_parser)
  stop("Missing required arguments. Please provide --input_info, --ref_prefix, and --ld_blocks.")
}

# 3. Initialization / 
# ------------------------------------------------------------------------------
out_dir <- opt$out_dir
if (!dir.exists(out_dir)) {
  dir.create(out_dir, recursive = TRUE)
}

# Session Info Logging
sink(file.path(out_dir, "session_info.txt"))
sessionInfo()
sink()

cat("====================================================\n")
cat("Starting LAVA Analysis\n")
cat("Working Directory:", out_dir, "\n")
cat("Chromosome:", opt$chr, "Region:", opt$start_bp, "-", opt$end_bp, "\n")
cat("====================================================\n")

# 4. Prepare LD Blocks & Reference SNPs / LDSNP
# ------------------------------------------------------------------------------
cat("[1] Reading LD blocks...\n")
ld_blocks <- fread(opt$ld_blocks)

# Filter blocks for the target region
target_blocks <- ld_blocks[
  (chr == as.character(opt$chr) | chr == paste0("chr", opt$chr) | chr == opt$chr) & 
  (start < opt$end_bp) & 
  (stop > opt$start_bp)
]

cat(sprintf("Found %d LD blocks overlapping the target region.\n", nrow(target_blocks)))

# Load reference info
chr_info_file <- paste0(opt$ref_prefix, "_chr", opt$chr, ".info")
if (!file.exists(chr_info_file)) {
  stop(paste("Reference info file not found:", chr_info_file))
}
info_dt <- fread(chr_info_file)
info_dt <- info_dt[CHR == opt$chr]

cat(sprintf("Total SNPs in reference info (Chr%d): %d\n", opt$chr, nrow(info_dt)))

# Filtering: MAF and Palindromic SNPs
info_dt[, MAF := pmin(FREQ, 1 - FREQ)]
n_low_maf <- nrow(info_dt[MAF < opt$maf_filter])
info_dt <- info_dt[MAF >= opt$maf_filter]
cat(sprintf("Removed %d SNPs with MAF < %s. Remaining: %d\n", n_low_maf, opt$maf_filter, nrow(info_dt)))

palindromic <- (info_dt$A1 == "A" & info_dt$A2 == "T") |
               (info_dt$A1 == "T" & info_dt$A2 == "A") |
               (info_dt$A1 == "C" & info_dt$A2 == "G") |
               (info_dt$A1 == "G" & info_dt$A2 == "C")
n_palin <- sum(palindromic)
info_dt <- info_dt[!palindromic, ]
cat(sprintf("Removed %d palindromic SNPs. Remaining: %d\n", n_palin, nrow(info_dt)))

# Construct strict rsID-driven loci definition
loci_list <- list()
for (i in 1:nrow(target_blocks)) {
  blk <- target_blocks[i, ]
  blk_snps <- info_dt[POS >= blk$start & POS <= blk$stop, SNP]
  
  if (length(blk_snps) > 0) {
    loc_id <- paste0("Target_Block_", i)
    loci_list[[i]] <- data.frame(
      LOC = loc_id,
      CHR = opt$chr,
      START = blk$start,
      STOP = blk$stop,
      SNPS = paste(blk_snps, collapse = ";"),
      stringsAsFactors = FALSE
    )
  }
}

if (length(loci_list) == 0) {
  stop("No valid loci found in the target region with overlapping SNPs after filtering.")
}

loci_df <- do.call(rbind, loci_list)
loci_out <- file.path(out_dir, "filtered_target.loci")
write.table(loci_df, loci_out, sep="\t", quote=FALSE, row.names=FALSE)
cat(sprintf("Created %s with %d blocks.\n", loci_out, nrow(loci_df)))

# 5. Core LAVA Analysis / LAVA
# ------------------------------------------------------------------------------
cat("\n[2] Initializing LAVA...\n")
input_info_dt <- fread(opt$input_info)
pheno_names <- input_info_dt$phenotype

lava_input <- process.input(input.info.file = opt$input_info,
                            sample.overlap.file = NULL,
                            ref.prefix = opt$ref_prefix,
                            phenos = pheno_names)

loci_data <- read.loci(loci_out)
results_list <- list()

cat(sprintf("Analyzing %d loci...\n", nrow(loci_data)))

for (i in 1:nrow(loci_data)) {
  loc_id <- loci_data$LOC[i]
  cat(sprintf(" -> Processing %s (%d/%d)...\n", loc_id, i, nrow(loci_data)))
  
  tryCatch({
    locus <- process.locus(loci_data[i,], lava_input)
    
    if (!is.null(locus)) {
      n_snps_used <- if(!is.null(locus$n.snps)) locus$n.snps else length(locus$snps)
      
      # run.bivar computes pairwise correlations for all input phenos if target is NULL
      res <- run.bivar(locus)
      
      if (!is.null(res)) {
        res$LOC <- loc_id
        res$CHR <- loci_data$CHR[i]
        res$START <- loci_data$START[i]
        res$STOP <- loci_data$STOP[i]
        res$n_snps <- n_snps_used
        results_list[[loc_id]] <- res
      }
    }
  }, error = function(e) {
    cat(sprintf("    [Error] %s: %s\n", loc_id, e$message))
  })
}

# 6. Save Results / 
# ------------------------------------------------------------------------------
cat("\n[3] Finalizing results...\n")
if (length(results_list) > 0) {
  final_res <- do.call(rbind, results_list)
  
  other_cols <- setdiff(names(final_res), c("LOC", "CHR", "START", "STOP", "n_snps", "phen1", "phen2"))
  cols <- c("LOC", "CHR", "START", "STOP", "n_snps", "phen1", "phen2", other_cols)
  final_res <- final_res[, cols]
  
  out_txt <- file.path(out_dir, "lava_bivar_results.txt")
  out_rds <- file.path(out_dir, "lava_bivar_results.rds")
  
  saveRDS(final_res, out_rds)
  write.table(final_res, out_txt, sep="\t", quote=FALSE, row.names=FALSE)
  
  cat(sprintf("Analysis complete. Results successfully saved to:\n  - %s\n  - %s\n", out_txt, out_rds))
} else {
  cat("Warning: No results generated. Check your inputs or regional constraints.\n")
}
