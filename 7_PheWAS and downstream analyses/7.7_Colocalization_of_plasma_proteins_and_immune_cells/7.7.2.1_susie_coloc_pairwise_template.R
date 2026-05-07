#!/usr/bin/env Rscript
# ==============================================================================
# [Script]: 7.7.2.1_susie_coloc_pairwise_template.R
# [Method]: SuSiE Colocalization Analysis
# [Step]: Run pairwise SuSiE coloc
# 
# [Function]:
# Run SuSiE-based colocalization analysis to identify shared causal variants.
# 
# [Input]:
#   --dataset1_dir : Directory containing formatted dataset 1
#   --dataset2_dir : Directory containing formatted dataset 2
#   --out_dir      : Output directory path
#   --ld_matrix    : Path to local LD matrix or reference panel
# 
# [Output]:
# SuSiE colocalization results (H4 posteriors, shared variants).
# ==============================================================================

rm(list = ls())
gc()

# 1. Load Dependencies
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(coloc)
  library(susieR)
})

# 2. Parse Command Line Arguments
option_list <- list(
  make_option(c("--dataset1_dir"), type="character", default=NULL, help="Directory for Dataset 1"),
  make_option(c("--dataset2_dir"), type="character", default=NULL, help="Directory for Dataset 2"),
  make_option(c("--out_dir"), type="character", default="./Coloc_Results", help="Output directory path"),
  make_option(c("--ld_matrix"), type="character", default=NULL, help="Path to LD reference panel")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$dataset1_dir) || is.null(opt$dataset2_dir)) {
  print_help(opt_parser)
  stop("Missing required input arguments.")
}

DATASET1_DIR <- opt$dataset1_dir
DATASET2_DIR <- opt$dataset2_dir
OUTPUT_DIR <- opt$out_dir
LD_MATRIX <- opt$ld_matrix

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

}))

## -----------------------------------------------------------------------------
## 1. Global config and parameter parsing / 
## -----------------------------------------------------------------------------
args <- commandArgs(trailingOnly = TRUE)
parse_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 0) return(default)
  if (idx + 1 <= length(args)) return(args[idx + 1]) else return(default)
}

WORK_DIR <- parse_arg("--work_dir", "./")
MR_RESULTS <- parse_arg("--mr_results", "./")
PQTL_DIR <- parse_arg("--pqtl_dir", "./")
IMMUNE_DIR <- parse_arg("--immune_dir", "./")
GENE_INFO <- parse_arg("--gene_info", "./")
LD_DIR <- parse_arg("--ld_dir", "./")
PLINK_BIN <- parse_arg("--plink_bin", "./")

TASKS <- as.integer(parse_arg("--tasks", "10"))
THREADS_PER_TASK <- as.integer(parse_arg("--threads_per_task", "2"))
RESUME <- as.logical(parse_arg("--resume", "TRUE"))
TEST_N <- as.integer(parse_arg("--test_n", "0"))
SPECIFIC_TASKS <- parse_arg("--specific_tasks", NULL)

WINDOW_SIZE <- 250000 # +/- 250kb / 250kb

## Verify Input Files / 
if (!file.exists(MR_RESULTS)) stop(paste("MR results file not found:", MR_RESULTS))
if (!dir.exists(PQTL_DIR)) stop(paste("pQTL directory not found:", PQTL_DIR))
if (!dir.exists(IMMUNE_DIR)) stop(paste("Immune directory not found:", IMMUNE_DIR))
if (!file.exists(GENE_INFO)) stop(paste("Gene info file not found:", GENE_INFO))

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
BASE_OUT <- file.path(WORK_DIR, paste0("susie_coloc_pairwise_", timestamp))
DIR_RESULTS <- file.path(BASE_OUT, "results")
DIR_LOGS <- file.path(BASE_OUT, "logs")
DIR_README <- file.path(BASE_OUT, "readme")
DIR_SCRIPTS <- file.path(BASE_OUT, "scripts")

for (d in c(BASE_OUT, DIR_RESULTS, DIR_LOGS, DIR_README, DIR_SCRIPTS)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

## Unified logger / 
LOG_FILE <- file.path(DIR_LOGS, paste0("run_log_", timestamp, ".txt"))
ERR_FILE <- file.path(DIR_LOGS, paste0("error_log_", timestamp, ".txt"))
PROGRESS_FILE <- file.path(DIR_LOGS, paste0("progress_", timestamp, ".csv"))
SUMMARY_FILE <- file.path(BASE_OUT, "summary_coloc_results.csv")

log_msg <- function(..., file = LOG_FILE) {
  msg <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", paste(..., collapse = " "))
  cat(msg, "\n"); write(msg, file = file, append = TRUE)
}
log_err <- function(..., file = ERR_FILE) {
  msg <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", paste(..., collapse = " "))
  cat(msg, "\n"); write(msg, file = file, append = TRUE)
}

log_msg("Start Pairwise Susie Coloc Pipeline. Work dir:", BASE_OUT)
log_msg("Window size:", WINDOW_SIZE, "bp up/downstream")

## Options
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(future.globals.maxSize = 48 * 1024^3)

safe_pkg <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) install.packages(pkg)
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}
for (p in c("data.table", "future", "future.apply", "RhpcBLASctl", "stringr", "tools", "ieugwasr", "susieR", "coloc", "ggplot2", "patchwork", "httr")) safe_pkg(p)

if (!requireNamespace("locuscomparer", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
  devtools::install_github("caleblare/locuscomparer")
}
library(locuscomparer)

head_preview <- function(dt, name) {
  log_msg(paste0("Head of ", name, ":"))
  capture.output(print(utils::head(dt, n = 5)), file = LOG_FILE, append = TRUE)
}

## -----------------------------------------------------------------------------
## 2. Load Data and Construct Tasks / 
## -----------------------------------------------------------------------------
log_msg("Loading MR results...")
mr_dt <- data.table::fread(MR_RESULTS)
head_preview(mr_dt, "MR Results")

# Filter for IVW + Wald / IVWWald
mr_filtered <- mr_dt[grepl("(?i)Inverse variance weighted|Wald ratio", method)]
if(nrow(mr_filtered) == 0) {
  log_msg("No IVW/Wald ratio results found in MR data. Exiting.")
  quit(status = 0)
}

# Unique pairs / -
tasks_dt <- unique(mr_filtered[, .(id.exposure, id.outcome, exposure, outcome)])
log_msg("Found", nrow(tasks_dt), "unique colocalization tasks.")

log_msg("Loading Gene Info...")
gene_info <- data.table::fread(GENE_INFO)
head_preview(gene_info, "Gene Info")

# Merge to get coordinates / 
# Note: id.exposure matches STUDY column in gene_info
gene_info_sub <- unique(gene_info[, .(STUDY, chromosome, start_position, end_position, resolved_symbol)])
tasks_dt <- merge(tasks_dt, gene_info_sub, by.x = "id.exposure", by.y = "STUDY", all.x = TRUE)

tasks_missing <- tasks_dt[is.na(chromosome)]
if (nrow(tasks_missing) > 0) {
  log_msg("Warning:", nrow(tasks_missing), "tasks missing gene coordinate info, these will be skipped.")
}

tasks_dt <- tasks_dt[!is.na(chromosome)]
tasks_dt[, start_max := pmax(1, start_position - WINDOW_SIZE)]
tasks_dt[, end_max := end_position + WINDOW_SIZE]
tasks_dt[, task_id := paste0(id.exposure, "_vs_", id.outcome)]
tasks_dt[, exposure_clean := trimws(gsub("(?i)\\\\s*\\(FinnGen\\)|\\\\s*\\(UKB\\)", "", exposure))]
tasks_dt[, outcome_clean := trimws(gsub("(?i)\\\\s*\\(FinnGen\\)|\\\\s*\\(UKB\\)", "", outcome))]

if (!is.null(SPECIFIC_TASKS)) {
  task_list <- unlist(strsplit(SPECIFIC_TASKS, ","))
  tasks_dt <- tasks_dt[task_id %in% task_list]
  log_msg("Specific tasks requested. Running only:", length(task_list), "tasks.")
} else if (TEST_N > 0 && TEST_N < nrow(tasks_dt)) {
  log_msg("Test mode enabled. Running only the first", TEST_N, "tasks.")
  tasks_dt <- tasks_dt[1:TEST_N]
}

## Load LD BIM into memory for allele alignment / LDBIM
log_msg("Loading LD panel BIM file for allele alignment...")
bim_file <- paste0(LD_DIR, ".bim")
if (!file.exists(bim_file)) stop("LD BIM file not found:", bim_file)
ld_bim <- data.table::fread(bim_file, select = c(2, 5, 6), col.names = c("snp", "a1", "a2"))
data.table::setkey(ld_bim, snp)

## -----------------------------------------------------------------------------
## 3. Helper Functions / 
## -----------------------------------------------------------------------------

## Build coloc dataset
build_coloc_dataset <- function(dt, type = "quant", s = NULL) {
  maf <- dt$eaf
  maf <- pmin(maf, 1 - maf)
  l <- list(beta = dt$beta, varbeta = dt$se^2, MAF = maf, N = dt$N, type = type, snp = dt$snp)
  if (type == "cc") l$s <- s
  return(l)
}

## Run SuSiE
run_susie <- function(dataset, R) {
  z <- dataset$beta / sqrt(dataset$varbeta)
  z[!is.finite(z)] <- 0
  n <- suppressWarnings(as.integer(stats::median(dataset$N, na.rm = TRUE)))
  if (is.na(n) || n <= 0) n <- length(z)

  snp_order <- intersect(dataset$snp, rownames(R))
  if (length(snp_order) == 0) stop("No intersection between dataset SNPs and LD matrix")

  R <- R[snp_order, snp_order, drop = FALSE]
  z <- z[match(snp_order, dataset$snp)]

  sus <- tryCatch({
    susieR::susie_rss(z = z, R = R, n = n, L = 10, coverage = 0.95, max_iter = 200)
  }, error = function(e) {
    tryCatch({
       susieR::susie_rss(z = z, LD = R, n = n, L = 10, coverage = 0.95, max_iter = 200)
    }, error = function(e2) { stop(paste("susie_rss failed:", e$message)) })
  })

  sus <- coloc::annotate_susie(sus, snp_order, LD = R)
  sus$snp_names <- snp_order
  return(sus)
}

## Visualization
create_locuscompare_plot <- function(merged_dt, title1, title2, results_dir, plot_filename, lead_snp, ld_mat) {
  tryCatch({
    # Prepare merged dataframe for locuscompare
    merged_data <- data.frame(
      rsid = merged_dt$snp,
      logp1 = -log10(pmax(merged_dt$pval_out, 1e-300)), # Outcome (Title1)
      logp2 = -log10(pmax(merged_dt$pval_exp, 1e-300)), # Exposure (Title2)
      chr = merged_dt$chr,
      pos = merged_dt$pos,
      stringsAsFactors = FALSE
    )
    
    # CRITICAL FIX: locuscomparer internal functions have a bug where merge() reorders 
    # the rsids alphabetically, causing the lead SNP to be miscolored and LD colors 
    # mapped to the wrong SNPs. Pre-sorting by rsid neutralizes this bug.
    merged_data <- merged_data[order(merged_data$rsid), ]
    
    # Check if lead_snp is provided and valid
    if (is.null(lead_snp) || is.na(lead_snp) || !(lead_snp %in% merged_data$rsid)) {
      # Fallback to the most significant SNP in Exposure
      lead_snp <- merged_data$rsid[which.max(merged_data$logp2)]
    }
    
    # Prepare LD dataframe from the local LD matrix
    if (lead_snp %in% rownames(ld_mat)) {
      r_values <- ld_mat[lead_snp, ]
      r2_values <- r_values^2
      ld_data <- data.frame(
        SNP_A = lead_snp,
        SNP_B = names(r2_values),
        R2 = r2_values,
        stringsAsFactors = FALSE
      )
    } else {
      # Fallback dummy LD if lead_snp is somehow not in the LD matrix
      ld_data <- data.frame(SNP_A = lead_snp, SNP_B = merged_data$rsid, R2 = 0, stringsAsFactors = FALSE)
      ld_data$R2[ld_data$SNP_B == lead_snp] <- 1
    }
    
    plot_file <- file.path(results_dir, plot_filename)
    
    # Call the internal locuscomparer function to bypass buggy API calls
    res <- tryCatch({
      locuscomparer:::make_combined_plot(
        merged = merged_data,
        title1 = title1,
        title2 = title2,
        ld = ld_data,
        chr = unique(merged_data$chr)[1],
        snp = lead_snp,
        combine = TRUE,
        legend = TRUE,
        legend_position = "bottomright",
        lz_ylab_linebreak = FALSE
      )
    }, error = function(e) {
      log_msg("Locuscompare plotting failed with lead_snp, falling back...", e$message)
      # In case of failure, try again without custom LD logic or default
      return(NULL)
    })
    
    if (inherits(res, "ggplot")) {
      ggplot2::ggsave(filename = plot_file, plot = res, width = 12, height = 10, dpi = 300, bg = "white")
      return(plot_file)
    } else if (is.list(res) && length(res) >= 1 && inherits(res[[1]], "ggplot")) {
      ggplot2::ggsave(filename = plot_file, plot = res[[1]], width = 12, height = 10, dpi = 300, bg = "white")
      return(plot_file)
    } else {
      log_err("Visualization failed: generated object is not a valid ggplot.")
      return(NA)
    }
  }, error = function(e) { log_err("Visualization overall failed:", e$message); return(NA) })
}

## -----------------------------------------------------------------------------
## 4. Single Task Execution Logic / 
## -----------------------------------------------------------------------------
execute_task <- function(task_row) {
  RhpcBLASctl::blas_set_num_threads(THREADS_PER_TASK)
  
  t_id <- task_row$task_id
  exp_id <- task_row$id.exposure
  out_id <- task_row$id.outcome
  chr_val <- task_row$chromosome
  s_max <- task_row$start_max
  e_max <- task_row$end_max
  exp_name <- task_row$exposure_clean
  out_name <- task_row$outcome_clean
  gene_name <- task_row$resolved_symbol
  
  out_folder <- file.path(DIR_RESULTS, t_id)
  out_folder_gene <- file.path(out_folder, "Gene_Level")
  out_folder_snp <- file.path(out_folder, "SNP_Level")
  
  if (!dir.exists(out_folder)) dir.create(out_folder, recursive = TRUE)
  if (!dir.exists(out_folder_gene)) dir.create(out_folder_gene, recursive = TRUE)
  if (!dir.exists(out_folder_snp)) dir.create(out_folder_snp, recursive = TRUE)
  
  out_prefix_gene <- file.path(out_folder_gene, t_id)
  out_prefix_snp <- file.path(out_folder_snp, t_id)
  
  if (RESUME && file.exists(paste0(out_prefix_gene, "_coloc_abf_global.csv"))) {
    return(data.table::data.table(task_id = t_id, status = "skipped"))
  }
  
  ## Read Data
  exp_file <- file.path(PQTL_DIR, paste0(exp_id, ".tsv.gz"))
  out_file <- file.path(IMMUNE_DIR, paste0(out_id, ".tsv.gz"))
  
  if (!file.exists(exp_file) || !file.exists(out_file)) {
    return(data.table::data.table(task_id = t_id, status = "missing_data"))
  }
  
  # Read subset directly using awk for extreme performance (avoid loading full GWAS into memory)
  # awk condition: column 2 is chr, column 3 is pos
  cmd_exp <- sprintf("gzcat '%s' | awk 'NR==1 || ($2==%d && $3>=%d && $3<=%d)'", exp_file, chr_val, s_max, e_max)
  cmd_out <- sprintf("gzcat '%s' | awk 'NR==1 || ($2==%d && $3>=%d && $3<=%d)'", out_file, chr_val, s_max, e_max)
  
  dt_exp <- tryCatch({ data.table::fread(cmd = cmd_exp) }, error = function(e) NULL)
  dt_out <- tryCatch({ data.table::fread(cmd = cmd_out) }, error = function(e) NULL)
  
  if (is.null(dt_exp) || is.null(dt_out) || nrow(dt_exp) < 10 || nrow(dt_out) < 10) {
    return(data.table::data.table(task_id = t_id, status = "few_snps"))
  }
  
  # Inner join by SNP
  merged <- merge(dt_exp, dt_out, by = c("snp", "chr", "pos"), suffixes = c("_exp", "_out"))
  if (nrow(merged) < 10) return(data.table::data.table(task_id = t_id, status = "few_snps_after_merge"))
  
  # Global Alignment Rule (Anchor to LD Panel) / （LD）
  # Extract LD BIM info
  merged <- merge(merged, ld_bim, by = "snp")
  if (nrow(merged) < 10) return(data.table::data.table(task_id = t_id, status = "few_snps_after_ld_match"))
  
  # Align Exposure
  status_exp <- ifelse(merged$A1_exp == merged$a1 & merged$A2_exp == merged$a2, "consistent",
                ifelse(merged$A1_exp == merged$a2 & merged$A2_exp == merged$a1, "flipped", "inconsistent"))
  
  # Align Outcome
  status_out <- ifelse(merged$A1_out == merged$a1 & merged$A2_out == merged$a2, "consistent",
                ifelse(merged$A1_out == merged$a2 & merged$A2_out == merged$a1, "flipped", "inconsistent"))
  
  # Keep only valid
  valid_idx <- which(status_exp != "inconsistent" & status_out != "inconsistent")
  if (length(valid_idx) < 10) return(data.table::data.table(task_id = t_id, status = "few_snps_after_alignment"))
  merged <- merged[valid_idx]
  status_exp <- status_exp[valid_idx]
  status_out <- status_out[valid_idx]
  
  # Adjust Beta & EAF
  merged[, final_beta_exp := ifelse(status_exp == "flipped", -beta_exp, beta_exp)]
  merged[, final_eaf_exp := ifelse(status_exp == "flipped", 1 - eaf_exp, eaf_exp)]
  merged[, final_beta_out := ifelse(status_out == "flipped", -beta_out, beta_out)]
  merged[, final_eaf_out := ifelse(status_out == "flipped", 1 - eaf_out, eaf_out)]
  
  # Set final A1/A2 strictly to LD
  merged[, A1 := a1]
  merged[, A2 := a2]
  
  # Compute LD Matrix
  snps_for_ld <- unique(merged$snp)
  ld_mat <- tryCatch({
    ieugwasr::ld_matrix_local(variants = snps_for_ld, bfile = LD_DIR, plink_bin = PLINK_BIN, with_alleles = FALSE)
  }, error = function(e) NULL)
  
  if (is.null(ld_mat)) return(data.table::data.table(task_id = t_id, status = "ld_fail"))
  
  # Clean LD matrix
  m <- if(is.data.frame(ld_mat)) as.matrix(ld_mat) else ld_mat
  common <- intersect(merged$snp, colnames(m))
  if (length(common) < 5) return(data.table::data.table(task_id = t_id, status = "ld_few_snps"))
  m <- m[common, common, drop = FALSE]
  m[is.na(m)] <- 0; diag(m) <- 1; m[m > 1] <- 1; m[m < -1] <- -1
  m[lower.tri(m)] <- t(m)[lower.tri(m)]
  
  merged <- merged[snp %in% common]
  ord <- match(rownames(m), merged$snp)
  merged <- merged[ord]
  
  # Build Datasets
  d_exp <- build_coloc_dataset(merged[, .(snp, beta = final_beta_exp, se = se_exp, eaf = final_eaf_exp, N = N_exp)], type = "quant")
  d_out <- build_coloc_dataset(merged[, .(snp, beta = final_beta_out, se = se_out, eaf = final_eaf_out, N = N_out)], type = "quant")
  
  # Traditional ABF
  coloc_abf <- tryCatch({ coloc::coloc.abf(dataset1 = d_exp, dataset2 = d_out) }, error = function(e) NULL)
  pp_h4_abf <- NA; lead_snp_abf <- NA; max_snp_pp_h4 <- NA
  if (!is.null(coloc_abf) && !is.null(coloc_abf$summary)) {
      abf_sum <- as.list(coloc_abf$summary)
      pp_h4_abf <- abf_sum$PP.H4.abf
      
      abf_res <- data.table::as.data.table(coloc_abf$results)
      if ("SNP.PP.H4" %in% colnames(abf_res)) {
          idx_max <- which.max(abf_res$SNP.PP.H4)
          if (length(idx_max) > 0) lead_snp_abf <- as.character(abf_res$snp[idx_max])
          max_snp_pp_h4 <- max(abf_res$SNP.PP.H4, na.rm = TRUE)
      }
      
      abf_dt <- data.table::data.table(task_id = t_id, exposure = exp_name, outcome = out_name, gene = gene_name, chr = chr_val, start = s_max, end = e_max, n_snps = nrow(merged), PP.H0.abf = abf_sum$PP.H0.abf, PP.H1.abf = abf_sum$PP.H1.abf, PP.H2.abf = abf_sum$PP.H2.abf, PP.H3.abf = abf_sum$PP.H3.abf, PP.H4.abf = abf_sum$PP.H4.abf, lead_snp = lead_snp_abf, method = "Traditional_ABF", level = "Gene_Region")
      data.table::fwrite(abf_dt, paste0(out_prefix_gene, "_coloc_abf_global.csv"))
      
      abf_res[, `:=`(task_id = t_id, exposure = exp_name, outcome = out_name, gene = gene_name, chr = chr_val, start = s_max, end = e_max, n_snps = nrow(merged), PP.H0.abf = abf_sum$PP.H0.abf, PP.H1.abf = abf_sum$PP.H1.abf, PP.H2.abf = abf_sum$PP.H2.abf, PP.H3.abf = abf_sum$PP.H3.abf, PP.H4.abf = abf_sum$PP.H4.abf, lead_snp = lead_snp_abf, method = "Traditional_ABF", level = "SNP")]
      data.table::fwrite(abf_res, paste0(out_prefix_snp, "_coloc_abf_summary.csv"))
      
      # Save strong SNP level signals per task
      if (!is.na(max_snp_pp_h4) && max_snp_pp_h4 > 0.8) {
          data.table::fwrite(abf_res[SNP.PP.H4 > 0.8], paste0(out_prefix_snp, "_strong_snps_abf.csv"))
      }
  }
  
  # SuSiE
  sus_exp <- tryCatch({ run_susie(d_exp, m) }, error = function(e) { log_err(paste("SuSiE error (Exposure) in task", t_id, ":", e$message)); NULL })
  sus_out <- tryCatch({ run_susie(d_out, m) }, error = function(e) { log_err(paste("SuSiE error (Outcome) in task", t_id, ":", e$message)); NULL })
  
  pp_h4_susie <- NA; max_susie_pip <- NA
  coloc_res <- NULL
  if (!is.null(sus_exp) && !is.null(sus_out)) {
      coloc_res <- tryCatch({ coloc::coloc.susie(sus_exp, sus_out) }, error = function(e) { log_err(paste("coloc.susie error in task", t_id, ":", e$message)); NULL })
      if (!is.null(coloc_res) && !is.null(coloc_res$summary)) {
          summary_dt <- data.table::as.data.table(coloc_res$summary)
          summary_dt[, `:=`(task_id = t_id, exposure = exp_name, outcome = out_name, gene = gene_name, chr = chr_val, method = "SuSiE", level = "Credible_Set_Pair")]
          data.table::fwrite(summary_dt, paste0(out_prefix_gene, "_coloc_susie_summary.csv"))
          pp_h4_susie <- max(summary_dt$PP.H4.abf, na.rm = TRUE)
          
          if (!is.null(sus_exp$sets$cs)) {
              cs_exp_dt <- data.table::as.data.table(sus_exp$sets$cs)
              cs_exp_dt[, `:=`(method = "SuSiE", level = "SNP_in_CS_Exposure")]
              data.table::fwrite(cs_exp_dt, paste0(out_prefix_snp, "_susie_exp_cs.csv"))
          }
          if (!is.null(sus_out$sets$cs)) {
              cs_out_dt <- data.table::as.data.table(sus_out$sets$cs)
              cs_out_dt[, `:=`(method = "SuSiE", level = "SNP_in_CS_Outcome")]
              data.table::fwrite(cs_out_dt, paste0(out_prefix_snp, "_susie_out_cs.csv"))
          }
      }
  } else {
      if (is.null(sus_exp)) log_msg(paste("Task", t_id, "- SuSiE skipped: Failed to build credible sets for Exposure (likely due to non-positive definite LD matrix or lack of strong signals)."))
      if (is.null(sus_out)) log_msg(paste("Task", t_id, "- SuSiE skipped: Failed to build credible sets for Outcome (likely due to non-positive definite LD matrix or lack of strong signals)."))
  }
  
  # Visualization Trigger / 
  viz_file <- NA
  snp_level_strong <- (!is.na(max_snp_pp_h4) && max_snp_pp_h4 > 0.8)
  susie_level_strong <- (!is.na(pp_h4_susie) && pp_h4_susie > 0.8)
  
  if (snp_level_strong || susie_level_strong) {
      log_msg(t_id, "Strong colocalization detected at SNP/CS level. Generating LocusCompare plot...")
      title1 <- paste0(out_name, " (Outcome)")
      title2 <- paste0(exp_name, " (Exposure)")
      plot_name <- paste0(t_id, "_locuscompare.png")
      viz_file <- create_locuscompare_plot(
          merged_dt = merged,
          title1 = title1, title2 = title2,
          results_dir = out_folder, plot_filename = plot_name, lead_snp = lead_snp_abf, ld_mat = m
      )
      if (is.na(viz_file) || is.null(viz_file)) {
          log_err(paste("Task", t_id, "- Visualization triggered but create_locuscompare_plot returned NA. Please check if there are sufficient matching SNPs or valid P-values."))
      }
  }
  
  # Save harmonized data
  data.table::fwrite(merged[, .(snp, chr, pos, A1, A2, beta_exp = final_beta_exp, se_exp, pval_exp, eaf_exp = final_eaf_exp, beta_out = final_beta_out, se_out, pval_out, eaf_out = final_eaf_out)], paste0(out_folder, "/", t_id, "_harmonized_data.csv"))
  
  return(data.table::data.table(task_id = t_id, status = "ok", pp_h4_abf = pp_h4_abf, pp_h4_susie = pp_h4_susie, viz_file = viz_file))
}

## -----------------------------------------------------------------------------
## 5. Parallel Execution / 
## -----------------------------------------------------------------------------
start_time <- Sys.time()
future::plan(future::multisession, workers = TASKS)

tasks_list <- split(tasks_dt, seq_len(nrow(tasks_dt)))
log_msg("Executing", length(tasks_list), "tasks in parallel...")

results_list <- future.apply::future_lapply(tasks_list, function(tsk) {
  out <- NULL; tries <- 0
  while (tries < 3) {
    out <- tryCatch({ execute_task(tsk) }, error = function(e) { Sys.sleep(1); NULL })
    if (!is.null(out)) break
    tries <- tries + 1
  }
  if (is.null(out)) out <- data.table::data.table(task_id = tsk$task_id, status = "error")
  data.table::fwrite(out, PROGRESS_FILE, append = file.exists(PROGRESS_FILE))
  return(out)
}, future.seed = TRUE, future.scheduling = 10)

res_dt <- data.table::rbindlist(results_list, fill = TRUE)
data.table::fwrite(res_dt, SUMMARY_FILE)

## -----------------------------------------------------------------------------
## 6. Wrap Up & Summary / 
## -----------------------------------------------------------------------------

## Detailed Summaries / 
tryCatch({
  # Traditional Coloc Summary
  trad_files <- list.files(DIR_RESULTS, pattern = "_coloc_abf_global\\.csv$", recursive = TRUE, full.names = TRUE)
  if (length(trad_files) > 0) {
    trad_dt <- data.table::rbindlist(lapply(trad_files, data.table::fread), fill = TRUE)
    data.table::fwrite(trad_dt, file.path(BASE_OUT, "summary_traditional_coloc_global.csv"))
  }

  # SuSiE Coloc Summary
  susie_files <- list.files(DIR_RESULTS, pattern = "_coloc_susie_summary\\.csv$", recursive = TRUE, full.names = TRUE)
  if (length(susie_files) > 0) {
    susie_dt <- data.table::rbindlist(lapply(susie_files, data.table::fread), fill = TRUE)
    data.table::fwrite(susie_dt, file.path(BASE_OUT, "summary_susie_detailed.csv"))
  }

  # Traditional SNP-level Summary (SNP.PP.H4 > 0.8)
  abf_snp_files <- list.files(DIR_RESULTS, pattern = "_coloc_abf_summary\\.csv$", recursive = TRUE, full.names = TRUE)
  if (length(abf_snp_files) > 0) {
    abf_snp_list <- lapply(abf_snp_files, function(f) {
      dt <- data.table::fread(f)
      if ("SNP.PP.H4" %in% names(dt)) {
        return(dt[SNP.PP.H4 > 0.8])
      }
      return(NULL)
    })
    strong_snps <- data.table::rbindlist(abf_snp_list, fill = TRUE)
    if (nrow(strong_snps) > 0) {
      # Sort by PP.H4 descending
      strong_snps <- strong_snps[order(-SNP.PP.H4)]
      data.table::fwrite(strong_snps, file.path(BASE_OUT, "summary_strong_snps_H4_08.csv"))
      log_msg("Extracted", nrow(strong_snps), "strong SNP signals (SNP.PP.H4 > 0.8)")
    }
  }
}, error = function(e) {
  log_err("Failed to generate detailed summaries:", e$message)
})

end_time <- Sys.time()
elapsed <- difftime(end_time, start_time, units = "mins")

stat_lines <- c(
  paste("Start:", start_time),
  paste("End:", end_time),
  paste("Elapsed(mins):", round(as.numeric(elapsed), 2)),
  paste("Total Tasks:", nrow(res_dt)),
  paste("OK:", sum(res_dt$status == "ok", na.rm = TRUE)),
  paste("Skipped:", sum(res_dt$status == "skipped", na.rm = TRUE)),
  paste("Failed:", sum(!(res_dt$status %in% c("ok","skipped")), na.rm = TRUE))
)
write(stat_lines, file = LOG_FILE, append = TRUE)

## Copy script to output
script_copy <- file.path(DIR_SCRIPTS, paste0("susie_coloc_pairwise_pQTL_Immune_", timestamp, ".R"))
file.copy("./", script_copy, overwrite = TRUE)

## Generate README
readme_lines <- c(
  "# Pairwise SuSiE Colocalization Pipeline (hg19) - pQTL vs Immune Cells",
  "# SuSiE (hg19)",
  "",
  paste0("**Created at/**: ", timestamp),
  "",
  "## Overview/",
  "This pipeline performs pairwise colocalization analysis between plasma proteins (eQTL/pQTL) and immune cell phenotypes.",
  "(pQTL). ",
  "Tasks are defined by MR results, and analysis is strictly aligned to the LD reference panel.",
  "MR, LD. ",
  "",
  "## Analysis Parameters/",
  "- Window Size/: +/- 250kb",
  "- Strong Signal Threshold/: PP.H4 > 0.8",
  "",
  "## Output Summary/",
  "- `summary_coloc_results.csv`: Final status and H4 probabilities for all tasks.",
  "- `summary_traditional_coloc_global.csv`: Traditional ABF global summary for all tasks.",
  "- `summary_susie_detailed.csv`: SuSiE colocalization detailed summary.",
  "- `summary_strong_snps_H4_08.csv`: Extracted strong SNP-level colocalization signals (SNP.PP.H4 > 0.8).",
  "- Individual task results and locuscompare plots are in the `results/` folder."
)
writeLines(readme_lines, file.path(DIR_README, "README.md"))

log_msg("Pipeline completed successfully.")
