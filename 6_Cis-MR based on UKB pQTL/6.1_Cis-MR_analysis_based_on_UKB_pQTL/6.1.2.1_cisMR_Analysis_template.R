#!/usr/bin/env Rscript
# ==============================================================================
# [Script]: 6.1.2.1_cisMR_Analysis_template.R
# [Method]: cis-MR Analysis (Two-Sample MR)
# [Step]: Perform Mendelian Randomization using cis-IVs
# 
# [Function]:
# Batch execute Two-Sample MR analysis of exposures (e.g., proteins) on an outcome.
# 
# [Input]:
#   --dir_iv       : Directory containing instrumental variable CSV files
#   --file_outcome : Path to the outcome GWAS summary statistics file
#   --file_pheno   : Path to the phenotype mapping/info CSV file
#   --out_dir      : Base output directory path
# 
# [Output]:
# MR results (OR, 95% CI), heterogeneity/pleiotropy tests, scatter/forest/funnel/LOO plots.
# ==============================================================================

rm(list = ls())
gc()

# 1. Load Dependencies
suppressPackageStartupMessages({
  library(optparse)
  library(TwoSampleMR)
  library(data.table)
  library(dplyr)
  library(parallel)
  library(ggplot2)
  library(readr)
})

# 2. Parse Command Line Arguments
option_list <- list(
  make_option(c("--dir_iv"), type="character", default=NULL, help="Directory containing IV CSV files"),
  make_option(c("--file_outcome"), type="character", default=NULL, help="Path to outcome GWAS summary stats"),
  make_option(c("--file_pheno"), type="character", default=NULL, help="Path to phenotype mapping/info CSV"),
  make_option(c("--out_dir"), type="character", default="./MR_Results", help="Output directory path")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$dir_iv) || is.null(opt$file_outcome) || is.null(opt$file_pheno)) {
  print_help(opt_parser)
  stop("Missing required input arguments.")
}

# Assign Parameters
dir_iv <- opt$dir_iv
file_outcome <- opt$file_outcome
file_pheno <- opt$file_pheno
dir_out_base <- opt$out_dir

# Create Output Directory
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
dir_out <- file.path(dir_out_base, paste0("MR_Analysis_Result_", timestamp))
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

dir_results <- file.path(dir_out, "results")
dir.create(dir_results, showWarnings = FALSE, recursive = TRUE)

#  (Log output file)
dir_logs <- file.path(dir_out, "logs")
dir.create(dir_logs, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(dir_logs, paste0("MR_run_log_", timestamp, ".txt"))
sink(log_file, split = TRUE)

cat("====================================================\n")
cat("Starting MR Analysis...\n")
cat("Time:", as.character(Sys.time()), "\n")
cat("====================================================\n")

# ==============================================================================
#  (Load Metadata and Outcome Data)
# ==============================================================================
cat("\n[1/3] Loading metadata and exposure data...\n")
pheno_info <- fread(file_pheno)

iv_files <- list.files(dir_iv, pattern = "\\.csv$", full.names = TRUE)
cat("  - Found", length(iv_files), "exposure IV files to process.\n")

# ,  SNP  (Pre-read all exposures to get unique SNPs)
cat("  - Pre-reading exposure IVs to extract all required SNPs...\n")
all_exp_list <- lapply(iv_files, function(f) {
  dat <- fread(f)
  if(nrow(dat) > 0 && "SNP" %in% names(dat)) {
    return(dat)
  }
  return(NULL)
})
names(all_exp_list) <- basename(iv_files)

#  SNP (Get all unique SNPs)
all_snps <- unique(unlist(lapply(all_exp_list, function(x) if(!is.null(x)) x$SNP else character(0))))
cat("  - Total unique SNPs across all exposures:", length(all_snps), "\n")

rm(all_exp_list)
gc(verbose = FALSE)

# ,  (Load outcome data into memory and filter immediately)
cat("  - Loading and filtering FinnGen outcome data (this may take a moment)...\n")
outcome_raw <- fread(file_outcome)
outcome_subset_global <- outcome_raw[rsids %in% all_snps]
rm(outcome_raw) # 
gc(verbose = FALSE)
cat("  - Outcome data loaded and filtered. Dimensions:", dim(outcome_subset_global)[1], "rows x", dim(outcome_subset_global)[2], "cols.\n")

# ==============================================================================
#  MR  (Define MR processing function for a single exposure)
# ==============================================================================
process_mr <- function(file_path) {
  
  #  (Regularly clean memory in worker processes)
  gc(verbose = FALSE)
  
  file_name <- basename(file_path)
  gc_id <- gsub("\\.csv$", "", file_name)
  
  #  (Extract phenotype name)
  pheno_name <- gc_id
  if ("STUDY" %in% names(pheno_info) && "DISEASE/TRAIT" %in% names(pheno_info)) {
    match_idx <- which(pheno_info$STUDY == gc_id)
    if (length(match_idx) > 0) {
      pheno_name <- pheno_info$`DISEASE/TRAIT`[match_idx[1]]
    }
  }
  
  # 1.  (Read exposure data directly to save memory instead of keeping all in RAM)
  exp_dat <- fread(file_path)
  if (nrow(exp_dat) == 0) return(NULL)
  
  #  TwoSampleMR  (Format to TwoSampleMR exposure format)
  exp_dat_fmt <- format_data(
    as.data.frame(exp_dat),
    type = "exposure",
    snp_col = "SNP",
    beta_col = "beta.exposure",
    se_col = "se.exposure",
    effect_allele_col = "effect_allele.exposure",
    other_allele_col = "other_allele.exposure",
    pval_col = "pval.exposure",
    eaf_col = "eaf.exposure"
  )
  exp_dat_fmt$exposure <- pheno_name
  exp_dat_fmt$id.exposure <- gc_id
  
  # 2.  SNP  (Extract matched outcome SNPs from memory)
  snps_to_keep <- exp_dat_fmt$SNP
  out_subset <- outcome_subset_global[rsids %in% snps_to_keep]
  if (nrow(out_subset) == 0) return(NULL)
  
  out_dat_fmt <- format_data(
    as.data.frame(out_subset),
    type = "outcome",
    snps = snps_to_keep,
    snp_col = "rsids",
    beta_col = "beta",
    se_col = "sebeta",
    effect_allele_col = "alt",
    other_allele_col = "ref",
    pval_col = "pval",
    eaf_col = "af_alt",
    chr_col = "#chrom",
    pos_col = "pos"
  )
  out_dat_fmt$outcome <- "HZ (FinnGen)"
  out_dat_fmt$id.outcome <- "HZ_FinnGen"
  out_dat_fmt$samplesize.outcome <- 487448
  
  # 3. : action = 2  EAF  0.5  SNP
  harm_dat <- suppressMessages(harmonise_data(
    exposure_dat = exp_dat_fmt, 
    outcome_dat = out_dat_fmt, 
    action = 2
  ))
  
  #  SNP (Keep only valid SNPs)
  harm_dat <- harm_dat[harm_dat$mr_keep == TRUE, ]
  
  #  EAF  0.5  SNP (Explicitly filter out ambiguous palindromic SNPs if they somehow slip through)
  #  SNP (Palindromic SNPs: A/T or C/G)
  is_palindromic <- function(a1, a2) {
    (a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") |
    (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C")
  }
  
  if ("eaf.exposure" %in% names(harm_dat)) {
    #  SNP  EAF  0.42  0.58 
    ambiguous_palindromic <- is_palindromic(harm_dat$effect_allele.exposure, harm_dat$other_allele.exposure) & 
                             harm_dat$eaf.exposure > 0.42 & harm_dat$eaf.exposure < 0.58
    
    #  SNP
    harm_dat <- harm_dat[!ambiguous_palindromic, ]
  }
  
  # : SNP  P  > 5e-8,  (Filter: Outcome P-value > 5e-8)
  if ("pval.outcome" %in% names(harm_dat)) {
    harm_dat <- harm_dat[!is.na(harm_dat$pval.outcome) & harm_dat$pval.outcome > 5e-8, ]
  }
  
  if (nrow(harm_dat) == 0) return(NULL)
  
  #  GCST  (Create subfolder for this GCST)
  dir_gcst <- file.path(dir_results, gc_id)
  dir.create(dir_gcst, showWarnings = FALSE, recursive = TRUE)
  
  #  (Save harmonized data)
  fwrite(harm_dat, file.path(dir_gcst, paste0(gc_id, "_harmonized_data.csv")))
  
  # 4.  MR  (Perform MR analysis)
  #  5  (Use 5 common methods)
  methods_to_run <- c("mr_wald_ratio", "mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_simple_mode", "mr_weighted_mode")
  mr_res <- suppressMessages(mr(harm_dat, method_list = methods_to_run))
  
  if (nrow(mr_res) == 0) return(NULL)
  
  #  OR  (Calculate OR and 95% CI)
  mr_res <- generate_odds_ratios(mr_res)
  fwrite(mr_res, file.path(dir_gcst, paste0(gc_id, "_MR_results.csv")))
  
  # 5.  (Heterogeneity and pleiotropy tests)
  if (nrow(harm_dat) >= 3) {
    plei <- suppressMessages(mr_pleiotropy_test(harm_dat))
    het <- suppressMessages(mr_heterogeneity(harm_dat))
    fwrite(plei, file.path(dir_gcst, paste0(gc_id, "_pleiotropy.csv")))
    fwrite(het, file.path(dir_gcst, paste0(gc_id, "_heterogeneity.csv")))
  }
  
  # 6.  (Visualizations)
  try({
    #  (Scatter plot)
    p1 <- mr_scatter_plot(mr_res, harm_dat)
    if (length(p1) > 0) ggsave(file = file.path(dir_gcst, paste0(gc_id, "_scatter.pdf")), plot = p1[[1]], width = 7, height = 7)
    
    #  (Forest plot)
    res_single <- mr_singlesnp(harm_dat)
    res_single <- res_single[!duplicated(res_single$SNP), ] # SNP, 
    if (nrow(res_single) > 0) {
      p2 <- mr_forest_plot(res_single)
      if (length(p2) > 0) ggsave(file = file.path(dir_gcst, paste0(gc_id, "_forest.pdf")), plot = p2[[1]], width = 7, height = 7)
      
      #  (Funnel plot)
      p4 <- mr_funnel_plot(res_single)
      if (length(p4) > 0) ggsave(file = file.path(dir_gcst, paste0(gc_id, "_funnel.pdf")), plot = p4[[1]], width = 7, height = 7)
    }
    
    #  (Leave-one-out plot)
    if (nrow(harm_dat) >= 3) {
      res_loo <- mr_leaveoneout(harm_dat)
      p3 <- mr_leaveoneout_plot(res_loo)
      if (length(p3) > 0) ggsave(file = file.path(dir_gcst, paste0(gc_id, "_loo.pdf")), plot = p3[[1]], width = 7, height = 7)
    }
  }, silent = TRUE)
  
  return(mr_res)
}

# ==============================================================================
#  (Parallel processing of all exposures)
# ==============================================================================
cat("\n[2/3] Running MR analysis in parallel...\n")
#  data.table  mclapply  (Disable data.table threads to prevent mclapply conflicts)
setDTthreads(1)
# ,  2  (Get available cores, reserve 2)
num_cores <- max(1, detectCores() - 2)
cat("  - Using", num_cores, "cores for parallel computation.\n")

#  mclapply  (Use mclapply with error handling)
all_results_list <- mclapply(iv_files, function(f) {
  tryCatch({
    process_mr(f)
  }, error = function(e) {
    cat("    [Error] Failed to process", basename(f), ":", e$message, "\n")
    return(NULL)
  })
}, mc.cores = num_cores, mc.preschedule = FALSE)

#  data.table  (Restore data.table threads)
setDTthreads(0)

# ==============================================================================
#  FDR  (Summarize results and FDR correction)
# ==============================================================================
cat("\n[3/3] Summarizing results and applying FDR correction...\n")
all_results <- bind_rows(all_results_list)

if (nrow(all_results) > 0) {
  #  FDR  (Extract main methods for FDR: IVW or Wald ratio)
  main_res <- all_results %>%
    filter(method %in% c("Inverse variance weighted", "Wald ratio"))
  
  #  FDR  (Apply FDR correction)
  main_res$fdr <- p.adjust(main_res$pval, method = "fdr")
  
  #  FDR  1 (Ensure FDR does not exceed 1)
  main_res$fdr <- ifelse(main_res$fdr > 1, 1, main_res$fdr)
  
  #  FDR  (Merge FDR back to all_results)
  all_results <- all_results %>%
    left_join(main_res %>% select(id.exposure, id.outcome, method, fdr), 
              by = c("id.exposure", "id.outcome", "method"))
  
  #  FDR  pval  (Sort by FDR or pval)
  all_results <- all_results %>%
    arrange(ifelse(is.na(fdr), 1, 0), fdr, pval)
  
  #  (Save global summary)
  summary_file <- file.path(dir_out, paste0("All_MR_Results_Summary_FDR_", timestamp, ".csv"))
  fwrite(all_results, summary_file)
  
  cat("  - Processed", length(all_results_list), "files.\n")
  cat("  - Global summary saved to:", summary_file, "\n")
} else {
  cat("  - No significant MR results generated across all exposures.\n")
}

# ==============================================================================
#  (Generate statistics info)
# ==============================================================================

end_time <- Sys.time()
run_time <- difftime(end_time, start_time, units = "mins")

cat("====================================================\n")
cat("MR Analysis completed in", round(as.numeric(run_time), 2), "minutes.\n")
cat("====================================================\n")

sink()
cat("MR Analysis fully completed. Logs saved to:", log_file, "\n")
