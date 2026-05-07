#!/usr/bin/env Rscript
# ==============================================================================
# [Script]: 7.6.2.1_cisMR_Analysis_template.R
# [Method]: cis-MR Analysis (Two-Sample MR)
# [Step]: Perform Mendelian Randomization using cis-IVs against Immune Cells
# 
# [Function]:
# Batch execute Two-Sample MR analysis of exposures on outcomes (e.g., 731 immune cells).
# 
# [Input]:
#   --dir_iv       : Directory containing instrumental variable CSV files
#   --dir_outcome  : Directory containing outcome GWAS files
#   --file_pheno   : Path to the phenotype mapping/info CSV file
#   --out_dir      : Base output directory path
# 
# [Output]:
# MR results, heterogeneity/pleiotropy tests, and visual plots.
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
  make_option(c("--dir_outcome"), type="character", default=NULL, help="Path to outcome GWAS summary stats directory"),
  make_option(c("--file_pheno"), type="character", default=NULL, help="Path to phenotype mapping/info CSV"),
  make_option(c("--out_dir"), type="character", default="./MR_Results", help="Output directory path")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$dir_iv) || is.null(opt$dir_outcome)) {
  print_help(opt_parser)
  stop("Missing required input arguments.")
}

# Assign Parameters
dir_iv <- opt$dir_iv
dir_outcome <- opt$dir_outcome
file_pheno <- opt$file_pheno
dir_out_base <- opt$out_dir

# Create Output Directory
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
dir_out <- file.path(dir_out_base, paste0("cisMR_Analysis_Result_", timestamp))
dir.create(dir_out, showWarnings = FALSE, recursive = TRUE)

dir_results <- file.path(dir_out, "results")
dir.create(dir_results, showWarnings = FALSE, recursive = TRUE)

#  (Log output file)
dir_logs <- file.path(dir_out, "logs")
dir.create(dir_logs, showWarnings = FALSE, recursive = TRUE)
log_file <- file.path(dir_logs, paste0("cisMR_run_log_", timestamp, ".txt"))
sink(log_file, split = TRUE)

cat("====================================================\n")
cat("Starting cis-MR Analysis...\n")
cat("Time:", as.character(Sys.time()), "\n")
cat("Test Mode:", is_test_mode, "\n")
cat("====================================================\n")

# ==============================================================================
# 2. IV (Get Exposure List and Read IV Data)
# ==============================================================================
cat("\n[1/4] Loading exposure list and extracting required SNPs...\n")

#  (Read significant protein list)
exp_info <- fread(file_exp_sig)
#  GCST  Phenotype (Extract unique GCST and Phenotype)
exp_to_run <- unique(exp_info[, .(GCST, Phenotype)])
exp_gcst_list <- unique(exp_to_run$GCST)

cat("  - Total exposures from list:", length(exp_gcst_list), "\n")

if (is_test_mode) {
  exp_gcst_list <- exp_gcst_list[1:min(20, length(exp_gcst_list))]
  cat("  - TEST MODE ON: Selected", length(exp_gcst_list), "exposures for testing.\n")
}

#  SNP
all_exp_data <- lapply(exp_gcst_list, function(gcst) {
  f <- file.path(dir_iv, paste0(gcst, ".csv"))
  if (file.exists(f)) {
    dat <- fread(f)
    if(nrow(dat) > 0 && "SNP" %in% names(dat)) {
      return(dat)
    }
  }
  return(NULL)
})
names(all_exp_data) <- exp_gcst_list
#  (Remove NULLs)
all_exp_data <- all_exp_data[!sapply(all_exp_data, is.null)]

#  unique SNPs 
all_custom_snps <- unique(unlist(lapply(all_exp_data, function(x) x$SNP)))
cat("  - Successfully loaded", length(all_exp_data), "exposure IV files.\n")
cat("  - Total unique SNPs (rsid) to extract from outcomes:", length(all_custom_snps), "\n")

# ==============================================================================
# 3.  (Prepare Outcome Data List)
# ==============================================================================
cat("\n[2/4] Preparing outcome data list...\n")
outcome_pheno_info <- fread(file_outcome_pheno)
out_files <- list.files(dir_outcome, pattern = "\\.tsv\\.gz$", full.names = TRUE)

if (is_test_mode) {
  out_files <- out_files[1:min(3, length(out_files))]
  cat("  - TEST MODE ON: Selected", length(out_files), "outcomes for testing.\n")
} else {
  cat("  - Total outcomes found:", length(out_files), "\n")
}

# ==============================================================================
# 4. : MR (Batch Processing: Serial by Outcome, Parallel by Exposure)
# ==============================================================================
cat("\n[3/4] Running MR analysis (Batch processing by Outcome)...\n")

#  data.table  mclapply 
setDTthreads(1)
num_cores <- 10 #  10  worker (Force 10 workers)
cat("  - Using", num_cores, "cores for parallel computation within each batch.\n")

all_results_global <- list()

for (out_idx in seq_along(out_files)) {
  out_file <- out_files[out_idx]
  out_gcst <- gsub("\\.tsv\\.gz$", "", basename(out_file))
  
  #  (Get outcome phenotype name)
  out_pheno_name <- out_gcst
  match_idx <- which(outcome_pheno_info$STUDY == out_gcst)
  if (length(match_idx) > 0) {
    out_pheno_name <- outcome_pheno_info$`DISEASE/TRAIT`[match_idx[1]]
  }
  
  out_samplesize <- NA
  if (length(match_idx) > 0 && "samplesize" %in% names(outcome_pheno_info)) {
    out_samplesize <- as.numeric(outcome_pheno_info$samplesize[match_idx[1]])
  }
  
  cat(sprintf("\n  -> Processing Outcome [%d/%d]: %s (%s)\n", out_idx, length(out_files), out_gcst, out_pheno_name))
  
  # ,  (Load and filter outcome data immediately to optimize memory)
  out_raw <- tryCatch({ fread(out_file) }, error = function(e) NULL)
  if (is.null(out_raw) || nrow(out_raw) == 0) {
    cat("     [Warning] Failed to read or empty outcome file. Skipping.\n")
    next
  }
  
  #  rsid  (Filter by rsid to match exposure)
  out_filtered <- out_raw[rsid %in% all_custom_snps]
  
  # , 
  if (!is.na(out_samplesize)) {
    out_filtered$N <- out_samplesize
  }
  
  #  (Clean up large raw data)
  rm(out_raw)
  gc(verbose = FALSE)
  
  if (nrow(out_filtered) == 0) {
    cat("     [Info] No matching SNPs found for this outcome. Skipping.\n")
    next
  }
  
  #  mclapply  (Parallel processing within batch using mclapply)
  # : mc.preschedule = TRUE  FD 
  res_list <- mclapply(names(all_exp_data), function(exp_gcst) {
    tryCatch({
      exp_dat <- all_exp_data[[exp_gcst]]
      exp_pheno <- exp_to_run$Phenotype[exp_to_run$GCST == exp_gcst][1]
      
      # ,  rsid  (Format exposure using rsid)
      exp_fmt <- format_data(
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
      exp_fmt$exposure <- exp_pheno
      exp_fmt$id.exposure <- exp_gcst
      
      #  SNP (Extract matching SNPs for current exposure)
      out_sub <- out_filtered[rsid %in% exp_fmt$SNP]
      if (nrow(out_sub) == 0) return(NULL)
      
      #  (Format outcome data)
      out_fmt <- format_data(
        as.data.frame(out_sub),
        type = "outcome",
        snps = exp_fmt$SNP,
        snp_col = "rsid",
        beta_col = "beta",
        se_col = "standard_error",
        effect_allele_col = "effect_allele",
        other_allele_col = "other_allele",
        pval_col = "p_value",
        eaf_col = "effect_allele_frequency",
        samplesize_col = "N",
        chr_col = "chromosome",
        pos_col = "base_pair_location"
      )
      out_fmt$outcome <- out_pheno_name
      out_fmt$id.outcome <- out_gcst
      
      #  (Harmonise data)
      harm_dat <- suppressMessages(harmonise_data(exp_fmt, out_fmt, action = 2))
      harm_dat <- harm_dat[harm_dat$mr_keep == TRUE, ]
      
      #  SNP (Filter out ambiguous palindromic SNPs)
      is_palindromic <- function(a1, a2) {
        (a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") |
        (a1 == "C" & a2 == "G") | (a1 == "G" & a2 == "C")
      }
      if ("eaf.exposure" %in% names(harm_dat)) {
        ambig <- is_palindromic(harm_dat$effect_allele.exposure, harm_dat$other_allele.exposure) & 
                 harm_dat$eaf.exposure > 0.42 & harm_dat$eaf.exposure < 0.58
        harm_dat <- harm_dat[!ambig, ]
      }
      
      if (nrow(harm_dat) == 0) return(NULL)
      
      #  Outcome-Exposure  (Create subfolder for output)
      dir_res_out_exp <- file.path(dir_results, out_gcst, exp_gcst)
      dir.create(dir_res_out_exp, showWarnings = FALSE, recursive = TRUE)
      
      #  (Save harmonized data)
      fwrite(harm_dat, file.path(dir_res_out_exp, paste0(exp_gcst, "_vs_", out_gcst, "_harmonized_data.csv")))
      
      #  MR  (Perform MR analysis)
      methods_to_run <- c("mr_wald_ratio", "mr_ivw", "mr_egger_regression", "mr_weighted_median", "mr_simple_mode", "mr_weighted_mode")
      mr_res <- suppressMessages(mr(harm_dat, method_list = methods_to_run))
      
      if (nrow(mr_res) == 0) return(NULL)
      
      #  OR  (Calculate OR and 95% CI)
      mr_res <- generate_odds_ratios(mr_res)
      
      #  (Heterogeneity and pleiotropy)
      plei_pval <- NA
      het_pval <- NA
      het_Q <- NA
      
      if (nrow(harm_dat) >= 3) {
        plei <- suppressMessages(mr_pleiotropy_test(harm_dat))
        if (nrow(plei) > 0) {
          plei_pval <- plei$pval[1]
          fwrite(plei, file.path(dir_res_out_exp, paste0(exp_gcst, "_vs_", out_gcst, "_pleiotropy.csv")))
        }
        
        het <- suppressMessages(mr_heterogeneity(harm_dat))
        if (nrow(het) > 0) {
          #  IVW 
          het_ivw <- het[het$method == "Inverse variance weighted", ]
          if (nrow(het_ivw) > 0) {
            het_pval <- het_ivw$Q_pval[1]
            het_Q <- het_ivw$Q[1]
          } else {
            het_pval <- het$Q_pval[1]
            het_Q <- het$Q[1]
          }
          fwrite(het, file.path(dir_res_out_exp, paste0(exp_gcst, "_vs_", out_gcst, "_heterogeneity.csv")))
        }
      }
      
      #  MR 
      mr_res$pleiotropy_pval <- plei_pval
      mr_res$heterogeneity_pval <- het_pval
      mr_res$heterogeneity_Q <- het_Q
      
      fwrite(mr_res, file.path(dir_res_out_exp, paste0(exp_gcst, "_vs_", out_gcst, "_MR_results.csv")))
      
      #  (Plotting)
      try({
        p1 <- mr_scatter_plot(mr_res, harm_dat)
        if (length(p1) > 0) ggsave(file = file.path(dir_res_out_exp, paste0(exp_gcst, "_vs_", out_gcst, "_scatter.pdf")), plot = p1[[1]], width = 7, height = 7)
        
        res_single <- mr_singlesnp(harm_dat)
        res_single <- res_single[!duplicated(res_single$SNP), ]
        if (nrow(res_single) > 0) {
          p2 <- mr_forest_plot(res_single)
          if (length(p2) > 0) ggsave(file = file.path(dir_res_out_exp, paste0(exp_gcst, "_vs_", out_gcst, "_forest.pdf")), plot = p2[[1]], width = 7, height = 7)
          p4 <- mr_funnel_plot(res_single)
          if (length(p4) > 0) ggsave(file = file.path(dir_res_out_exp, paste0(exp_gcst, "_vs_", out_gcst, "_funnel.pdf")), plot = p4[[1]], width = 7, height = 7)
        }
        if (nrow(harm_dat) >= 3) {
          res_loo <- mr_leaveoneout(harm_dat)
          p3 <- mr_leaveoneout_plot(res_loo)
          if (length(p3) > 0) ggsave(file = file.path(dir_res_out_exp, paste0(exp_gcst, "_vs_", out_gcst, "_loo.pdf")), plot = p3[[1]], width = 7, height = 7)
        }
      }, silent = TRUE)
      
      return(mr_res)
      
    }, error = function(e) {
      # , 
      return(NULL)
    })
  }, mc.cores = num_cores, mc.preschedule = TRUE)
  
  #  (Combine batch results)
  valid_res <- res_list[!sapply(res_list, is.null)]
  if (length(valid_res) > 0) {
    all_results_global[[out_gcst]] <- bind_rows(valid_res)
  }
  
  #  (Clean memory after each batch)
  rm(out_filtered)
  gc(verbose = FALSE)
}

#  data.table 
setDTthreads(0)

# ==============================================================================
# 5.  FDR  (Summarize results and FDR correction)
# ==============================================================================
cat("\n[4/4] Summarizing global results and applying FDR correction...\n")

if (length(all_results_global) > 0) {
  all_results_df <- bind_rows(all_results_global)
  
  #  FDR 
  all_results_df$fdr <- NA
  
  #  FDR  (Extract IVW and Wald ratio for FDR)
  main_res_idx <- which(all_results_df$method %in% c("Inverse variance weighted", "Wald ratio"))
  
  if (length(main_res_idx) > 0) {
    #  P  FDR 
    fdr_values <- p.adjust(all_results_df$pval[main_res_idx], method = "fdr")
    fdr_values <- ifelse(fdr_values > 1, 1, fdr_values)
    all_results_df$fdr[main_res_idx] <- fdr_values
  }
  
  #  ( FDR ,  FDR )
  all_results_df <- all_results_df %>%
    arrange(ifelse(is.na(fdr), 1, 0), fdr, pval)
  
  #  (Save all results)
  summary_file <- file.path(dir_out, paste0("All_MR_Results_Summary_FDR_", timestamp, ".csv"))
  fwrite(all_results_df, summary_file)
  cat("  - Global summary saved to:", summary_file, "\n")
  
  # : FDR < 0.05  (id.exposure, id.outcome) 
  sig_pairs <- all_results_df %>% 
    filter(!is.na(fdr) & fdr < 0.05) %>% 
    select(id.exposure, id.outcome) %>% 
    distinct()
    
  if (nrow(sig_pairs) > 0) {
    #  (Extract all methods for significant pairs)
    sig_results <- all_results_df %>%
      inner_join(sig_pairs, by = c("id.exposure", "id.outcome"))
    
    sig_file <- file.path(dir_out, paste0("Significant_MR_Results_FDR0.05_", timestamp, ".csv"))
    fwrite(sig_results, sig_file)
    cat("  - Significant summary (FDR < 0.05) saved to:", sig_file, "\n")
  } else {
    cat("  - No significant results found with FDR < 0.05.\n")
  }
} else {
  cat("  - No successful MR results generated across all batches.\n")
}

# ==============================================================================
# 6.  README  (Generate README and statistics info)
# ==============================================================================
dir_readme <- file.path(dir_out, "readme")
dir.create(dir_readme, showWarnings = FALSE, recursive = TRUE)

end_time <- Sys.time()
run_time <- difftime(end_time, start_time, units = "mins")

readme_content_en <- paste0(
  "Project: cis-MR Analysis of rs1800628-associated UKB Plasma Proteins vs 731 Immune Cells\n",
  "Date: ", Sys.Date(), "\n",
  "Time: ", format(Sys.time(), "%H:%M:%S"), "\n",
  "Output Folder: ", dir_out, "\n\n",
  "Description:\n",
  "This directory contains the results of cis-Mendelian Randomization (MR) analysis between rs1800628-associated UKB plasma proteins (exposures) and 731 immune cell traits (outcomes).\n\n",
  "Workflow:\n",
  "1. Exposures were extracted and IVs matched by GCST.\n",
  "2. Outcome data was loaded in batches, caching in memory to optimize read speed.\n",
  "3. Data harmonization aligned SNPs via rsid and dropped ambiguous palindromic SNPs (EAF near 0.5).\n",
  "4. MR analysis used up to 5 methods. Parallel execution used 10 workers.\n",
  "5. Pleiotropy/Heterogeneity tests and visualizations (scatter, forest, funnel, LOO) were generated.\n",
  "6. FDR multiple testing correction was applied on primary methods (IVW / Wald ratio).\n\n",
  "Run Time: ", round(as.numeric(run_time), 2), " minutes.\n",
  "Total Outcomes Processed: ", length(out_files), "\n"
)

readme_content_zh <- paste0(
  ": rs1800628UKB731cis-MR\n",
  ": ", Sys.Date, "\n",
  ": ", format(Sys.time, "%H:%M:%S"), "\n",
  ": ", dir_out, "\n\n",
  ":\n",
  "rs1800628UKB731(cis-MR). \n\n",
  ":\n",
  "1. GCST. \n",
  "2. , . \n",
  "3. rsid, EAF0.5SNP. \n",
  "4. 5MR, (OR/95% CI). \n",
  "5. , 、、. \n",
  "6. (IVW / Wald)PFDR. \n\n",
  ": ", round(as.numeric(run_time), 2), " . \n",
  ": ", length(out_files), "\n"
)

readme_file_en <- file.path(dir_readme, paste0("cisMR_Analysis_README_EN_", timestamp, ".txt"))
readme_file_zh <- file.path(dir_readme, paste0("cisMR_Analysis_README_ZH_", timestamp, ".txt"))

writeLines(readme_content_en, readme_file_en)
writeLines(readme_content_zh, readme_file_zh)

cat("====================================================\n")
cat("MR Analysis completed in", round(as.numeric(run_time), 2), "minutes.\n")
cat("====================================================\n")

sink()
cat("MR Analysis fully completed. Logs saved to:", log_file, "\n")
