#!/usr/bin/env Rscript
# ==============================================================================
# [Script]: 2.2.1.1_MR_analysis_template.R
# [Method]: Mendelian Randomization 
# [Step]: MR Analysis based on Colocalization SNPs (/SNP MR)
# 
# [Function]:
# Conduct MR analysis using SNPs identified from colocalization.
#       Computes F-statistics, harmonizes data, runs MR, MR-PRESSO, and generates plots.
#       F, , MR、MR-PRESSO, . 
# 
# [Data Availability / ]:
# Requires directory of exposure SNPs (CSV) and outcome GWAS file.
# ==============================================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(TwoSampleMR)
  library(MRInstruments)
  library(plyr)
  library(dplyr)
  library(ggplot2)
  library(forestplot)
  library(MRPRESSO)
})

option_list <- list(
  make_option(c("--snps_dir"), type="character", default=NULL,
              help="Directory containing exposure SNPs CSV files", metavar="character"),
  make_option(c("--outcome_path"), type="character", default=NULL,
              help="Path to outcome GWAS file (e.g. .tsv.gz)", metavar="character"),
  make_option(c("--outcome_name"), type="character", default="Outcome",
              help="Name of the outcome phenotype [default: %default]", metavar="character"),
  make_option(c("--out_dir"), type="character", default="./MR_Results",
              help="Output directory path [default: %default]", metavar="character"),
  make_option(c("--test_n"), type="integer", default=0,
              help="Test mode: process only N files (0 for all) [default: %default]", metavar="integer")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$snps_dir) || is.null(opt$outcome_path)) {
  print_help(opt_parser)
  stop("Missing required arguments: --snps_dir, --outcome_path")
}

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
result_base_dir <- file.path(opt$out_dir, paste0("MR_Analysis_", timestamp))
dir.create(result_base_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(result_base_dir, "run_log.txt")
sink(log_file, append = TRUE, split = TRUE)

cat("========================================\n")
cat("Mendelian Randomization Analysis\n")
cat("========================================\n")
cat("Start Time:", as.character(Sys.time()), "\n")
cat("SNPs Dir:", opt$snps_dir, "\n")
cat("Outcome:", opt$outcome_path, "\n")

snp_files <- list.files(opt$snps_dir, pattern = "\\.csv$", full.names = TRUE)

if (opt$test_n > 0) {
  snp_files <- snp_files[1:min(opt$test_n, length(snp_files))]
  cat(sprintf("Test mode enabled: Processing %d files\n", length(snp_files)))
}

if (length(snp_files) == 0) {
  stop("No CSV SNP files found in the specified directory.")
}
cat(sprintf("Found %d exposure SNP files\n", length(snp_files)))

# ==============================================================================
# Load Outcome Data
# ==============================================================================
cat("Loading outcome data...\n")
outcome_data <- fread(opt$outcome_path)
outcome_data$PHENO <- opt$outcome_name
outcome_data <- as.data.frame(outcome_data)

# Try to detect columns
cols <- colnames(outcome_data)
out_snp <- intersect(cols, c("rsid", "SNP", "variant_id"))[1]
out_beta <- intersect(cols, c("beta", "BETA", "effect_size"))[1]
out_se <- intersect(cols, c("standard_error", "se", "SE"))[1]
out_eaf <- intersect(cols, c("effect_allele_frequency", "eaf", "EAF", "freq"))[1]
out_pval <- intersect(cols, c("p_value", "pval", "P", "p"))[1]
out_ea <- intersect(cols, c("effect_allele", "A1", "EA"))[1]
out_oa <- intersect(cols, c("other_allele", "A2", "OA", "NEA"))[1]
out_chr <- intersect(cols, c("chromosome", "chr", "CHR"))[1]
out_pos <- intersect(cols, c("base_pair_location", "pos", "BP", "position"))[1]
out_n <- intersect(cols, c("n", "N", "sample_size"))[1]

out_data <- format_data(outcome_data,
                       type = "outcome",
                       phenotype_col = "PHENO",
                       snp_col = out_snp,
                       beta_col = out_beta,
                       se_col = out_se,
                       eaf_col = out_eaf,
                       pval_col = out_pval,
                       effect_allele_col = out_ea,
                       other_allele_col = out_oa,
                       chr_col = out_chr,
                       pos_col = out_pos,
                       samplesize_col = out_n)

palindromic_stats <- data.frame()

# ==============================================================================
# Process Exposures
# ==============================================================================
for (file_idx in 1:length(snp_files)) {
  snp_file <- snp_files[file_idx]
  exposure_name <- tools::file_path_sans_ext(basename(snp_file))
  
  cat(sprintf("\n[%d/%d] Processing: %s\n", file_idx, length(snp_files), exposure_name))
  
  tryCatch({
    result_dir <- file.path(result_base_dir, exposure_name)
    dir.create(result_dir, recursive = TRUE, showWarnings = FALSE)
    
    exp_data <- as.data.frame(fread(snp_file))
    
    exp_data_formatted <- format_data(
      exp_data,
      type = "exposure",
      snp_col = "SNP",
      beta_col = "beta.exposure",
      se_col = "se.exposure",
      eaf_col = "eaf.exposure",
      effect_allele_col = "effect_allele.exposure",
      other_allele_col = "other_allele.exposure",
      pval_col = "pval.exposure",
      samplesize_col = "samplesize.exposure",
      chr_col = "chr.exposure",
      pos_col = "pos.exposure"
    )
    exp_data_formatted$exposure <- exposure_name
    
    if ("eaf.exposure" %in% colnames(exp_data_formatted)) {
      invalid_eaf <- is.na(exp_data_formatted$eaf.exposure) | exp_data_formatted$eaf.exposure < 0 | exp_data_formatted$eaf.exposure > 1
      if (any(invalid_eaf)) exp_data_formatted <- exp_data_formatted[!invalid_eaf, ]
    }
    
    dat <- harmonise_data(exposure_dat = exp_data_formatted, outcome_dat = out_data, action = 2)
    
    empty_res <- data.frame(
        id.exposure = exposure_name, id.outcome = opt$outcome_name, outcome = opt$outcome_name,
        exposure = exposure_name, method = "Wald ratio", nsnp = 0, b = NA, se = NA, pval = NA,
        lo_ci = NA, up_ci = NA, or = NA, or_lci95 = NA, or_uci95 = NA, SNP = NA, F_statistic = NA, R2 = NA, mr_presso_global_test_pval = NA, stringsAsFactors = FALSE
    )
      
    if (is.null(dat) || nrow(dat) == 0) {
      write.csv(empty_res, file.path(result_dir, "mr_results.csv"), row.names = FALSE)
      next
    }
    
    write.csv(dat, file.path(result_dir, "harmonised_data.csv"), row.names = FALSE)
    
    # Calculate F-statistic
    if (all(c("beta.exposure", "se.exposure", "eaf.exposure", "samplesize.exposure") %in% colnames(dat))) {
      beta_val <- dat$beta.exposure; se_val <- dat$se.exposure; eaf_val <- dat$eaf.exposure; N_val <- dat$samplesize.exposure
      numerator <- 2 * beta_val^2 * eaf_val * (1 - eaf_val)
      denominator <- numerator + 2 * se_val^2 * N_val * eaf_val * (1 - eaf_val)
      dat$R2 <- numerator / denominator
      dat$F <- (dat$R2 * (N_val - 2)) / (1 - dat$R2)
      dat$R2[is.infinite(dat$R2) | is.na(dat$R2) | dat$R2 < 0] <- 0
      dat$F[is.infinite(dat$F) | is.na(dat$F) | dat$F < 0] <- 0
    } else {
      dat$R2 <- NA; dat$F <- NA
    }
    
    res <- mr(dat)
    het <- mr_heterogeneity(dat)
    pleio <- mr_pleiotropy_test(dat)
    steiger <- directionality_test(dat)
    
    res$exposure <- exposure_name; res$outcome <- opt$outcome_name
    res$het_Q <- if(nrow(het)>0) het$Q[1] else NA
    res$het_Q_pval <- if(nrow(het)>0) het$Q_pval[1] else NA
    res$pleio_pval <- if(nrow(pleio)>0) pleio$pval else NA
    res$egger_intercept <- if(nrow(pleio)>0) pleio$egger_intercept else NA
    res$steiger_dir <- if(nrow(steiger)>0) steiger$correct_causal_direction[1] else NA
    res$steiger_pval <- if(nrow(steiger)>0) steiger$steiger_pval[1] else NA
    
    res_or <- generate_odds_ratios(res)
    res_or$SNP <- if(nrow(dat)==1) dat$SNP[1] else paste(unique(dat$SNP), collapse=";")
    res_or$F_statistic <- if("F" %in% colnames(dat)) mean(dat$F, na.rm=TRUE) else NA
    res_or$R2 <- if("R2" %in% colnames(dat)) mean(dat$R2, na.rm=TRUE) else NA
    res_or$mr_presso_global_test_pval <- NA
    
    # MR-PRESSO
    if (nrow(dat) >= 4) {
      presso_data <- dat %>% dplyr::select(SNP, beta.exposure, beta.outcome, se.exposure, se.outcome) %>%
        dplyr::rename(rsid=SNP, beta_exposure=beta.exposure, beta_outcome=beta.outcome, se_exposure=se.exposure, se_outcome=se.outcome) %>%
        dplyr::filter(!is.na(beta_exposure) & !is.na(beta_outcome) & !is.na(se_exposure) & !is.na(se_outcome) & se_exposure>0 & se_outcome>0) %>% as.data.frame()
      
      tryCatch({
        presso_results <- mr_presso(BetaOutcome="beta_outcome", BetaExposure="beta_exposure", SdOutcome="se_outcome", SdExposure="se_exposure", OUTLIERtest=TRUE, DISTORTIONtest=TRUE, data=presso_data, NbDistribution=5000, SignifThreshold=0.05)
        global_p <- tryCatch(presso_results$`MR-PRESSO results`$`Global Test`$Pvalue, error = function(e) NA)
        outlier_tbl <- tryCatch(presso_results$`MR-PRESSO results`$`Outlier Test`, error = function(e) NULL)
        
        presso_summary <- data.frame(exposure=exposure_name, outcome=opt$outcome_name, global_test_pval=global_p, snp_count_used=nrow(presso_data), outlier_count=if(is.null(outlier_tbl)) 0 else nrow(outlier_tbl))
        write.csv(presso_summary, file.path(result_dir, "mr_presso_result.csv"), row.names = FALSE)
        res_or$mr_presso_global_test_pval <- global_p
      }, error = function(e) cat("MR-PRESSO failed:", e$message, "\n"))
    }
    
    write.csv(res_or, file.path(result_dir, "mr_results.csv"), row.names = FALSE)
    
    # Plots
    if (nrow(dat) > 2) {
      tryCatch({
        loo <- mr_leaveoneout(dat)
        p_loo <- mr_leaveoneout_plot(loo)
        if(is.list(p_loo)) ggsave(file.path(result_dir, "leaveoneout_plot.png"), p_loo[[1]], width=8, height=min(15, max(6, nrow(loo)*0.3)), limitsize=FALSE)
      }, error=function(e) cat("Leave-one-out failed\n"))
    }
    
    tryCatch({
      p_scatter <- mr_scatter_plot(res, dat)
      for(i in seq_along(p_scatter)) ggsave(file.path(result_dir, paste0("scatter_plot_", i, ".png")), p_scatter[[i]], width=8, height=6)
    }, error=function(e) cat("Scatter plot failed\n"))
    
    tryCatch({
      p_forest <- mr_forest_plot(mr_singlesnp(dat))
      if(length(p_forest)>0) ggsave(file.path(result_dir, "forest_plot.png"), p_forest[[1]], width=8, height=6)
    }, error=function(e) cat("Forest plot failed\n"))
    
  }, error = function(e) {
    cat(sprintf("Error processing %s: %s\n", exposure_name, e$message))
  })
}

# Combine Results
cat("Combining results...\n")
all_files <- list.files(result_base_dir, pattern = "mr_results\\.csv$", recursive = TRUE, full.names = TRUE)
if (length(all_files) > 0) {
  results_all <- unique(rbindlist(lapply(all_files, fread), fill = TRUE))
  results_all$or <- exp(results_all$b)
  results_all$or_lci95 <- exp(results_all$b - 1.96 * results_all$se)
  results_all$or_uci95 <- exp(results_all$b + 1.96 * results_all$se)
  results_all$fdr_pval <- NA
  target_idx <- which(results_all$method %in% c("Inverse variance weighted", "Wald ratio"))
  if (length(target_idx) > 0) results_all$fdr_pval[target_idx] <- p.adjust(results_all$pval[target_idx], method="fdr")
  write.csv(results_all, file.path(result_base_dir, "all_mr_results.csv"), row.names = FALSE)
}

cat("Analysis complete.\n")
sink()
