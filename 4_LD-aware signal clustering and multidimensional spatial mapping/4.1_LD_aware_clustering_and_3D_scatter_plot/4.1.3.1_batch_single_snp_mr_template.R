# ==============================================================================
# [Script]: batch_single_snp_mr.R
# [Method]: LD+3d
# [Step]: 
# 
# [Function]:
# Execute the '' step within the 'LD+3d' analytical framework.
# 
# [Parameters / ]:
# Standard predefined thresholds. LD r2 > 0.6 for clustering.
# 
# [Steps / ]:
#   1. Data loading and initialization / 
#   2. Core analytical execution / 
#   3. Results formatting and output / 
# ==============================================================================

#!/usr/bin/env Rscript

# ==============================================================================
# Script Information / 
# ==============================================================================
# Script Name: batch_single_snp_mr.R
# Date: 2026-02-28
# Version: 1.5
# Description: 
#   Performs batch Single-SNP Mendelian Randomization (MR) analysis.
#   Based on pre-defined Signal Sets (Clusters) with R2 > 0.6.
#   Optimized for memory efficiency by processing outcomes sequentially.
#   Key Features:
#     - Reads Signal Set Summary to guide analysis
#     - Generates tasks based on Signal-Outcome pairs
#     - Loads Outcome GWAS data sequentially (one at a time) to save memory
#     - Processes relevant Signals for the current Outcome in parallel
#     - Results saved in cell/cluster specific subdirectories
#     - Clears memory after each Outcome batch
#
#   SNP(MR).
#   (R2 > 0.6).
#   , .
#   :
#     - 
#     - -
#     - GWAS
#     - 
#     - /
#     - 
#
# Usage: 
#   Rscript batch_single_snp_mr.R
# ==============================================================================

# ==============================================================================
# Configuration / 
# ==============================================================================

# ==============================================================================
# Configuration / 
# ==============================================================================

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(TwoSampleMR)
  library(dplyr)
  library(parallel)
  library(glue)
  library(tools)
})

# ==============================================================================
# 0. Command Line Arguments / 
# ==============================================================================
option_list <- list(
  make_option(c("--exposure_summary"), type="character", default=NULL,
              help="Path to Exposure Cluster Summary CSV file / CSV"),
  make_option(c("--exposure_dir"), type="character", default=NULL,
              help="Directory containing Exposure SNP sets / SNP"),
  make_option(c("--outcome_aging"), type="character", default=NULL,
              help="Path to Aging Outcome GWAS file / AgingGWAS"),
  make_option(c("--outcome_ra"), type="character", default=NULL,
              help="Path to RA Outcome GWAS file / RAGWAS"),
  make_option(c("--outcome_hz"), type="character", default=NULL,
              help="Path to HZ Outcome GWAS file / HZGWAS"),
  make_option(c("--out_dir"), type="character", default="./Single_SNP_MR_Results",
              help="Output directory path / ")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$exposure_summary) || is.null(opt$exposure_dir) || is.null(opt$outcome_aging) || is.null(opt$outcome_ra) || is.null(opt$outcome_hz)) {
  print_help(opt_parser)
  stop("Required arguments missing. Please provide all input paths.")
}

# Paths / 
exposure_summary_file <- opt$exposure_summary
exposure_root_dir <- opt$exposure_dir

outcome_files <- list(
  list(
    name = "Aging",
    path = opt$outcome_aging,
    type = tools::file_ext(opt$outcome_aging),
    phenotype = "aging",
    samplesize = 1958774,
    cols = list(snp="SNP", beta="beta", se="se", pval="Pvalue", effect_allele="effect_allele", other_allele="other_allele", eaf="eaf")
  ),
  list(
    name = "RA",
    path = opt$outcome_ra,
    type = tools::file_ext(opt$outcome_ra),
    phenotype = "Rheumatoid arthritis",
    samplesize = 315668,
    cols = list(snp="rsid", beta="beta", se="se", pval="p_value", effect_allele="effect_allele", other_allele="other_allele", eaf="effect_allele_frequency")
  ),
  list(
    name = "HZ",
    path = opt$outcome_hz,
    type = tools::file_ext(opt$outcome_hz),
    phenotype = "herpes zoster",
    samplesize = 458440,
    cols = list(snp="rsid", beta="beta", se="standard_error", pval="p_value", effect_allele="effect_allele", other_allele="other_allele", eaf="effect_allele_frequency")
  )
)

base_output_dir <- opt$out_dir

# ==============================================================================
# Setup / 
# ==============================================================================

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- file.path(base_output_dir, timestamp)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

dir_logs <- file.path(output_dir, "logs")
dir_readme <- file.path(output_dir, "readme")
dir.create(dir_logs, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_readme, recursive = TRUE, showWarnings = FALSE)

# Logger
log_file <- file.path(dir_logs, "execution.log")
log_msg <- function(msg) {
  timestamp_log <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  message(paste(timestamp_log, msg))
  cat(paste(timestamp_log, msg, "\n"), file = log_file, append = TRUE)
}

log_msg("Starting Batch Single-SNP MR Analysis (R2 > 0.6) - Optimized Flow")
log_msg(glue("Output Directory: {output_dir}"))
log_msg(glue("Exposure Summary File: {exposure_summary_file}"))
log_msg(glue("Exposure SNP Sets Root: {exposure_root_dir}"))

# Check Input Data Existence / 
log_msg("Checking input data existence...")
missing_files <- FALSE

if (!file.exists(exposure_summary_file)) {
  log_msg(glue("Error: Exposure summary file not found: {exposure_summary_file}"))
  missing_files <- TRUE
}

if (!dir.exists(exposure_root_dir)) {
  log_msg(glue("Error: Exposure root directory not found: {exposure_root_dir}"))
  missing_files <- TRUE
}

for (info in outcome_files) {
  if (!file.exists(info$path)) {
    log_msg(glue("Error: Outcome file not found for {info$name}: {info$path}"))
    missing_files <- TRUE
  }
}

if (missing_files) {
  log_msg("Critical input files missing. Aborting.")
  quit(save = "no", status = 1)
} else {
  log_msg("All input files verified.")
}

# Copy Script
script_path <- commandArgs(trailingOnly = FALSE)[4]
if (grepl("--file=", script_path)) {
  script_path <- sub("--file=", "", script_path)
  file.copy(script_path, file.path(output_dir, basename(script_path)))
} else {
  script_path <- "./4.1.3.1_batch_single_snp_mr_template.R"
  file.copy(script_path, file.path(output_dir, basename(script_path)))
}

# ==============================================================================
# Helper Functions / 
# ==============================================================================

# Prepare Exposure Data / 
prepare_exposure <- function(exposure_path, cluster_id) {
  if (!file.exists(exposure_path)) return(NULL)
  
  exp_dat <- fread(exposure_path)
  exp_dat <- as.data.frame(exp_dat)
  
  if (!"SNP" %in% names(exp_dat)) {
    # log_msg(glue("Error: SNP column missing in {exposure_path}")) # Avoid logging inside parallel worker if possible
    return(NULL)
  }
  
  if (nrow(exp_dat) == 0) return(NULL)
  
  # Remove duplicate rows if any
  exp_dat <- unique(exp_dat)
  
  # Create unique exposure ID: Cell_Gene_SNP
  # This ensures each row is treated as a unique instrument/exposure
  exp_dat$unique_exposure_id <- paste(exp_dat$CELL_ID, exp_dat$GENE, exp_dat$SNP, sep = "_")
  
  # Format data for TwoSampleMR
  exp_formatted <- format_data(
    exp_dat,
    type = "exposure",
    snp_col = "SNP",
    beta_col = "BETA",
    se_col = "SE",
    eaf_col = "eaf",
    effect_allele_col = "effect_allele",
    other_allele_col = "other_allele",
    pval_col = "P_VALUE",
    samplesize_col = "N",
    phenotype_col = "unique_exposure_id" # Use unique ID as exposure name
  )
  
  # Ensure id.exposure matches exposure
  exp_formatted$id.exposure <- exp_formatted$exposure
  
  # Add cell and gene columns for later retrieval
  meta_info <- unique(exp_dat[, c("unique_exposure_id", "CELL_ID", "GENE")])
  names(meta_info) <- c("exposure", "cell", "gene")
  
  # Merge metadata into exp_formatted
  exp_formatted <- merge(exp_formatted, meta_info, by.x = "exposure", by.y = "exposure", all.x = TRUE)
  
  return(exp_formatted)
}

# Process Single Signal-Outcome Task / -
process_single_mr <- function(task, outcome_data_obj, output_base_dir) {
  cluster_id <- task$cluster_id
  cell <- task$cell
  exp_file <- task$exp_file
  outcome_name <- task$outcome_name
  
  # Create output directory: output/cell/cluster_id
  unit_out_dir <- file.path(output_base_dir, cell, cluster_id)
  if (!dir.exists(unit_out_dir)) {
    dir.create(unit_out_dir, recursive = TRUE, showWarnings = FALSE)
  }
  
  # Read Exposure
  exp_data <- prepare_exposure(exp_file, cluster_id)
  if (is.null(exp_data)) {
    return(list(status = "error", msg = "Failed to prepare exposure data"))
  }
  
  out_raw <- outcome_data_obj$data
  info <- outcome_data_obj$info
  
  # Format Outcome (filter by SNPs in exposure)
  # Pre-filter out_raw using data.table logic if possible for speed, but format_data is safer
  # If out_raw is data.table, subsetting is fast.
  
  # Try simple subset first to reduce size before format_data
  snps_needed <- unique(exp_data$SNP)
  
  # Ensure out_raw is data frame or data table
  if (inherits(out_raw, "data.table")) {
     out_subset <- out_raw[get(info$cols$snp) %in% snps_needed]
     out_subset <- as.data.frame(out_subset)
  } else {
     out_subset <- out_raw[out_raw[[info$cols$snp]] %in% snps_needed, ]
  }
  
  if (nrow(out_subset) == 0) {
    return(list(status = "skipped", msg = "No matching SNPs in Outcome"))
  }

  out_formatted <- format_data(
    out_subset,
    type = "outcome",
    snps = snps_needed,
    snp_col = info$cols$snp,
    beta_col = info$cols$beta,
    se_col = info$cols$se,
    eaf_col = info$cols$eaf,
    effect_allele_col = info$cols$effect_allele,
    other_allele_col = info$cols$other_allele,
    pval_col = info$cols$pval
  )
  
  if (nrow(out_formatted) == 0) {
    return(list(status = "skipped", msg = "No matching SNPs after format"))
  }
  
  out_formatted$samplesize.outcome <- info$samplesize
  out_formatted$outcome <- info$phenotype
  out_formatted$id.outcome <- outcome_name
  
  # Harmonise
  dat <- harmonise_data(exp_data, out_formatted, action = 2)
  
  if (nrow(dat) == 0) {
    return(list(status = "skipped", msg = "No SNPs after harmonisation"))
  }
  
  # Perform Single SNP MR
  res <- mr_singlesnp(dat, single_method = "mr_wald_ratio")
  
  # Combine with data to see SNP details
  res <- merge(res, dat, by = c("id.exposure", "id.outcome", "SNP"))
  
  # Restore metadata
  if (!"cell" %in% names(res) || !"gene" %in% names(res)) {
     meta_cols <- unique(exp_data[, c("exposure", "cell", "gene")])
     res <- merge(res, meta_cols, by.x = "id.exposure", by.y = "exposure", all.x = TRUE)
  }
  
  # Calculate OR and 95% CI
  res$OR <- exp(res$b)
  res$OR_lower <- exp(res$b - 1.96 * res$se)
  res$OR_upper <- exp(res$b + 1.96 * res$se)
  
  # Reorder columns
  desired_order <- c("cell", "gene", "SNP")
  existing_cols <- names(res)
  final_cols <- c(desired_order, setdiff(existing_cols, desired_order))
  final_cols <- intersect(final_cols, existing_cols)
  res <- res[, final_cols, drop = FALSE]
  
  # Save Results
  safe_phenotype <- gsub(" ", "_", info$phenotype)
  out_file_name <- glue("{cluster_id}-{safe_phenotype}.csv")
  fwrite(res, file.path(unit_out_dir, out_file_name))
  
  return(list(status = "success", msg = glue("Success ({nrow(res)} SNPs)")))
}

# ==============================================================================
# Task Generation / 
# ==============================================================================

log_msg("Generating task list from summary file...")
summary_dt <- fread(exposure_summary_file)
log_msg(glue("Read {nrow(summary_dt)} rows from summary file."))

# Build a registry of tasks: Cluster -> Outcome pairs
# : Cluster -> Outcome 
task_registry <- list()

for (i in 1:nrow(summary_dt)) {
  row <- summary_dt[i, ]
  cluster_id <- row$Exposure_Cluster_ID
  cell <- row$Exposure_Cell
  outcomes_str <- row$Outcomes
  
  # Construct Path to SNP file
  exp_file <- file.path(exposure_root_dir, cell, paste0(cluster_id, ".csv"))
  
  if (!file.exists(exp_file)) {
    log_msg(glue("Warning: SNP file not found for {cluster_id}: {exp_file}"))
    next
  }
  
  # Identify Target Outcomes
  outcome_names_raw <- unlist(strsplit(outcomes_str, ";\\s*"))
  outcome_names_raw <- trimws(outcome_names_raw)
  outcome_names_raw <- outcome_names_raw[outcome_names_raw != ""]
  
  # Normalize outcome names to match our config
  # Note: The summary file might use full names or short names. 
  # We need to map them to our 'outcome_files' keys (Aging, RA, HZ).
  
  valid_outcomes <- c()
  for (nm in outcome_names_raw) {
    nm_lower <- tolower(nm)
    found <- FALSE
    
    # Check against our config
    if (nm_lower == "aging") { valid_outcomes <- c(valid_outcomes, "Aging"); found <- TRUE }
    else if (nm_lower == "ra" || nm_lower == "rheumatoid arthritis") { valid_outcomes <- c(valid_outcomes, "RA"); found <- TRUE }
    else if (nm_lower == "hz" || nm_lower == "herpes zoster") { valid_outcomes <- c(valid_outcomes, "HZ"); found <- TRUE }
    
    if (!found) {
      # log_msg(glue("Warning: Unknown outcome '{nm}' for cluster {cluster_id}"))
    }
  }
  
  if (length(valid_outcomes) == 0) {
    log_msg(glue("Warning: No valid outcomes mapped for {cluster_id} (Raw: {outcomes_str})"))
    next
  }
  
  # Register tasks
  for (outcome_name in unique(valid_outcomes)) {
    task_registry[[length(task_registry) + 1]] <- list(
      cluster_id = cluster_id,
      cell = cell,
      exp_file = exp_file,
      outcome_name = outcome_name
    )
  }
}

log_msg(glue("Total individual tasks (Cluster-Outcome pairs) generated: {length(task_registry)}"))

# ==============================================================================
# Execution / 
# ==============================================================================

# Determine number of cores
n_cores <- parallel::detectCores(logical = FALSE)
if (is.na(n_cores) || n_cores < 4) n_cores <- 4
n_cores <- max(1, n_cores - 2)
log_msg(glue("Using {n_cores} cores for parallel processing."))

# Initialize results container
# : stored as a list of lists
all_results_log <- list()

# Loop through each defined outcome type
for (outcome_info in outcome_files) {
  current_outcome_name <- outcome_info$name
  
  # 1. Filter tasks for this outcome
  # 1. 
  current_tasks <- Filter(function(t) t$outcome_name == current_outcome_name, task_registry)
  
  if (length(current_tasks) == 0) {
    log_msg(glue("No tasks found for outcome: {current_outcome_name}. Skipping."))
    next
  }
  
  log_msg(glue("Processing Outcome: {current_outcome_name} ({length(current_tasks)} tasks)..."))
  
  # 2. Load Outcome Data (Once per outcome batch)
  # 2.  
  log_msg(glue("  Loading data from {outcome_info$path}..."))
  outcome_data_obj <- NULL
  tryCatch({
    dt <- fread(outcome_info$path)
    # Keep as data.table for speed if possible, but process_single_mr expects obj with info
    outcome_data_obj <- list(data = dt, info = outcome_info)
    log_msg(glue("  Loaded {nrow(dt)} rows."))
  }, error = function(e) {
    log_msg(glue("  Error loading {current_outcome_name}: {e$message}"))
  })
  
  if (is.null(outcome_data_obj)) {
    log_msg(glue("  Skipping {current_outcome_name} due to load error."))
    next
  }
  
  # 3. Parallel Execution
  # 3. 
  # We pass 'outcome_data_obj' to the function. 
  # Note: In 'mclapply' (forking), large objects in parent env are shared (copy-on-write).
  # Since we only read from outcome_data_obj, this is memory efficient.
  
  batch_results <- mclapply(current_tasks, function(task) {
    res <- tryCatch({
      process_single_mr(task, outcome_data_obj, output_dir)
    }, error = function(e) {
      return(list(status = "error", msg = e$message))
    })
    
    # Return structured result
    return(list(
      cluster_id = task$cluster_id,
      cell = task$cell,
      outcome = task$outcome_name,
      status = res$status,
      msg = as.character(res$msg)
    ))
  }, mc.cores = n_cores)
  
  # 4. Collect Results
  # 4. 
  all_results_log <- c(all_results_log, batch_results)
  
  # 5. Cleanup
  # 5. 
  log_msg(glue("  Finished {current_outcome_name}. Cleaning memory..."))
  rm(outcome_data_obj)
  rm(batch_results)
  gc()
}

# ==============================================================================
# Summary / 
# ==============================================================================

log_msg("Generating final summary...")

# Convert results log to data frame for easier aggregation
if (length(all_results_log) > 0) {
  results_df <- rbindlist(all_results_log)
} else {
  results_df <- data.table(cluster_id=character(), cell=character(), outcome=character(), status=character(), msg=character())
}

# We need to pivot/aggregate to match the original format: One row per Cluster
# Format: Cell, Cluster_ID, Status, Details (Outcome1: Status; Outcome2: Status)

if (nrow(results_df) > 0) {
  # Create a summary string per row
  results_df[, detail_str := paste0(outcome, ": ", msg)]
  
  # Aggregate by Cluster
  summary_agg <- results_df[, .(
    Cell = first(cell),
    Details = paste(detail_str, collapse = "; "),
    Overall_Status = ifelse(all(status == "success"), "success", "mixed") # Simple status logic
  ), by = cluster_id]
  
  # Rename for output
  setnames(summary_agg, "cluster_id", "Cluster_ID")
  setnames(summary_agg, "Overall_Status", "Status")
  
  # Write CSV
  fwrite(summary_agg, file.path(dir_readme, "mr_analysis_summary.csv"))
  
  total_processed <- nrow(summary_agg)
} else {
  total_processed <- 0
  writeLines("No results generated.", file.path(dir_readme, "mr_analysis_summary.csv"))
}

# Readme
readme_content <- glue("
# Single-SNP MR Analysis Summary / SNP MR

## Project Information
- **Project**: Multi-omics and Herpes Zoster (Single SNP MR)
- **Date**: {Sys.Date()}
- **Script**: batch_single_snp_mr.R
- **Parameters**: LD R2 >= 0.6

## Overview
This analysis performs Single-SNP Mendelian Randomization for each SNP identified in the LD analysis.
Analysis is driven by Signal Sets (Clusters) defined in `Exposure_Cluster_Summary.csv`.

## Input Data
- **Exposure Summary**: {basename(exposure_summary_file)}
- **SNP Sets Root**: {exposure_root_dir}
- **Outcomes**:
  1. Aging (N=1,958,774)
  2. Rheumatoid Arthritis (N=315,668)
  3. Herpes Zoster (N=458,440)

## Methodology
1. **Signal Set Processing**: Iterate through each signal set in the summary file.
2. **Exposure Loading**: Load corresponding SNP file from `cell/cluster_id.csv`.
3. **Outcome Selection**: Only run MR for outcomes specified in the summary file for that cluster.
4. **MR Analysis**: Perform Single SNP MR (`mr_singlesnp`), calculating Wald Ratio.
5. **Statistics**: Calculate OR and 95% CI.

## Results
- Results are saved in: `[Cell]/[Cluster_ID]/[Cluster_ID]-[Phenotype].csv`

## Execution Status
- Total Clusters Processed: {total_processed}
- See `mr_analysis_summary.csv` for details.
")

writeLines(readme_content, file.path(dir_readme, "README.md"))

log_msg("Analysis Complete.")
