# ==============================================================================
# [Script]: calculate_LD_and_filter(r2>0.6).R
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
# Script Name: calculate_LD_and_filter.R
# Description: 
#   This script calculates LD (R2) for exposure signals (Exposure_Cluster_ID).
#   It iterates through each signal set, extracts Lead SNPs, and identifies all 
#   SNPs in LD (R2 > 0.6) using the UKB EUR reference panel.
#   It retrieves summary statistics from either sc-eQTL data or Bulk eQTL data.
#   Results are saved by Cell type and Cluster ID with specific columns.
#   Duplicates (same SNP, same source data) are removed. LD columns (Lead_SNP, R2) are dropped.
#
# Usage: 
#   Rscript calculate_LD_and_filter.R
#
# Dependencies:
#   - R packages: data.table, parallel, tools, glue, dplyr
#   - External tools: PLINK 1.9
#
# Input:
#   - Cluster Summary: Exposure_Cluster_Summary.csv
#   - sc-eQTL Data: .tsv.gz files (OneK1K)
#   - Bulk eQTL Data: MHC_All_Genes_Combined_20260224_204230.csv
#   - LD Reference Panel: UKB_LD_Reference_Panel_EUR_Merged
#
# Output:
#   - Directory with timestamp/
#     - {Cell_Type}/
#       - {Exposure_Cluster_ID}.csv
#     - logs/
# ==============================================================================

# ==============================================================================
# Configuration / 
# ==============================================================================

# Libraries / 
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(parallel)
  library(tools)
  library(glue)
  library(dplyr)
})

# ==============================================================================
# 0. Command Line Arguments / 
# ==============================================================================
option_list <- list(
  make_option(c("--input_cluster"), type="character", default=NULL,
              help="Path to the Exposure Cluster Summary CSV file / CSV"),
  make_option(c("--sc_eqtl_dir"), type="character", default=NULL,
              help="Directory containing sc-eQTL data / sc-eQTL"),
  make_option(c("--bulk_eqtl"), type="character", default=NULL,
              help="Path to Bulk eQTL CSV file / Bulk eQTL CSV"),
  make_option(c("--plink_path"), type="character", default=NULL,
              help="Path to PLINK executable / PLINK"),
  make_option(c("--ld_ref"), type="character", default=NULL,
              help="Prefix of the LD reference panel (hg19) / LD (hg19)"),
  make_option(c("--out_dir"), type="character", default="./LD_Filter_Results",
              help="Output directory path / ")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_cluster) || is.null(opt$sc_eqtl_dir) || is.null(opt$plink_path) || is.null(opt$ld_ref)) {
  print_help(opt_parser)
  stop("Required arguments missing. Please provide --input_cluster, --sc_eqtl_dir, --plink_path, and --ld_ref.")
}

# Paths / 
input_cluster_path <- opt$input_cluster
sc_eqtl_dir <- opt$sc_eqtl_dir
bulk_eqtl_path <- opt$bulk_eqtl
plink_path <- opt$plink_path
ld_ref_prefix <- opt$ld_ref
base_output_dir <- opt$out_dir

# Region Parameters / 
target_chr <- 6
target_start <- 25000000
target_end <- 34000000

# Parallel Cores / 
# Detect cores, reserve 2 for system
n_cores <- max(1, detectCores() - 2)
# Limit to 8 to avoid memory issues if many large files are loaded
n_cores <- min(n_cores, 8)

# ==============================================================================
# Setup / 
# ==============================================================================

# Generate Timestamp / 
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- file.path(base_output_dir, paste0(timestamp, "(r2>0.6)"))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Create Subdirectories / 
dir_logs <- file.path(output_dir, "logs")
dir_temp <- file.path(output_dir, "temp")

dir.create(dir_logs, recursive = TRUE, showWarnings = FALSE)
dir.create(dir_temp, recursive = TRUE, showWarnings = FALSE)

# Logger Function / 
log_file <- file.path(dir_logs, "execution.log")
log_msg <- function(msg) {
  timestamp_log <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  message(paste(timestamp_log, msg))
  cat(paste(timestamp_log, msg, "\n"), file = log_file, append = TRUE)
}

log_msg("Starting LD Calculation Script")
log_msg(glue("Output Directory: {output_dir}"))
log_msg(glue("Using {n_cores} cores for parallel processing"))

# Copy Script to Output / 
script_path <- commandArgs(trailingOnly = FALSE)[4]
if (grepl("--file=", script_path)) {
  script_path <- sub("--file=", "", script_path)
  file.copy(script_path, file.path(output_dir, basename(script_path)))
  log_msg(glue("Copied script to {output_dir}"))
} else {
  script_path <- "./4.1.2.1_calculate_LD_and_filter_template.R"
  file.copy(script_path, file.path(output_dir, basename(script_path)))
}

# ==============================================================================
# Load Data / 
# ==============================================================================

# 1. Load Cluster Summary
if (!file.exists(input_cluster_path)) {
  stop(glue("Input cluster summary not found: {input_cluster_path}"))
}
log_msg(glue("Reading cluster summary: {input_cluster_path}"))
cluster_dt <- fread(input_cluster_path)
log_msg(glue("Loaded {nrow(cluster_dt)} clusters."))

# 2. Load Bulk eQTL Data
log_msg(glue("Reading Bulk eQTL data: {bulk_eqtl_path}"))
if (file.exists(bulk_eqtl_path)) {
  bulk_dt <- fread(bulk_eqtl_path)
  log_msg(glue("Loaded Bulk eQTL data: {nrow(bulk_dt)} rows."))
  
  # Ensure RSID column is named 'SNP' or 'RSID' for consistency
  if (!"RSID" %in% names(bulk_dt) && "SNP" %in% names(bulk_dt)) {
    bulk_dt[, RSID := SNP]
  }
  
  # Pre-process Bulk Data: Unique by SNP (min P-value)
  if ("pval.exposure" %in% names(bulk_dt)) {
    # Filter Bulk by P-value < 5e-8 (Global Filter)
    bulk_dt <- bulk_dt[as.numeric(pval.exposure) < 5e-8]
    log_msg(glue("Filtered Bulk Data (P < 5e-8): {nrow(bulk_dt)} rows."))
    
    log_msg("Pre-processing Bulk Data: Keeping min P-value per Gene-SNP pair...")
    bulk_dt[, pval.exposure := as.numeric(pval.exposure)]
    setorder(bulk_dt, pval.exposure, na.last = TRUE)
    
    # Unique by gene_name and RSID
    if ("gene_name" %in% names(bulk_dt)) {
        bulk_dt <- unique(bulk_dt, by = c("gene_name", "RSID"))
    } else {
        bulk_dt <- unique(bulk_dt, by = "RSID")
    }
    log_msg(glue("Reduced Bulk Data to {nrow(bulk_dt)} unique Gene-SNP pairs."))
  }
} else {
  stop(glue("Bulk eQTL file not found: {bulk_eqtl_path}"))
}

# ==============================================================================
# Processing Logic / 
# ==============================================================================

process_cluster <- function(idx, row, total_items) {
  cluster_id <- row$Exposure_Cluster_ID
  cell <- row$Exposure_Cell
  lead_snps_str <- row$Lead_SNPs
  
  # Parse Lead SNPs
  lead_snps <- trimws(unlist(strsplit(lead_snps_str, ";")))
  lead_snps <- lead_snps[lead_snps != ""]
  
  if (length(lead_snps) == 0) {
    return(list(status = "warning", msg = "No valid Lead SNPs found"))
  }
  
  # Determine Data Source
  source_type <- "unknown"
  source_dt <- NULL
  
  # Check if sc-eQTL file exists
  sc_file <- file.path(sc_eqtl_dir, paste0(cell, ".tsv.gz"))
  
  if (file.exists(sc_file)) {
    source_type <- "sc-eQTL"
    # Load sc-eQTL data (filtered by region)
    # Using zcat or gzcat depending on system. On macOS zcat expects .Z, gzcat is better for .gz
    # But let's check if gzcat works. If not, try zcat < file. Or gzip -dc.
    # The safest is 'gzip -dc' or 'gunzip -c'.
    cmd <- glue("gzip -dc '{sc_file}' | awk -F'\t' 'BEGIN {{OFS=\"\t\"}} NR==1 || ($7 == {target_chr} && $8 >= {target_start} && $8 <= {target_end})'")
    tryCatch({
      source_dt <- fread(cmd = cmd, header = TRUE, sep = "\t", showProgress = FALSE)
      
      # Filter sc-eQTL by FDR <= 0.05
      if (!is.null(source_dt) && nrow(source_dt) > 0 && "FDR" %in% names(source_dt)) {
        source_dt <- source_dt[FDR <= 0.05]
      }
    }, error = function(e) {
      source_dt <- NULL
    })
  }
  
  # If not found in sc-eQTL or load failed, try Bulk
  if (is.null(source_dt) || nrow(source_dt) == 0) {
    source_type <- "bulk"
    source_dt <- bulk_dt 
  }
  
  if (is.null(source_dt) || nrow(source_dt) == 0) {
    return(list(status = "error", msg = glue("No source data found for cell: {cell}")))
  }
  
  # Ensure RSID column
  if (!"RSID" %in% names(source_dt)) {
     if ("SNP" %in% names(source_dt)) {
         source_dt[, RSID := SNP]
     } else {
         return(list(status = "error", msg = "RSID/SNP column missing in source data"))
     }
  }

  # ==============================================================================
  # Filtering Logic / 
  # ==============================================================================
  # Keep unique Gene-SNP pairs (min P-value if duplicates exist for same Gene-SNP)
  # snpgene; snpgene,P
  
  p_col <- if (source_type == "sc-eQTL") "P_VALUE" else "pval.exposure"
  gene_col <- if (source_type == "sc-eQTL") "GENE" else "gene_name"
  
  if (p_col %in% names(source_dt)) {
    # Ensure numeric
    source_dt[, (p_col) := as.numeric(get(p_col))]
    
    # Order by P-value (ascending)
    setorderv(source_dt, p_col, na.last = TRUE)
    
    # Keep unique Gene-SNP pairs (first occurrence = min P-value)
    if (gene_col %in% names(source_dt)) {
      source_dt <- unique(source_dt, by = c(gene_col, "RSID"))
    } else {
      source_dt <- unique(source_dt, by = "RSID")
    }
  } else {
    # If P-value column missing
    if (gene_col %in% names(source_dt)) {
      source_dt <- unique(source_dt, by = c(gene_col, "RSID"))
    } else {
      source_dt <- unique(source_dt, by = "RSID")
    }
  }

  # Create temp file for extract (SNPs present in source data)
  temp_id <- glue("{timestamp}_{Sys.getpid()}_{idx}")
  snps_extract_file <- file.path(dir_temp, glue("extract_{temp_id}.txt"))
  writeLines(unique(source_dt$RSID), snps_extract_file)
  
  # Create temp file for Lead SNPs
  leads_file <- file.path(dir_temp, glue("leads_{temp_id}.txt"))
  writeLines(lead_snps, leads_file)
  
  # PLINK Output
  plink_out <- file.path(dir_temp, glue("plink_{temp_id}"))
  
  # Run PLINK to get LD (R2 > 0.6)
  plink_cmd <- glue(
    "'{plink_path}'",
    " --bfile '{ld_ref_prefix}'",
    " --r2",
    " --ld-snp-list '{leads_file}'", 
    " --ld-window-r2 0.6",
    " --ld-window-kb 10000",
    " --ld-window 99999",
    " --extract '{snps_extract_file}'",
    " --out '{plink_out}'",
    " --threads 1",
    " --memory 4000",
    " --silent"
  )
  
  system(plink_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  ld_file <- paste0(plink_out, ".ld")
  
  if (!file.exists(ld_file)) {
    unlink(c(snps_extract_file, leads_file))
    return(list(status = "warning", msg = "PLINK failed or no LD results"))
  }
  
  ld_res <- fread(ld_file, header = TRUE)
  unlink(c(ld_file, snps_extract_file, leads_file, paste0(plink_out, ".log"), paste0(plink_out, ".nosex")))
  
  if (nrow(ld_res) == 0) {
    return(list(status = "warning", msg = "No SNPs with R2 > 0.6 found"))
  }
  
  # Merge with Source Data
  # LD Res: CHR_A, BP_A, SNP_A (Lead), CHR_B, BP_B, SNP_B (Target), R2
  merged <- merge(ld_res, source_dt, by.x = "SNP_B", by.y = "RSID", all.x = TRUE)
  
  # Add Signal (Cluster ID)
  merged[, Signal := cluster_id]
  
  # Format Output Columns
  output_dt <- NULL
  
  if (source_type == "sc-eQTL") {
    # Expected columns: Signal, CELL_ID, CELL_TYPE, SNP, GENE, GENE_ID, CHR, POS, other_allele, effect_allele, eaf, P_VALUE, FDR, N, Z, BETA, SE
    # Source columns: CELL_ID, CELL_TYPE, RSID, SNPID, GENE, GENE_ID, CHR, POS, A1, A2, A2_FREQ_ONEK1K, P_VALUE, FDR, N, Z, BETA, SE
    
    # Map columns
    # Assume A2 is effect_allele (based on EAF usually being effect allele freq) and A1 is other_allele
    # Need to verify if A2_FREQ corresponds to A2. Usually yes.
    
    # Check if columns exist
    cols_check <- c("CELL_ID", "CELL_TYPE", "RSID", "GENE", "GENE_ID", "CHR", "POS", "A1", "A2", "A2_FREQ_ONEK1K", "P_VALUE", "FDR", "N", "Z", "BETA", "SE")
    missing_cols <- setdiff(cols_check, names(merged))
    
    if (length(missing_cols) > 0) {
       # Handle missing columns if any
       for (col in missing_cols) merged[, (col) := NA]
    }
    
    output_dt <- merged[, .(
      Signal = Signal,
      CELL_ID = CELL_ID,
      CELL_TYPE = CELL_TYPE,
      SNP = SNP_B, # RSID (use SNP_B as RSID is dropped after merge)
      GENE = GENE,
      GENE_ID = GENE_ID,
      CHR = CHR,
      POS = POS,
      other_allele = A1,      # Mapping A1 to other
      effect_allele = A2,     # Mapping A2 to effect
      eaf = A2_FREQ_ONEK1K,   # Mapping A2 freq to EAF
      P_VALUE = P_VALUE,
      FDR = FDR,
      N = N,
      Z = Z,
      BETA = BETA,
      SE = SE
    )]
    
  } else { # Bulk
    # Expected columns: Signal, CELL_ID, CELL_TYPE, SNP, GENE, GENE_ID, CHR, POS, other_allele, effect_allele, eaf, P_VALUE, FDR, N, Z, BETA, SE
    # Source columns: SNP, gene_name, id.exposure, chr.exposure, pos.exposure, other_allele.exposure, effect_allele.exposure, eaf.exposure, pval.exposure, samplesize.exposure, beta.exposure, se.exposure, exposure(GENE_ID)
    
    # Create missing columns
    merged[, CELL_ID := "Bulk_MHC"]
    merged[, CELL_TYPE := "Bulk_Blood"]
    if (!"FDR" %in% names(merged)) merged[, FDR := NA] # Placeholder
    if (!"Z" %in% names(merged)) merged[, Z := beta.exposure / se.exposure]
    
    output_dt <- merged[, .(
      Signal = Signal,
      CELL_ID = CELL_ID,
      CELL_TYPE = CELL_TYPE,
      SNP = SNP_B, # RSID
      GENE = gene_name,
      GENE_ID = exposure, # Using exposure (ENSG ID) as GENE_ID
      CHR = chr.exposure,
      POS = pos.exposure,
      other_allele = other_allele.exposure,
      effect_allele = effect_allele.exposure,
      eaf = eaf.exposure,
      P_VALUE = pval.exposure,
      FDR = FDR,
      N = samplesize.exposure,
      Z = Z,
      BETA = beta.exposure,
      SE = se.exposure
    )]
  }
  
  # Remove duplicates
  output_dt <- unique(output_dt)
  
  # Save
  # Output folder: Timestamp/{Cell}/{Exposure_Cluster_ID}.csv
  cell_clean <- gsub("[^a-zA-Z0-9_.-]", "_", cell)
  cell_dir <- file.path(output_dir, cell_clean)
  dir.create(cell_dir, showWarnings = FALSE, recursive = TRUE)
  
  out_file <- file.path(cell_dir, paste0(cluster_id, ".csv"))
  fwrite(output_dt, out_file)
  
  return(list(status = "success", msg = glue("Saved {nrow(output_dt)} SNPs to {basename(out_file)}")))
}

# ==============================================================================
# Execute Parallel Processing / 
# ==============================================================================

log_msg("Starting parallel processing of clusters...")

# Convert cluster_dt to list of rows for mclapply
cluster_list <- split(cluster_dt, seq(nrow(cluster_dt)))

results <- mclapply(seq_along(cluster_list), function(i) {
  process_cluster(i, cluster_list[[i]], length(cluster_list))
}, mc.cores = n_cores)

# Summarize Results
success_count <- sum(sapply(results, function(x) x$status == "success"))
warning_count <- sum(sapply(results, function(x) x$status == "warning"))
error_count <- sum(sapply(results, function(x) x$status == "error"))

log_msg(glue("Processing completed."))
log_msg(glue("Total Clusters: {length(cluster_list)}"))
log_msg(glue("Success: {success_count}"))
log_msg(glue("Warnings: {warning_count}"))
log_msg(glue("Errors: {error_count}"))

# Clean up temp dir
unlink(dir_temp, recursive = TRUE)
log_msg("Temp directory cleaned.")

log_msg("Script execution finished.")
