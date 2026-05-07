# ==============================================================================
# [Script]: 4.1.4.1_batch_intersection_3d_prep_template.R
# [Method]: Multi-outcome Signal Intersection & FDR Control
# [Step]: 4.1.4.1_batch_intersection_3d_prep
#
# [Function]:
# Processes MR analysis results to identify common signals across three outcomes
# (Aging, Rheumatoid Arthritis, and Herpes Zoster), calculates FDR, and prepares data for 3D visualization.
#
# [Usage]: 
#   Rscript 4.1.4.1_batch_intersection_3d_prep_template.R \
#     --mr_dir <path> \
#     --summary <path> \
#     --out_dir <path>
# ==============================================================================

# ==============================================================================
# 1. Setup and Configuration
# ==============================================================================


# Clear environment
rm(list = ls())

# Load libraries
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(tools)
})

# Define command line arguments
option_list <- list(
  make_option(c("--mr_dir"), type="character", default=NULL,
              help="Directory containing MR results / MR"),
  make_option(c("--summary"), type="character", default=NULL,
              help="Exposure Cluster Summary CSV file / CSV"),
  make_option(c("--out_dir"), type="character", default="./Intersection_3D_Data",
              help="Output directory path / ")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$mr_dir) || is.null(opt$summary)) {
  print_help(opt_parser)
  stop("Missing required arguments. / . ", call.=FALSE)
}

# Define Paths
MR_RESULTS_DIR <- opt$mr_dir
SUMMARY_FILE <- opt$summary
OUTPUT_BASE_DIR <- opt$out_dir


# Create timestamped output directory
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTPUT_DIR <- file.path(OUTPUT_BASE_DIR, timestamp)

# Create subdirectories
DATA_DIR <- file.path(OUTPUT_DIR, "data")
STATS_DIR <- file.path(OUTPUT_DIR, "stats")
LOGS_DIR <- file.path(OUTPUT_DIR, "logs")

for (dir in c(OUTPUT_DIR, DATA_DIR, STATS_DIR, LOGS_DIR)) {
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
    message(paste("Created directory:", dir))
  }
}

# Define Outcome Mappings (Summary Name -> File Suffix/Folder Name)
# Note: Summary uses "HZ", "RA", "aging". Files use "herpes_zoster", "Rheumatoid_arthritis", "aging".
OUTCOME_MAP <- list(
  "HZ" = "herpes_zoster",
  "RA" = "Rheumatoid_arthritis",
  "aging" = "aging"
)

# Start logging
log_file <- file.path(LOGS_DIR, "analysis_log.txt")
sink(log_file, split = TRUE)

cat("==============================================================================\n")
cat("Analysis Started at:", as.character(Sys.time()), "\n")
cat("MR Results Directory:", MR_RESULTS_DIR, "\n")
cat("Summary File:", SUMMARY_FILE, "\n")
cat("Output Directory:", OUTPUT_DIR, "\n")
cat("==============================================================================\n\n")

# ==============================================================================
# 2. Read and Filter Summary Data
# ==============================================================================

cat("Reading summary file...\n")
if (!file.exists(SUMMARY_FILE)) {
  stop(paste("Summary file not found:", SUMMARY_FILE))
}

summary_df <- fread(SUMMARY_FILE)

# Function to check if a row has all three outcomes
has_all_outcomes <- function(outcomes_str) {
  if (is.na(outcomes_str) || outcomes_str == "") return(FALSE)
  
  # Check for presence of all three keywords
  has_hz <- str_detect(outcomes_str, "HZ")
  has_ra <- str_detect(outcomes_str, "RA")
  has_aging <- str_detect(outcomes_str, "aging")
  
  return(has_hz && has_ra && has_aging)
}

# Filter for intersection signals
intersection_signals <- summary_df %>%
  filter(sapply(Outcomes, has_all_outcomes))

cat(paste("Found", nrow(intersection_signals), "intersection signals (Aging & RA & HZ).\n"))

if (nrow(intersection_signals) == 0) {
  warning("No intersection signals found. Exiting.")
  sink()
  q(save = "no")
}

# ==============================================================================
# 3. Process Each Signal
# ==============================================================================

# Initialize list to store statistics
stats_list <- list()

process_signal <- function(row_idx) {
  signal_info <- intersection_signals[row_idx, ]
  cluster_id <- signal_info$Exposure_Cluster_ID
  cell_type <- signal_info$Exposure_Cell
  associations_str <- signal_info$Associations
  
  cat(paste0("\nProcessing Signal ", row_idx, "/", nrow(intersection_signals), ": ", cluster_id, " (", cell_type, ")\n"))
  
  # Initialize stats for this signal
  signal_stats <- data.table(
    Signal_ID = cluster_id,
    Cell_Type = cell_type,
    Initial_Intersection_SNPs = 0,
    RA_Matched = 0,
    RA_Swapped = 0,
    RA_Removed = 0,
    HZ_Matched = 0,
    HZ_Swapped = 0,
    HZ_Removed = 0,
    Final_SNPs = 0,
    Status = "Failed"
  )
  
  # Construct paths for the three outcome files
  # Structure: MR_RESULTS_DIR / cell_type / cluster_id / cluster_id-outcome.csv
  
  # Function to read and process a single outcome file
  read_and_process_outcome <- function(outcome_key) {
    outcome_name <- OUTCOME_MAP[[outcome_key]]
    file_path <- file.path(MR_RESULTS_DIR, cell_type, cluster_id, paste0(cluster_id, "-", outcome_name, ".csv"))
    
    if (!file.exists(file_path)) {
      warning(paste("File not found:", file_path))
      return(NULL)
    }
    
    # Read CSV
    df <- fread(file_path)
    
    # Check if empty
    if (nrow(df) == 0) return(NULL)
    
    # Calculate FDR
    if (!"p" %in% names(df)) {
      if ("pval" %in% names(df)) {
        df$p <- df$pval
      } else {
        warning(paste("P-value column 'p' not found in", file_path))
        return(NULL)
      }
    }
    
    df$fdr <- p.adjust(df$p, method = "fdr")
    
    # Filter FDR < 0.05
    df_filtered <- df %>% filter(fdr < 0.05)
    
    if (nrow(df_filtered) == 0) return(NULL)

    # Check for required columns
    required_cols <- c("cell", "gene", "SNP", "b", "se", "p", "fdr", 
                       "effect_allele.exposure", "other_allele.exposure", 
                       "beta.exposure", "eaf.exposure")
    
    missing_cols <- setdiff(required_cols, names(df_filtered))
    if (length(missing_cols) > 0) {
      # Try alternative names for beta/eaf/alleles if needed, but assuming standard TwoSampleMR format for now based on file inspection
      warning(paste("Missing columns in", file_path, ":", paste(missing_cols, collapse=", ")))
      return(NULL)
    }
    
    # Select relevant columns
    df_subset <- df_filtered %>% select(all_of(required_cols))
    
    # Rename value columns with outcome suffix to distinguish them after merge
    # Keep join keys (cell, gene, SNP) as is
    suffix <- paste0("_", outcome_key)
    cols_to_rename <- setdiff(names(df_subset), c("cell", "gene", "SNP"))
    
    # Rename using a loop or vectorized approach
    for (col in cols_to_rename) {
      names(df_subset)[names(df_subset) == col] <- paste0(col, suffix)
    }
    
    return(df_subset)
  }
  
  # Read all three
  df_aging <- read_and_process_outcome("aging")
  df_ra <- read_and_process_outcome("RA")
  df_hz <- read_and_process_outcome("HZ")
  
  if (is.null(df_aging) || is.null(df_ra) || is.null(df_hz)) {
    cat("Skipping signal due to missing files or empty data for one or more outcomes.\n")
    signal_stats$Status <- "Missing Data"
    return(signal_stats)
  }
  
  if (nrow(df_aging) == 0 || nrow(df_ra) == 0 || nrow(df_hz) == 0) {
    cat("Skipping signal due to no significant results (FDR < 0.05) in one or more outcomes.\n")
    signal_stats$Status <- "No Sig Results"
    return(signal_stats)
  }
  
  # Merge (Inner Join)
  # First merge Aging and RA
  merged_df <- inner_join(df_aging, df_ra, by = c("cell", "gene", "SNP"))
  
  # Then merge with HZ
  merged_df <- inner_join(merged_df, df_hz, by = c("cell", "gene", "SNP"))
  
  if (nrow(merged_df) == 0) {
    cat("No common significant SNPs found across all three outcomes.\n")
    signal_stats$Status <- "No Intersection"
    return(signal_stats)
  }

  # Record initial intersection count
  signal_stats$Initial_Intersection_SNPs <- nrow(merged_df)

  
  # --------------------------------------------------------------------------
  # Check Effect Allele Alignment
  # --------------------------------------------------------------------------
  # We use Aging as the reference outcome for effect allele, beta, and eaf.
  # We check if RA and HZ align with Aging.
  
  # Check RA alignment
  align_ra <- merged_df$effect_allele.exposure_aging == merged_df$effect_allele.exposure_RA
  swapped_ra <- merged_df$effect_allele.exposure_aging == merged_df$other_allele.exposure_RA
  problematic_ra <- !align_ra & !swapped_ra
  
  signal_stats$RA_Matched <- sum(align_ra)
  signal_stats$RA_Swapped <- sum(swapped_ra)
  signal_stats$RA_Removed <- sum(problematic_ra)
  
  if (any(problematic_ra)) {
    cat(paste("Critical Warning:", sum(problematic_ra), "SNPs have alleles that do not match or swap in RA. Removing them.\n"))
    merged_df <- merged_df[!problematic_ra, ]
  }
  
  if (any(!align_ra & !problematic_ra)) {
      cat(paste("Warning: Found", sum(!align_ra & !problematic_ra), "SNPs with swapped alleles between Aging and RA.\n"))
  }

  
  # Check HZ alignment
  align_hz <- merged_df$effect_allele.exposure_aging == merged_df$effect_allele.exposure_HZ
  swapped_hz <- merged_df$effect_allele.exposure_aging == merged_df$other_allele.exposure_HZ
  problematic_hz <- !align_hz & !swapped_hz

  signal_stats$HZ_Matched <- sum(align_hz)
  signal_stats$HZ_Swapped <- sum(swapped_hz)
  signal_stats$HZ_Removed <- sum(problematic_hz)

  if (any(problematic_hz)) {
    cat(paste("Critical Warning:", sum(problematic_hz), "SNPs have alleles that do not match or swap in HZ. Removing them.\n"))
    merged_df <- merged_df[!problematic_hz, ]
  }

  if (any(!align_hz & !problematic_hz)) {
      cat(paste("Warning: Found", sum(!align_hz & !problematic_hz), "SNPs with swapped alleles between Aging and HZ.\n"))
  }
  
  if (nrow(merged_df) == 0) {
    cat("No valid SNPs left after allele alignment check.\n")
    signal_stats$Status <- "Failed - Allele Mismatch"
    return(signal_stats)
  }
  
  # Update final count
  signal_stats$Final_SNPs <- nrow(merged_df)
  signal_stats$Status <- "Success"
  
  # Construct Final DataFrame
  # We use the Aging columns as the representative ones for the SNP info
  merged_df$effect_allele <- merged_df$effect_allele.exposure_aging
  merged_df$other_allele <- merged_df$other_allele.exposure_aging
  merged_df$beta_exposure <- merged_df$beta.exposure_aging
  merged_df$eaf_exposure <- merged_df$eaf.exposure_aging
  
  # Calculate Z-scores (b/se)
  merged_df$z_aging <- merged_df$b_aging / merged_df$se_aging
  merged_df$z_RA <- merged_df$b_RA / merged_df$se_RA
  merged_df$z_HZ <- merged_df$b_HZ / merged_df$se_HZ
  
  # Select final columns
  # Format: cell, gene, SNP, effect_allele, other_allele, beta_exposure, eaf_exposure, 
  #         b_aging, se_aging, z_aging, p_aging, fdr_aging, ...
  
  final_cols <- c("cell", "gene", "SNP", 
                  "effect_allele", "other_allele", "beta_exposure", "eaf_exposure",
                  "b_aging", "se_aging", "z_aging", "p_aging", "fdr_aging",
                  "b_RA", "se_RA", "z_RA", "p_RA", "fdr_RA",
                  "b_HZ", "se_HZ", "z_HZ", "p_HZ", "fdr_HZ")
  
  # Ensure we have these columns (created above)
  # We need to subset merged_df but keep the newly created columns and the b/p/fdr columns
  
  # Let's clean up merged_df to only have what we want
  # We can just select directly since we added the non-suffixed columns
  
  merged_df <- merged_df %>% select(all_of(final_cols))


  
  # Annotate Associations
  # Parse Associations string: "bin-BTN3A2-rs9358932-HZ; bin-BTN3A2-rs9379859-HZ"
  # We want to check if the current row's cell-gene-snp is in this list.
  
  # Split by "; "
  assoc_list <- str_split(associations_str, "; ")[[1]]
  
  # For each association string, remove the last part (outcome) to get the key
  # Example: "bin-BTN3A2-rs9358932-HZ" -> "bin-BTN3A2-rs9358932"
  # Assuming the format is always cell-gene-snp-outcome, and outcome is the last part separated by hyphen.
  # However, gene names might have hyphens.
  # Safer: "outcome" is strictly HZ, RA, aging.
  # Let's remove the suffix "-HZ", "-RA", "-aging".
  
  clean_assoc_keys <- gsub("-(HZ|RA|aging)$", "", assoc_list)
  
  # Create a key for the current row
    # Use cell_type from summary (Exposure_Cell) to match the Associations string format
    # The raw data 'cell' column (e.g., "Bulk_MHC") might differ from the summary 'cell' name (e.g., "bulk")
    merged_df$row_key <- paste(cell_type, merged_df$gene, merged_df$SNP, sep = "-")
  
  # Check match
  merged_df$is_association_signal <- merged_df$row_key %in% clean_assoc_keys
  
  # Drop the temporary key column
  merged_df$row_key <- NULL
  
  # Standardize cell column
  # Ensure output cell column matches the summary file cell type (e.g., change "Bulk_MHC" to "bulk")
  merged_df$cell <- cell_type
  
  # Save result
  output_filename <- paste0(cluster_id, "_intersection_3d_data.csv")
  output_path <- file.path(DATA_DIR, output_filename)
  fwrite(merged_df, output_path)
  
  cat(paste("Saved merged data to:", output_path, "\n"))
  cat(paste("Rows:", nrow(merged_df), "\n"))
  
  return(signal_stats)
}

# Run loop
stats_results <- list()

for (i in 1:nrow(intersection_signals)) {
  tryCatch({
    stats_results[[i]] <- process_signal(i)
  }, error = function(e) {
    cat(paste("Error processing signal", i, ":", e$message, "\n"))
    # Try to capture basic info even on error if possible, but for now just skip or add a failed placeholder
    # Ideally, process_signal should handle errors internally and return a failed status row.
    # Since we have tryCatch inside the loop, we can create a dummy failed row here if needed.
  })
}

# Combine stats
final_stats <- rbindlist(stats_results, fill = TRUE)

# Save stats
stats_filename <- paste0("Intersection_Processing_Stats_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".csv")
stats_path <- file.path(STATS_DIR, stats_filename)
fwrite(final_stats, stats_path)
cat(paste("Saved processing stats to:", stats_path, "\n"))

# ==============================================================================
# 4. Finalize
# ==============================================================================

cat("\n==============================================================================\n")
cat("Analysis Completed at:", as.character(Sys.time()), "\n")
cat("Output saved to:", OUTPUT_DIR, "\n")
cat("==============================================================================\n")

sink()
