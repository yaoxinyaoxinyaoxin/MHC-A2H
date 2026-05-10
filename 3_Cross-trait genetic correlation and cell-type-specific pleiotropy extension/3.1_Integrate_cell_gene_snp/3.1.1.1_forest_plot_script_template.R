# ==============================================================================
# [Script]: 0.2_forest_plot_script.R
# [Method]: cell-gene-snp
# [Step]: +
# 
# [Function]:
# Execute the '+' step within the 'cell-gene-snp' analytical framework.
# 
# [Parameters / ]:
# Standard predefined thresholds.
# 
# [Steps / ]:
#   1. Data loading and initialization / 
#   2. Core analytical execution / 
#   3. Results formatting and output / 
# ==============================================================================

# R script for generating Forest Plot with Table (Hierarchical / )
# Description: 
#   This script reads multiple CSV files containing MR results from a specified directory.
#   It standardizes column names and merges data from Bulk and Single-Cell sources.
#   It generates hierarchical forest plots:
#     - Level 1 (Group): Gene Name
#     - Level 2 (Subgroup): Cell Type (indented)
#   It now supports enhanced directory structure, specific outcome pair analysis, single outcome plots, and Gene-Level shared outcome plots.
#   It also generates comprehensive statistics (including Gene-Cell level), GO/KEGG Enrichment Analysis.
#   Saves the results as PDF/PNG and copies the script to the output directory.
#   Note: Beta correction/alignment features have been removed to preserve original directionality.

# Steps / :
#   1. Load required libraries 
#   2. Define input directory and output directories 
#   3. Scan and load data files 
#   4. Standardize and merge data 
#   5. Define hierarchical plotting and enrichment functions 
#   6. Process data (Sort by Gene -> Cell) (: ->)
#   7. Generate plot 
#   8. Generate Single Outcome Plots 
#   9. Generate Gene-Level Shared Outcome Plots 
#   10. Generate Specific Outcome Pairs (Same Cell + Same SNP) 
#   11. Generate Statistics (Statistics on Genes, SNPs, and Overlaps) 
#   12. Perform Enrichment Analysis 
#   13. Generate log
#   14. Save outputs and backup the script 

# Load required libraries / 
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(forestploter)
  library(grid)
  library(scales)
  library(EnsDb.Hsapiens.v75)
  library(AnnotationDbi)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggrepel)
  library(ggsci)
  library(patchwork)
  library(cowplot)
  library(future)
  library(future.apply)
})

# ==============================================================================
# 0. Command Line Arguments / 
# ==============================================================================
option_list <- list(
  make_option(c("--input_main"), type="character", default=NULL,
              help="Directory containing main analysis CSV files / CSV"),
  make_option(c("--input_ra"), type="character", default=NULL,
              help="Directory containing RA analysis CSV files / RACSV"),
  make_option(c("--input_val"), type="character", default=NULL,
              help="Directory containing validation CSV files / CSV"),
  make_option(c("--out_dir"), type="character", default="./ForestPlot_Results",
              help="Output directory path / ")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_main)) {
  print_help(opt_parser)
  stop("Main input directory must be provided. / . ")
}

# ==============================================================================
# 1. Settings and Paths / 
# ==============================================================================

# Input directories / 
input_dir <- opt$input_main
input_dir_ra <- opt$input_ra
validation_dir <- opt$input_val
output_base <- opt$out_dir

# Visualization Parameters / 
col_padding_spaces <- 1 
forest_plot_width_chars <- 60

# Output directory with Dynamic Timestamp / 
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- file.path(output_base, paste0("ForestPlot_Analysis__", timestamp))

# Subdirectories / 
dir_single_outcome <- file.path(output_dir, "1_Single_Outcome")
dir_joint_outcome <- file.path(output_dir, "2_Joint_Outcome")
dir_gene_shared <- file.path(output_dir, "3_Intersection_Gene")
dir_cell_gene_shared <- file.path(output_dir, "4_Intersection_Cell_Gene")
dir_cell_gene_snp_shared <- file.path(output_dir, "5_Intersection_Cell_Gene_SNP")
dir_statistics <- file.path(output_dir, "6_Statistics")
dir_enrichment <- file.path(output_dir, "7_Enrichment_Analysis")
dir_added_cell_gene_snp <- file.path(output_dir, "8_Added_Cell_Gene_SNP_Forest_Correlation")
dir_validation <- file.path(output_dir, "Validation_Analysis") 
dir_logs <- file.path(output_dir, "Logs")
dir_scripts <- file.path(output_dir, "Scripts")

# Create directories / 
for (d in c(dir_single_outcome, dir_joint_outcome, dir_gene_shared, dir_cell_gene_shared, dir_cell_gene_snp_shared, dir_statistics, dir_enrichment, dir_added_cell_gene_snp, dir_validation, dir_logs, dir_scripts)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE)
}

# Setup Logging / 
log_file <- file.path(dir_logs, "analysis_log.txt")
sink(log_file, split = TRUE) # Redirect output to file and console / 
options(device = function(...) grDevices::pdf(file = file.path(dir_logs, "Rplots.pdf"), ...))

cat("Analysis started at:", as.character(Sys.time()), "\n")
cat("Input directory:", input_dir, "\n")
cat("Output directory:", output_dir, "\n\n")

# ==============================================================================
# 2. Data Loading and Processing / 
# ==============================================================================

# List all CSV files in the input directory / CSV
files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)

if (length(files) == 0) {
  stop("No CSV files found in the input directory! / CSV！")
}

cat("Found", length(files), "CSV files. / ", length(files), "CSV. \n")

# Function to determine source type from filename / 
get_source_type <- function(filename) {
  if (grepl("OneK1K", filename, ignore.case = TRUE)) {
    return("Single-Cell")
  } else {
    return("Bulk")
  }
}

# Function to standardize columns / 
standardize_cols <- function(df, type, filename) {
  cat("Processing file:", basename(filename), "(Type:", type, ")\n")
  
  # Check Gene Name / 
  if (!"gene_name" %in% colnames(df)) {
    if ("gene_id" %in% colnames(df)) {
      cat("  Mapping Gene IDs to Symbols...\n")
      tryCatch({
        df$gene_name <- mapIds(EnsDb.Hsapiens.v75, keys = df$gene_id, column = "SYMBOL", keytype = "GENEID", multiVals = "first")
        df$gene_name[is.na(df$gene_name)] <- df$gene_id[is.na(df$gene_name)]
      }, error = function(e) {
        cat("  Error mapping gene IDs: ", e$message, "\n")
        df$gene_name <- df$gene_id 
      })
    } else if ("exposure" %in% colnames(df)) {
      df$gene_name <- df$exposure
    } else {
      warning("  No 'gene_name', 'gene_id', or 'exposure' column found in ", basename(filename))
      df$gene_name <- NA
    }
  }
  
  # Ensure numeric columns / 
  numeric_cols <- c("b", "se", "pval", "or", "or_lci95", "or_uci95", "nsnp", "fdr_pval", "F_statistic", "beta.exposure", "eaf.exposure")
  for (col in numeric_cols) {
    if (col %in% colnames(df)) {
      df[[col]] <- as.numeric(df[[col]])
    } else {
      df[[col]] <- NA
    }
  }

  # Ensure exposure columns exist
  exposure_cols <- c("effect_allele.exposure", "other_allele.exposure")
  for (col in exposure_cols) {
    if (!col %in% colnames(df)) {
      df[[col]] <- NA
    }
  }
  
  # Map 'b' to 'beta.exposure' if missing / ,  'b'  'beta.exposure'
  if (all(is.na(df$beta.exposure)) && "b" %in% colnames(df) && !all(is.na(df$b))) {
      df$beta.exposure <- df$b
  } else {
      # Row-wise fill
      df$beta.exposure <- ifelse(is.na(df$beta.exposure), df$b, df$beta.exposure)
  }

  # Standardize 'lead_snp' /  'lead_snp'
  if (!"lead_snp" %in% colnames(df)) {
    if ("SNP" %in% colnames(df)) {
      df$lead_snp <- df$SNP
    } else if ("exposure" %in% colnames(df)) {
      df$lead_snp <- df$exposure
    } else {
      df$lead_snp <- NA
    }
  }
  
  # Ensure 'cell' column exists /  'cell' 
  if (!"cell" %in% colnames(df)) {
    if (type == "Bulk") {
      df$cell <- "Bulk"
    } else {
      warning("  'cell' column missing in Single-Cell file: ", basename(filename))
      df$cell <- "Unknown"
    }
  }

  # Ensure 'outcome' column exists /  'outcome' 
  if (!"outcome" %in% colnames(df)) {
    warning("  'outcome' column missing in file: ", basename(filename))
    df$outcome <- NA
  }
  
  # Method cleaning / 
  if ("method" %in% colnames(df)) {
    df$method <- gsub("Inverse variance weighted", "IVW", df$method)
  } else {
    df$method <- "Wald ratio" 
  }
  
  # Ensure 'steiger_dir' is character /  'steiger_dir' 
  if ("steiger_dir" %in% colnames(df)) {
    df$steiger_dir <- as.character(df$steiger_dir)
  }
  
  # Add source type / 
  df$source_type <- type
  
  return(df)
}

# Loop through files and load data / 
data_list <- list()

# Load Main Data / 
cat("Loading Main Analysis Data...\n")
files <- list.files(input_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(files) == 0) {
  stop("No CSV files found in the main input directory!")
}

for (f in files) {
  temp_df <- read.csv(f, stringsAsFactors = FALSE)
  sType <- get_source_type(basename(f))
  temp_df <- standardize_cols(temp_df, sType, f)
  temp_df$Analysis_Group <- "Main"
  data_list[[length(data_list) + 1]] <- temp_df
}

# Load RA Data / RA
cat("Loading RA Analysis Data...\n")
files_ra <- list.files(input_dir_ra, pattern = "\\.csv$", full.names = TRUE)
if (length(files_ra) > 0) {
  for (f in files_ra) {
    temp_df <- read.csv(f, stringsAsFactors = FALSE)
    sType <- get_source_type(basename(f))
    temp_df <- standardize_cols(temp_df, sType, f)
    temp_df$Analysis_Group <- "Main" # RA is also part of Main analysis now
    data_list[[length(data_list) + 1]] <- temp_df
  }
} else {
  warning("No RA data found in: ", input_dir_ra)
}

# Load Validation Data / 
cat("Loading Validation Analysis Data...\n")
files_val <- list.files(validation_dir, pattern = "\\.csv$", full.names = TRUE)
if (length(files_val) > 0) {
  for (f in files_val) {
    temp_df <- read.csv(f, stringsAsFactors = FALSE)
    sType <- get_source_type(basename(f))
    
    # Pre-process Outcome Names for Validation
    if ("outcome" %in% colnames(temp_df)) {
        # Rename FinnGen HZ to distinguish
        temp_df$outcome[temp_df$outcome == "Herpes Zoster (Finngen)"] <- "Herpes Zoster (FinnGen)"
        # Ensure standard HZ name matches if needed (though standardization happens in standardize_cols)
    }
    
    temp_df <- standardize_cols(temp_df, sType, f)
    temp_df$Analysis_Group <- "Validation"
    data_list[[length(data_list) + 1]] <- temp_df
  }
} else {
  warning("No validation data found in: ", validation_dir)
}

# Bind all dataframes / 
cat("Merging all datasets... / ...\n")
data_combined <- bind_rows(data_list)

# Handle Colocalization Columns / 
coloc_cols <- c("SuSiE_block", "SuSiE_100kb", "SuSiE_50kb", "SuSiE_10kb")

for (col in coloc_cols) {
  if (!col %in% colnames(data_combined)) {
    data_combined[[col]] <- NA
  }
}

# Select relevant columns / 
base_cols <- c("gene_name", "cell", "source_type", "method", "nsnp", "b", "se", "pval", "or", "or_lci95", "or_uci95", "fdr_pval", "lead_snp", "F_statistic", "gene_pos", "outcome", "effect_allele.exposure", "other_allele.exposure", "eaf.exposure", "beta.exposure", "Analysis_Group")
keep_cols <- c(base_cols, coloc_cols)

data_all <- data_combined %>% dplyr::select(any_of(keep_cols))

# Global Renaming of Outcomes as requested / 
data_all$outcome[data_all$outcome == "Herpes Zoster"] <- "HZ (UKB)"
data_all$outcome[data_all$outcome == "HZ"] <- "HZ (UKB)"
data_all$outcome[data_all$outcome == "Rheumatoid arthritis"] <- "RA"
data_all$outcome[data_all$outcome == "Herpes Zoster (FinnGen)"] <- "HZ (FinnGen)"
data_all$outcome[data_all$outcome == "aging"] <- "mvAge"
# Also ensure consistency for mvAge if needed, but keeping as 'mvAge' based on existing script usage


# Sort Data: Position -> Gene -> Main (Bulk/SC) -> Validation
data_all$sort_order <- case_when(
  data_all$Analysis_Group == "Validation" ~ 3,
  data_all$source_type == "Bulk" ~ 1,
  TRUE ~ 2
)

data_all <- data_all %>%
  mutate(
    chr_str = str_extract(gene_pos, "^[^:]+"),
    start_pos = as.numeric(str_extract(gene_pos, "(?<=:)\\d+")),
    chr_num = case_when(
      str_detect(chr_str, "chr[0-9]+$") ~ as.numeric(str_remove(chr_str, "chr")),
      chr_str == "chrX" ~ 23,
      chr_str == "chrY" ~ 24,
      chr_str == "chrM" | chr_str == "chrMT" ~ 25,
      TRUE ~ 99
    ),
    outcome_priority_main = case_when(
      outcome == "mvAge" ~ 1,
      outcome == "HZ (UKB)" ~ 2,
      outcome == "RA" ~ 3,
      outcome == "HZ (FinnGen)" ~ 4,
      TRUE ~ 99
    )
  ) %>%
  arrange(chr_num, start_pos, gene_name, sort_order, cell, lead_snp, outcome_priority_main, outcome)

cat("Total records loaded:", nrow(data_all), "\n")

# Separate Main and Validation Data
data <- data_all %>% filter(Analysis_Group == "Main")
data_validation <- data_all %>% filter(Analysis_Group == "Validation")

cat("Main Data Records:", nrow(data), "\n")
cat("Validation Data Records:", nrow(data_validation), "\n")

# ==============================================================================
# 3. Helper Functions / 
# ==============================================================================

# Helper functions for formatting
format_f <- function(f) {
  if (is.na(f) || f == "") return("-")
  sprintf("%.1f", as.numeric(f))
}

format_pval <- function(p) {
  if (is.na(p) || p == "") return("-")
  format(as.numeric(p), scientific = TRUE, digits = 3)
}

format_est <- function(or, lci, uci) {
  if (is.na(or) || is.na(lci) || is.na(uci)) return("-")
  sprintf("%.2f (%.2f-%.2f)", or, lci, uci)
}

format_coloc_pct <- function(val) {
  if (is.na(val) || val == "" || val == "-") return("-")
  num_val <- suppressWarnings(as.numeric(val))
  if (is.na(num_val)) return("-")
  return(paste0(round(num_val * 100, 1), "%"))
}

# Function to generate hierarchical plot data / 
prepare_hierarchical_data <- function(df_sorted) {
  plot_data_list <- list()
  
  genes <- unique(df_sorted$gene_name)
  
  for (g in genes) {
    if (is.na(g)) next
    
    g_pos <- df_sorted$gene_pos[df_sorted$gene_name == g][1]
    if (is.na(g_pos)) g_pos <- ""
    
    display_name_str <- if (g_pos != "") paste0(g, " (", g_pos, ")") else g
    
    # 1. Create Header Row (Gene Name)
    header_row <- data.frame(
      Display_Name = display_name_str,
      Outcome = "",
      Lead_SNP = "",
      Effect_Allele = "",
      Other_Allele = "",
      EAF = "",
      Beta_Exposure = "",
      F_stat = "",
      Method = "",
      N_SNP = "",
      `OR(95%CI)` = "",
      P_FDR = "",
      or = NA,
      low = NA,
      high = NA,
      is_summary = TRUE,
      is_bulk = FALSE,
      cell = "Summary", # Added for coloring
      stringsAsFactors = FALSE
    )
    
    for (col in coloc_cols) header_row[[col]] <- ""
    
    plot_data_list[[length(plot_data_list) + 1]] <- header_row
    
    # 2. Get subset for this gene
    gene_subset <- df_sorted[df_sorted$gene_name == g, ]
    
    for (i in 1:nrow(gene_subset)) {
      row <- gene_subset[i, ]
      
      display_name <- paste0("        ", row$cell)
      
      data_row <- data.frame(
        Display_Name = display_name,
        Outcome = ifelse(is.na(row$outcome), "", row$outcome),
        Lead_SNP = ifelse(is.na(row$lead_snp), "", row$lead_snp),
        Effect_Allele = ifelse(is.na(row$effect_allele.exposure), "-", row$effect_allele.exposure),
        Other_Allele = ifelse(is.na(row$other_allele.exposure), "-", row$other_allele.exposure),
        EAF = ifelse(is.na(row$eaf.exposure), "-", sprintf("%.3f", as.numeric(row$eaf.exposure))),
        Beta_Exposure = ifelse(is.na(row$beta.exposure), "-", sprintf("%.3f", as.numeric(row$beta.exposure))),
        F_stat = format_f(if ("F_statistic" %in% names(row)) row$F_statistic else NA),
        Method = row$method,
        N_SNP = row$nsnp,
        `OR(95%CI)` = format_est(row$or, row$or_lci95, row$or_uci95),
        P_FDR = format_pval(row$fdr_pval),
        or = row$or,
        low = row$or_lci95,
        high = row$or_uci95,
        is_summary = FALSE,
        is_bulk = (row$source_type == "Bulk"),
        cell = row$cell, # Added for coloring
        stringsAsFactors = FALSE
      )
      
      for (col in coloc_cols) {
        val <- row[[col]]
        data_row[[col]] <- format_coloc_pct(val)
      }
      
      plot_data_list[[length(plot_data_list) + 1]] <- data_row
    }
  }
  
  plot_df <- do.call(rbind, plot_data_list)
  colnames(plot_df)[colnames(plot_df) == "OR.95.CI."] <- "OR(95%CI)"
  return(plot_df)
}

# Function to generate and save hierarchical forest plot / 
# Updated to accept target directory and hide beta/allele cols / beta/allele
generate_hierarchical_plot <- function(plot_df, output_prefix, target_dir, show_beta_allele = TRUE, color_cell_only = FALSE) {
  
  cat(paste0("Generating forest plot: ", output_prefix, "... / : ", output_prefix, "...\n"))
  
  display_coloc_cols <- c("SuSiE (Block)", "SuSiE (100kb)", "SuSiE (50kb)", "SuSiE (10kb)")
  
  names(plot_df)[names(plot_df) == "SuSiE_block"] <- "SuSiE (Block)"
  names(plot_df)[names(plot_df) == "SuSiE_100kb"] <- "SuSiE (100kb)"
  names(plot_df)[names(plot_df) == "SuSiE_50kb"] <- "SuSiE (50kb)"
  names(plot_df)[names(plot_df) == "SuSiE_10kb"] <- "SuSiE (10kb)"
  
  # Ensure character padding columns exist
  forest_plot_width_chars <- 30 # Default or global
  plot_df$` ` <- paste(rep(" ", forest_plot_width_chars), collapse = " ")
  plot_df$`  ` <- "   " 
  
  # Select columns based on show_beta_allele
  base_cols_to_show <- c("Display_Name", "Outcome", "Lead_SNP")
  if (show_beta_allele) {
    base_cols_to_show <- c(base_cols_to_show, "Effect_Allele", "Other_Allele", "EAF", "Beta_Exposure")
  }
  base_cols_to_show <- c(base_cols_to_show, "F_stat", "Method", "P_FDR", "OR(95%CI)", " ", "  ")
  
  # Check if columns exist before selecting
  existing_cols <- intersect(c(base_cols_to_show, display_coloc_cols), colnames(plot_df))
  final_df <- plot_df[, existing_cols]
  
  # Rename columns
  colnames(final_df)[colnames(final_df) == "Display_Name"] <- "Gene (Location)/Cell Type"
  colnames(final_df)[colnames(final_df) == "Lead_SNP"] <- "Lead SNP" 
  
  if (show_beta_allele) {
    colnames(final_df)[colnames(final_df) == "Effect_Allele"] <- "EA" # Abbreviation
    colnames(final_df)[colnames(final_df) == "Other_Allele"] <- "OA" # Abbreviation
    colnames(final_df)[colnames(final_df) == "Beta_Exposure"] <- "Beta" # Abbreviation
    # EAF is already EAF
  }
  
  colnames(final_df)[colnames(final_df) == "P_FDR"] <- "FDR"
  colnames(final_df)[colnames(final_df) == "F_stat"] <- "F-statistic"
  
  # Padding logic (using global col_padding_spaces if available, else default)
  col_padding_spaces <- if(exists("col_padding_spaces")) col_padding_spaces else 2
  if (col_padding_spaces > 0) {
    pad_str <- strrep(" ", col_padding_spaces)
    cols_to_pad <- setdiff(colnames(final_df), c(" ", "  ", "OR(95%CI)"))
    
    for (col in cols_to_pad) {
      final_df[[col]] <- as.character(final_df[[col]])
      final_df[[col]] <- paste0(pad_str, final_df[[col]], pad_str)
    }
    
    if ("OR(95%CI)" %in% colnames(final_df)) {
      col <- "OR(95%CI)"
      final_df[[col]] <- as.character(final_df[[col]])
      final_df[[col]] <- paste0(pad_str, final_df[[col]]) 
    }
  }
  
  # Row fills - Simple style from backup script (White/Grey/LightBlue)
  row_fills <- ifelse(plot_df$is_summary, "#D9D9D9", 
                      ifelse(plot_df$is_bulk, "#E6F3FF", "white")) 
  
  # Adjust ci_column index
  # If show_beta_allele is TRUE:
  # 1: Display_Name, 2: Outcome, 3: Lead_SNP, 4: EA, 5: OA, 6: EAF, 7: Beta, 8: F_stat, 9: Method, 10: FDR, 11: OR, 12: " ", 13: "  "
  # CI column is 12 (" ")
  # If show_beta_allele is FALSE:
  # 1: Display_Name, 2: Outcome, 3: Lead_SNP, 4: F_stat, 5: Method, 6: FDR, 7: OR, 8: " ", 9: "  "
  # CI column is 8 (" ")
  
  ci_col_idx <- if (show_beta_allele) 12 else 8

  tm <- forest_theme(
    base_size = 10,
    ci_pch = 16,
    ci_col = "#377eb8",
    refline_col = "red",
    refline_lty = "dashed",
    vertline_lty = "dashed",
    vertline_col = "grey90",
    summary_col = "#000000",
    core = list(bg_params=list(fill = row_fills),
                fg_params=list(hjust = 0, vjust = 0.5)), 
    colhead = list(fg_params=list(hjust = 0, vjust = 0.5, fontface = "bold")) 
  )
  
  n_rows <- nrow(final_df)
  height_in <- max(4, n_rows * 0.25 + 1.5) 
  width_in <- 22.5 # Reduced width by half (45 -> 22.5) / 
  
  out_pdf <- file.path(target_dir, paste0(output_prefix, ".pdf"))
  out_png <- file.path(target_dir, paste0(output_prefix, ".png"))
  out_csv <- file.path(target_dir, paste0(output_prefix, ".csv"))
  
  # Save data to CSV
  write.csv(plot_df, out_csv, row.names = FALSE, na = "")
  
  # X-axis limits
  all_vals <- c(plot_df$low, plot_df$high, 1)
  all_vals <- all_vals[!is.na(all_vals) & all_vals > 0]
  
  if (length(all_vals) > 0) {
    min_val <- min(all_vals)
    max_val <- max(all_vals)
    max_log_dist <- max(abs(log10(min_val)), abs(log10(max_val)))
    padded_log_dist <- max_log_dist * 1.1
    final_lower <- 10^(-padded_log_dist)
    final_upper <- 10^(padded_log_dist)
    xlim_range <- c(final_lower, final_upper)
  } else {
    xlim_range <- c(0.1, 10) 
  }
  
  possible_ticks <- c(0.01, 0.05, 0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100)
  ticks_vec <- possible_ticks[possible_ticks >= xlim_range[1] & possible_ticks <= xlim_range[2]]
  if(!1 %in% ticks_vec) ticks_vec <- sort(unique(c(ticks_vec, 1)))
  
  align_vec <- rep("l", ncol(final_df))
  align_vec[which(colnames(final_df) == " ")] <- "c"
  align_vec[which(colnames(final_df) == "  ")] <- "c"
  
  p <- forest(
    data = final_df,
    est = plot_df$or,
    lower = plot_df$low,
    upper = plot_df$high,
    sizes = 0.25,
    ci_column = ci_col_idx, 
    ref_line = 1,
    vertline = ticks_vec,
    arrow_lab = c("Protective", "Risk"),
    xlim = xlim_range,
    ticks_at = ticks_vec,
    x_trans = "log10",
    theme = tm,
    is_summary = plot_df$is_summary, 
    align = align_vec
  )
  
  # Headers
  header_cols <- if (show_beta_allele) 1:11 else 1:7
  coloc_cols_indices <- if (show_beta_allele) 14:17 else 10:13
  
  p <- insert_text(p, 
                   text = "Mendelian Randomization Analysis", 
                   part = "header", 
                   col = header_cols, 
                   gp = gpar(fontface = "bold", cex = 1.1, hjust = 0.5))
  
  p <- add_text(p, 
                text = "Colocalization Analysis (PP.H4)", 
                part = "header", 
                row = 1, 
                col = coloc_cols_indices, 
                gp = gpar(fontface = "bold", cex = 1.1, hjust = 0.5))
  
  p <- add_border(p, part = "header", row = 1, col = header_cols, where = "bottom", gp = gpar(lwd = 1.5))
  p <- add_border(p, part = "header", row = 1, col = coloc_cols_indices, where = "bottom", gp = gpar(lwd = 1.5))
  
  pdf(out_pdf, width = width_in, height = height_in)
  print(p)
  dev.off()
  
  png(out_png, width = width_in * 300, height = height_in * 300, res = 300)
  print(p)
  dev.off()
  
  cat("Saved plots and data to: / : \n")
  cat("  PDF:", out_pdf, "\n")
  cat("  PNG:", out_png, "\n")
  cat("  CSV:", out_csv, "\n")
  
  return(p)
}

# Function to perform Enrichment Analysis (GO & KEGG) / 
perform_enrichment <- function(gene_list, output_subdir, title_prefix) {
  
  if (length(gene_list) == 0) {
    cat(paste0("  No genes provided for enrichment: ", title_prefix, "\n"))
    return(NULL)
  }
  
  if (!dir.exists(output_subdir)) dir.create(output_subdir, recursive = TRUE)
  
  cat(paste0("  Performing Enrichment Analysis for: ", title_prefix, " (", length(gene_list), " genes)...\n"))
  
  # Save Gene List for Traceability / 
  gene_df <- data.frame(Gene = gene_list, stringsAsFactors = FALSE)
  write.csv(gene_df, file.path(output_subdir, paste0("Gene_List_", title_prefix, ".csv")), row.names = FALSE)
  
  # 1. Map Gene Symbols to Entrez IDs
  entrez_ids <- mapIds(org.Hs.eg.db, keys = gene_list, column = "ENTREZID", keytype = "SYMBOL", multiVals = "first")
  entrez_ids <- entrez_ids[!is.na(entrez_ids)]
  
  if (length(entrez_ids) == 0) {
    cat("    No valid Entrez IDs mapped. Skipping enrichment.\n")
    return(NULL)
  }
  
  # 2. GO Enrichment (BP, CC, MF)
  for (ont in c("BP", "CC", "MF")) {
    cat(paste0("    Running GO-", ont, "...\n"))
    tryCatch({
      ego <- enrichGO(gene = entrez_ids, OrgDb = org.Hs.eg.db, ont = ont, pAdjustMethod = "BH", pvalueCutoff = 1, qvalueCutoff = 1, readable = TRUE)
      
      if (!is.null(ego)) {
        # Save raw results
        write.csv(as.data.frame(ego), file.path(output_subdir, paste0("GO_", ont, "_", title_prefix, "_All.csv")), row.names = FALSE)
        
        # Filter for significant results (p.adjust < 0.05)
        ego_sig <- filter(ego, p.adjust < 0.05)
        
        if (nrow(ego_sig) > 0) {
          write.csv(as.data.frame(ego_sig), file.path(output_subdir, paste0("GO_", ont, "_", title_prefix, "_Significant.csv")), row.names = FALSE)
          
          # Plots
          p_dot <- dotplot(ego_sig, showCategory = 20) + 
            ggtitle(paste0("GO-", ont, " Enrichment: ", title_prefix)) +
            scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 60)) +
            theme(axis.text.y = element_text(size = 10))
            
          ggsave(file.path(output_subdir, paste0("Dotplot_GO_", ont, "_", title_prefix, ".pdf")), p_dot, width = 10, height = 8)
          
          p_bar <- barplot(ego_sig, showCategory = 20) + 
            ggtitle(paste0("GO-", ont, " Enrichment: ", title_prefix)) +
            scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 60)) +
            theme(axis.text.y = element_text(size = 10))
            
          ggsave(file.path(output_subdir, paste0("Barplot_GO_", ont, "_", title_prefix, ".pdf")), p_bar, width = 10, height = 8)
        }
      }
    }, error = function(e) {
      cat(paste0("    Error in GO-", ont, ": ", e$message, "\n"))
    })
  }
  
  # 3. KEGG Enrichment
  cat("    Running KEGG...\n")
  tryCatch({
    kk <- enrichKEGG(gene = entrez_ids, organism = 'hsa', pvalueCutoff = 1)
    
    if (!is.null(kk)) {
      # Save raw results
      write.csv(as.data.frame(kk), file.path(output_subdir, paste0("KEGG_", title_prefix, "_All.csv")), row.names = FALSE)
      
      # Filter
      kk_sig <- filter(kk, p.adjust < 0.05)
      
      if (nrow(kk_sig) > 0) {
        write.csv(as.data.frame(kk_sig), file.path(output_subdir, paste0("KEGG_", title_prefix, "_Significant.csv")), row.names = FALSE)
        
        p_dot <- dotplot(kk_sig, showCategory = 20) + 
          ggtitle(paste0("KEGG Enrichment: ", title_prefix)) +
          scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 60)) +
          theme(axis.text.y = element_text(size = 10))
          
        ggsave(file.path(output_subdir, paste0("Dotplot_KEGG_", title_prefix, ".pdf")), p_dot, width = 10, height = 8)
      }
    }
  }, error = function(e) {
     cat(paste0("    Error in KEGG: ", e$message, "\n"))
  })
}

# ==============================================================================
# 4. Main Execution / 
# ==============================================================================

# Function to generate correlation plot / 
# Correlation Method Assessment / :
# We evaluate both Pearson and Spearman. 
# Pearson assesses linear relationship (y = mx + c), suitable for direct replication.
# Spearman assesses monotonic relationship (rank order), more robust to outliers and suitable for checking consistency of direction across different scales/traits.
# Decision: We will calculate both, but prioritize Spearman for the subtitle and interpretation of "consistency".
# : , Spearman“”. 
# Update: Use OR values (MR Effect Size) and optimized Spearman correlation plot with Cell Type coloring and Gene Labels.
# : OR（MR）Spearman, . 
# Added filter_type to support Positive/Negative versions
generate_correlation_plot <- function(data_df, outcome1, outcome2, output_dir, filter_type = "all") {
  
  cat(paste0("    Generating optimized correlation plot for ", outcome1, " vs ", outcome2, " (using OR, Filter: ", filter_type, ")...\n"))
  
  # Filter data for the two outcomes
  df_subset <- data_df %>% dplyr::filter(outcome %in% c(outcome1, outcome2))
  
  # Check if we have both outcomes
  if (n_distinct(df_subset$outcome) < 2) {
    cat("      Not enough outcomes for correlation analysis. Skipping.\n")
    return(NULL)
  }
  
  # Ensure OR column exists or calculate from b
  if (!"or" %in% colnames(df_subset) || all(is.na(df_subset$or))) {
    if ("b" %in% colnames(df_subset)) {
       df_subset$or <- exp(df_subset$b)
    } else {
       cat("      Error: No 'or' or 'b' column found for correlation analysis.\n")
       return(NULL)
    }
  }
  
  # Pivot to wide format to get ORs side-by-side
  # OR
  # We need to keep 'cell' for coloring and 'gene_name' for labeling
  df_wide <- df_subset %>%
    select(gene_name, cell, lead_snp, outcome, or) %>%
    pivot_wider(names_from = outcome, values_from = or)
  
  # Remove rows with NA (where one outcome is missing for the pair)
  df_wide <- na.omit(df_wide)

  # Apply Filter based on filter_type
  #  filter_type 
  out1_col <- outcome1
  out2_col <- outcome2
  
  if (filter_type == "positive") {
    # Positive correlation: Both OR > 1 OR Both OR < 1? 
    # Usually "Positive Correlation" means the points lie on the x=y diagonal.
    # User request: "," (One positive correlation version, one negative correlation version).
    # This likely implies filtering for points that contribute to a positive correlation (slope > 0) vs negative correlation (slope < 0).
    # Slope > 0 means (x > 1 and y > 1) OR (x < 1 and y < 1).
    # Slope < 0 means (x > 1 and y < 1) OR (x < 1 and y > 1).
    
    # Let's use log(OR) to determine direction more easily.
    log_x <- log(df_wide[[out1_col]])
    log_y <- log(df_wide[[out2_col]])
    
    # Consistent direction (Positive Correlation subset)
    keep_idx <- (log_x * log_y) > 0
    df_wide <- df_wide[keep_idx, ]
    subtitle_prefix <- " (Positive Trend Subset)"
    
  } else if (filter_type == "negative") {
    # Negative Correlation subset (Inconsistent direction)
    log_x <- log(df_wide[[out1_col]])
    log_y <- log(df_wide[[out2_col]])
    
    # Inconsistent direction (Negative Correlation subset)
    keep_idx <- (log_x * log_y) < 0
    df_wide <- df_wide[keep_idx, ]
    subtitle_prefix <- " (Negative Trend Subset)"
    
  } else {
    subtitle_prefix <- ""
  }
  
  if (nrow(df_wide) < 3) {
    cat(paste0("      Not enough points (<3) for correlation analysis (Filter: ", filter_type, "). Skipping.\n"))
    return(NULL)
  }
  
  # Rename columns for easier access (handle spaces in outcome names)
  # (Already done above via variables, but for ggplot we need the names)
  
  # Calculate Spearman correlation ONLY
  # Spearman
  cor_spearman <- cor(df_wide[[out1_col]], df_wide[[out2_col]], method = "spearman")
  
  # Determine statistical significance of correlation
  test_res <- cor.test(df_wide[[out1_col]], df_wide[[out2_col]], method = "spearman", exact = FALSE)
  p_val <- test_res$p.value
  
  p <- ggplot(df_wide, aes_string(x = paste0("`", out1_col, "`"), y = paste0("`", out2_col, "`"))) +
    # Add regression line (linear fit to show trend)
    geom_smooth(method = "lm", color = "black", linetype = "dashed", se = TRUE, alpha = 0.15, size = 0.5) +
    # Add points colored by Cell Type
    geom_point(aes(fill = cell), shape = 21, size = 3.5, color = "white", alpha = 0.85, stroke = 0.5) +
    # Add labels for Gene Name with repulsion to avoid overlap
    geom_text_repel(aes(label = gene_name), 
                    size = 2.2, # Reduced font size as requested
                    box.padding = 0.5, 
                    point.padding = 0.3,
                    segment.color = "grey60",
                    segment.size = 0.3,
                    max.overlaps = 30,
                    force = 2) +
    scale_fill_igv() + 
    # Theme optimization
    theme_classic(base_size = 14) +
    labs(# title = paste0(outcome1, " vs ", outcome2, subtitle_prefix), # Title removed as requested
         subtitle = paste0("Spearman rho = ", round(cor_spearman, 3), 
                           ", P = ", formatC(p_val, format = "e", digits = 2)),
         x = paste0("Odds Ratio (", outcome1, ")"),
         y = paste0("Odds Ratio (", outcome2, ")"),
         fill = "Cell Type") +
    theme(
      plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
      plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey30"),
      axis.title = element_text(face = "bold", size = 14),
      axis.text = element_text(size = 12, color = "black"),
      legend.position = "right",
      legend.title = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 10),
      axis.line = element_line(color = "black", size = 0.8)
    )
    
  # Log scale axes for OR plots
  p <- p + scale_x_log10() + scale_y_log10() + annotation_logticks()
  
  # Construct filename suffix based on filter
  filter_suffix <- ""
  if (filter_type == "positive") filter_suffix <- "_Positive"
  if (filter_type == "negative") filter_suffix <- "_Negative"
  
  # Save plots
  ggsave(file.path(output_dir, paste0("Correlation_OR_", gsub("[^[:alnum:]]", "_", outcome1), "_vs_", gsub("[^[:alnum:]]", "_", outcome2), filter_suffix, ".pdf")), p, width = 8, height = 7)
  ggsave(file.path(output_dir, paste0("Correlation_OR_", gsub("[^[:alnum:]]", "_", outcome1), "_vs_", gsub("[^[:alnum:]]", "_", outcome2), filter_suffix, ".png")), p, width = 8, height = 7, dpi = 300)
  
  cat("      Optimized correlation plot saved.\n")
  return(p)
}

# Function to perform Validation Analysis Module / 
perform_validation_analysis_module <- function(df_main, df_val, outcome_main, outcome_val, module_name, output_dir, skip_correlation_plot = FALSE, use_or_for_correlation = FALSE) {
  
  cat(paste0("  Running Validation Module: ", module_name, " (", outcome_main, " vs ", outcome_val, ")...\n"))
  
  # Ensure sub-directory exists
  mod_dir <- file.path(output_dir, module_name)
  if (!dir.exists(mod_dir)) dir.create(mod_dir, recursive = TRUE)
  
  # 1. Join Data (Inner Join on Gene + Cell + Lead SNP)
  # 1. （++Lead SNP）
  
  # Rename columns in validation to avoid collision during join
  df_val_renamed <- df_val %>%
    select(gene_name, cell, lead_snp, b, se, pval, or, or_lci95, or_uci95, source_type) %>%
    rename(
      b_val = b,
      se_val = se,
      pval_val = pval,
      or_val = or,
      or_lci_val = or_lci95,
      or_uci_val = or_uci95,
      source_val = source_type
    )
  
  combined_df <- df_main %>%
    inner_join(df_val_renamed, by = c("gene_name", "cell", "lead_snp"))
  
  if (nrow(combined_df) == 0) {
    cat("    No matching Gene-Cell-SNP triplets found. Trying Gene-SNP only...\n")
    # Fallback: Join on Gene + SNP (ignoring Cell if mismatched naming)
    combined_df <- df_main %>%
      inner_join(df_val_renamed %>% select(-cell), by = c("gene_name", "lead_snp"))
    
    if (nrow(combined_df) == 0) {
      cat("    No matching Gene-SNP pairs found. Skipping module.\n")
      return(NULL)
    } else {
      cat(paste0("    Found ", nrow(combined_df), " matches based on Gene-SNP (ignoring Cell).\n"))
    }
  } else {
    cat(paste0("    Found ", nrow(combined_df), " matches based on Gene-Cell-SNP.\n"))
  }
  
  # 2. Statistics: Direction Consistency / : 
  # Consistent if b_main * b_val > 0 (Comparing MR Effect Directions)
  combined_df$direction_consistent <- (combined_df$b * combined_df$b_val) > 0
  
  # Deduplicate Logic (Optimized for Validation Statistics) / 
  # Strategy: Deduplicate by Gene + Cell, keeping the best SNP (lowest pval in Main outcome)
  # : +, PSNP
  
  # 1. Filter to unique Gene-Cell pairs based on best SNP
  unique_gene_cell_pairs <- combined_df %>%
    arrange(pval) %>% # Sort by Main Outcome P-value (ascending) / P
    group_by(gene_name, cell) %>%
    slice(1) %>% # Keep top 1 / 
    ungroup()
    
  # 2. Calculate Consistency on these unique pairs
  gene_cell_stats <- unique_gene_cell_pairs %>%
    select(gene_name, cell, direction_consistent) %>%
    rename(is_consistent_gene_cell = direction_consistent)

  n_consistent <- sum(gene_cell_stats$is_consistent_gene_cell, na.rm = TRUE)
  n_inconsistent <- sum(!gene_cell_stats$is_consistent_gene_cell, na.rm = TRUE)
  total <- nrow(gene_cell_stats)
  pct_consistent <- (n_consistent / total) * 100
  
  cat(paste0("    Direction Consistency (Gene-Cell Level): ", n_consistent, "/", total, " (", round(pct_consistent, 1), "%)\n"))
  
  # Save Statistics
  stats_df <- data.frame(
    Module = module_name,
    Outcome_Main = outcome_main,
    Outcome_Validation = outcome_val,
    N_Matches_Gene_Cell = total,
    N_Consistent = n_consistent,
    N_Inconsistent = n_inconsistent,
    Pct_Consistent = pct_consistent
  )
  write.csv(stats_df, file.path(mod_dir, paste0("Statistics_Consistency_", module_name, ".csv")), row.names = FALSE)
  
  # 3. Correlation Analysis / 
  # Select Data for Plotting based on Module / 
  if (module_name == "Validation_mvAge_vs_Frailty") {
     # Use ALL intersection pairs for mvAge vs Frailty (User Request) / mvAge vs Frailty 
     plot_df <- combined_df
     cat("    Using ALL cell-gene-snp pairs for mvAge vs Frailty Correlation Plot.\n")
  } else {
     # Use Unique pairs for others (e.g. Shingles Validation) / 
     plot_df <- unique_gene_cell_pairs
     cat("    Using UNIQUE gene-cell pairs for Validation Correlation Plot.\n")
  }
  
  if (!skip_correlation_plot) {
    if (nrow(plot_df) >= 3) {
      
      if (use_or_for_correlation) {
         # Use OR values / OR
         val_main <- plot_df$or
         val_validation <- plot_df$or_val
         label_main <- paste0("Odds Ratio (", outcome_main, ")")
         label_val <- paste0("Odds Ratio (", outcome_val, ")")
      } else {
         # Use Beta values / Beta
         val_main <- plot_df$b
         val_validation <- plot_df$b_val
         label_main <- paste0("Beta (", outcome_main, ")")
         label_val <- paste0("Beta (", outcome_val, ")")
      }
    
      cor_pearson <- cor(val_main, val_validation, method = "pearson", use = "complete.obs")
      cor_spearman <- cor(val_main, val_validation, method = "spearman", use = "complete.obs")
      
      # Significance test (Spearman)
      test_res <- cor.test(val_main, val_validation, method = "spearman", exact = FALSE)
      p_val <- test_res$p.value
      
      # Determine Subtitle based on User Request (mvAge vs Frailty -> Only Spearman)
      if (outcome_main == "mvAge" && outcome_val == "Frailty Index") {
          subtitle_text <- paste0("Spearman rho = ", round(cor_spearman, 3), " (P = ", formatC(p_val, format = "e", digits = 2), ")",
                                  "\nConsistency (Gene-Cell): ", n_consistent, "/", total, " (", round(pct_consistent, 1), "%)")
      } else {
          subtitle_text <- paste0("Spearman rho = ", round(cor_spearman, 3), " (P = ", formatC(p_val, format = "e", digits = 2), ")",
                                  "\nPearson r = ", round(cor_pearson, 3), 
                                  "\nConsistency (Gene-Cell): ", n_consistent, "/", total, " (", round(pct_consistent, 1), "%)")
      }

      # Plot with Optimization (IGV colors)
      # : IGV, 
      p <- ggplot(plot_df, aes(x = val_main, y = val_validation)) +
        geom_point(aes(fill = cell), shape = 21, size = 3, color = "black", alpha = 0.8) +
        geom_smooth(method = "lm", color = "black", linetype = "dashed", se = TRUE) +
        scale_fill_igv() + 
        labs(
          title = paste0("Validation: ", outcome_main, " vs ", outcome_val),
          subtitle = subtitle_text,
          x = label_main,
          y = label_val,
          fill = "Cell Type"
        ) +
        theme_bw + # Use theme_bw for border / theme_bw
        theme(
          plot.title = element_text(face = "bold", hjust = 0.5, size = 16),
          plot.subtitle = element_text(hjust = 0.5, size = 12, color = "grey30"),
          axis.title = element_text(face = "bold", size = 14),
          axis.text = element_text(size = 12, color = "black"),
          legend.position = "right",
          legend.title = element_text(face = "bold", size = 12),
          legend.text = element_text(size = 10),
          axis.line = element_line(color = "black", size = 0.8)
        )
      
      # Logic for Gene Labels (No labels for mvAge vs Frailty)
      if (outcome_main != "mvAge" || outcome_val != "Frailty Index") {
         # Add labels for non-mvAge-Frailty plots
         p <- p + geom_text_repel(aes(label = gene_name), 
                        size = 2.2,
                        box.padding = 0.5, 
                        point.padding = 0.3,
                        segment.color = "grey60",
                        segment.size = 0.3,
                        max.overlaps = 30,
                        force = 2)
      }

      if (use_or_for_correlation) {
          # Use log scales for OR plots / OR
          p <- p + scale_x_log10() + scale_y_log10() + annotation_logticks()
      }

      # Save Plot Data (Using the data actually used for plotting) / 
      write.csv(plot_df, file.path(mod_dir, paste0("Data_for_Correlation_Plot_", module_name, ".csv")), row.names = FALSE)

      ggsave(file.path(mod_dir, paste0("Correlation_Plot_", module_name, ".pdf")), p, width = 8, height = 7)
      ggsave(file.path(mod_dir, paste0("Correlation_Plot_", module_name, ".png")), p, width = 8, height = 7)
    } else {
      cat("    Not enough points for correlation plot.\n")
    }
  } else {
    cat("    Skipping correlation plot as requested.\n")
  }
  
  # 4. Save Inconsistent Data for Review / 
  if (n_inconsistent > 0) {
    inconsistent_df <- unique_gene_cell_pairs %>% filter(!direction_consistent)
    write.csv(inconsistent_df, file.path(mod_dir, paste0("Data_Inconsistent_", module_name, ".csv")), row.names = FALSE)
  }
  
  # 5. Forest Plot for Comparison (Optional: Side-by-side or just listing)
  # For forest plot, we use the same dataset as plotting
  
  keys <- plot_df %>% select(gene_name, cell, lead_snp)
  
  df_main_subset <- df_main %>% semi_join(keys, by = c("gene_name", "cell", "lead_snp"))
  df_val_subset <- df_val %>% semi_join(keys, by = c("gene_name", "cell", "lead_snp"))
  
  df_for_plot <- bind_rows(df_main_subset, df_val_subset)
  
  run_analysis(suffix = paste0("Forest_Plot_", module_name), 
               custom_data = df_for_plot, 
               subdir = mod_dir,
               show_beta_allele = TRUE) # Show beta to visualize direction

  return(stats_df)
}

# Helper function to run analysis / 
# Updated for Directory Redesign and Enrichment / 
run_analysis <- function(gene_subset = NULL, suffix, custom_data = NULL, subdir = NULL, show_beta_allele = TRUE, color_cell_only = FALSE) {
  cat("DEBUG: New run_analysis function loaded and running for suffix:", suffix, "\n")
  
  # Determine base target directory / 
  if (!is.null(subdir)) {
    base_target <- subdir
  } else {
    base_target <- output_dir
  }
  
  # Create a specific subfolder for this version / 
  target_dir <- file.path(base_target, suffix)
  if (!dir.exists(target_dir)) dir.create(target_dir, recursive = TRUE)

  if (!is.null(custom_data)) {
    data_subset <- custom_data
    if (nrow(data_subset) == 0) {
       cat(paste0("No data provided for: ", suffix, ". Skipping.\n"))
       # Cleanup empty dir if created
       if (length(list.files(target_dir)) == 0) unlink(target_dir, recursive = TRUE)
       return(NULL)
    }
    gene_count <- length(unique(data_subset$gene_name))
    cat(paste0("\nProcessing version: ", suffix, " (", gene_count, " genes, filtered by specific rows)...\n"))
    
  } else {
    if (length(gene_subset) == 0) {
      cat(paste0("No genes found for: ", suffix, ". Skipping.\n"))
      # Cleanup empty dir if created
      if (length(list.files(target_dir)) == 0) unlink(target_dir, recursive = TRUE)
      return(NULL)
    }
    
    cat(paste0("\nProcessing version: ", suffix, " (", length(gene_subset), " genes)...\n"))
    data_subset <- data[data$gene_name %in% gene_subset, ]
  }
  
  if (nrow(data_subset) == 0) {
     cat("  No data rows found after filtering. Skipping.\n")
     if (length(list.files(target_dir)) == 0) unlink(target_dir, recursive = TRUE)
     return(NULL)
  }
  
  # 1. Save Data / 
  write.csv(data_subset, file.path(target_dir, paste0("Data_", suffix, ".csv")), row.names = FALSE)
  
  # 2. Save Gene List / 
  genes <- unique(data_subset$gene_name)
  genes <- genes[!is.na(genes)]
  write.table(genes, file.path(target_dir, paste0("Gene_List_", suffix, ".txt")), 
              row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # 3. Prepare plot data / 
  plot_data_hierarchical <- prepare_hierarchical_data(data_subset)
  
  # 4. Generate and Save Plot / 
  p <- generate_hierarchical_plot(plot_data_hierarchical, paste0("Forest_Plot_Hierarchical_", suffix), target_dir, show_beta_allele = show_beta_allele, color_cell_only = color_cell_only)
  
  # 5. Perform Enrichment Analysis / 
  perform_enrichment(genes, file.path(target_dir, "Enrichment"), suffix)
  
  return(p)
}

# Helper function to run FinnGen validation / FinnGen
run_finngen_validation <- function(base_data, base_subdir, suffix_name, match_snp = FALSE) {
  # Check if FinnGen data is available in global env / FinnGen
  if (!exists("data_validation") || nrow(data_validation) == 0) return()
  
  finngen_hz <- data_validation %>% dplyr::filter(outcome == "HZ (FinnGen)")
  if (nrow(finngen_hz) == 0) return()
  
  # Identify keys / 
  if (match_snp) {
    keys <- base_data %>% select(gene_name, cell, lead_snp) %>% distinct()
    finngen_matched <- finngen_hz %>% inner_join(keys, by = c("gene_name", "cell", "lead_snp"))
  } else {
    keys <- base_data %>% select(gene_name, cell) %>% distinct()
    finngen_matched <- finngen_hz %>% inner_join(keys, by = c("gene_name", "cell"))
  }
  
  if (nrow(finngen_matched) > 0) {
    cat(paste0("    Adding FinnGen Validation for: ", suffix_name, " (", nrow(finngen_matched), " rows)...\n"))
    
    # Combine / 
    combined_data <- bind_rows(base_data, finngen_matched)
    
    # Sort Data to ensure consistency / 
    # User Request: Integrate FinnGen data sorting with Main data
    # Sort by: Chr -> Pos -> Gene -> Cell -> Lead SNP -> Outcome Priority
    # This ensures that for the same Gene-Cell-SNP, the three outcomes are listed together.
    
    # Define outcome priority: Main Outcomes first, FinnGen last
     # User Request: aging -> Rheumatoid arthritis -> Herpes Zoster -> Herpes Zoster (FinnGen)
     combined_data <- combined_data %>%
       mutate(
         outcome_priority = case_when(
           outcome == "mvAge" ~ 1,
           outcome == "RA" ~ 1.5,
           outcome == "HZ (UKB)" ~ 2,
           outcome == "HZ (FinnGen)" ~ 3,
           TRUE ~ 99
         )
       )
 
     if (all(c("chr_num", "start_pos") %in% colnames(combined_data))) {
        combined_data <- combined_data %>%
          arrange(chr_num, start_pos, gene_name, cell, lead_snp, outcome_priority, outcome)
     } else {
        combined_data <- combined_data %>%
          arrange(gene_name, cell, lead_snp, outcome_priority, outcome)
     }
    
    # Create subdir / 
    subdir_val <- file.path(base_subdir, "With_FinnGen_Validation")
    if (!dir.exists(subdir_val)) dir.create(subdir_val, recursive = TRUE)

    # Statistics: Direction Consistency with FinnGen / : FinnGen
    cat(paste0("    Calculating consistency statistics for FinnGen Validation...\n"))
    
    # Filter for Main HZ and FinnGen HZ
    hz_main <- combined_data %>% dplyr::filter(outcome == "HZ (UKB)")
    hz_finngen <- combined_data %>% dplyr::filter(outcome == "HZ (FinnGen)")
    
    if (nrow(hz_main) > 0 && nrow(hz_finngen) > 0) {
       # Join to compare
       # Ensure we use the correct keys based on match_snp
       join_keys <- c("gene_name", "cell")
       if (match_snp) join_keys <- c(join_keys, "lead_snp")
       
       # Select minimal columns for join to avoid duplication
       hz_finngen_clean <- hz_finngen %>% 
          select(all_of(c(join_keys, "b", "or", "pval"))) %>%
          rename(b_finngen = b, or_finngen = or, pval_finngen = pval)
       
       consistency_df <- hz_main %>%
          inner_join(hz_finngen_clean, by = join_keys)
          
       if (nrow(consistency_df) > 0) {
          # Check direction (b * b_finngen > 0)
          # Ensure b is numeric
          consistency_df$direction_consistent <- (as.numeric(consistency_df$b) * as.numeric(consistency_df$b_finngen)) > 0
          
          n_consistent <- sum(consistency_df$direction_consistent, na.rm = TRUE)
          total_pairs <- nrow(consistency_df)
          pct_consistent <- (n_consistent / total_pairs) * 100
          
          stats_out <- data.frame(
             Dataset = suffix_name,
             Outcome_Main = "HZ (UKB)",
             Outcome_Validation = "HZ (FinnGen)",
             N_Pairs = total_pairs,
             N_Consistent = n_consistent,
             Pct_Consistent = round(pct_consistent, 2)
          )
          
          write.csv(stats_out, file.path(subdir_val, paste0("Statistics_Consistency_FinnGen_", suffix_name, ".csv")), row.names = FALSE)
          
          # Also save the detailed comparison
          write.csv(consistency_df, file.path(subdir_val, paste0("Details_Consistency_FinnGen_", suffix_name, ".csv")), row.names = FALSE)
       }
    }

    run_analysis(suffix = paste0(suffix_name, "_with_FinnGen"),
                 custom_data = combined_data,
                 subdir = subdir_val,
                 show_beta_allele = TRUE)
  } else {
    cat(paste0("    No matching FinnGen data found for: ", suffix_name, "\n"))
  }
}

# ==============================================================================
# 1. Single Outcome Plots / 
# ==============================================================================
cat("\nProcessing 1. Single Outcome Plots (Multi-threaded)...\n")
unique_outcomes <- unique(data$outcome)
unique_outcomes <- unique_outcomes[!is.na(unique_outcomes)]

for (out in unique_outcomes) {
  tryCatch({
    out_safe <- gsub("[^[:alnum:]]", "_", out)
    single_data <- data %>% dplyr::filter(outcome == out)
    if (nrow(single_data) > 0) {
      run_analysis(
        suffix = paste0("Single_Outcome_", out_safe),
        custom_data = single_data,
        subdir = dir_single_outcome
      )
    }
  }, error = function(e) {
    cat(paste0("Error in Single Outcome Plot for ", out, ": ", e$message, "\n"))
  })
}

# ==============================================================================
# 2. Joint Outcome Plot (mvAge + Herpes Zoster) / 
# ==============================================================================
cat("\nProcessing 2. Joint Outcome Plot (mvAge + Herpes Zoster)...\n")
joint_outcomes <- c("mvAge", "HZ (UKB)")
joint_data <- data %>% dplyr::filter(outcome %in% joint_outcomes)

if (nrow(joint_data) > 0) {
  run_analysis(suffix = "Joint_mvAge_Herpes_Zoster", 
               custom_data = joint_data, 
               subdir = dir_joint_outcome)
}

# ==============================================================================
# 2.1 Joint Outcome Plot (mvAge + RA + HZ) /  
# ==============================================================================
cat("\nProcessing 2.1 Joint Outcome Plot (mvAge + RA + HZ) - Union...\n")
union_outcomes <- c("mvAge", "RA", "HZ (UKB)")
union_data <- data %>% dplyr::filter(outcome %in% union_outcomes)

if (nrow(union_data) > 0) {
  run_analysis(suffix = "Joint_mvAge_RA_HZ_Union", 
               custom_data = union_data, 
               subdir = dir_joint_outcome)
}

# ==============================================================================
# 3. Gene-Level Shared Outcome Plots / 
# ==============================================================================
cat("\nProcessing 3. Gene-Level Shared Outcome Plots (Same Gene, Multiple Outcomes, Any Cell)...\n")

gene_shared_pairs_list <- list(
  c("mvAge", "HZ (UKB)"),
  c("RA", "HZ (UKB)")
)

gene_shared_pair_names <- c(
  "mvAge_vs_Herpes_Zoster",
  "RA_vs_Herpes_Zoster"
)

for (i in 1:length(gene_shared_pairs_list)) {
  pair <- gene_shared_pairs_list[[i]]
  p_name <- gene_shared_pair_names[i]
  
  cat(paste0("  Checking Shared Genes for pair: ", pair[1], " vs ", pair[2], "...\n"))
  
  # 1. Identify Shared Genes
  genes1 <- unique(data$gene_name[data$outcome == pair[1]])
  genes2 <- unique(data$gene_name[data$outcome == pair[2]])
  shared_genes <- intersect(genes1, genes2)
  shared_genes <- shared_genes[!is.na(shared_genes)]
  
  if (length(shared_genes) > 0) {
    cat(paste0("    Found ", length(shared_genes), " shared genes.\n"))
    
    # 2. Extract data for these genes AND these outcomes
    shared_gene_data <- data %>%
      dplyr::filter(gene_name %in% shared_genes) %>%
      dplyr::filter(outcome %in% pair)
    
    # 3. Generate Plot
    subdir_path <- file.path(dir_gene_shared, p_name)
    run_analysis(suffix = paste0("Gene_Shared_", p_name), 
                 custom_data = shared_gene_data, 
                 subdir = subdir_path)
    
  } else {
    cat(paste0("    No shared genes found for pair: ", p_name, "\n"))
  }
}

# ==============================================================================
# 4. Intersection Cell-Gene Plots / -
# ==============================================================================
# Version 1: Complete Version (Gene + Cell same, regardless of SNP/Direction)
# Version 2: Clean Version (Gene + Cell same, AND direction consistent within each outcome)
cat("\nProcessing 4. Intersection Cell-Gene Plots...\n")

cell_gene_pairs_list <- list(
  c("mvAge", "HZ (UKB)"),
  c("RA", "HZ (UKB)")
)
cell_gene_pair_names <- c(
  "mvAge_vs_Herpes_Zoster",
  "RA_vs_Herpes_Zoster"
)

for (i in 1:length(cell_gene_pairs_list)) {
  pair <- cell_gene_pairs_list[[i]]
  p_name <- cell_gene_pair_names[i]
  
  cat(paste0("  Checking Shared Cell-Gene for pair: ", pair[1], " vs ", pair[2], "...\n"))
  
  # Filter data for this pair
  pair_data <- data %>% dplyr::filter(outcome %in% pair)
  
  # Group by Gene, Cell to find those present in both outcomes
  pair_counts <- pair_data %>%
    group_by(gene_name, cell) %>%
    summarize(n_out = n_distinct(outcome), .groups = "drop") %>%
    dplyr::filter(n_out == 2) # Must have both
  
  if (nrow(pair_counts) > 0) {
    # Extract original rows for these valid gene-cell pairs
    final_pair_data <- pair_data %>%
      inner_join(pair_counts[, c("gene_name", "cell")], by = c("gene_name", "cell"))
      
    # 4.1 Complete Version
    subdir_complete <- file.path(dir_cell_gene_shared, p_name, "Complete")
    run_analysis(suffix = paste0("Cell_Gene_Shared_", p_name, "_Complete"), 
                 custom_data = final_pair_data, 
                 subdir = subdir_complete)
                 
    # Correlation Plot for Complete Set
    generate_correlation_plot(final_pair_data, pair[1], pair[2], file.path(subdir_complete, paste0("Cell_Gene_Shared_", p_name, "_Complete")))
    # Add Positive/Negative versions for mvAge/RA vs HZ
    if (p_name %in% c("mvAge_vs_Herpes_Zoster", "RA_vs_Herpes_Zoster")) {
       generate_correlation_plot(final_pair_data, pair[1], pair[2], file.path(subdir_complete, paste0("Cell_Gene_Shared_", p_name, "_Complete")), filter_type = "positive")
       generate_correlation_plot(final_pair_data, pair[1], pair[2], file.path(subdir_complete, paste0("Cell_Gene_Shared_", p_name, "_Complete")), filter_type = "negative")
       
       # Add FinnGen Validation Version
       run_finngen_validation(final_pair_data, subdir_complete, paste0("Cell_Gene_Shared_", p_name, "_Complete"), match_snp = FALSE)
    }

    # 4.2 Clean Version (Consistent Direction per Outcome)
    cat("    Filtering for Clean Version (Consistent Direction within Outcome)...\n")
    
    # Function to check consistency for a gene-cell-outcome triplet
    # We need to know if all SNPs for (Gene, Cell, Outcome) have the same OR direction.
    # Check if all OR > 1 or all OR < 1 for each (Gene, Cell, Outcome)
    
    # Ensure 'or' is available
    if (!"or" %in% colnames(final_pair_data)) {
       if ("beta.exposure" %in% colnames(final_pair_data)) {
           final_pair_data$or <- exp(final_pair_data$beta.exposure)
       } else if ("b" %in% colnames(final_pair_data)) {
           final_pair_data$or <- exp(final_pair_data$b)
       }
    }
    
    # Calculate consistency per Gene-Cell-Outcome using OR values
    consistency_check <- final_pair_data %>%
      group_by(gene_name, cell, outcome) %>%
      summarize(
        all_risk = all(or > 1, na.rm = TRUE),
        all_prot = all(or < 1, na.rm = TRUE),
        is_consistent = (all(or > 1, na.rm = TRUE) | all(or < 1, na.rm = TRUE)),
        .groups = "drop"
      )
    
    # Now check if a Gene-Cell pair is consistent for BOTH outcomes
    valid_gene_cells <- consistency_check %>%
      group_by(gene_name, cell) %>%
      summarize(
        both_consistent = all(is_consistent),
        n_outcomes_checked = n(),
        .groups = "drop"
      ) %>%
      dplyr::filter(both_consistent == TRUE & n_outcomes_checked == 2)
      
    if (nrow(valid_gene_cells) > 0) {
      clean_data <- final_pair_data %>%
        inner_join(valid_gene_cells[, c("gene_name", "cell")], by = c("gene_name", "cell"))
        
      subdir_clean <- file.path(dir_cell_gene_shared, p_name, "Clean")
      run_analysis(suffix = paste0("Cell_Gene_Shared_", p_name, "_Clean"), 
                   custom_data = clean_data, 
                   subdir = subdir_clean)
                   
      # Correlation Plot for Clean Set
      generate_correlation_plot(clean_data, pair[1], pair[2], file.path(subdir_clean, paste0("Cell_Gene_Shared_", p_name, "_Clean")))
      
      if (p_name %in% c("mvAge_vs_Herpes_Zoster", "RA_vs_Herpes_Zoster")) {
         generate_correlation_plot(clean_data, pair[1], pair[2], file.path(subdir_clean, paste0("Cell_Gene_Shared_", p_name, "_Clean")), filter_type = "positive")
         generate_correlation_plot(clean_data, pair[1], pair[2], file.path(subdir_clean, paste0("Cell_Gene_Shared_", p_name, "_Clean")), filter_type = "negative")
         
         # Add FinnGen Validation Version
         run_finngen_validation(clean_data, subdir_clean, paste0("Cell_Gene_Shared_", p_name, "_Clean"), match_snp = FALSE)
      }
      
    } else {
      cat("    No consistent gene-cell pairs found for Clean Version.\n")
    }

  } else {
    cat(paste0("  No shared cell-gene pairs found for: ", p_name, "\n"))
  }
}

# ==============================================================================
# 4.3 Three-way Intersection Cell-Gene Plots (mvAge, RA, HZ) / -
# ==============================================================================
cat("\nProcessing 4.3 Three-way Intersection Cell-Gene Plots (mvAge, RA, HZ)...\n")

triplet_outcomes <- c("mvAge", "RA", "HZ (UKB)")
triplet_name <- "mvAge_vs_RA_vs_HZ"

cat(paste0("  Checking Shared Cell-Gene for triplet: ", paste(triplet_outcomes, collapse=", "), "...\n"))

# Filter data for this triplet
triplet_data <- data %>% dplyr::filter(outcome %in% triplet_outcomes)

# Group by Gene, Cell to find those present in ALL 3 outcomes
triplet_counts <- triplet_data %>%
  group_by(gene_name, cell) %>%
  summarize(n_out = n_distinct(outcome), .groups = "drop") %>%
  dplyr::filter(n_out == 3) # Must have all 3

if (nrow(triplet_counts) > 0) {
  # Extract original rows
  final_triplet_data <- triplet_data %>%
    inner_join(triplet_counts[, c("gene_name", "cell")], by = c("gene_name", "cell"))
    
  # Complete Version
  subdir_triplet <- file.path(dir_cell_gene_shared, triplet_name, "Complete")
  run_analysis(suffix = paste0("Cell_Gene_Shared_", triplet_name, "_Complete"), 
               custom_data = final_triplet_data, 
               subdir = subdir_triplet)

  # Generate Pairwise Correlation Plots for the Triplet Data
  pairs_to_plot <- list(
      c("mvAge", "RA"),
      c("mvAge", "HZ (UKB)"),
      c("RA", "HZ (UKB)")
  )
  
  for (pair in pairs_to_plot) {
      # Use tryCatch to avoid stopping if one pair fails
      tryCatch({
          generate_correlation_plot(final_triplet_data, pair[1], pair[2], file.path(subdir_triplet, paste0("Cell_Gene_Shared_", triplet_name, "_Complete")))
      }, error = function(e) {
          cat(paste0("    Error generating correlation plot for ", pair[1], " vs ", pair[2], ": ", e$message, "\n"))
      })
  }

} else {
  cat(paste0("    No shared cell-gene triplets found for: ", triplet_name, "\n"))
}

# ==============================================================================
# 5. Intersection Cell-Gene-SNP Plots / --SNP
# ==============================================================================
cat("\nProcessing 5. Intersection Cell-Gene-SNP Plots (Same Gene, Cell, SNP)...\n")

snp_pair_names <- c(
  "mvAge_vs_Herpes_Zoster",
  "RA_vs_Herpes_Zoster",
  "mvAge_vs_RA"
)

snp_pairs_list <- list(
  c("mvAge", "HZ (UKB)"),
  c("RA", "HZ (UKB)"),
  c("mvAge", "RA")
)

for (i in 1:length(snp_pairs_list)) {
  pair <- snp_pairs_list[[i]]
  p_name <- snp_pair_names[i]
  
  # Filter data for this pair
  pair_data <- data %>% dplyr::filter(outcome %in% pair)
  
  # Group by Gene, Cell, SNP
  pair_counts <- pair_data %>%
    group_by(gene_name, cell, lead_snp) %>%
    summarize(n_out = n_distinct(outcome), .groups = "drop") %>%
    dplyr::filter(n_out == 2)
    
  if (nrow(pair_counts) > 0) {
    final_pair_data <- pair_data %>%
      inner_join(pair_counts[, c("gene_name", "cell", "lead_snp")], by = c("gene_name", "cell", "lead_snp"))
      
    # Create subdir / 
    subdir_path <- file.path(dir_cell_gene_snp_shared, p_name)
    
    # 5.1 Standard Version (Intersection)
    p_forest <- run_analysis(suffix = paste0("Cell_Gene_SNP_Shared_", p_name), 
                 custom_data = final_pair_data, 
                 subdir = subdir_path,
                 show_beta_allele = TRUE)
    
    # Correlation Plot
    generate_correlation_plot(final_pair_data, pair[1], pair[2], file.path(subdir_path, paste0("Cell_Gene_SNP_Shared_", p_name)))
    
    if (p_name %in% c("mvAge_vs_Herpes_Zoster", "RA_vs_Herpes_Zoster", "mvAge_vs_RA")) {
       p_cor_pos <- generate_correlation_plot(final_pair_data, pair[1], pair[2], file.path(subdir_path, paste0("Cell_Gene_SNP_Shared_", p_name)), filter_type = "positive")
       p_cor_neg <- generate_correlation_plot(final_pair_data, pair[1], pair[2], file.path(subdir_path, paste0("Cell_Gene_SNP_Shared_", p_name)), filter_type = "negative")
       
       # Combine plots if all exist
       if (!is.null(p_forest) && !is.null(p_cor_pos) && !is.null(p_cor_neg)) {
          cat(paste0("    Generating Combined Plot for ", pair[1], " vs ", pair[2], "...\n"))
          
          # Convert forest plot (gtable) to ggplot/grob for patchwork
          # patchwork works with grobs using wrap_elements
          p_forest_grob <- wrap_elements(p_forest)
          
          # Combine: Forest Plot (Left) | (Pos (Top Right) / Neg (Bottom Right))
          # Use shared legend (guides = "collect")
          # Adjusted layout for wider forest plot and better legend positioning
          combined_plot <- (p_forest_grob | (p_cor_pos / p_cor_neg)) + 
            plot_layout(widths = c(3, 1), guides = "collect") + 
            plot_annotation(tag_levels = 'A') & 
            theme(legend.position = "right", legend.box = "vertical") # Legend at right
            
          # Save combined plot
          combined_dir <- file.path(subdir_path, paste0("Combined_Plot_", p_name))
          if (!dir.exists(combined_dir)) dir.create(combined_dir, recursive = TRUE)
          
          ggsave(file.path(combined_dir, paste0("Combined_Plot_", p_name, ".pdf")), combined_plot, width = 32, height = 15) # Significantly increased width
          ggsave(file.path(combined_dir, paste0("Combined_Plot_", p_name, ".png")), combined_plot, width = 32, height = 15, dpi = 300)
          
          # Generate Legend Markdown
          legend_content <- paste0(
            "## Figure Legend\n\n",
            "**Figure 1. Causal Associations between ", pair[1], " and ", pair[2], " driven by Shared Genetic Variants.**\n\n",
            "**(A)** Hierarchical forest plot showing the Mendelian Randomization (MR) estimates for the causal effect of molecular traits (gene expression/protein levels) on both ", pair[1], " and ", pair[2], ". Results are grouped by Gene and indented by Cell Type. Only traits with significant effects on both outcomes driven by the same Lead SNP are shown. Columns display the Lead SNP, F-statistic, MR Method, FDR-adjusted P-value, and Odds Ratio (95% CI). Colocalization probabilities (SuSiE and ABF methods) are shown on the right.\n\n",
            "**(B)** Scatter plot showing the correlation of MR effect sizes (Odds Ratios) between ", pair[1], " and ", pair[2], " for the subset of traits showing a **positive correlation** (consistent direction of effect). Each point represents a Gene-Cell-SNP triplet, colored by cell type. The regression line (dashed) indicates the overall trend. Spearman correlation coefficient and P-value are provided.\n\n",
            "**(C)** Scatter plot showing the correlation of MR effect sizes for the subset of traits showing a **negative correlation** (opposite direction of effect). This highlights potential antagonistic pleiotropy or tissue-specific regulatory mechanisms.\n"
          )
          
          writeLines(legend_content, file.path(combined_dir, "Figure_Legend.md"))
          cat("    Combined plot and legend saved.\n")
       }
       
       # Add FinnGen Validation Version
       run_finngen_validation(final_pair_data, subdir_path, paste0("Cell_Gene_SNP_Shared_", p_name), match_snp = TRUE)
    }
    
    # 5.2 Sub-version: Inconsistent Directions / : 
    # (Checking if effect on mvAge is opposite to effect on HZ)
    if (p_name %in% c("mvAge_vs_Herpes_Zoster", "RA_vs_Herpes_Zoster", "mvAge_vs_RA")) {
       cat(paste0("    Analyzing Inconsistent Directions for ", pair[1], " vs ", pair[2], " (SNP Level)...\n"))
       
       if (!"b" %in% colnames(final_pair_data)) {
          if ("or" %in% colnames(final_pair_data)) {
             final_pair_data$b <- log(final_pair_data$or)
          } else {
             final_pair_data$b <- NA
          }
       }

       # Pivot to check directions
       wide_check <- final_pair_data %>%
         select(gene_name, cell, lead_snp, outcome, b) %>%
         pivot_wider(names_from = outcome, values_from = b) %>%
         na.omit()
       
       # Filter where signs are different (Product of betas < 0)
       inconsistent_pairs <- wide_check %>%
         filter((.data[[pair[1]]] * .data[[pair[2]]]) < 0)
       
       if (nrow(inconsistent_pairs) > 0) {
          inconsistent_data <- final_pair_data %>%
            semi_join(inconsistent_pairs, by = c("gene_name", "cell", "lead_snp"))
            
          run_analysis(suffix = paste0("Cell_Gene_SNP_Shared_", p_name, "_Inconsistent"), 
                       custom_data = inconsistent_data, 
                       subdir = subdir_path)
                       
          # Correlation for Inconsistent
          generate_correlation_plot(inconsistent_data, pair[1], pair[2], file.path(subdir_path, paste0("Cell_Gene_SNP_Shared_", p_name, "_Inconsistent")))
          
          if (p_name %in% c("mvAge_vs_Herpes_Zoster", "RA_vs_Herpes_Zoster", "mvAge_vs_RA")) {
             generate_correlation_plot(inconsistent_data, pair[1], pair[2], file.path(subdir_path, paste0("Cell_Gene_SNP_Shared_", p_name, "_Inconsistent")), filter_type = "positive")
             generate_correlation_plot(inconsistent_data, pair[1], pair[2], file.path(subdir_path, paste0("Cell_Gene_SNP_Shared_", p_name, "_Inconsistent")), filter_type = "negative")
             
             # Add FinnGen Validation Version
             run_finngen_validation(inconsistent_data, subdir_path, paste0("Cell_Gene_SNP_Shared_", p_name, "_Inconsistent"), match_snp = TRUE)
          }
          
       } else {
          cat("      No inconsistent pairs found.\n")
       }
    }
  }
}

# ==============================================================================
# 5.3 Three-way Intersection Cell-Gene-SNP Plots (mvAge, RA, HZ) / --SNP
# ==============================================================================
cat("\nProcessing 5.3 Three-way Intersection Cell-Gene-SNP Plots (mvAge, RA, HZ)...\n")

triplet_outcomes <- c("mvAge", "RA", "HZ (UKB)")
triplet_name <- "mvAge_vs_RA_vs_HZ"

cat(paste0("  Checking Shared Cell-Gene-SNP for triplet: ", paste(triplet_outcomes, collapse=", "), "...\n"))

# Filter data for this triplet
triplet_data <- data %>% dplyr::filter(outcome %in% triplet_outcomes)

# Group by Gene, Cell, SNP to find those present in ALL 3 outcomes
triplet_counts <- triplet_data %>%
  group_by(gene_name, cell, lead_snp) %>%
  summarize(n_out = n_distinct(outcome), .groups = "drop") %>%
  dplyr::filter(n_out == 3) # Must have all 3

if (nrow(triplet_counts) > 0) {
  # Extract original rows
  final_triplet_data <- triplet_data %>%
    inner_join(triplet_counts[, c("gene_name", "cell", "lead_snp")], by = c("gene_name", "cell", "lead_snp"))
    
  # Create subdir
  subdir_triplet <- file.path(dir_cell_gene_snp_shared, triplet_name)
  
  run_analysis(suffix = paste0("Cell_Gene_SNP_Shared_", triplet_name), 
               custom_data = final_triplet_data, 
               subdir = subdir_triplet,
               show_beta_allele = TRUE)
               
  # Correlation Plots (Pairwise)
  pairs_to_plot <- list(
      c("mvAge", "RA"),
      c("mvAge", "HZ (UKB)"),
      c("RA", "HZ (UKB)")
  )
  
  for (pair in pairs_to_plot) {
      tryCatch({
          generate_correlation_plot(final_triplet_data, pair[1], pair[2], file.path(subdir_triplet, paste0("Cell_Gene_SNP_Shared_", triplet_name)))
      }, error = function(e) {
          cat(paste0("    Error generating correlation plot for ", pair[1], " vs ", pair[2], ": ", e$message, "\n"))
      })
  }

} else {
  cat(paste0("    No shared cell-gene-snp triplets found for: ", triplet_name, "\n"))
}

# ==============================================================================
# 5.4 Added Output Version (Integrated Forest + Specific Spearman) / （ + ）
# ==============================================================================
cat("\nProcessing 5.4 Added Output Version for Cell-Gene-SNP (Integrated mvAge/RA/HZ)...\n")

# Define pairs and their filters
pair_specs <- list(
  list(pair = c("mvAge", "RA"), name = "mvAge_vs_RA", corr_filters = c("negative")),
  list(pair = c("mvAge", "HZ (UKB)"), name = "mvAge_vs_HZ", corr_filters = c("positive", "negative")),
  list(pair = c("RA", "HZ (UKB)"), name = "RA_vs_HZ", corr_filters = c("positive", "negative"))
)

# 1. Identify all triplets (Gene-Cell-SNP) involved in ANY of the specified intersections
triplet_keys <- data.frame(gene_name = character(), cell = character(), lead_snp = character(), stringsAsFactors = FALSE)

for (spec in pair_specs) {
  pair <- spec$pair
  # Filter data for this pair
  pair_data <- data %>% dplyr::filter(outcome %in% pair)
  
  # Find shared Cell-Gene-SNP
  pair_counts <- pair_data %>%
    group_by(gene_name, cell, lead_snp) %>%
    summarize(n_out = n_distinct(outcome), .groups = "drop") %>%
    dplyr::filter(n_out == 2)
  
  if (nrow(pair_counts) > 0) {
    triplet_keys <- bind_rows(triplet_keys, pair_counts[, c("gene_name", "cell", "lead_snp")])
  }
}

# Deduplicate keys
triplet_keys <- unique(triplet_keys)

if (nrow(triplet_keys) == 0) {
  cat("    No shared cell-gene-snp triplets found for any pair. Skipping Section 5.4.\n")
} else {
  cat(paste0("    Found ", nrow(triplet_keys), " unique triplets across all intersections.\n"))

  # 2. Extract data for these triplets (All outcomes: mvAge, RA, HZ)
  target_outcomes <- c("mvAge", "RA", "HZ (UKB)")
  
  final_data <- data %>%
    dplyr::filter(outcome %in% target_outcomes) %>%
    inner_join(triplet_keys, by = c("gene_name", "cell", "lead_snp"))
    
  # Deduplicate
  final_data <- final_data %>% distinct(gene_name, cell, lead_snp, outcome, .keep_all = TRUE)

  # 3. Generate Integrated Forest Plot
  cat("  Generating Integrated Forest Plot...\n")
  
  p_forest <- run_analysis(
    suffix = "Integrated_Cell_Gene_SNP_Forest",
    custom_data = final_data,
    subdir = dir_added_cell_gene_snp,
    show_beta_allele = TRUE,
    color_cell_only = TRUE
  )
  
  # 4. Generate Specific Correlation Plots
  correlation_plots <- list()
  
  for (spec in pair_specs) {
    pair <- spec$pair
    p_name <- spec$name
    corr_filters <- spec$corr_filters
    
    # Re-calculate intersection for this specific pair
    pair_data <- data %>% dplyr::filter(outcome %in% pair)
    pair_counts <- pair_data %>%
      group_by(gene_name, cell, lead_snp) %>%
      summarize(n_out = n_distinct(outcome), .groups = "drop") %>%
      dplyr::filter(n_out == 2)
      
    if (nrow(pair_counts) > 0) {
      pair_intersection_data <- pair_data %>%
        inner_join(pair_counts[, c("gene_name", "cell", "lead_snp")], by = c("gene_name", "cell", "lead_snp"))
        
      for (ftype in corr_filters) {
        p_cor <- generate_correlation_plot(
          pair_intersection_data, pair[1], pair[2],
          dir_added_cell_gene_snp, 
          filter_type = ftype
        )
        
        if (!is.null(p_cor)) {
          key <- paste0(p_name, "_", ftype)
          correlation_plots[[key]] <- p_cor
        }
      }
    }
  }

  # 5. Combine Plots
  if (!is.null(p_forest) && length(correlation_plots) > 0) {
    cat("  Combining Plots...\n")
    
    p_forest_grob <- wrap_elements(p_forest)
    
    # Order: mvAge-RA(Neg), mvAge-HZ(Pos), mvAge-HZ(Neg), RA-HZ(Pos), RA-HZ(Neg)
    plot_list <- list()
    if (!is.null(correlation_plots[["mvAge_vs_RA_negative"]])) plot_list[[length(plot_list)+1]] <- correlation_plots[["mvAge_vs_RA_negative"]]
    if (!is.null(correlation_plots[["mvAge_vs_HZ_positive"]])) plot_list[[length(plot_list)+1]] <- correlation_plots[["mvAge_vs_HZ_positive"]]
    if (!is.null(correlation_plots[["mvAge_vs_HZ_negative"]])) plot_list[[length(plot_list)+1]] <- correlation_plots[["mvAge_vs_HZ_negative"]]
    if (!is.null(correlation_plots[["RA_vs_HZ_positive"]])) plot_list[[length(plot_list)+1]] <- correlation_plots[["RA_vs_HZ_positive"]]
    if (!is.null(correlation_plots[["RA_vs_HZ_negative"]])) plot_list[[length(plot_list)+1]] <- correlation_plots[["RA_vs_HZ_negative"]]
    
    if (length(plot_list) > 0) {
      # Create right panel with collected guides
      right_panel <- wrap_plots(plot_list, ncol = 1) + 
        plot_layout(guides = "collect") & 
        theme(legend.position = "right")
      
      # Calculate dynamic height
      n_genes <- length(unique(final_data$gene_name))
      n_rows_data <- nrow(final_data)
      estimated_height <- (n_genes + n_rows_data) * 0.3 + 3
      total_height <- max(25, estimated_height) # Ensure at least 25
      
      # Tagged Version
      combined_plot_tagged <- (p_forest_grob | right_panel) +
        plot_layout(widths = c(2.5, 1)) + 
        plot_annotation(tag_levels = "A") &
        theme(legend.position = "right")
        
      # Untagged Version
      combined_plot_untagged <- (p_forest_grob | right_panel) +
        plot_layout(widths = c(2.5, 1)) &
        theme(legend.position = "right")
      
      # Save
      ggsave(file.path(dir_added_cell_gene_snp, "Integrated_Combined_Plot_Tagged.pdf"), combined_plot_tagged, width = 40, height = total_height, limitsize = FALSE)
      ggsave(file.path(dir_added_cell_gene_snp, "Integrated_Combined_Plot_Tagged.png"), combined_plot_tagged, width = 40, height = total_height, dpi = 300, limitsize = FALSE)
      ggsave(file.path(dir_added_cell_gene_snp, "Integrated_Combined_Plot_Untagged.pdf"), combined_plot_untagged, width = 40, height = total_height, limitsize = FALSE)
      ggsave(file.path(dir_added_cell_gene_snp, "Integrated_Combined_Plot_Untagged.png"), combined_plot_untagged, width = 40, height = total_height, dpi = 300, limitsize = FALSE)
      
      cat("  Integrated Combined plots saved.\n")
    }
  }
}


# ==============================================================================
# 6. Statistics Analysis / 
# ==============================================================================

cat("\nGenerating Statistics Analysis... / ...\n")

# 10.1 Basic Counts per Outcome / 
stats_counts <- data %>%
  group_by(outcome) %>%
  summarise(
    N_Genes = n_distinct(gene_name, na.rm = TRUE),
    N_SNPs = n_distinct(lead_snp, na.rm = TRUE),
    N_Gene_SNP_Pairs = n_distinct(paste(gene_name, lead_snp), na.rm = TRUE),
    N_Cell_Gene_Pairs = n_distinct(paste(cell, gene_name), na.rm = TRUE),
    N_Cell_Gene_SNP_Triplets = n_distinct(paste(cell, gene_name, lead_snp), na.rm = TRUE),
    N_Unique_Cells = n_distinct(cell, na.rm = TRUE)
  )

print("Basic Counts per Outcome:")
print(stats_counts)
write.csv(stats_counts, file.path(dir_statistics, "Statistics_Basic_Counts.csv"), row.names = FALSE)

# Visualization 1: Basic Counts Bar Plot
tryCatch({
  stats_counts_long <- stats_counts %>%
    pivot_longer(cols = c(N_Genes, N_SNPs, N_Gene_SNP_Pairs, N_Cell_Gene_Pairs, N_Cell_Gene_SNP_Triplets), names_to = "Type", values_to = "Count")
  
  custom_colors <- c(
    "N_Genes" = "#E64B35",                  # Red
    "N_SNPs" = "#3C5488",                   # Dark Blue
    "N_Gene_SNP_Pairs" = "#F39B7F",         # Orange
    "N_Cell_Gene_Pairs" = "#4DBBD5",        # Blue
    "N_Cell_Gene_SNP_Triplets" = "#00A087"  # Green
  )
  
  # Define Factor Levels for Order
  stats_counts_long$Type <- factor(stats_counts_long$Type, levels = c("N_Genes", "N_SNPs", "N_Gene_SNP_Pairs", "N_Cell_Gene_Pairs", "N_Cell_Gene_SNP_Triplets"))
  stats_counts_long$outcome <- factor(stats_counts_long$outcome, levels = c("mvAge", "RA", "HZ (UKB)"))
  
  p_counts <- ggplot(stats_counts_long, aes(x = outcome, y = Count, fill = Type)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.5) +
    geom_text(aes(label = Count), position = position_dodge(width = 0.9), vjust = -0.5, size = 3, fontface = "bold") +
    theme_classic(base_size = 14) + # Journal style theme
    labs(x = NULL, y = "Count") + # Remove Title and X-axis label (Outcome)
    scale_y_continuous(expand = expansion(mult = c(0, 0.1))) + # Reduce gap at bottom
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 12, color = "black", face = "bold"),
      axis.text.y = element_text(size = 12, color = "black"),
      axis.title.y = element_text(size = 14, face = "bold"),
      legend.position = "top",
      legend.title = element_blank(),
      legend.text = element_text(size = 10),
      panel.grid.major.y = element_line(color = "grey90", linetype = "dashed"),
      axis.line = element_line(size = 0.8)
    ) +
    scale_fill_manual(values = custom_colors, 
                      labels = c("Genes", "SNPs", "Gene-SNP Pairs", "Cell-Gene Pairs", "Cell-Gene-SNP Triplets"))
  
  ggsave(file.path(dir_statistics, "BarPlot_Counts_Per_Outcome.pdf"), p_counts, width = 10, height = 6)
  ggsave(file.path(dir_statistics, "BarPlot_Counts_Per_Outcome.png"), p_counts, width = 10, height = 6, dpi = 300)
}, error = function(e) {
  cat("Error generating counts plot: ", e$message, "\n")
})

# 10.2 Gene-Cell Level Statistics / -
cat("Generating Gene-Cell Level Statistics...\n")
stats_gene_cell <- data %>%
  group_by(gene_name, cell) %>%
  summarise(
    N_Outcomes = n_distinct(outcome, na.rm = TRUE),
    Outcomes = paste(unique(outcome), collapse = "; "),
    .groups = "drop"
  )

write.csv(stats_gene_cell, file.path(dir_statistics, "Statistics_Gene_Cell_Level.csv"), row.names = FALSE)

# 10.3 Pairwise Overlaps / 
cat("Generating Pairwise Overlaps...\n")
outcomes <- unique(data$outcome)
outcomes <- outcomes[!is.na(outcomes)]
overlap_results <- list()

if (length(outcomes) >= 2) {
  pairs <- combn(outcomes, 2, simplify = FALSE)
  
  for (pair in pairs) {
    out1 <- pair[1]
    out2 <- pair[2]
    
    # Gene Overlap
    genes1 <- unique(data$gene_name[data$outcome == out1])
    genes2 <- unique(data$gene_name[data$outcome == out2])
    genes_shared <- intersect(genes1, genes2)
    
    # SNP Overlap
    snps1 <- unique(data$lead_snp[data$outcome == out1])
    snps2 <- unique(data$lead_snp[data$outcome == out2])
    snps_shared <- intersect(snps1, snps2)
    
    # Gene-SNP Pair Overlap
    gs1 <- unique(paste(data$gene_name[data$outcome == out1], data$lead_snp[data$outcome == out1]))
    gs2 <- unique(paste(data$gene_name[data$outcome == out2], data$lead_snp[data$outcome == out2]))
    gs_shared <- intersect(gs1, gs2)
    
    overlap_results[[length(overlap_results) + 1]] <- data.frame(
      Outcome_1 = out1,
      Outcome_2 = out2,
      N_Shared_Genes = length(genes_shared),
      N_Shared_SNPs = length(snps_shared),
      N_Shared_Gene_SNP_Pairs = length(gs_shared),
      stringsAsFactors = FALSE
    )
  }
  
  overlap_df <- do.call(rbind, overlap_results)
  write.csv(overlap_df, file.path(dir_statistics, "Statistics_Pairwise_Overlaps.csv"), row.names = FALSE)
  
  cat("Pairwise Overlaps Generated.\n")
  
  # Visualization 2: Overlap Heatmap
  tryCatch({
    # Prepare data for heatmap (Gene-SNP Pairs)
    # To make a full symmetric matrix for heatmap, we need both (A,B) and (B,A) and (A,A)
    
    # Initialize matrix
    mat_data <- matrix(0, nrow = length(outcomes), ncol = length(outcomes))
    rownames(mat_data) <- outcomes
    colnames(mat_data) <- outcomes
    
    # Fill diagonal
    for (out in outcomes) {
      val <- stats_counts$N_Gene_SNP_Pairs[stats_counts$outcome == out]
      if(length(val) > 0) mat_data[out, out] <- val
    }
    
    # Fill off-diagonal
    for (i in 1:nrow(overlap_df)) {
      o1 <- overlap_df$Outcome_1[i]
      o2 <- overlap_df$Outcome_2[i]
      val <- overlap_df$N_Shared_Gene_SNP_Pairs[i]
      mat_data[o1, o2] <- val
      mat_data[o2, o1] <- val
    }
    
    # Convert to long format for ggplot
    mat_long <- as.data.frame(as.table(mat_data))
    colnames(mat_long) <- c("Outcome1", "Outcome2", "Count")
    
    p_heatmap <- ggplot(mat_long, aes(x = Outcome1, y = Outcome2, fill = Count)) +
      geom_tile(color = "white") +
      geom_text(aes(label = Count), color = "black", size = 4) +
      scale_fill_gradient(low = "#ebf3fb", high = "#377eb8") +
      theme_minimal() +
      labs(title = "Overlap of Gene-SNP Pairs between Outcomes",
           subtitle = "Diagonal: Total pairs; Off-diagonal: Shared pairs",
           x = "", y = "") +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
        axis.text.y = element_text(size = 10),
        plot.title = element_text(face = "bold", hjust = 0.5)
      )
    
    ggsave(file.path(dir_statistics, "Heatmap_Shared_Gene_SNP_Pairs.pdf"), p_heatmap, width = 8, height = 7)
    ggsave(file.path(dir_statistics, "Heatmap_Shared_Gene_SNP_Pairs.png"), p_heatmap, width = 8, height = 7, dpi = 300)
    
  }, error = function(e) {
    cat("Error generating heatmap: ", e$message, "\n")
  })
}

# 8.3 Multi-Outcome Venn Diagram (Try ggVennDiagram) / 
target_outcomes <- c("mvAge", "HZ (UKB)", "RA") # Updated to use abbreviations
present_targets <- intersect(target_outcomes, outcomes)

if (length(present_targets) >= 2) {
  cat("  Generating Venn Diagrams for: ", paste(present_targets, collapse=", "), "...\n")
  
  if (requireNamespace("ggVennDiagram", quietly = TRUE)) {
    tryCatch({
      library(ggVennDiagram)
      
      # Colors for Venn Levels (matching BarPlot)
      color_g <- "#E64B35"    # Genes
      color_s <- "#3C5488"    # SNPs
      color_gs <- "#F39B7F"   # Gene-SNP Pairs
      color_cg <- "#4DBBD5"   # Cell-Gene Pairs
      color_cgs <- "#00A087"  # Triplets
      
      # 1. Gene Level (N_Genes)
      venn_list_g <- list()
      # Explicitly set order: mvAge, RA, HZ (UKB)
      ordered_targets <- c("mvAge", "RA", "HZ (UKB)")
      ordered_targets <- ordered_targets[ordered_targets %in% present_targets]
      
      for(out in ordered_targets) venn_list_g[[out]] <- unique(data$gene_name[data$outcome == out])
      p_venn_g <- ggVennDiagram(venn_list_g, label_alpha = 0, category.names = ordered_targets) + 
        scale_fill_gradient(low = "white", high = color_g) + labs(title = NULL, subtitle = NULL) + theme(legend.position = "none")
      ggsave(file.path(dir_statistics, "Venn_Genes.pdf"), p_venn_g, width = 6, height = 6)
      ggsave(file.path(dir_statistics, "Venn_Genes.png"), p_venn_g, width = 6, height = 6, dpi = 300)

      # 2. SNP Level (N_SNPs)
      venn_list_s <- list()
      for(out in ordered_targets) venn_list_s[[out]] <- unique(data$lead_snp[data$outcome == out])
      p_venn_s <- ggVennDiagram(venn_list_s, label_alpha = 0, category.names = ordered_targets) + 
        scale_fill_gradient(low = "white", high = color_s) + labs(title = NULL, subtitle = NULL) + theme(legend.position = "none")
      ggsave(file.path(dir_statistics, "Venn_SNPs.pdf"), p_venn_s, width = 6, height = 6)
      ggsave(file.path(dir_statistics, "Venn_SNPs.png"), p_venn_s, width = 6, height = 6, dpi = 300)

      # 3. Gene-SNP Pair Level (N_Gene_SNP_Pairs)
      venn_list_gs <- list()
      for(out in ordered_targets) venn_list_gs[[out]] <- unique(paste(data$gene_name[data$outcome == out], data$lead_snp[data$outcome == out]))
      p_venn_gs <- ggVennDiagram(venn_list_gs, label_alpha = 0, category.names = ordered_targets) + 
        scale_fill_gradient(low = "white", high = color_gs) + labs(title = NULL, subtitle = NULL) + theme(legend.position = "none")
      ggsave(file.path(dir_statistics, "Venn_Gene_SNP_Pairs.pdf"), p_venn_gs, width = 6, height = 6)
      ggsave(file.path(dir_statistics, "Venn_Gene_SNP_Pairs.png"), p_venn_gs, width = 6, height = 6, dpi = 300)

      # 4. Cell-Gene Level (N_Cell_Gene_Pairs)
      venn_list_cg <- list()
      for(out in ordered_targets) venn_list_cg[[out]] <- unique(paste(data$cell[data$outcome == out], data$gene_name[data$outcome == out]))
      p_venn_cg <- ggVennDiagram(venn_list_cg, label_alpha = 0, category.names = ordered_targets) + 
        scale_fill_gradient(low = "white", high = color_cg) + labs(title = NULL, subtitle = NULL) + theme(legend.position = "none")
      ggsave(file.path(dir_statistics, "Venn_Cell_Gene_Pairs.pdf"), p_venn_cg, width = 6, height = 6)
      ggsave(file.path(dir_statistics, "Venn_Cell_Gene_Pairs.png"), p_venn_cg, width = 6, height = 6, dpi = 300)

      # 5. Cell-Gene-SNP Triplet Level (N_Cell_Gene_SNP_Triplets)
      venn_list_cgs <- list()
      for(out in ordered_targets) venn_list_cgs[[out]] <- unique(paste(data$cell[data$outcome == out], data$gene_name[data$outcome == out], data$lead_snp[data$outcome == out]))
      p_venn_cgs <- ggVennDiagram(venn_list_cgs, label_alpha = 0, category.names = ordered_targets) + 
        scale_fill_gradient(low = "white", high = color_cgs) + labs(title = NULL, subtitle = NULL) + theme(legend.position = "none")
      ggsave(file.path(dir_statistics, "Venn_Cell_Gene_SNP_Triplets.pdf"), p_venn_cgs, width = 6, height = 6)
      ggsave(file.path(dir_statistics, "Venn_Cell_Gene_SNP_Triplets.png"), p_venn_cgs, width = 6, height = 6, dpi = 300)

      cat("  Venn diagrams (5 Levels) generated.\n")
      
      # --- Combined Plot Generation ---
      cat("  Generating Combined Plot...\n")
      if (requireNamespace("patchwork", quietly = TRUE) && exists("p_counts")) {
        library(patchwork)
        
        # Combine: BarPlot on top (A), Venns row below (B)
        # Wrap Venns to be a single element so they share the tag "B"
        p_venns_all <- wrap_elements(p_venn_g | p_venn_s | p_venn_gs | p_venn_cg | p_venn_cgs)
        
        # Design to narrow the top plot (BarPlot)
        # Using area() to be robust against string formatting issues
        # Row 1: A in middle (cols 2-6 of 7)
        # Row 2: B full width (cols 1-7 of 7)
        design <- c(area(1, 2, 1, 6), area(2, 1, 2, 7))
        
        # Version 1: No Labels (Maintain style and layout, just hide labels)
        p_combined_nolabels <- p_counts + p_venns_all +
          plot_layout(design = design, heights = c(1.2, 1))
          # No plot_annotation(tag_levels = 'A')
        
        ggsave(file.path(dir_statistics, "Combined_Statistics_Plot_NoLabels.pdf"), p_combined_nolabels, width = 20, height = 10)
        ggsave(file.path(dir_statistics, "Combined_Statistics_Plot_NoLabels.png"), p_combined_nolabels, width = 20, height = 10, dpi = 300)
        cat("  Combined plot (No Labels) saved successfully.\n")

        # Version 2: Aligned Labels (A and B aligned left)
        # To align A and B while keeping the visual "narrowness" of A:
        # We wrap A in a row with spacers, so the whole row is treated as "A" starting at x=0
        # Layout: Spacer + Plot + Spacer (1:5:1 ratio roughly matches 2-6/7)
        p_counts_centered <- wrap_elements(plot_spacer() + p_counts + plot_spacer() + plot_layout(widths = c(1, 5, 1)))
        
        p_combined_aligned <- p_counts_centered / p_venns_all +
          plot_layout(heights = c(1.2, 1)) +
          plot_annotation(tag_levels = 'A') & 
          theme(plot.tag = element_text(size = 18, face = "bold"))
        
        ggsave(file.path(dir_statistics, "Combined_Statistics_Plot_Aligned.pdf"), p_combined_aligned, width = 20, height = 10)
        ggsave(file.path(dir_statistics, "Combined_Statistics_Plot_Aligned.png"), p_combined_aligned, width = 20, height = 10, dpi = 300)
        cat("  Combined plot (Aligned Labels) saved successfully.\n")
        
      } else {
        cat("  'patchwork' package not found or p_counts missing. Skipping combined plot.\n")
      }

    }, error = function(e) {
      cat("  Error creating Venn diagrams or combined plot: ", e$message, "\n")
    })
  } else {
    cat("  'ggVennDiagram' package not found. Skipping Venn diagram generation.\n")
  }
}

# ==============================================================================
# 11. Perform Enrichment Analysis / 
# ==============================================================================
cat("\nPerforming Enrichment Analysis...\n")

# 11.1 Single Outcome Enrichment / 
cat("Processing Single Outcome Enrichment...\n")
for (out in unique_outcomes) {
  out_safe <- gsub("[^[:alnum:]]", "_", out)
  genes <- unique(data$gene_name[data$outcome == out])
  genes <- genes[!is.na(genes)]
  
  if (length(genes) > 0) {
    perform_enrichment(genes, file.path(dir_enrichment, "Single_Outcome", out_safe), paste0("Single_Outcome_", out_safe))
  }
}

# 11.2 Shared Enrichment (Pairwise and 3-Way)
cat("Processing Shared Enrichment (Pairwise and 3-Way)...\n")

# Define pairs and triplets
combinations_list <- list(
  c("mvAge", "HZ (UKB)"),
  c("mvAge", "RA"),
  c("HZ (UKB)", "RA"),
  c("mvAge", "HZ (UKB)", "RA")
)

for (combo in combinations_list) {
  combo_name <- paste(combo, collapse = "_vs_")
  combo_name_safe <- gsub("[^[:alnum:]]", "_", combo_name)
  combo_name_safe <- gsub(" ", "", combo_name_safe) # Remove spaces
  
  # Check if all outcomes in combo exist in data
  if (all(combo %in% outcomes)) {
    cat(paste0("  Enrichment for: ", combo_name, "...\n"))
    
    # Filter data for each outcome
    gene_sets <- list()
    for (out in combo) {
      gene_sets[[out]] <- unique(data$gene_name[data$outcome == out])
    }
    
    # Find intersection
    shared_genes <- Reduce(intersect, gene_sets)
    shared_genes <- shared_genes[!is.na(shared_genes)]
    
    if (length(shared_genes) > 0) {
      perform_enrichment(shared_genes, file.path(dir_enrichment, "Shared", combo_name_safe), paste0("Shared_", combo_name_safe))
    } else {
        cat(paste0("    No shared genes found for: ", combo_name, "\n"))
    }
  }
}

# ==============================================================================
# 12. Generate Scientific log
# ==============================================================================
cat("\nGenerating Scientific log...\n")

# Prepare Interpretation Variables
n_mvAge_genes <- stats_counts$N_Genes[stats_counts$outcome == "mvAge"]
n_hz_genes <- stats_counts$N_Genes[stats_counts$outcome == "HZ (UKB)"]
n_ra_genes <- stats_counts$N_Genes[stats_counts$outcome == "RA"]

if (length(n_mvAge_genes) == 0) n_mvAge_genes <- 0
if (length(n_hz_genes) == 0) n_hz_genes <- 0
if (length(n_ra_genes) == 0) n_ra_genes <- 0

# Shared Genes (mvAge vs HZ vs RA)
mvAge_hz_ra_shared_genes <- intersect(
  intersect(data$gene_name[data$outcome == "mvAge"], data$gene_name[data$outcome == "HZ (UKB)"]),
  data$gene_name[data$outcome == "RA"]
)
n_mvAge_hz_ra_shared <- length(mvAge_hz_ra_shared_genes)

# ==============================================================================
# 14. Validation Analysis / 
# ==============================================================================
cat("\nProcessing Validation Analysis...\n")

if (nrow(data_validation) > 0) {
  
  # 14.1 mvAge (Main) vs Frailty Index (Validation)
  cat("  14.1 Checking mvAge (Main) vs Frailty Index (Validation)...\n")
  
  # Extract mvAge from Main Data
  mvAge_main <- data %>% filter(outcome == "mvAge")
  
  # Extract Frailty from Validation Data
  frailty_val <- data_validation %>% filter(outcome == "Frailty index")
  
  if (nrow(mvAge_main) > 0 && nrow(frailty_val) > 0) {
    perform_validation_analysis_module(
      df_main = mvAge_main,
      df_val = frailty_val,
      outcome_main = "mvAge",
      outcome_val = "Frailty Index",
      module_name = "Validation_mvAge_vs_Frailty",
      output_dir = dir_validation,
      use_or_for_correlation = TRUE
    )
  } else {
    cat("    Missing mvAge (Main) or Frailty Index (Validation) data.\n")
  }
  
  # 14.2 Herpes Zoster (Main) vs Herpes Zoster (FinnGen/Validation) - REMOVED AS REQUESTED
  # cat("  14.2 Checking Herpes Zoster (Main) vs Herpes Zoster (FinnGen)...\n")
  # Module removed.

  
} else {
  cat("  No validation data available. Skipping validation analysis.\n")
}

# ==============================================================================
# 15. Save Script Backup / 
# ==============================================================================
cat("\nBacking up script...\n")
script_path <- commandArgs(trailingOnly = FALSE)[4]
if (grepl("--file=", script_path)) {
  script_path <- sub("--file=", "", script_path)
} else {
  script_path <- "./3.1.1.1_forest_plot_script_template.R"
}

if (file.exists(script_path)) {
  file.copy(script_path, file.path(dir_scripts, basename(script_path)))
  cat("Script backed up to output directory. / . \n")
}

# End logging / 
sink()
