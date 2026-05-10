# ==============================================================================
# [Script]: calculate_LD_clustering_all(MHC).R
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
# Script Name: calculate_LD_clustering_all.R
# Description: 
#   This script performs LD-based clustering for Exposure (Cell-specific), Outcome (Phenotype-specific),
#   and Global (Outcome-agnostic) signals.
#   ,    LD.
#
#   Definition of "Same Signal" / :
#   -----------------------------------------------------------------------------------------
#   1. Exposure Side (Exposure Clusters):
#      - Criteria: Same Cell Type + SNPs in LD (R2 >= 0.6).
#      - Logic: Associations from the same cell type, where lead SNPs are in high LD.
#      - ID Format: Exposure_Cluster_X
#
#   2. Outcome Side (Global Clusters / Global Signals):
#      - Criteria: SNPs in LD (R2 >= 0.6), regardless of Cell Type or Outcome Phenotype.
#      - Logic: Groups all SNPs based purely on genetic correlation (LD).
#      - ID Format: Global_Cluster_Y
#   -----------------------------------------------------------------------------------------
#
#   Workflow / :
#     1. Data Loading & Filtering: Load joint Forest Plot data, filter for MHC region (chr6:25-34Mb).
#     2. Global LD Calculation: Calculate pairwise R2 for all unique SNPs using PLINK.
#     3. Exposure Analysis: Group SNPs by Cell Type, cluster by R2 >= 0.6.
#     4. Outcome Analysis: Group SNPs by Outcome Phenotype, cluster by R2 >= 0.6.
#     5. Global Clustering: Cluster all SNPs together (Outcome-agnostic) by R2 >= 0.6.
#     6. Visualization: Generate Heatmaps, Venn/Bar charts, and updated Sankey Diagram.
#
# Usage: 
#   Rscript calculate_LD_clustering_all.R
#
# Dependencies:
#   - R packages: data.table, dplyr, pheatmap, RColorBrewer, grid, gridExtra, tools, glue, SNPlocs.Hsapiens.dbSNP155.GRCh37, ggplot2, VennDiagram, ggalluvial
#   - External tools: PLINK 1.9
# ==============================================================================

# ==============================================================================
# Configuration / 
# ==============================================================================

# Clear environment / 
rm(list = ls())
gc()

# Libraries / 
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(pheatmap)
  library(RColorBrewer)
  library(grid)
  library(gridExtra)
  library(tools)
  library(glue)
  library(SNPlocs.Hsapiens.dbSNP155.GRCh37)
  library(ggplot2)
  library(VennDiagram)
  if(require(ggalluvial)) {
    library(ggalluvial)
  } else {
    warning("Package 'ggalluvial' is not installed. Sankey diagram might be skipped or simplified.")
  }
})

# ==============================================================================
# 0. Command Line Arguments / 
# ==============================================================================
option_list <- list(
  make_option(c("--input_csv"), type="character", default=NULL,
              help="Path to the input joint Forest Plot CSV file / CSV"),
  make_option(c("--plink_path"), type="character", default=NULL,
              help="Path to PLINK executable / PLINK"),
  make_option(c("--ld_ref"), type="character", default=NULL,
              help="Prefix of the LD reference panel (hg19) / LD (hg19)"),
  make_option(c("--out_dir"), type="character", default="./LD_Clustering_Results",
              help="Output directory path / ")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_csv) || is.null(opt$plink_path) || is.null(opt$ld_ref)) {
  print_help(opt_parser)
  stop("Input CSV, PLINK path, and LD Reference Panel must be provided. / CSV, PLINKLD. ")
}

# Paths / 
input_file <- opt$input_csv
plink_path <- opt$plink_path
ld_ref_prefix <- opt$ld_ref
base_root_dir <- opt$out_dir

# Timestamped Output Directory / 
timestamp_folder <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- file.path(base_root_dir, timestamp_folder)

# Ensure output directory exists
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Logger Function / 
log_file <- file.path(output_dir, "execution.log")
log_msg <- function(msg) {
  timestamp_log <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  message(paste(timestamp_log, msg))
  cat(paste(timestamp_log, msg, "\n"), file = log_file, append = TRUE)
}

log_msg(glue("Starting Analysis Script. Output directory: {output_dir}"))

# ==============================================================================
# 1. Data Processing & Filtering / 
# ==============================================================================

log_msg("Reading and processing input data...")

# Read data
if (!file.exists(input_file)) {
  log_msg(glue("Error: Input file not found: {input_file}"))
  stop("Input file missing.")
}
dt <- fread(input_file, header = TRUE, fill = TRUE)

# Standardize Outcome Names: Rename 'aging' to 'mvAge'
# ------------------------------------------------------------------------------
if ("Outcome" %in% names(dt)) {
  dt[Outcome == "aging", Outcome := "mvAge"]
  log_msg("Renamed 'aging' to 'mvAge' in Outcome column.")
}
if ("Outcome_Phenotype" %in% names(dt)) {
  dt[Outcome_Phenotype == "aging", Outcome_Phenotype := "mvAge"]
}


# Process Gene Information from Summary Rows (Hierarchical structure)
dt[, Gene := NA_character_]
current_gene <- NA
for (i in 1:nrow(dt)) {
  is_summary_val <- dt$is_summary[i]
  if (!is.na(is_summary_val) && (is_summary_val == TRUE || is_summary_val == "TRUE")) {
    raw_name <- dt$Display_Name[i]
    gene_name <- sub(" \\(.*", "", raw_name)
    current_gene <- gene_name
  }
  dt$Gene[i] <- current_gene
}

# Filter Data: Keep rows with valid Lead_SNP
dt_filtered <- dt[Lead_SNP != "" & !is.na(Lead_SNP)]
dt_clean <- dt_filtered

log_msg(glue("Total data rows with SNPs before location filtering: {nrow(dt_clean)}"))
log_msg(glue("Unique SNPs before location filtering: {length(unique(dt_clean$Lead_SNP))}"))

# 1.5. SNP Location Filtering (MHC Region: chr6:25Mb-34Mb) / SNP
# ==============================================================================
unique_snps_initial <- unique(dt_clean$Lead_SNP)
log_msg(glue("Filtering SNPs by location (chr6: 25Mb-34Mb) using SNPlocs.Hsapiens.dbSNP155.GRCh37..."))

tryCatch({
  # Query SNP locations
  snps_gpos <- snpsById(SNPlocs.Hsapiens.dbSNP155.GRCh37, ids = unique_snps_initial, ifnotfound = "drop")
  snps_df <- as.data.frame(snps_gpos)
  
  # Filter for chr6 and position 25Mb-34Mb
  # SNPlocs usually uses "ch" prefix or just numbers. We check seqnames.
  mhc_snps_df <- snps_df %>%
    filter(seqnames == "ch6" | seqnames == "6") %>%
    filter(pos >= 25000000 & pos <= 34000000)
  
  valid_mhc_snps <- mhc_snps_df$RefSNP_id
  
  log_msg(glue("  Found {length(valid_mhc_snps)} SNPs in MHC region (chr6:25-34Mb)."))
  log_msg(glue("  Dropped {length(unique_snps_initial) - length(valid_mhc_snps)} SNPs outside region or not found."))
  
  # Filter main data
  dt_clean <- dt_clean[Lead_SNP %in% valid_mhc_snps]
  
  log_msg(glue("Total data rows after location filtering: {nrow(dt_clean)}"))
  log_msg(glue("Unique SNPs after location filtering: {length(unique(dt_clean$Lead_SNP))}"))
  
}, error = function(e) {
  log_msg(glue("Error during SNP location filtering: {e$message}"))
  stop("SNP filtering failed.")
})

if (nrow(dt_clean) == 0) {
  log_msg("No data remaining after filtering. Exiting.")
  quit(save = "no")
}

# ==============================================================================
# 2. Global LD Calculation / LD
# ==============================================================================

unique_snps <- unique(dt_clean$Lead_SNP)
unique_snps <- unique_snps[unique_snps != ""]

log_msg(glue("Calculating Global LD matrix for {length(unique_snps)} unique SNPs..."))

# Write SNP list for PLINK
snp_list_file <- file.path(output_dir, "all_snps_list.txt")
writeLines(unique_snps, snp_list_file)

plink_out_prefix <- file.path(output_dir, "global_ld")

# Run PLINK
# --r2 to get pairwise R2
# --ld-window-r2 0 to get all pairs (even low LD)
cmd_plink <- glue(
  "'{plink_path}'",
  " --bfile '{ld_ref_prefix}'",
  " --extract '{snp_list_file}'",
  " --r2",
  " --ld-window-r2 0", 
  " --ld-window 99999",
  " --ld-window-kb 99999",
  " --out '{plink_out_prefix}'",
  " --threads 4",
  " --silent"
)

log_msg("Running PLINK...")
system(cmd_plink, ignore.stdout = TRUE)

ld_file <- paste0(plink_out_prefix, ".ld")
if (!file.exists(ld_file)) {
  stop("PLINK failed to generate LD file.")
}

log_msg("Reading LD results...")
ld_data <- fread(ld_file)

# Construct Full Matrix
all_snps_in_ld <- unique(c(ld_data$SNP_A, ld_data$SNP_B))
missing_snps <- setdiff(unique_snps, all_snps_in_ld)
if (length(missing_snps) > 0) {
  log_msg(glue("Warning: {length(missing_snps)} SNPs not found in LD reference panel."))
}
all_snps_final <- c(all_snps_in_ld, missing_snps)
n_snps <- length(all_snps_final)

log_msg(glue("Constructing {n_snps} x {n_snps} LD Matrix..."))

mat <- matrix(0, nrow = n_snps, ncol = n_snps)
rownames(mat) <- all_snps_final
colnames(mat) <- all_snps_final
diag(mat) <- 1

# Fill matrix
idx_a <- match(ld_data$SNP_A, all_snps_final)
idx_b <- match(ld_data$SNP_B, all_snps_final)
mat[cbind(idx_a, idx_b)] <- ld_data$R2
mat[cbind(idx_b, idx_a)] <- ld_data$R2

# Save Global Matrix
matrix_file <- file.path(output_dir, "Global_LD_Matrix.csv")
write.csv(mat, matrix_file, quote = FALSE)
log_msg(glue("Global LD Matrix saved to {matrix_file}"))

# Clean up PLINK files
unlink(c(snp_list_file, paste0(plink_out_prefix, ".ld"), paste0(plink_out_prefix, ".log"), paste0(plink_out_prefix, ".nosex")))

# ==============================================================================
# 3. Clustering Function Definition / 
# ==============================================================================

perform_clustering <- function(data, group_col, analysis_name, output_subdir, global_mat, r2_threshold = 0.6) {
  # Create output subdirectory
  analysis_dir <- file.path(output_dir, output_subdir)
  if (!dir.exists(analysis_dir)) dir.create(analysis_dir, recursive = TRUE)
  
  log_msg(glue("Starting {analysis_name} Analysis..."))
  if (!is.null(group_col)) {
    log_msg(glue("  Grouping by: {group_col}"))
  } else {
    log_msg("  Grouping by: Global (None)")
  }
  log_msg(glue("  Clustering Threshold: R2 >= {r2_threshold}"))
  
  # Prepare result container
  cluster_results <- data.table()
  
  # Determine groups
  if (!is.null(group_col)) {
    groups <- unique(data[[group_col]])
    groups <- groups[!is.na(groups) & groups != ""]
  } else {
    groups <- "Global"
  }
  
  combined_plots_list <- list()
  
  for (grp in groups) {
    log_msg(glue("  Processing group: {grp}"))
    
    # Get SNPs for this group
    if (!is.null(group_col)) {
      grp_data <- data[get(group_col) == grp]
      grp_snps <- unique(grp_data$Lead_SNP)
    } else {
      # Global analysis: Use ALL unique SNPs in the dataset
      grp_snps <- unique(data$Lead_SNP)
    }
    
    # Filter SNPs in LD matrix
    valid_snps <- intersect(grp_snps, rownames(global_mat))
    n_grp_snps <- length(valid_snps)
    
    # Initialize cluster assignment for this group
    grp_clusters <- data.table(
      Lead_SNP = grp_snps,
      Group = grp,
      Cluster_ID = NA_character_
    )
    
    # Handle missing SNPs (assign unique cluster ID)
    missing <- setdiff(grp_snps, valid_snps)
    cluster_counter <- 0
    if (length(missing) > 0) {
      for (s in missing) {
        cluster_counter <- cluster_counter + 1
        cid <- paste0(grp, "_Signal_M", cluster_counter) # M for missing LD
        grp_clusters[Lead_SNP == s, Cluster_ID := cid]
      }
    }
    
    # Perform Clustering if we have valid SNPs
    if (n_grp_snps > 0) {
      sub_mat <- global_mat[valid_snps, valid_snps, drop = FALSE]
      
      if (nrow(sub_mat) > 1) {
        # Distance = 1 - R2
        dist_mat <- as.dist(1 - sub_mat)
        hc <- hclust(dist_mat, method = "single")
        # Cut tree at height 1 - threshold
        cut_height <- 1 - r2_threshold
        cut_groups <- cutree(hc, h = cut_height)
      } else {
        cut_groups <- setNames(1, valid_snps)
      }
      
      # Assign Cluster IDs
      unique_cids <- unique(cut_groups)
      for (g in unique_cids) {
        snps_in_c <- names(cut_groups)[cut_groups == g]
        cluster_name <- paste0(grp, "_Signal_", g)
        grp_clusters[Lead_SNP %in% snps_in_c, Cluster_ID := cluster_name]
      }
      
      # Visualization (Heatmap)
      if (n_grp_snps >= 2) {
        grp_out_dir <- file.path(analysis_dir, make.names(grp)) # Ensure valid folder name
        if (!dir.exists(grp_out_dir)) dir.create(grp_out_dir, recursive = TRUE)
        
        # Save Sub-Matrix
        write.csv(sub_mat, file.path(grp_out_dir, "LD_Matrix.csv"), quote = FALSE)
        
        # Heatmap settings
        plot_dim <- min(30, max(10, n_grp_snps * 0.2))
        colors <- colorRampPalette(brewer.pal(9, "Reds"))(100)
        
        # PDF Heatmap
        pdf_path <- file.path(grp_out_dir, "LD_Heatmap.pdf")
        tryCatch({
          pdf(pdf_path, width = plot_dim, height = plot_dim)
          pheatmap(sub_mat,
                   color = colors,
                   cluster_rows = TRUE,
                   cluster_cols = TRUE,
                   display_numbers = FALSE,
                   border_color = NA,
                   fontsize_row = max(5, 400/n_grp_snps),
                   fontsize_col = max(5, 400/n_grp_snps),
                   main = paste0("LD (R2) Heatmap - ", grp))
          dev.off()
        }, error = function(e) log_msg(glue("    Error heatmap PDF {grp}: {e$message}")))
        
        # Store for Combined Plot
        combined_plots_list[[grp]] <- list(mat = sub_mat, name = grp, n_snps = n_grp_snps)
      }
    }
    
    # Bind to main results
    cluster_results <- rbind(cluster_results, grp_clusters)
  }
  
  # Save Clustering Results
  fwrite(cluster_results, file.path(analysis_dir, "Clustering_Results.csv"))
  
  # Generate Combined Heatmap if applicable
  if (length(combined_plots_list) > 0) {
    # Sort by number of SNPs
    n_snps_vec <- sapply(combined_plots_list, function(x) x$n_snps)
    sorted_list <- combined_plots_list[order(n_snps_vec, decreasing = TRUE)]
    
    pdf_combined <- file.path(analysis_dir, "Combined_Heatmaps.pdf")
    png_combined <- file.path(analysis_dir, "Combined_Heatmaps.png")
    
    # Layout Logic: Use 'Main + Grid' for Exposure, simple grid for others
    if (analysis_name == "Exposure") {
      log_msg("  Generating Exposure Combined Heatmap (Main + Grid layout)...")
      
      grobs_list <- list()
      grobs_list_nolabels <- list()
      labels <- LETTERS[1:length(sorted_list)]
      colors <- colorRampPalette(brewer.pal(9, "Reds"))(100)
      
      for (i in seq_along(sorted_list)) {
        item <- sorted_list[[i]]
        lbl <- labels[i]
        cl_name <- item$name
        mat_data <- item$mat
        
        # Main plot vs Others
        if (i == 1) {
          # Main: Legend, Labels, No Grid
          font_size <- max(1, 200/nrow(mat_data))
          p_obj <- pheatmap(mat_data, color = colors, cluster_rows = TRUE, cluster_cols = TRUE,
                            display_numbers = FALSE, border_color = NA, legend = TRUE,
                            show_rownames = TRUE, show_colnames = TRUE,
                            fontsize_row = font_size, fontsize_col = font_size,
                            main = "", silent = TRUE)
        } else {
          # Others: No Legend, No Labels
          p_obj <- pheatmap(mat_data, color = colors, cluster_rows = TRUE, cluster_cols = TRUE,
                            display_numbers = FALSE, border_color = NA, legend = FALSE,
                            show_rownames = FALSE, show_colnames = FALSE,
                            main = "", silent = TRUE)
        }
        
        heatmap_grob <- p_obj$gtable
        
        # Add labels
        title_grob <- textGrob(label = lbl, x = unit(0.05, "npc"), y = unit(0.5, "npc"), 
                               just = "left", gp = gpar(fontsize = 20, fontface = "bold"))
        cell_grob <- textGrob(label = cl_name, x = unit(0.5, "npc"), y = unit(0.5, "npc"),
                              just = "center", gp = gpar(fontsize = 16, fontface = "bold"))
        
        # Labeled version: Title (A, B...) + Heatmap + Cell Name
        labeled_grob <- arrangeGrob(
          title_grob, heatmap_grob, cell_grob,
          ncol = 1, heights = unit(c(1, 10, 1), c("null", "null", "null"))
        )
        grobs_list[[i]] <- labeled_grob

        # Unlabeled version: Heatmap + Cell Name (Remove Title)
        unlabeled_grob <- arrangeGrob(
          heatmap_grob, cell_grob,
          ncol = 1, heights = unit(c(10, 1), c("null", "null"))
        )
        grobs_list_nolabels[[i]] <- unlabeled_grob
      }
      
      # Helper to arrange and save
      save_exposure_layout <- function(g_list, p_pdf, p_png, right_cols = NULL) {
        n_plots <- length(g_list)
        if (n_plots == 1) {
          final_plot <- g_list[[1]]
        } else {
          rest_plots <- g_list[-1]
          n_rest <- length(rest_plots)
          
          if (!is.null(right_cols)) {
             n_cols_rest <- right_cols
          } else {
             n_cols_rest <- 2
             if (n_rest > 6) n_cols_rest <- 3
          }

          n_rows_rest <- ceiling(n_rest / n_cols_rest)
          
          if (n_rest < n_cols_rest * n_rows_rest) {
            for (k in (n_rest + 1):(n_cols_rest * n_rows_rest)) rest_plots[[k]] <- nullGrob()
          }
          rest_grid <- arrangeGrob(grobs = rest_plots, ncol = n_cols_rest)
          
          final_plot <- arrangeGrob(g_list[[1]], rest_grid, ncol = 2, widths = c(1.5, 1))
        }
        
        tryCatch({
          pdf(p_pdf, width = 20, height = 15)
          grid.draw(final_plot)
          dev.off()
          
          png(p_png, width = 20*300, height = 15*300, res = 300)
          grid.draw(final_plot)
          dev.off()
          
          log_msg(glue("  Saved: {p_pdf} & {p_png}"))
        }, error = function(e) log_msg(glue("  Error saving combined plots: {e$message}")))
      }
      
      # Save Labeled Version
      save_exposure_layout(grobs_list, pdf_combined, png_combined)
      
      # Save Unlabeled Version
      pdf_combined_no <- file.path(analysis_dir, "Combined_Heatmaps_NoLabels.pdf")
      png_combined_no <- file.path(analysis_dir, "Combined_Heatmaps_NoLabels.png")
      save_exposure_layout(grobs_list_nolabels, pdf_combined_no, png_combined_no)

      # Save Unlabeled Version (4 columns) - Requested by User
      pdf_combined_no_4col <- file.path(analysis_dir, "Combined_Heatmaps_NoLabels_4col.pdf")
      png_combined_no_4col <- file.path(analysis_dir, "Combined_Heatmaps_NoLabels_4col.png")
      save_exposure_layout(grobs_list_nolabels, pdf_combined_no_4col, png_combined_no_4col, right_cols = 4)
      
    } else {
      # Default Grid Layout for Outcome or Global
      tryCatch({
        grobs <- list()
        for (i in seq_along(sorted_list)) {
          item <- sorted_list[[i]]
          p <- pheatmap(item$mat, color = colorRampPalette(brewer.pal(9, "Reds"))(100),
                        cluster_rows = TRUE, cluster_cols = TRUE, border_color = NA,
                        main = paste0(item$name, " (N=", item$n_snps, ")"), silent = TRUE)
          grobs[[i]] <- p$gtable
        }
        
        # Save PDF
        pdf(pdf_combined, width = 20, height = 15)
        marrangeGrob(grobs, nrow = 2, ncol = 2) %>% grid.draw()
        dev.off()
        
        # Save PNG (only first page if multiple, or use %d pattern if needed, but simple for now)
        # marrangeGrob produces a list of grobs (pages) if it spans multiple pages.
        # grid.draw on arrangelist draws the first page.
        # To be safe for single page common case:
        png(png_combined, width = 20*300, height = 15*300, res = 300)
        marrangeGrob(grobs, nrow = 2, ncol = 2) %>% grid.draw()
        dev.off()
        
        log_msg(glue("  Combined Heatmaps saved to {pdf_combined} and {png_combined}"))
      }, error = function(e) log_msg(glue("  Error combining heatmaps: {e$message}")))
    }
  }
  
  return(cluster_results)
}

# ==============================================================================
# 4. Run Analyses / 
# ==============================================================================

# 4.1 Exposure Analysis (Group by Cell Type)
# ------------------------------------------------------------------------------
log_msg("--- Running Exposure Analysis (Cell-grouped) ---")
exposure_clusters <- perform_clustering(
  data = dt_clean,
  group_col = "cell",
  analysis_name = "Exposure",
  output_subdir = "Exposure_Analysis",
  global_mat = mat,
  r2_threshold = 0.6
)
setnames(exposure_clusters, "Cluster_ID", "Exposure_Cluster_ID")
setnames(exposure_clusters, "Group", "Exposure_Cell")

# 4.1.1 Exposure Intersection Analysis (Venn & Stacked Bar Chart)
# ... (Kept as is) ...
# For brevity, reusing the logic from previous version, inserted here.
# NOTE: The user requested to keep functionality.
log_msg("--- Running Exposure Intersection Analysis (Venn & Bar Chart) ---")
tryCatch({
  cluster_outcome_dt <- merge(
    exposure_clusters, 
    dt_clean[, .(Lead_SNP, cell, Outcome)], 
    by.x = c("Lead_SNP", "Exposure_Cell"), 
    by.y = c("Lead_SNP", "cell"),
    all.x = TRUE
  )
  cluster_summary <- cluster_outcome_dt[, .(
    Outcomes = list(sort(unique(Outcome[!is.na(Outcome)])))
  ), by = .(Exposure_Cell, Exposure_Cluster_ID)]
  
  all_outcomes <- sort(unique(dt_clean$Outcome))
  all_outcomes <- all_outcomes[!is.na(all_outcomes) & all_outcomes != ""]
  
  if (length(all_outcomes) > 0) {
    venn_sets <- list()
    for (out in all_outcomes) {
      ids <- cluster_summary[sapply(Outcomes, function(x) out %in% x), Exposure_Cluster_ID]
      # Rename HZ to HZ (UKB) for the plot label
      label_name <- ifelse(out == "HZ", "HZ (UKB)", out)
      venn_sets[[label_name]] <- ids
    }
    
    # Reorder specific outcomes for Venn Diagram layout: mvAge (Top-Left), RA (Top-Right), HZ (Bottom)
    desired_order <- c("mvAge", "RA", "HZ (UKB)")
    if (all(desired_order %in% names(venn_sets)) && length(venn_sets) == 3) {
       venn_sets <- venn_sets[desired_order]
    }

    venn_dir <- file.path(output_dir, "Exposure_Analysis")
    
    if (length(all_outcomes) <= 5) {
      venn_pdf <- file.path(venn_dir, "Cluster_Intersection_Venn.pdf")
      venn_png <- file.path(venn_dir, "Cluster_Intersection_Venn.png")
      my_colors <- brewer.pal(max(3, length(all_outcomes)), "Set2")[1:length(all_outcomes)]
      
      venn.plot <- venn.diagram(x = venn_sets, filename = NULL, output = TRUE, imagetype = "png", 
                                height = 3000, width = 3000, resolution = 300, compression = "lzw", 
                                lwd = 2, col = my_colors, fill = alpha(my_colors, 0.3), 
                                cex = 1.5, fontfamily = "sans", cat.cex = 1.5, cat.fontfamily = "sans", 
                                rotation.degree = 0,
                                main = "Intersection of Exposure Signals by Outcome")
      pdf(venn_pdf, width = 10, height = 10); grid.draw(venn.plot); dev.off()
      png(venn_png, width = 3000, height = 3000, res = 300); grid.draw(venn.plot); dev.off()
    }
    
    all_combinations <- list()
    for (k in 1:length(all_outcomes)) {
      combs <- combn(all_outcomes, k, simplify = FALSE)
      all_combinations <- c(all_combinations, combs)
    }
    bar_data_list <- list()
    for (comb in all_combinations) {
      comb_name <- paste(comb, collapse = " & ")
      has_comb <- sapply(cluster_summary$Outcomes, function(x) all(comb %in% x))
      subset_dt <- cluster_summary[has_comb, ]
      if (nrow(subset_dt) > 0) {
        counts <- subset_dt[, .(Count = .N), by = Exposure_Cell]
        counts[, Intersection := comb_name]
        counts[, N_Outcomes := length(comb)]
        bar_data_list[[length(bar_data_list) + 1]] <- counts
      }
    }
    if (length(bar_data_list) > 0) {
      bar_data <- rbindlist(bar_data_list)
      get_len <- function(s) length(strsplit(s, " & ")[[1]])
      inter_levels <- unique(bar_data$Intersection)
      inter_order_idx <- order(sapply(inter_levels, get_len), inter_levels)
      bar_data$Intersection <- factor(bar_data$Intersection, levels = inter_levels[inter_order_idx])
      n_cells <- length(unique(bar_data$Exposure_Cell))
      cell_colors <- colorRampPalette(brewer.pal(12, "Paired"))(n_cells)
      p_bar <- ggplot(bar_data, aes(x = Intersection, y = Count, fill = Exposure_Cell)) + geom_bar(stat = "identity", position = "stack", width = 0.6, color = "white", size = 0.2) + scale_fill_manual(values = cell_colors) + theme_classic(base_size = 14) + labs(title = "Inclusive Exposure Signal Counts by Outcome Intersection", subtitle = "Cumulative counts: Clusters containing at least the specified outcomes", x = "Intersection Type", y = "Number of Signals (Clusters)", fill = "Cell Type") + theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, color = "black", size = 12), axis.text.y = element_text(color = "black", size = 12), axis.title = element_text(face = "bold"), legend.position = "right", plot.title = element_text(hjust = 0.5, face = "bold", size = 16), panel.grid.major.y = element_line(color = "grey90", linetype = "dashed")) + geom_text(aes(label = after_stat(y), group = Intersection), stat = 'summary', fun = sum, vjust = -0.5, size = 4)
      bar_pdf <- file.path(venn_dir, "Intersection_Bar_Chart_Inclusive.pdf")
      bar_png <- file.path(venn_dir, "Intersection_Bar_Chart_Inclusive.png")
      ggsave(bar_pdf, p_bar, width = 12, height = 8)
      ggsave(bar_png, p_bar, width = 12, height = 8, dpi = 300)
    }
  }
}, error = function(e) {
  log_msg(glue("  Error in Intersection Analysis: {e$message}"))
})

# ==============================================================================
# 3.1. Exposure Cluster Summary
# ==============================================================================
# Generate summary CSV for Exposure Analysis
# Format: Cluster_ID, Lead_SNPs, Associations, Outcomes
# Sort by: Exposure_Cell

log_msg("--- Generating Exposure Cluster Summary ---")

# Define associations locally for this summary (cluster_associations is defined later in script)
#  (cluster_associations )
cluster_associations_temp <- merge(exposure_clusters, dt_clean, 
                                   by.x = c("Lead_SNP", "Exposure_Cell"), 
                                   by.y = c("Lead_SNP", "cell"), 
                                   all.x = TRUE)

# Generate summary directly from valid associations to avoid Cartesian product errors
# ,  (Fix: Avoid fictitious associations)
exposure_summary_final <- cluster_associations_temp[, .(
  Exposure_Cell = unique(Exposure_Cell),
  Lead_SNPs = paste(sort(unique(Lead_SNP)), collapse = "; "),
  Associations = paste(sort(unique(paste(Exposure_Cell, Gene, Lead_SNP, Outcome, sep = "-"))), collapse = "; "),
  Outcomes = paste(sort(unique(Outcome)), collapse = "; ")
), by = Exposure_Cluster_ID]

# Verification Step / 
log_msg("Verifying Associations in Exposure Cluster Summary... /  Exposure Cluster Summary ...")

# Construct valid associations set from source data
valid_associations_set <- unique(paste(dt_clean$cell, dt_clean$Gene, dt_clean$Lead_SNP, dt_clean$Outcome, sep = "-"))

# Check each association in the summary
all_summary_associations <- unlist(strsplit(exposure_summary_final$Associations, "; "))
# Remove any empty strings if present / 
all_summary_associations <- all_summary_associations[all_summary_associations != ""]

invalid_assocs <- setdiff(all_summary_associations, valid_associations_set)

if (length(invalid_assocs) > 0) {
  error_msg <- glue("CRITICAL ERROR: Found {length(invalid_assocs)} invalid associations in summary! / :  {length(invalid_assocs)} ! Examples: {paste(head(invalid_assocs, 3), collapse=', ')}")
  log_msg(error_msg)
  stop(error_msg)
} else {
  log_msg("Verification passed: All associations in summary match source data. / : .")
}

# Sort by Exposure_Cell
exposure_summary_final <- exposure_summary_final[order(Exposure_Cell, Exposure_Cluster_ID)]

# Save
exp_analysis_dir <- file.path(output_dir, "Exposure_Analysis")
if (!dir.exists(exp_analysis_dir)) dir.create(exp_analysis_dir, recursive = TRUE)

fwrite(exposure_summary_final, file.path(exp_analysis_dir, "Exposure_Cluster_Summary.csv"))
log_msg(glue("Exposure Cluster Summary saved to {file.path(exp_analysis_dir, 'Exposure_Cluster_Summary.csv')}"))

# 4.2 Outcome Analysis (Group by Phenotype/Outcome)
# ------------------------------------------------------------------------------
log_msg("--- Running Outcome Analysis (Phenotype-grouped) ---")
outcome_clusters <- perform_clustering(
  data = dt_clean,
  group_col = "Outcome",
  analysis_name = "Outcome",
  output_subdir = "Outcome_Analysis",
  global_mat = mat,
  r2_threshold = 0.6
)
setnames(outcome_clusters, "Cluster_ID", "Outcome_Cluster_ID")
setnames(outcome_clusters, "Group", "Outcome_Phenotype")

# 4.3 Global Clustering (Outcome-agnostic) for Right Side of Sankey
# ------------------------------------------------------------------------------
log_msg("--- Running Global SNP Clustering (Outcome-agnostic) ---")
# Group by NULL to treat all SNPs as one group
global_clustering <- perform_clustering(
  data = dt_clean,
  group_col = NULL, # No grouping column -> All SNPs
  analysis_name = "Global_SNP",
  output_subdir = "Global_SNP_Analysis",
  global_mat = mat,
  r2_threshold = 0.6
)
# Rename columns
setnames(global_clustering, "Cluster_ID", "Global_Cluster_ID")
setnames(global_clustering, "Group", "Global_Group") # Will be "Global"

# ==============================================================================
# 4.3.5 Outcome Signal Intersection Analysis / 
# ==============================================================================
log_msg("--- Running Outcome Signal Intersection Analysis ---")
log_msg("Logic: Identify which Outcome Clusters (Phenotype-specific) map to the same Global Cluster (R2>=0.6).")

# 1. Map Outcome Clusters to Global Clusters
# outcome_clusters: Lead_SNP, Outcome_Phenotype, Outcome_Cluster_ID
# global_clustering: Lead_SNP, Global_Cluster_ID

# Merge on Lead_SNP
outcome_global_map <- merge(outcome_clusters, global_clustering[, .(Lead_SNP, Global_Cluster_ID)], 
                            by = "Lead_SNP", all.x = TRUE)

# 2. Summarize by Global_Cluster_ID
# We want to see for each Global Cluster, which Outcome Clusters are present.
# Identify mapping: Global_ID -> [Aging_Signal_X, RA_Signal_Y, HZ_Signal_Z]
outcome_intersection_summary <- outcome_global_map[, .(
  Outcomes_Present = paste(sort(unique(Outcome_Phenotype)), collapse = ";"),
  N_Outcomes = uniqueN(Outcome_Phenotype),
  mvAge_Signals = paste(unique(Outcome_Cluster_ID[Outcome_Phenotype == "mvAge"]), collapse = "; "),
  RA_Signals = paste(unique(Outcome_Cluster_ID[Outcome_Phenotype == "RA"]), collapse = "; "),
  HZ_Signals = paste(unique(Outcome_Cluster_ID[Outcome_Phenotype == "HZ"]), collapse = "; ")
), by = Global_Cluster_ID]

# Sort by N_Outcomes (descending)
outcome_intersection_summary <- outcome_intersection_summary[order(-N_Outcomes)]

# Save Intersection Summary
intersection_file <- file.path(output_dir, "Outcome_Signal_Intersection.csv")
fwrite(outcome_intersection_summary, intersection_file)
log_msg(glue("Outcome Intersection Summary saved to {intersection_file}"))

# 3. Venn Diagram of Global Clusters by Phenotype
# Create a list of Global_Cluster_IDs for each phenotype
venn_list_outcome <- list(
  mvAge = unique(outcome_global_map[Outcome_Phenotype == "mvAge", Global_Cluster_ID]),
  RA = unique(outcome_global_map[Outcome_Phenotype == "RA", Global_Cluster_ID]),
  HZ = unique(outcome_global_map[Outcome_Phenotype == "HZ", Global_Cluster_ID])
)

# Plot Venn using VennDiagram package (ggvenn not available)
library(VennDiagram)
library(grid)

# Generate Venn diagram (returns a grob)
venn_plot <- venn.diagram(
  x = venn_list_outcome,
  category.names = names(venn_list_outcome),
  filename = NULL,
  output = TRUE,
  
  # Circles
  lwd = 2,
  lty = 'blank',
  fill = c("#E41A1C", "#377EB8", "#4DAF4A"),
  
  # Numbers
  cex = 1.5,
  fontface = "bold",
  fontfamily = "sans",
  
  # Set names
  cat.cex = 1.5,
  cat.fontface = "bold",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135),
  cat.dist = c(0.055, 0.055, 0.085),
  cat.fontfamily = "sans",
  rotation = 1
)

# Add Title
title_grob <- textGrob(label = "Intersection of Global Signals (R2>=0.6) by Phenotype",
                       x = 0.5, y = 1.05, gp = gpar(fontsize = 20, fontface = "bold"))
combined_grob <- gTree(children = gList(venn_plot, title_grob))

# Save PDF
venn_outcome_pdf <- file.path(output_dir, "Outcome_Signal_Venn.pdf")
pdf(venn_outcome_pdf, width = 8, height = 8)
grid.draw(combined_grob)
dev.off()
log_msg(glue("Outcome Venn Diagram saved to {venn_outcome_pdf}"))

# Save PNG
venn_outcome_png <- file.path(output_dir, "Outcome_Signal_Venn.png")
png(venn_outcome_png, width = 2400, height = 2400, res = 300)
grid.draw(combined_grob)
dev.off()

# ==============================================================================
  # 4.3.6. Global Cluster Summary
  # ==============================================================================
  # Generate summary CSV for Global SNP Analysis
  # Format: Global_Cluster_ID, Lead_SNPs, Outcomes, N_Outcomes, Aging_Signals, RA_Signals, HZ_Signals, Exposure_Signals, N_Exposure_Signals
  # Sort by: Global_Cluster_ID
  # Includes ALL Global Signals (not just shared ones)
  
  log_msg("--- Generating Global Cluster Summary ---")
  
  # 1. Base: Global Clusters (Lead_SNP, Global_Cluster_ID)
  #    Join with dt_clean to get Outcomes (raw)
  global_base <- merge(global_clustering, dt_clean[, .(Lead_SNP, Outcome)], by = "Lead_SNP", allow.cartesian = TRUE)
  
  global_summary_final <- global_base[, .(
    Lead_SNPs = paste(sort(unique(Lead_SNP)), collapse = "; "),
    Outcomes = paste(sort(unique(Outcome)), collapse = "; ")
  ), by = Global_Cluster_ID]
  
  # 2. Add Outcome Signals (mvAge/RA/HZ)
  #    Reuse outcome_intersection_summary from Section 4.3.5 if available
  if (exists("outcome_intersection_summary")) {
    global_summary_final <- merge(global_summary_final, 
                                  outcome_intersection_summary[, .(Global_Cluster_ID, N_Outcomes, mvAge_Signals, RA_Signals, HZ_Signals)], 
                                  by = "Global_Cluster_ID", all.x = TRUE)
  } else {
    log_msg("  Warning: outcome_intersection_summary not found. Skipping Outcome Signal details.")
    global_summary_final[, N_Outcomes := 0]
    global_summary_final[, mvAge_Signals := ""]
    global_summary_final[, RA_Signals := ""]
    global_summary_final[, HZ_Signals := ""]
  }
  
  # 3. Add Exposure Signals (Exposure Clusters from different cells)
  #    Logic: Identify all Exposure Clusters that contain any SNP in the Global Cluster
  
  # Join global_clustering with exposure_clusters
  # global_clustering: Lead_SNP, Global_Cluster_ID
  # exposure_clusters: Lead_SNP, Exposure_Cluster_ID, Exposure_Cell
  
  global_exposure_map <- merge(global_clustering, exposure_clusters[, .(Lead_SNP, Exposure_Cell, Exposure_Cluster_ID)], 
                               by = "Lead_SNP", allow.cartesian = TRUE)
  
  global_exposure_agg <- global_exposure_map[, .(
    Exposure_Signals = paste(sort(unique(paste(Exposure_Cell, Exposure_Cluster_ID, sep = "-"))), collapse = "; "),
    N_Exposure_Signals = uniqueN(Exposure_Cluster_ID)
  ), by = Global_Cluster_ID]
  
  # Merge Exposure Info
  global_summary_final <- merge(global_summary_final, global_exposure_agg, by = "Global_Cluster_ID", all.x = TRUE)
  
  # Fill NAs
  cols_to_fill <- c("mvAge_Signals", "RA_Signals", "HZ_Signals", "Exposure_Signals")
  for (col in cols_to_fill) {
    if (col %in% names(global_summary_final)) {
       set(global_summary_final, i = which(is.na(global_summary_final[[col]])), j = col, value = "")
    }
  }
  
  cols_to_zero <- c("N_Outcomes", "N_Exposure_Signals")
  for (col in cols_to_zero) {
     if (col %in% names(global_summary_final)) {
       set(global_summary_final, i = which(is.na(global_summary_final[[col]])), j = col, value = 0)
    }
  }
  
  # Sort by Global_Cluster_ID (extract numeric part for correct sorting)
  global_summary_final[, ID_Num := as.numeric(sub("Global_Cluster_", "", Global_Cluster_ID))]
  global_summary_final <- global_summary_final[order(ID_Num)]
  global_summary_final[, ID_Num := NULL]
  
  # Save
  global_analysis_dir <- file.path(output_dir, "Global_SNP_Analysis")
  if (!dir.exists(global_analysis_dir)) dir.create(global_analysis_dir, recursive = TRUE)
  
  fwrite(global_summary_final, file.path(global_analysis_dir, "Global_Cluster_Summary.csv"))
  log_msg(glue("Global Cluster Summary saved to {file.path(global_analysis_dir, 'Global_Cluster_Summary.csv')}"))

# ==============================================================================
# 4.4 Exposure Cluster Sharing Classification / 
# ==============================================================================
log_msg("--- Running Exposure Cluster Sharing Classification ---")
log_msg("Logic: Classify Exposure Clusters based on the number of unique Outcomes they are associated with.")

# 1. Join exposure_clusters with original data to get Outcomes for each SNP in the cluster
# exposure_clusters: Lead_SNP, Exposure_Cluster_ID, Exposure_Cell
# dt_clean: Lead_SNP, cell, Outcome, Display_Name (Gene) -> Use 'Gene' column populated in data processing step
# Note: exposure_clusters$Exposure_Cell matches dt_clean$cell

cluster_associations <- merge(exposure_clusters, dt_clean, 
                              by.x = c("Lead_SNP", "Exposure_Cell"), 
                              by.y = c("Lead_SNP", "cell"), 
                              all.x = TRUE)

# 2. Count unique outcomes per Exposure_Cluster_ID
cluster_sharing_status <- cluster_associations[, .(
  N_Outcomes = uniqueN(Outcome),
  Outcomes_List = paste(sort(unique(Outcome)), collapse = ";")
), by = .(Exposure_Cluster_ID, Exposure_Cell)]

# 3. Classify Sharing Type
cluster_sharing_status[, Sharing_Type := fcase(
  N_Outcomes == 3, "3-way Shared",
  N_Outcomes == 2, "2-way Shared",
  N_Outcomes == 1, "Unique",
  default = "Unknown"
)]

# Save Classification Result
sharing_file <- file.path(output_dir, "Exposure_Cluster_Sharing_Status.csv")
fwrite(cluster_sharing_status, sharing_file)
log_msg(glue("Sharing classification saved to {sharing_file}"))

# Log counts
sharing_counts <- cluster_sharing_status[, .N, by = Sharing_Type]
print(sharing_counts)
log_msg(paste("Sharing Counts:\n", paste(capture.output(print(sharing_counts)), collapse = "\n")))

# ==============================================================================
# 5. Signal Mapping Visualization (Sankey Diagram) / 
# ==============================================================================
log_msg("--- Generating Signal Mapping Visualization (Exposure -> Global SNP Clusters) ---")
log_msg("Logic: Exposure Clusters (Cell-specific, R2>=0.6) -> Global Clusters (Outcome-agnostic, R2>=0.6)")

mapping_base_dir <- file.path(output_dir, "Signal_Mapping")
if (!dir.exists(mapping_base_dir)) dir.create(mapping_base_dir, recursive = TRUE)

# Function to generate signal mapping for a specific set of target outcomes
generate_signal_mapping <- function(target_outcomes, subfolder_name) {
  
  log_msg(glue("Running Signal Mapping for: {subfolder_name}"))
  log_msg(glue("Target Outcomes: {paste(target_outcomes, collapse = ', ')}"))
  
  # Create subfolder
  current_mapping_dir <- file.path(mapping_base_dir, subfolder_name)
  if (!dir.exists(current_mapping_dir)) dir.create(current_mapping_dir, recursive = TRUE)
  
  # --- Outcome-Anchored Filtering ---
  # Goal: Identify Global Clusters that are associated with ALL target outcomes.
  #       Then retrieve ALL Exposure Clusters (Cell-specific) linked to these Global Clusters.
  #       This ensures we capture signals mediated by DIFFERENT cells (e.g. Global_Signal_91).
  
  # 1. Identify "Shared" Global Clusters (Outcome-centric)
  #    Use global_clustering (Lead_SNP, Global_Cluster_ID) and dt_clean (Lead_SNP, Outcome)
  
  # Create a mapping of Global_ID -> Outcome
  global_outcome_map <- merge(global_clustering, dt_clean[, .(Lead_SNP, Outcome)], by="Lead_SNP", allow.cartesian=TRUE)
  
  # Filter for Global Clusters that contain associations for ALL target_outcomes
  valid_global_dt <- global_outcome_map[, .(Has_All_Targets = all(target_outcomes %in% Outcome)), by = Global_Cluster_ID]
  shared_global_ids <- valid_global_dt[Has_All_Targets == TRUE, Global_Cluster_ID]
  
  if (length(shared_global_ids) == 0) {
    log_msg(glue("  No shared Global Clusters found for: {paste(target_outcomes, collapse = ', ')}"))
    return(NULL)
  }
  
  log_msg(glue("  Found {length(shared_global_ids)} shared Global Clusters (Outcome-centric)."))
  
  # 2. Retrieve Exposure Clusters linked to these Shared Global Clusters
  #    Condition: The Exposure Cluster must be linked to a Shared Global Cluster via SNPs
  #    AND must contain at least one association relevant to the target outcomes (to avoid noise).
  
  # Get all SNPs belonging to the shared global clusters
  shared_snps <- global_clustering[Global_Cluster_ID %in% shared_global_ids, Lead_SNP]
  
  # Filter Exposure Clusters that contain these SNPs
  # exposure_clusters: Exposure_Cluster_ID, Lead_SNP, Exposure_Cell
  exp_candidate <- exposure_clusters[Lead_SNP %in% shared_snps]
  
  # Further filter: Ensure the Exposure Cluster has associations with at least one target outcome
  # This prevents showing Exposure Clusters that might be in LD but affect a totally different outcome (if any exist).
  # But since we are only analyzing 3 outcomes, this is less of an issue, but good for safety.
  
  # Join with cluster_associations to check outcomes
  exp_assoc_check <- merge(exp_candidate, cluster_associations[, .(Exposure_Cluster_ID, Outcome)], 
                           by = "Exposure_Cluster_ID", allow.cartesian = TRUE)
  
  # Keep Exposure_Cluster_IDs that have at least one target outcome
  valid_exp_ids <- unique(exp_assoc_check[Outcome %in% target_outcomes, Exposure_Cluster_ID])
  
  # Final Filtered Exposure Clusters
  exp_filtered <- exposure_clusters[Exposure_Cluster_ID %in% valid_exp_ids]
  
  # Update shared_snps to only those present in valid Exposure Clusters (for clean mapping)
  shared_snps_final <- unique(exp_filtered$Lead_SNP)
  
  # 3. Filter Global Clusters (for Sankey Flow) based on final SNPs
  glob_filtered <- global_clustering[Lead_SNP %in% shared_snps_final]
  glob_filtered_unique <- unique(glob_filtered[, .(Lead_SNP, Global_Cluster_ID)])
  
  # 4. Create Sankey Data
  # Merge Exposure (Left) and Global (Right) via Lead_SNP
  sankey_data <- merge(exp_filtered, glob_filtered_unique, by = "Lead_SNP", suffixes = c("_Exp", "_Glob"), allow.cartesian = TRUE)
  
  # Handle missing Exposure_Cell if any (should not happen given data structure)
  if (!"Exposure_Cell" %in% names(sankey_data)) {
    # Try to recover Exposure_Cell from exposure_clusters if missing
    sankey_data <- merge(sankey_data, unique(exposure_clusters[, .(Exposure_Cluster_ID, Exposure_Cell)]), 
                         by="Exposure_Cluster_ID", all.x=TRUE)
  }

  # Aggregate Flow: Count unique SNPs per flow (Exposure_Cluster -> Global_Cluster)
  flow_data <- sankey_data[, .(Weight = uniqueN(Lead_SNP)), by = .(Exposure_Cluster_ID, Exposure_Cell, Global_Cluster_ID)]
  
  # 5. Filter Associations for Labels (Restrict to target outcomes)
  # We only want to list the target outcomes in the labels
  target_associations <- cluster_associations[Exposure_Cluster_ID %in% valid_exp_ids & Outcome %in% target_outcomes]
  
  # 6. Generate Exposure Labels (Restricted to Target Outcomes)
  exp_labels <- sapply(valid_exp_ids, function(id) {
    assocs <- unique(target_associations[Exposure_Cluster_ID == id, .(Exposure_Cell, Gene, Lead_SNP, Outcome)])
    # Order for consistency
    assocs <- assocs[order(Outcome, Gene, Lead_SNP)]
    
    label_lines <- apply(assocs, 1, function(row) {
      gene_name <- if (is.na(row['Gene']) || row['Gene'] == "") "Unknown" else row['Gene']
      paste(row['Exposure_Cell'], gene_name, row['Lead_SNP'], row['Outcome'], sep = "-")
    })
    paste(label_lines, collapse = "\n")
  })
  
  # 7. Generate Global Labels
  # For Global Labels, list the SNPs that are part of this flow context
  unique_glob_ids <- unique(flow_data$Global_Cluster_ID)
  glob_labels <- sapply(unique_glob_ids, function(id) {
    snps <- unique(global_clustering[Global_Cluster_ID == id, Lead_SNP])
    relevant_snps <- intersect(snps, shared_snps_final)
    paste(relevant_snps, collapse = "\n")
  })
  
  # 8. Save Reference Info and Mapping Table
  # Add labels to flow_data for reference
  # Note: labels are named vectors, need to match by ID
  flow_data$Exposure_Label <- exp_labels[flow_data$Exposure_Cluster_ID]
  flow_data$Global_Label <- glob_labels[flow_data$Global_Cluster_ID]
  
  ref_info <- flow_data[, .(Exposure_Cluster_ID, Global_Cluster_ID, Exposure_Cell, Exposure_Label, Global_Label)]
  fwrite(ref_info, file.path(current_mapping_dir, "Cluster_ID_Reference_Info.csv"))
  fwrite(flow_data, file.path(current_mapping_dir, "Filtered_Exposure_to_Global_Mapping.csv"))
  
  # 9. Plot Sankey
  if (require(ggalluvial)) {
    library(ggalluvial)
    library(ggplot2)
  
  # Colors (Cell Type)
  cell_types <- unique(flow_data$Exposure_Cell)
  if (length(cell_types) <= 8) {
    pal <- RColorBrewer::brewer.pal(max(3, length(cell_types)), "Set2")
  } else {
    pal <- rainbow(length(cell_types))
  }
  names(pal) <- cell_types
  
  # Fill mapping
  left_ids <- unique(flow_data$Exposure_Cluster_ID)
  right_ids <- unique(flow_data$Global_Cluster_ID)
  
  # Right = Gray
  fill_values <- rep("grey80", length(right_ids))
  names(fill_values) <- right_ids
  
  # Left = Cell Color
  id_to_cell <- unique(flow_data[, .(Exposure_Cluster_ID, Exposure_Cell)])
  id_to_cell <- id_to_cell[!duplicated(Exposure_Cluster_ID)]
  left_fills <- pal[id_to_cell$Exposure_Cell]
  names(left_fills) <- id_to_cell$Exposure_Cluster_ID
  
  final_fill_values <- c(fill_values, left_fills)
  
  # Remove NA names if any
  final_fill_values <- final_fill_values[!is.na(names(final_fill_values))]
  
  p_sankey <- ggplot(flow_data,
                     aes(y = Weight, axis1 = Exposure_Cluster_ID, axis2 = Global_Cluster_ID)) +
    geom_alluvium(aes(fill = Exposure_Cell), width = 1/8, alpha = 0.7) +
    geom_stratum(width = 1/8, aes(fill = after_stat(stratum)), color = "grey30", size = 0.2) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3.0, min.y = 1) +
    scale_x_discrete(limits = c("Exposure Clusters", "Global Clusters"), expand = c(.1, .1)) +
    scale_fill_manual(values = final_fill_values, name = "Cell Type") +
    labs(title = glue("Signal Mapping: {subfolder_name}"),
         subtitle = glue("Shared by: {paste(target_outcomes, collapse = ', ')}"),
         y = "Number of Shared SNPs") +
    theme_void() +
    theme(
      legend.position = "none",
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      plot.subtitle = element_text(hjust = 0.5, color = "grey40", size = 12),
      axis.text.x = element_text(size = 14, face = "bold", vjust = 1),
      plot.margin = margin(20, 20, 20, 20)
    )
  
  ggsave(file.path(current_mapping_dir, "Exposure_Global_Sankey.pdf"), p_sankey, width = 16, height = 12)
  ggsave(file.path(current_mapping_dir, "Exposure_Global_Sankey.png"), p_sankey, width = 16, height = 12, dpi = 300)
  
  log_msg(glue("  Saved results to {current_mapping_dir}"))
  } else {
    log_msg("  Package 'ggalluvial' not installed. Skipping Sankey diagram.")
  }
}

# --- Execute Analysis for Different Combinations ---

# 1. 3-Way Shared (mvAge, RA, HZ)
# Note: Ensure outcome names match your data exactly. 
# Based on logs: "mvAge", "RA", "HZ"
outcomes_3way <- c("mvAge", "RA", "HZ")
generate_signal_mapping(outcomes_3way, "3_Way_Shared_mvAge_RA_HZ")

# 2. 2-Way Shared (Inclusive)
# Pairs: mvAge-RA, mvAge-HZ, RA-HZ
pairs <- list(
  c("mvAge", "RA"),
  c("mvAge", "HZ"),
  c("RA", "HZ")
)

for (pair in pairs) {
  folder_name <- glue("2_Way_Shared_{paste(pair, collapse = '_')}")
  generate_signal_mapping(pair, folder_name)
}


# ==============================================================================
# 6. Final Summary / 
# ==============================================================================

# Update Master Summary
# Join Exposure, Outcome, AND Global clusters
# Start with dt_clean
master_summary <- merge(dt_clean, exposure_clusters[, .(Lead_SNP, Exposure_Cell, Exposure_Cluster_ID)], 
                        by.x = c("Lead_SNP", "cell"), by.y = c("Lead_SNP", "Exposure_Cell"), all.x = TRUE)

master_summary <- merge(master_summary, outcome_clusters[, .(Lead_SNP, Outcome_Phenotype, Outcome_Cluster_ID)],
                        by.x = c("Lead_SNP", "Outcome"), by.y = c("Lead_SNP", "Outcome_Phenotype"), all.x = TRUE)

# Global clustering is unique per SNP (no group key needed)
global_unique <- unique(global_clustering[, .(Lead_SNP, Global_Cluster_ID)])
master_summary <- merge(master_summary, global_unique, by = "Lead_SNP", all.x = TRUE)

fwrite(master_summary, file.path(output_dir, "Master_Clustered_Summary.csv"))
log_msg(glue("Master summary saved to {file.path(output_dir, 'Master_Clustered_Summary.csv')}"))

log_msg("Analysis Complete.")
