# ==============================================================================
# [Script]: 4.1.5.1_visualize_3D_signals_template.R
# [Method]: 3D Association Signal Visualization
# [Step]: 4.1.5.1_visualize_3D_signals
#
# [Function]:
# Performs comprehensive visualization of 3D trait associations accounting for LD structure.
#
# [Usage]: 
#   Rscript 4.1.5.1_visualize_3D_signals_template.R \
#     --input <path_to_merged_intersection_3d_data.csv> \
#     --out_dir <path_to_output_directory>
# ==============================================================================

# ------------------------------------------------------------------------------
# 1. Setup & Dependencies 
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(readr)
  library(tidyr)
  library(stringr)
  library(ggrepel)
  library(pheatmap)
  library(gridExtra)
  library(RColorBrewer)
  library(scales)
  library(ggforce) # For geom_mark_ellipse/hull
  library(tibble) # For column_to_rownames
})


# Clear environment 
rm(list = ls())

# Define command line arguments
suppressPackageStartupMessages({
  library(optparse)
})

option_list <- list(
  make_option(c("--input"), type="character", default=NULL,
              help="Input CSV file containing merged intersection 3D data / 3DCSV"),
  make_option(c("--out_dir"), type="character", default="./3D_Signal_Visualization",
              help="Output directory path / ")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input)) {
  print_help(opt_parser)
  stop("Missing required argument: --input / : --input", call.=FALSE)
}

# Define paths 
INPUT_FILE <- opt$input
# Dynamic output directory with timestamp 
OUTPUT_DIR <- file.path(opt$out_dir, paste0(format(Sys.time(), "%Y%m%d_%H%M%S"), "_3D_signal_visualization"))


if (!dir.exists(OUTPUT_DIR)) {
  dir.create(OUTPUT_DIR, recursive = TRUE)
}

# ------------------------------------------------------------------------------
# 2. Data Loading & Preprocessing 
# ------------------------------------------------------------------------------
cat("Loading data... (...)\n")
df <- read_csv(INPUT_FILE, show_col_types = FALSE)

# Basic checks 
if (nrow(df) == 0) stop("Input data is empty. ")
required_cols <- c("b_aging", "b_RA", "b_HZ", "cell", "Cluster_ID", "SNP", "gene")
if (!all(required_cols %in% colnames(df))) {
  stop(paste("Missing columns :", paste(setdiff(required_cols, colnames(df)), collapse=", ")))
}

# Define Signal Source (cell + Cluster_ID) 
# Note: Distinct signal source even if LD overlaps, as per user rule.
# Fix: Avoid redundant cell name if Cluster_ID already starts with cell name
# Use regex to be more robust and ensure no whitespace issues
df <- df %>%
  mutate(
    cell = str_trim(cell),
    Cluster_ID = str_trim(Cluster_ID),
    # Check if Cluster_ID starts with cell name (case insensitive just in case, though usually consistent)
    starts_with_cell = str_detect(Cluster_ID, regex(paste0("^", cell), ignore_case = TRUE))
  ) %>%
  filter(Cluster_ID != "cd4nc_Signal_20") %>%
  mutate(signal_source = if_else(
    starts_with_cell, 
    Cluster_ID, 
    paste(cell, Cluster_ID, sep = "_")
  )) %>%
  select(-starts_with_cell) # Clean up temporary column

# Standardization (Z-score) 
# Critical due to scale differences (Aging ~0.03 vs RA/HZ ~0.5)
traits <- c("b_aging", "b_RA", "b_HZ")
df_scaled <- df %>%
  mutate(across(all_of(traits), ~ scale(.)[,1], .names = "z_{.col}"))

# ------------------------------------------------------------------------------
# 3. Dimensionality Reduction (PCA) 
# ------------------------------------------------------------------------------
cat("Performing PCA... (PCA...)\n")

# Use Z-scored data for PCA
pca_input <- df_scaled %>% select(starts_with("z_"))
pca_res <- prcomp(pca_input, center = FALSE, scale. = FALSE) # Already scaled

# Extract PC scores
df_pca <- df_scaled %>%
  bind_cols(as.data.frame(pca_res$x))

# Calculate variance explained 
var_exp <- pca_res$sdev^2 / sum(pca_res$sdev^2)
pc1_var <- percent(var_exp[1], accuracy=0.1)
pc2_var <- percent(var_exp[2], accuracy=0.1)

# Extract Loadings (for interpretation) 
loadings <- as.data.frame(pca_res$rotation)
loadings$Variable <- rownames(loadings)

cat("Variance Explained :\n")
print(var_exp)
cat("Loadings :\n")
print(loadings)

# Save PCA results (PCA)
write_csv(df_pca, file.path(OUTPUT_DIR, "pca_results_with_metadata.csv"))
write_csv(loadings, file.path(OUTPUT_DIR, "pca_loadings.csv"))

# ------------------------------------------------------------------------------
# 4. Visualization: Gene Level (Bubble Plot) (: -)
# ------------------------------------------------------------------------------
cat("Plotting Gene-Level Visualization (Bubble Plot)... (...)\n")

# Prepare matrix for bubble plot: Signal Source vs Gene (Value = PC1 score)
# CRITICAL UPDATE: Handling LD redundancy / Multi-SNP per Gene
# Instead of simple mean(), we should select the "Representative SNP" for each gene
# within each signal source to avoid averaging out strong signals with weak LD proxies.
# Strategy: Select the SNP with the largest absolute PC1 score (most extreme effect)
# or the one with the smallest p-value (if available, here we use PC1 magnitude as proxy for relevance).

snp_counts <- df_pca %>%
  group_by(gene, signal_source) %>%
  summarise(snp_n = dplyr::n_distinct(SNP), .groups = "drop")

gene_matrix_long <- df_pca %>%
  group_by(gene, signal_source) %>%
  slice_max(order_by = abs(PC1), n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(gene, signal_source, score = PC1) %>%
  left_join(snp_counts, by = c("gene", "signal_source"))

# Reorder axes for publication-style readability 
# Goal: Signal set on X axis, Gene on Y axis; use hierarchical clustering to order both axes.
mat_for_order <- gene_matrix_long %>%
  select(signal_source, gene, score) %>%
  pivot_wider(names_from = gene, values_from = score, values_fill = 0) %>%
  column_to_rownames("signal_source") %>%
  as.matrix()

row_order <- hclust(dist(mat_for_order))$order
col_order <- hclust(dist(t(mat_for_order)))$order
signal_levels <- rownames(mat_for_order)[row_order]
gene_levels <- colnames(mat_for_order)[col_order]

gene_matrix_long <- gene_matrix_long %>%
  mutate(
    signal_source = factor(signal_source, levels = signal_levels),
    gene = factor(gene, levels = gene_levels)
  )

# Create Bubble Plot using ggplot2 (ggplot2)
# Requirement: bubble size reflects SNP count per (signal_source, gene) in source data.
max_abs_score <- max(abs(gene_matrix_long$score), na.rm = TRUE)
n_genes <- n_distinct(gene_matrix_long$gene)
n_signals <- n_distinct(gene_matrix_long$signal_source)

size_cap <- as.numeric(stats::quantile(gene_matrix_long$snp_n, probs = 0.95, na.rm = TRUE, names = FALSE))
size_cap <- max(size_cap, 10)
gene_matrix_long <- gene_matrix_long %>%
  mutate(snp_n_capped = pmin(snp_n, size_cap))

plot_width <- max(9, min(18, 3.2 + 0.075 * n_signals))
plot_height <- max(10, min(36, 6 + 0.22 * n_genes))

p_bubble <- ggplot(gene_matrix_long, aes(x = signal_source, y = gene, fill = score, size = snp_n_capped)) +
  geom_point(
    shape = 21,
    color = "grey20",
    stroke = 0.25,
    alpha = 0.95
  ) +
  scale_size_continuous(
    range = c(2.2, 8.0),
    trans = "sqrt",
    breaks = scales::pretty_breaks(n = 4),
    name = "LD SNP Count"
  ) +
  scale_fill_gradient2(
    low = "#3B4CC0",
    mid = "white",
    high = "#B40426",
    midpoint = 0,
    limits = c(-max_abs_score, max_abs_score),
    oob = scales::squish,
    name = "PC1 Score"
  ) +
  scale_x_discrete(expand = expansion(mult = c(0, 0), add = c(0, 0))) +
  scale_y_discrete(
    expand = expansion(mult = c(0, 0), add = c(0.32, 0.32)),
    labels = function(x) paste0(x, "        ") # Add padding spaces to force spacing
  ) +
  guides(
    fill = guide_colorbar(
      barheight = unit(90, "pt"),
      barwidth = unit(10, "pt"),
      title.position = "top"
    )
  ) +
  theme_minimal(base_size = 12, base_family = "Helvetica") +
  theme(
    panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
    panel.grid.minor = element_blank(),
    axis.title = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1,
      size = 6.8,
      color = "grey10",
      margin = margin(t = 6)
    ),
    axis.text.y = element_text(hjust = 1, size = 7.2, color = "grey10", face = "italic", margin = margin(r = 0)), # Use padding spaces instead of margin
    axis.ticks.y = element_line(color = "transparent"), # Use transparent ticks to enforce spacing
    axis.ticks.length.y = unit(5, "pt"),
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    plot.margin = margin(10, 18, 38, 44)
  ) +
  labs(
    x = NULL,
    y = NULL
  ) +
  coord_cartesian(clip = "off")

# Save Bubble Plot 
ggsave(file.path(OUTPUT_DIR, "3_Gene_Signal_Source_BubblePlot.pdf"), p_bubble, width = plot_width, height = plot_height, limitsize = FALSE)
ggsave(file.path(OUTPUT_DIR, "3_Gene_Signal_Source_BubblePlot.png"), p_bubble, width = plot_width, height = plot_height, dpi = 400, limitsize = FALSE)

# Save the data used for Bubble Plot as supplementary data 
write_csv(gene_matrix_long, file.path(OUTPUT_DIR, "3_Gene_Signal_Source_BubblePlot_Data.csv"))

# ------------------------------------------------------------------------------
# 5. Summary & Logs 
# ------------------------------------------------------------------------------
cat("Analysis complete. Results saved to (, ):", OUTPUT_DIR, "\n")
cat("Summary of generated files :\n")
print(list.files(OUTPUT_DIR))

# Save session info 
sink(file.path(OUTPUT_DIR, "session_info.txt"))
sessionInfo()
sink()
