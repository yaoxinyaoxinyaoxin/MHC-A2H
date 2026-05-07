# ==============================================================================
# [Script]: 5.2.1.1_visualize_integrated_PPI_template.R
# [Description / ]: 
# Integrates and visualizes Protein-Protein Interaction (PPI) networks from datasets and their unions.
# [Date / ]: 2026-05-04
# [Version / ]: 1.0
# [Usage / ]: 
#   Rscript 5.2.1.1_visualize_integrated_PPI_template.R \
#     --file_aging <path> \
#     --file_hz <path> \
#     --file_ra <path> \
#     --file_union_aging_hz <path> \
#     --file_union_aging_ra <path> \
#     --file_union_hz_ra <path> \
#     --file_union_all <path> \
#     --out_dir <path>
# ==============================================================================
# ==============================================================================
# [Script]: visualize_integrated_PPI.R
# [Method]: (PPI)
# [Step]: (PPI)
# 
# [Function]:
# Conduct functional enrichment analysis (GO/KEGG) and construct Protein-Protein Interaction (PPI) networks.
# 
# [Parameters / ]:
# Standard predefined thresholds.
# 
# [Steps / ]:
#   1. Data loading and initialization / 
#   2. Core analytical execution / 
#   3. Results formatting and output / 
# ==============================================================================

#!/usr/bin/env Rscript

# ---------------------------------------------------------------------------------------------------------------------
# Script Name: visualize_integrated_PPI.R
# Author: AI Assistant
# Date: 2026-02-12
# Version: 2.0
# Description: 
#   This script integrates and visualizes Protein-Protein Interaction (PPI) networks from three datasets (Aging, RA, HZ)
#   and their unions. It generates a combined figure panel with 7 plots:
#   - 3 Single trait PPIs (Edge Bundling)
#   - 3 Pairwise union PPIs (Edge Bundling)
#   - 1 Three-way union PPI (Integrated Network)
#   （、、）-(PPI). 
#   7: 
#   - 3PPI
#   - 3PPI
#   - 1PPI
#
# Usage: 
#   Rscript visualize_integrated_PPI.R
#
# Dependencies:
#   - tidyverse, igraph, ggraph, RColorBrewer, readr, ggrepel, patchwork
#
# Steps:
#   1. Initialize environment and load libraries. . 
#   2. Load all 7 PPI datasets. 7PPI. 
#   3. Process and merge network data. . 
#   4. Define node categories (Specific vs Intersections). （ vs ）. 
#   5. Generate Edge Bundling Plots for the 6 subset networks. 6. 
#   6. Generate Integrated Network Plot for the 3-way union. . 
#   7. Combine plots into a specific layout. . 
#   8. Save outputs. . 
# ---------------------------------------------------------------------------------------------------------------------


# ---------------------------------------------------------------------------------------------------------------------
# 1. Environment Setup 
# ---------------------------------------------------------------------------------------------------------------------

# Clear environment
rm(list = ls())

# Load required libraries
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(igraph)
  library(ggraph)
  library(RColorBrewer)
  library(readr)
  library(ggrepel)
  library(patchwork)
})

# Define command line arguments
option_list <- list(
  make_option(c("--file_aging"), type="character", default=NULL, help="Path to aging PPI TSV"),
  make_option(c("--file_hz"), type="character", default=NULL, help="Path to HZ PPI TSV"),
  make_option(c("--file_ra"), type="character", default=NULL, help="Path to RA PPI TSV"),
  make_option(c("--file_union_aging_hz"), type="character", default=NULL, help="Path to union Aging+HZ PPI TSV"),
  make_option(c("--file_union_aging_ra"), type="character", default=NULL, help="Path to union Aging+RA PPI TSV"),
  make_option(c("--file_union_hz_ra"), type="character", default=NULL, help="Path to union HZ+RA PPI TSV"),
  make_option(c("--file_union_all"), type="character", default=NULL, help="Path to union all three PPI TSV"),
  make_option(c("--out_dir"), type="character", default="./Integrated_PPI_Result", help="Output directory path")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$file_aging) || is.null(opt$file_hz) || is.null(opt$file_ra) || 
    is.null(opt$file_union_aging_hz) || is.null(opt$file_union_aging_ra) || 
    is.null(opt$file_union_hz_ra) || is.null(opt$file_union_all)) {
  print_help(opt_parser)
  stop("Missing required input files. / . ", call.=FALSE)
}

FILE_AGING <- opt$file_aging
FILE_HZ <- opt$file_hz
FILE_RA <- opt$file_ra
FILE_UNION_AGING_HZ <- opt$file_union_aging_hz
FILE_UNION_AGING_RA <- opt$file_union_aging_ra
FILE_UNION_HZ_RA <- opt$file_union_hz_ra
FILE_UNION_ALL <- opt$file_union_all

# Output directory (Timestamped)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTPUT_BASE <- opt$out_dir
OUTPUT_DIR <- file.path(OUTPUT_BASE, paste0(timestamp, "_Integrated_PPI_Analysis"))

if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

# Logs directory
LOG_DIR <- file.path(OUTPUT_DIR, "logs")
if (!dir.exists(LOG_DIR)) dir.create(LOG_DIR, recursive = TRUE)

# Logging function
log_message <- function(message) {
  current_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  formatted_msg <- paste0("[", current_time, "] ", message)
  cat(formatted_msg, "
")
  write(formatted_msg, file = file.path(LOG_DIR, "run_log.txt"), append = TRUE)
}

log_message("Starting integrated PPI visualization (7 datasets). PPI（7）. ")
log_message(paste0("Output Directory: ", OUTPUT_DIR))


# ---------------------------------------------------------------------------------------------------------------------
# 2. Data Loading 
# ---------------------------------------------------------------------------------------------------------------------

log_message("Loading PPI datasets... PPI...")

# Helper to read and standardize
read_ppi <- function(file_path, source_name) {
  if (!file.exists(file_path)) {
    log_message(paste0("Warning: File not found: ", file_path))
    return(NULL)
  }
  
  df <- read_tsv(file_path, show_col_types = FALSE)
  
  # Standardize column names
  if ("#node1" %in% colnames(df)) df <- df %>% rename(node1 = `#node1`)
  
  # Select essential columns
  if(all(c("node1", "node2", "combined_score") %in% colnames(df))) {
    df <- df %>% select(node1, node2, combined_score) %>%
      mutate(source = source_name)
    return(df)
  } else {
    log_message(paste0("Warning: Invalid columns in ", source_name))
    return(NULL)
  }
}

# Load all 7 files
df_aging <- read_ppi(FILE_AGING, "mvAge")
df_hz    <- read_ppi(FILE_HZ, "HZ (UKB)")
df_ra    <- read_ppi(FILE_RA, "RA")

df_union_aging_hz <- read_ppi(FILE_UNION_AGING_HZ, "Union_mvAge_HZ (UKB)")
df_union_aging_ra <- read_ppi(FILE_UNION_AGING_RA, "Union_mvAge_RA")
df_union_hz_ra    <- read_ppi(FILE_UNION_HZ_RA, "Union_HZ (UKB)_RA")

df_union_all <- read_ppi(FILE_UNION_ALL, "Union_All")

log_message("Data loaded successfully. . ")

# ---------------------------------------------------------------------------------------------------------------------
# 3. Define Node Categories for Integrated Plot 
# ---------------------------------------------------------------------------------------------------------------------

log_message("Classifying nodes... ...")

# Extract unique nodes from base sets
nodes_aging <- if(!is.null(df_aging)) unique(c(df_aging$node1, df_aging$node2)) else c()
nodes_hz    <- if(!is.null(df_hz)) unique(c(df_hz$node1, df_hz$node2)) else c()
nodes_ra    <- if(!is.null(df_ra)) unique(c(df_ra$node1, df_ra$node2)) else c()

# Master node list from Union All
all_nodes <- if(!is.null(df_union_all)) unique(c(df_union_all$node1, df_union_all$node2)) else c()

# Classification Logic
node_metadata <- data.frame(name = all_nodes, stringsAsFactors = FALSE) %>%
  mutate(
    in_aging = name %in% nodes_aging,
    in_hz = name %in% nodes_hz,
    in_ra = name %in% nodes_ra
  ) %>%
  mutate(
    Category = case_when(
      in_aging & in_hz & in_ra ~ "Intersection (All 3)",
      in_aging & in_hz ~ "Intersection (mvAge & HZ (UKB))",
      in_aging & in_ra ~ "Intersection (mvAge & RA)",
      in_hz & in_ra    ~ "Intersection (HZ (UKB) & RA)",
      in_aging ~ "mvAge Specific",
      in_hz    ~ "HZ (UKB) Specific",
      in_ra    ~ "RA Specific",
      TRUE     ~ "Other"
    )
  )

log_message("Node Categories Summary:")
print(table(node_metadata$Category))

# Define Colors for Categories
category_colors <- c(
  "Intersection (All 3)"      = "#E41A1C", # Red (Most important)
  "Intersection (mvAge & HZ (UKB))" = "#984EA3", # Purple
  "Intersection (mvAge & RA)" = "#F781BF", # Pink
  "Intersection (HZ (UKB) & RA)"    = "#FF7F00", # Orange
  "mvAge Specific"            = "#377EB8", # Blue
  "HZ (UKB) Specific"               = "#4DAF4A", # Green
  "RA Specific"               = "#A65628", # Brown
  "Other"                     = "#999999"  # Grey
)

# ---------------------------------------------------------------------------------------------------------------------
# 4. Visualization Functions 
# ---------------------------------------------------------------------------------------------------------------------

log_message("Generating plots... ...")

# 4.1 Edge Bundling Plot (For Single Traits)
# -------------------------------------------
create_edge_bundling_plot <- function(edges_df, cluster_colors = NULL) {
  
  if (nrow(edges_df) == 0) return(NULL)

  # Create graph
  g_temp <- graph_from_data_frame(d = edges_df, directed = FALSE)
  g_temp <- simplify(g_temp, edge.attr.comb = "mean")
  
  # Calculate degree and clustering
  V(g_temp)$degree <- degree(g_temp)
  cluster_res <- cluster_louvain(g_temp)
  V(g_temp)$cluster <- as.factor(membership(cluster_res))
  
  # Hierarchy for bundling
  gene_cluster_info <- data.frame(
    gene = V(g_temp)$name,
    cluster = paste0("Cluster_", V(g_temp)$cluster),
    degree = V(g_temp)$degree,
    stringsAsFactors = FALSE
  ) %>% arrange(cluster)
  
  hierarchy_edges <- rbind(
    data.frame(from = "Origin", to = unique(gene_cluster_info$cluster)),
    data.frame(from = gene_cluster_info$cluster, to = gene_cluster_info$gene)
  )
  
  hierarchy_graph <- graph_from_data_frame(hierarchy_edges)
  V(hierarchy_graph)$degree <- 0
  V(hierarchy_graph)$is_leaf <- FALSE
  V(hierarchy_graph)$cluster <- NA
  
  leaf_indices <- match(gene_cluster_info$gene, V(hierarchy_graph)$name)
  V(hierarchy_graph)$degree[leaf_indices] <- gene_cluster_info$degree
  V(hierarchy_graph)$is_leaf[leaf_indices] <- TRUE
  V(hierarchy_graph)$cluster[leaf_indices] <- as.character(gene_cluster_info$cluster)
  
  # Connections
  edges_df_graph <- as_data_frame(g_temp, what = "edges")
  from_id <- match(edges_df_graph$from, V(hierarchy_graph)$name)
  to_id <- match(edges_df_graph$to, V(hierarchy_graph)$name)
  valid_edges <- !is.na(from_id) & !is.na(to_id)
  connect_df <- data.frame(from = from_id[valid_edges], to = to_id[valid_edges])
  
  # Colors
  nb_clusters <- length(unique(V(g_temp)$cluster))
  if (is.null(cluster_colors)) {
    if(nb_clusters > 12) {
      my_colors <- colorRampPalette(brewer.pal(12, "Set3"))(nb_clusters)
    } else {
      my_colors <- brewer.pal(max(3, nb_clusters), "Set3")[1:nb_clusters]
    }
  } else {
    my_colors <- cluster_colors
  }
  
  # Plot
  p <- ggraph(hierarchy_graph, layout = 'dendrogram', circular = TRUE) +
    geom_conn_bundle(data = get_con(from = connect_df$from, to = connect_df$to), 
                     alpha = 0.2, colour = "skyblue", tension = 0.9) +
    geom_node_point(aes(filter = is_leaf, x = x*1.05, y = y*1.05, color = cluster, size = degree), 
                    alpha = 0.8) +
    geom_node_text(aes(x = x*1.1, y = y*1.1, filter = is_leaf, 
                       label = name, angle = -((-node_angle(x, y) + 90) %% 180) + 90, hjust = 'outward'),
                   size = 1.5, fontface = "bold") +
    scale_size(range = c(1, 4)) +
    scale_color_manual(values = my_colors) + 
    coord_fixed() +
    theme_void() +
    theme(
      legend.position = "none", # Hide legend for compact panels
      plot.margin = margin(5, 5, 5, 5)
    )
  return(p)
}

# 4.2 Concentric Layout Plot (For Unions)
# ---------------------------------------
create_concentric_plot <- function(edges_df, title = "", node_metadata_full = NULL) {
  
  if (nrow(edges_df) == 0) return(NULL)
  
  # Create graph
  g_sub <- graph_from_data_frame(d = edges_df, directed = FALSE)
  g_sub <- simplify(g_sub, edge.attr.comb = "mean")
  
  # Add node attributes
  V(g_sub)$degree <- degree(g_sub)
  V(g_sub)$betweenness <- betweenness(g_sub)
  
  # Map categories if metadata provided
  if (!is.null(node_metadata_full)) {
    V(g_sub)$Category <- node_metadata_full$Category[match(V(g_sub)$name, node_metadata_full$name)]
  } else {
    V(g_sub)$Category <- "Other"
  }
  
  # Assign concentric layers based on Betweenness
  # Reverse rank: Higher betweenness -> Lower rank index (Center)
  V(g_sub)$rank_bet <- rank(-V(g_sub)$betweenness, ties.method = "first")
  n_nodes <- vcount(g_sub)
  
  # Define 5 layers based on percentiles
  cuts <- c(0, 0.05, 0.15, 0.35, 0.65, 1.0)
  V(g_sub)$layer <- cut(V(g_sub)$rank_bet / n_nodes, breaks = cuts, labels = FALSE, include.lowest = TRUE)
  
  # Custom layout function for concentric circles
  create_concentric_layout <- function(g) {
    layout_df <- data.frame(name = V(g)$name, layer = V(g)$layer)
    layout_matrix <- matrix(0, nrow = vcount(g), ncol = 2)
    
    for (l in 1:max(V(g)$layer)) {
      idx <- which(V(g)$layer == l)
      n_layer <- length(idx)
      if (n_layer > 0) {
        radius <- l  # Radius increases with layer index
        # Distribute evenly
        angles <- seq(0, 2*pi, length.out = n_layer + 1)[1:n_layer]
        layout_matrix[idx, 1] <- radius * cos(angles)
        layout_matrix[idx, 2] <- radius * sin(angles)
      }
    }
    return(layout_matrix)
  }
  
  layout_concentric <- create_concentric_layout(g_sub)
  
  # Determine top nodes for labeling
  threshold_degree <- quantile(V(g_sub)$degree, 0.90)
  V(g_sub)$is_hub <- V(g_sub)$degree >= threshold_degree
  
  p <- ggraph(g_sub, layout = layout_concentric) +
    geom_edge_arc(aes(alpha = combined_score), color = "grey80", strength = 0.1, width = 0.3) +
    geom_node_point(aes(fill = Category, size = degree), shape = 21, color = "white", stroke = 0.2) +
    geom_node_text(aes(label = ifelse(is_hub, name, "")),
                   repel = TRUE, size = 2.0, fontface = "bold", 
                   bg.color = "white", bg.r = 0.1, max.overlaps = 50,
                   color = "black") +
    scale_fill_manual(values = category_colors) +
    scale_size(range = c(1, 5)) +
    scale_edge_alpha(range = c(0.1, 0.6), guide = "none") +
    coord_fixed() +
    theme_void() +
    theme(
      legend.position = "none",
      plot.margin = margin(5, 5, 5, 5)
    )
    
  return(p)
}

# ---------------------------------------------------------------------------------------------------------------------
# 5. Generate Plots 
# ---------------------------------------------------------------------------------------------------------------------

log_message("Generating individual plots... ...")

# 1. Single Trait Plots (Edge Bundling) - Keep logic
p_aging <- create_edge_bundling_plot(df_aging) + ggtitle("mvAge")
p_hz <- create_edge_bundling_plot(df_hz) + ggtitle("HZ (UKB)")
p_ra <- create_edge_bundling_plot(df_ra) + ggtitle("RA")

# 2. Pairwise Union Plots (Concentric) - Changed from Edge Bundling to Concentric
p_union_aging_hz <- create_concentric_plot(df_union_aging_hz, "Union (mvAge & HZ (UKB))", node_metadata) + ggtitle("Union: mvAge & HZ (UKB)")
p_union_aging_ra <- create_concentric_plot(df_union_aging_ra, "Union (mvAge & RA)", node_metadata) + ggtitle("Union: mvAge & RA")
p_union_hz_ra <- create_concentric_plot(df_union_hz_ra, "Union (HZ (UKB) & RA)", node_metadata) + ggtitle("Union: HZ (UKB) & RA")

# 3. Integrated Three-Way Plot (Concentric) - Changed from FR to Concentric (Big Plot)
log_message("Generating integrated network plot (Concentric)... ...")

# Construct graph object for the full union
g_all <- graph_from_data_frame(d = df_union_all, vertices = node_metadata, directed = FALSE)
g_all <- simplify(g_all, edge.attr.comb = "mean") 
V(g_all)$degree <- degree(g_all)
V(g_all)$betweenness <- betweenness(g_all)

# Concentric layout for big plot
V(g_all)$rank_bet <- rank(-V(g_all)$betweenness, ties.method = "first")
n_nodes_all <- vcount(g_all)
cuts_all <- c(0, 0.05, 0.15, 0.35, 0.65, 1.0)
V(g_all)$layer <- cut(V(g_all)$rank_bet / n_nodes_all, breaks = cuts_all, labels = FALSE, include.lowest = TRUE)

create_concentric_layout_all <- function(g) {
    layout_df <- data.frame(name = V(g)$name, layer = V(g)$layer)
    layout_matrix <- matrix(0, nrow = vcount(g), ncol = 2)
    for (l in 1:max(V(g)$layer)) {
      idx <- which(V(g)$layer == l)
      n_layer <- length(idx)
      if (n_layer > 0) {
        radius <- l 
        angles <- seq(0, 2*pi, length.out = n_layer + 1)[1:n_layer]
        layout_matrix[idx, 1] <- radius * cos(angles)
        layout_matrix[idx, 2] <- radius * sin(angles)
      }
    }
    return(layout_matrix)
}
layout_concentric_all <- create_concentric_layout_all(g_all)

# Identify Hubs for labeling
threshold_degree_all <- quantile(V(g_all)$degree, 0.90)
V(g_all)$is_hub <- V(g_all)$degree >= threshold_degree_all

p_union_all <- ggraph(g_all, layout = layout_concentric_all) +
  geom_edge_arc(aes(alpha = combined_score), color = "grey80", strength = 0.1, width = 0.3) +
  geom_node_point(aes(fill = Category, size = degree), shape = 21, color = "white", stroke = 0.2) +
  geom_node_text(aes(label = ifelse(Category == "Intersection (All 3)" | is_hub, name, "")),
                 repel = TRUE, size = 3.0, fontface = "bold", 
                 bg.color = "white", bg.r = 0.1, max.overlaps = 50,
                 color = "black") +
  scale_fill_manual(values = category_colors) +
  scale_size(range = c(2, 8)) +
  scale_edge_alpha(range = c(0.1, 0.6), guide = "none") +
  coord_fixed() +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.box = "vertical",
    plot.margin = margin(10, 10, 10, 10)
  ) +
  labs(
    fill = "Gene Category",
    size = "Degree",
    title = "Three-Way Integrated Network (Concentric)"
  ) +
  guides(
    fill = guide_legend(override.aes = list(size = 4)),
    size = guide_legend(override.aes = list(fill = "grey50", color = "black"))
  )

# ---------------------------------------------------------------------------------------------------------------------
# 6. Combine and Save 
# ---------------------------------------------------------------------------------------------------------------------

log_message("Combining plots... ...")

# Use Patchwork to create the layout
# Layout design:
# Left Panel (2x3):
#   [Aging]         [RA]             [HZ]
#   [Union A+HZ]    [Union A+RA]     [Union HZ+RA]
# Right Panel (1x1):
#   [Union All 3 (Big)]

left_panel <- (p_aging | p_ra | p_hz) / (p_union_aging_hz | p_union_aging_ra | p_union_hz_ra)
combined_plot <- left_panel | p_union_all 

combined_plot <- combined_plot + 
  plot_layout(widths = c(1.8, 1)) # Left side wider because it has 3 columns

# Create vertical layout version
combined_plot_vertical <- wrap_elements(left_panel) / p_union_all + 
  plot_layout(heights = c(1, 1.2)) # Adjust height ratio for better visualization

log_message("Saving combined figure (Horizontal layout)... ...")
outfile_pdf <- file.path(OUTPUT_DIR, "Integrated_PPI_Panel_7_Plots_Horizontal.pdf")
outfile_png <- file.path(OUTPUT_DIR, "Integrated_PPI_Panel_7_Plots_Horizontal.png")

# Save Horizontal as PDF
ggsave(outfile_pdf, plot = combined_plot, width = 24, height = 16, device = "pdf")

# Save Horizontal as PNG
ggsave(outfile_png, plot = combined_plot, width = 24, height = 16, dpi = 300)

log_message("Saving combined figure (Vertical layout)... ...")
outfile_vertical_pdf <- file.path(OUTPUT_DIR, "Integrated_PPI_Panel_7_Plots_Vertical.pdf")
outfile_vertical_png <- file.path(OUTPUT_DIR, "Integrated_PPI_Panel_7_Plots_Vertical.png")

# Save Vertical as PDF
ggsave(outfile_vertical_pdf, plot = combined_plot_vertical, width = 16, height = 24, device = "pdf")

# Save Vertical as PNG
ggsave(outfile_vertical_png, plot = combined_plot_vertical, width = 16, height = 24, dpi = 300)

log_message("Saving individual plots (PNG and PDF)... （PNGPDF）...")
ggsave(file.path(OUTPUT_DIR, "1_mvAge_PPI.png"), plot = p_aging, width = 8, height = 8)
ggsave(file.path(OUTPUT_DIR, "1_mvAge_PPI.pdf"), plot = p_aging, width = 8, height = 8, device = "pdf")

ggsave(file.path(OUTPUT_DIR, "2_HZ (UKB)_PPI.png"), plot = p_hz, width = 8, height = 8)
ggsave(file.path(OUTPUT_DIR, "2_HZ (UKB)_PPI.pdf"), plot = p_hz, width = 8, height = 8, device = "pdf")

ggsave(file.path(OUTPUT_DIR, "3_RA_PPI.png"), plot = p_ra, width = 8, height = 8)
ggsave(file.path(OUTPUT_DIR, "3_RA_PPI.pdf"), plot = p_ra, width = 8, height = 8, device = "pdf")

ggsave(file.path(OUTPUT_DIR, "4_Union_mvAge_HZ (UKB)_PPI.png"), plot = p_union_aging_hz, width = 8, height = 8)
ggsave(file.path(OUTPUT_DIR, "4_Union_mvAge_HZ (UKB)_PPI.pdf"), plot = p_union_aging_hz, width = 8, height = 8, device = "pdf")

ggsave(file.path(OUTPUT_DIR, "5_Union_mvAge_RA_PPI.png"), plot = p_union_aging_ra, width = 8, height = 8)
ggsave(file.path(OUTPUT_DIR, "5_Union_mvAge_RA_PPI.pdf"), plot = p_union_aging_ra, width = 8, height = 8, device = "pdf")

ggsave(file.path(OUTPUT_DIR, "6_Union_HZ (UKB)_RA_PPI.png"), plot = p_union_hz_ra, width = 8, height = 8)
ggsave(file.path(OUTPUT_DIR, "6_Union_HZ (UKB)_RA_PPI.pdf"), plot = p_union_hz_ra, width = 8, height = 8, device = "pdf")

ggsave(file.path(OUTPUT_DIR, "7_Union_All_Network.png"), plot = p_union_all, width = 12, height = 12)
ggsave(file.path(OUTPUT_DIR, "7_Union_All_Network.pdf"), plot = p_union_all, width = 12, height = 12, device = "pdf")

log_message(paste0("All plots saved to: ", OUTPUT_DIR))

log_message("Analysis completed successfully. . ")
