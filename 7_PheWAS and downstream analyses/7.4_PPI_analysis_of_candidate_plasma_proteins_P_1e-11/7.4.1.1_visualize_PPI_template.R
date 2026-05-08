#!/usr/bin/env Rscript
# ==============================================================================
# [Script]: 7.4.1.1_visualize_PPI_template.R
# [Method]: PPI Network Analysis
# [Step]: Network Construction and Visualization
# 
# [Function]:
# Construct and visualize Protein-Protein Interaction (PPI) networks from STRING data.
# 
# [Input]:
#   --input_file   : Path to STRING TSV file
#   --out_dir      : Output directory path
# 
# [Output]:
# PPI network plots (PDF/PNG) and network statistics.
# ==============================================================================

rm(list = ls())
gc()

# 1. Load Dependencies
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(igraph)
  library(ggraph)
  library(RColorBrewer)
})

# 2. Parse Command Line Arguments
option_list <- list(
  make_option(c("--input_file"), type="character", default=NULL, help="Input STRING TSV file"),
  make_option(c("--out_dir"), type="character", default="./PPI_Visualization", help="Output directory path")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_file)) {
  print_help(opt_parser)
  stop("Missing required input file.")
}

INPUT_FILE <- opt$input_file
OUTPUT_DIR <- opt$out_dir

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

INPUT_FILE <- "./"

# Generate timestamped output directory
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTPUT_BASE <- file.path("./", paste0(timestamp, "_PPI_Visualization"))
OUTPUT_ANALYSIS <- file.path(OUTPUT_BASE, "analysis")
OUTPUT_LOGS <- file.path(OUTPUT_BASE, "logs")
OUTPUT_README <- file.path(OUTPUT_BASE, "readme")
OUTPUT_SCRIPTS <- file.path(OUTPUT_BASE, "scripts")

# Create directories
if (!dir.exists(OUTPUT_ANALYSIS)) dir.create(OUTPUT_ANALYSIS, recursive = TRUE)
if (!dir.exists(OUTPUT_LOGS)) dir.create(OUTPUT_LOGS, recursive = TRUE)
if (!dir.exists(OUTPUT_README)) dir.create(OUTPUT_README, recursive = TRUE)
if (!dir.exists(OUTPUT_SCRIPTS)) dir.create(OUTPUT_SCRIPTS, recursive = TRUE)

# Set CRAN mirror for stability
# CRAN
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# Function to log messages
log_message <- function(message) {
  current_time <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  formatted_msg <- paste0("[", current_time, "] ", message)
  cat(formatted_msg, "\n")
  write(formatted_msg, file = file.path(OUTPUT_LOGS, "run_log.txt"), append = TRUE)
}

# Start logging
log_message("Starting PPI visualization script. PPI. ")
log_message(paste0("Input File: ", INPUT_FILE))
log_message(paste0("Output Directory: ", OUTPUT_BASE))

# Backup current script to output scripts folder
current_script_path <- "./"
if (file.exists(current_script_path)) {
  file.copy(current_script_path, file.path(OUTPUT_SCRIPTS, paste0("visualize_PPI_", timestamp, ".R")))
  log_message("Script copied to output scripts directory. . ")
}

# Load required libraries with error handling
# , 
load_lib <- function(pkg) {
  if (!require(pkg, character.only = TRUE)) {
    log_message(paste0("Installing package: ", pkg))
    install.packages(pkg)
    if (!require(pkg, character.only = TRUE)) {
      stop(paste0("Failed to install/load package: ", pkg))
    }
  }
  log_message(paste0("Loaded package: ", pkg))
}

libs <- c("tidyverse", "igraph", "ggraph", "RColorBrewer", "readr", "ggrepel")
sapply(libs, load_lib)

# ---------------------------------------------------------------------------------------------------------------------
# 2. Data Loading 
# ---------------------------------------------------------------------------------------------------------------------

log_message("Reading input file... ...")

# Check if file exists
if (!file.exists(INPUT_FILE)) {
  stop(paste0("Input file not found: ", INPUT_FILE))
}

# Read data using read_tsv
# read_tsv
ppi_data <- read_tsv(INPUT_FILE, show_col_types = FALSE)

# Check column mapping and rename #node1 to node1 if present
# , #node1, node1
if ("#node1" %in% colnames(ppi_data)) {
  ppi_data <- ppi_data %>% rename(node1 = `#node1`)
  log_message("Renamed column '#node1' to 'node1'.  '#node1'  'node1'. ")
}

# Check for essential columns
required_cols <- c("node1", "node2", "combined_score")
missing_cols <- setdiff(required_cols, colnames(ppi_data))
if (length(missing_cols) > 0) {
  stop(paste0("Missing required columns: ", paste(missing_cols, collapse = ", ")))
}

# Generate column mapping record
column_mapping_file <- file.path(OUTPUT_README, "column_mapping_record.txt")
writeLines(paste("Original columns:", paste(colnames(ppi_data), collapse=", ")), column_mapping_file)
log_message("Column mapping verified and recorded. . ")

# Check data structure
log_message("Data loaded successfully. Head of data: . : ")
print(head(ppi_data))

# ---------------------------------------------------------------------------------------------------------------------
# 3. Data Processing & Network Construction 
# ---------------------------------------------------------------------------------------------------------------------

log_message("Constructing network object. . ")

# Select relevant columns: node1, node2, combined_score
# : node1, node2, combined_score
edges <- ppi_data %>%
  select(node1, node2, combined_score) %>%
  filter(combined_score > 0) # Ensure score is positive 0

# Create igraph object
# igraph
g <- graph_from_data_frame(d = edges, directed = FALSE)

# Simplify graph (remove loops and multiple edges)
g <- simplify(g, edge.attr.comb = "mean")

# Calculate topological properties
log_message("Calculating network topology and clusters. . ")

V(g)$degree <- degree(g)
V(g)$betweenness <- betweenness(g)
V(g)$closeness <- closeness(g)

# Read the input gene list to get the direction information earlier
gene_list_file <- "./"
if (file.exists(gene_list_file)) {
  gene_list_data <- read_csv(gene_list_file, show_col_types = FALSE)
  # Map direction to nodes
  node_names <- V(g)$name
  direction_map <- gene_list_data$direction[match(node_names, gene_list_data$gene)]
  V(g)$Plasma_Protein_Level <- ifelse(is.na(direction_map), "Unknown", direction_map)
  log_message("Mapped Plasma_Protein_Level to network nodes. . ")
} else {
  V(g)$Plasma_Protein_Level <- "Unknown"
  log_message(paste0("Warning: Gene list file not found: ", gene_list_file, " : . "))
}

# Community Detection (Clustering) for Cytoscape-like coloring
# cluster_louvain
#  combined_score, 
set.seed(123) # Ensure reproducibility for heuristic clustering / 
E(g)$weight <- E(g)$combined_score
cluster_res <- cluster_louvain(g)

# Calculate Modularity to evaluate clustering reliability
# (Modularity)
mod_score <- modularity(cluster_res)
log_message(paste0("Clustering Modularity Score : ", round(mod_score, 4)))
if (mod_score > 0.4) {
  log_message("Assessment: Modularity > 0.4. Very strong community structure (highly reliable). / : >0.4, , . ")
} else if (mod_score > 0.3) {
  log_message("Assessment: Modularity > 0.3. Good community structure (reliable). / : >0.3, , . ")
} else {
  log_message("Assessment: Modularity <= 0.3. Community structure might be weak. / : <=0.3, . ")
}

V(g)$cluster <- as.factor(membership(cluster_res))
log_message(paste0("Identified ", length(unique(V(g)$cluster)), " clusters."))

# Identify hub genes (e.g., top 10% by degree)
# Hub（, 10%）
threshold_degree <- quantile(V(g)$degree, 0.90)
V(g)$type <- ifelse(V(g)$degree >= threshold_degree, "Hub", "Non-Hub")

# Save network statistics to readme folder
# readme
node_stats <- data.frame(
  Gene = V(g)$name,
  Degree = V(g)$degree,
  Betweenness = V(g)$betweenness,
  Closeness = V(g)$closeness,
  Cluster = V(g)$cluster,
  Type = V(g)$type
) %>% arrange(desc(Degree))

# Read the input gene list to get the direction information
gene_list_file <- "./"
if (file.exists(gene_list_file)) {
  gene_list_data <- read_csv(gene_list_file, show_col_types = FALSE)
  # Merge direction into node_stats and rename it to Plasma_Protein_Level
  # node_stats, Plasma_Protein_Level
  node_stats <- node_stats %>%
    left_join(gene_list_data %>% select(gene, direction), by = c("Gene" = "gene")) %>%
    rename(Plasma_Protein_Level = direction)
}

stats_file <- file.path(OUTPUT_README, "network_statistics.csv")
write_csv(node_stats, stats_file)
log_message(paste0("Saved network statistics to: ", stats_file))

# ---------------------------------------------------------------------------------------------------------------------
# 4. Visualization 
# ---------------------------------------------------------------------------------------------------------------------

log_message("Generating network plots. . ")

# Set random seed for reproducibility
set.seed(123)

# Define layout (Fruchterman-Reingold is standard for PPI)
# （Fruchterman-ReingoldPPI）
layout <- create_layout(g, layout = "fr") 

# Define a custom palette with distinct colors (Cytoscape style often uses distinct colors for clusters)
# （Cytoscape）
nb_clusters <- length(unique(V(g)$cluster))
if(nb_clusters > 12) {
    cluster_colors <- colorRampPalette(brewer.pal(12, "Set3"))(nb_clusters)
} else {
    cluster_colors <- brewer.pal(max(3, nb_clusters), "Set3")[1:nb_clusters]
}

# Adjust colors for better visibility
# Replace #FFFFB3 (Light Yellow) with #FDBF6F (Light Orange)
if (any(cluster_colors == "#FFFFB3")) {
  cluster_colors[cluster_colors == "#FFFFB3"] <- "#FDBF6F"
}
# Replace #FFED6F (another yellow) with #CAB2D6 (Light Purple)
if (any(cluster_colors == "#FFED6F")) {
  cluster_colors[cluster_colors == "#FFED6F"] <- "#CAB2D6" 
}
# Adjust color conflict if necessary
if (length(cluster_colors) >= 6) {
    cluster_colors[6] <- "#6A3D9A" # Dark Purple for better contrast
}

# Plot 3: Cytoscape Style (Cluster Coloring + Node Borders)
# 3: Cytoscape（ + ）
p3 <- ggraph(layout) +
  # Edges: Grey, slightly transparent, width proportional to score? Or fixed.
  # : , 
  geom_edge_link(aes(alpha = combined_score), color = "grey60", width = 0.3) +
  
  # Nodes: Filled circle with border (shape=21), fill by cluster, size by degree
  # : （shape=21）, , 
  geom_node_point(aes(fill = cluster, size = degree), shape = 21, color = "black", stroke = 0.5) +
  
  # Color scale for fill
  scale_fill_manual(values = cluster_colors) +
  
  # Size scale
  scale_size(range = c(3, 10)) +
  
  # Labels for All genes
  geom_node_text(aes(label = name), 
                 repel = TRUE, size = 2.5, fontface = "bold", family = "sans", bg.color = "white", bg.r = 0.15, max.overlaps = 50) +
  
  theme_void(base_family = "sans") +
  theme(
    legend.position = "right",
    legend.title = element_text(face = "bold"),
    text = element_text(family = "sans")
  ) +
  labs(
    fill = "Cluster",
    size = "Degree",
    edge_alpha = "Confidence"
  )

# Create a bottom legend version
p3_bottom <- p3 + theme(legend.position = "bottom")

# Plot 4: Concentric Circles by Betweenness Centrality
# 4: 

# Assign concentric layers based on Betweenness
# Reverse rank: Higher betweenness -> Lower rank index (Center)
V(g)$rank_bet <- rank(-V(g)$betweenness, ties.method = "first")
n_nodes <- vcount(g)

# Define 5 layers based on percentiles
# 5
cuts <- c(0, 0.05, 0.15, 0.35, 0.65, 1.0)
V(g)$layer <- cut(V(g)$rank_bet / n_nodes, breaks = cuts, labels = FALSE, include.lowest = TRUE)

# Custom layout function for concentric circles
create_concentric_layout <- function(g) {
  layout_df <- data.frame(name = V(g)$name, layer = V(g)$layer, direction = V(g)$Plasma_Protein_Level)
  layout_matrix <- matrix(0, nrow = vcount(g), ncol = 2)
  
  for (l in 1:max(V(g)$layer)) {
    idx <- which(V(g)$layer == l)
    n_layer <- length(idx)
    if (n_layer > 0) {
      radius <- l  # Radius increases with layer index
      
      # UpDown
      # Separate Up and Down nodes
      idx_up <- idx[layout_df$direction[idx] == "Up"]
      idx_down <- idx[layout_df$direction[idx] == "Down"]
      idx_other <- idx[!layout_df$direction[idx] %in% c("Up", "Down")]
      
      # Assign angles for different groups
      # Up: 0  180  (0  pi) -> 
      # Down: 180  360  (pi  2*pi) -> 
      
      n_up <- length(idx_up)
      n_down <- length(idx_down)
      n_other <- length(idx_other)
      
      if (n_up > 0) {
        angles_up <- seq(0, pi, length.out = n_up + 2)[2:(n_up + 1)]
        layout_matrix[idx_up, 1] <- radius * cos(angles_up)
        layout_matrix[idx_up, 2] <- radius * sin(angles_up)
      }
      
      if (n_down > 0) {
        angles_down <- seq(pi, 2*pi, length.out = n_down + 2)[2:(n_down + 1)]
        layout_matrix[idx_down, 1] <- radius * cos(angles_down)
        layout_matrix[idx_down, 2] <- radius * sin(angles_down)
      }
      
      if (n_other > 0) {
        # , , 
        angles_other <- seq(0, 2*pi, length.out = n_other + 1)[1:n_other] + (pi/4)
        layout_matrix[idx_other, 1] <- radius * cos(angles_other)
        layout_matrix[idx_other, 2] <- radius * sin(angles_other)
      }
    }
  }
  return(layout_matrix)
}

layout_concentric <- create_concentric_layout(g)

p4 <- ggraph(g, layout = layout_concentric) +
  # Edges: curved for aesthetics in concentric layout
  # : 
  geom_edge_arc(aes(alpha = combined_score), color = "grey70", strength = 0.1, width = 0.2) +
  
  # Nodes
  geom_node_point(aes(fill = betweenness, size = degree), shape = 21, color = "black", stroke = 0.3) +
  
  # Color scale (Betweenness: High=Red, Low=Blue)
  # （: =, =）
  scale_fill_gradient(low = "#377EB8", high = "#E41A1C") +
  
  # Size scale
  scale_size(range = c(2, 8)) +
  
  # Labels for All genes
  geom_node_text(aes(label = name), 
                 repel = TRUE, size = 2.5, fontface = "bold", family = "sans", bg.color = "white", bg.r = 0.15, max.overlaps = 50) +
  
  coord_fixed() +
  theme_void(base_family = "sans") +
  theme(
    legend.position = "right",
    text = element_text(family = "sans")
  ) +
  labs(
    fill = "Betweenness",
    size = "Degree",
    edge_alpha = "Confidence"
  )

# Create a bottom legend version
p4_bottom <- p4 + theme(legend.position = "bottom")

# Plot 5: Hierarchical Edge Bundling (Chord-like)
# 5: 

log_message("Generating Hierarchical Edge Bundling plot. . ")

# Create a hierarchical structure: Origin -> Cluster -> Gene
# : ->  -> 
# This is required for hierarchical edge bundling in ggraph

# 1. Prepare hierarchy data frame
# 1. 
# Prepend "Cluster_" to cluster IDs to avoid name collisions with genes
# ID"Cluster_"
gene_cluster_info <- data.frame(
  gene = V(g)$name,
  cluster = paste0("Cluster_", V(g)$cluster),
  degree = V(g)$degree,
  stringsAsFactors = FALSE
) %>% arrange(cluster)

# Construct hierarchy edges
hierarchy_edges <- rbind(
  data.frame(from = "Origin", to = unique(gene_cluster_info$cluster)),
  data.frame(from = gene_cluster_info$cluster, to = gene_cluster_info$gene)
)

# Create hierarchy igraph object
# igraph
hierarchy_graph <- graph_from_data_frame(hierarchy_edges)

# Add node info (degree, is_leaf) to the hierarchy graph
# （, ）
V(hierarchy_graph)$degree <- 0
V(hierarchy_graph)$is_leaf <- FALSE
V(hierarchy_graph)$cluster <- NA

# Map info to leaf nodes (genes)
leaf_indices <- match(gene_cluster_info$gene, V(hierarchy_graph)$name)
V(hierarchy_graph)$degree[leaf_indices] <- gene_cluster_info$degree
V(hierarchy_graph)$is_leaf[leaf_indices] <- TRUE
V(hierarchy_graph)$cluster[leaf_indices] <- as.character(gene_cluster_info$cluster)

# 2. Prepare connections (the actual PPI edges)
# 2. （PPI）
edges_df <- as_data_frame(g, what = "edges")
from_id <- match(edges_df$from, V(hierarchy_graph)$name)
to_id <- match(edges_df$to, V(hierarchy_graph)$name)

# Filter out invalid edges
valid_edges <- !is.na(from_id) & !is.na(to_id)
connect_df <- data.frame(from = from_id[valid_edges], to = to_id[valid_edges])

# Create the Edge Bundling Plot
p5 <- ggraph(hierarchy_graph, layout = 'dendrogram', circular = TRUE) +
  # Bundled edges
  geom_conn_bundle(data = get_con(from = connect_df$from, to = connect_df$to), 
                   alpha = 0.2, colour = "skyblue", tension = 0.9) +
  
  # Nodes (only for genes/leaves)
  # （/）
  geom_node_point(aes(filter = is_leaf, x = x*1.05, y = y*1.05, color = cluster, size = degree), 
                  alpha = 0.8) +
  
  # Node Labels (All genes, oriented outwards)
  # （, ）
  geom_node_text(aes(x = x*1.1, y = y*1.1, filter = is_leaf, 
                     label = name, angle = -((-node_angle(x, y) + 90) %% 180) + 90, hjust = 'outward'),
                 size = 2.0, fontface = "bold", family = "sans") +
  
  scale_size(range = c(1, 5)) +
  scale_color_manual(values = cluster_colors) + 
  
  coord_fixed() +
  theme_void(base_family = "sans") +
  theme(
    legend.position = "right",
    plot.margin = margin(50, 50, 50, 50), # Extra margin for labels
    text = element_text(family = "sans")
  ) +
  labs(
    color = "Cluster",
    size = "Degree"
  )

# Create a bottom legend version
p5_bottom <- p5 + theme(legend.position = "bottom")

# ---------------------------------------------------------------------------------------------------------------------
# 5. Save Outputs 
# ---------------------------------------------------------------------------------------------------------------------

log_message("Saving plots... ...")

# Save plots
plot3_file <- file.path(OUTPUT_ANALYSIS, "PPI_Network_Cytoscape_Style.pdf")
plot3_png <- file.path(OUTPUT_ANALYSIS, "PPI_Network_Cytoscape_Style.png")
plot3_bottom_file <- file.path(OUTPUT_ANALYSIS, "PPI_Network_Cytoscape_Style_bottom_legend.pdf")
plot3_bottom_png <- file.path(OUTPUT_ANALYSIS, "PPI_Network_Cytoscape_Style_bottom_legend.png")

plot4_file <- file.path(OUTPUT_ANALYSIS, "PPI_Network_Concentric.pdf")
plot4_png <- file.path(OUTPUT_ANALYSIS, "PPI_Network_Concentric.png")
plot4_bottom_file <- file.path(OUTPUT_ANALYSIS, "PPI_Network_Concentric_bottom_legend.pdf")
plot4_bottom_png <- file.path(OUTPUT_ANALYSIS, "PPI_Network_Concentric_bottom_legend.png")

plot5_file <- file.path(OUTPUT_ANALYSIS, "PPI_Network_EdgeBundling.pdf")
plot5_png <- file.path(OUTPUT_ANALYSIS, "PPI_Network_EdgeBundling.png")
plot5_bottom_file <- file.path(OUTPUT_ANALYSIS, "PPI_Network_EdgeBundling_bottom_legend.pdf")
plot5_bottom_png <- file.path(OUTPUT_ANALYSIS, "PPI_Network_EdgeBundling_bottom_legend.png")

# Save p3
ggsave(plot3_file, plot = p3, width = 12, height = 10, device = "pdf")
ggsave(plot3_png, plot = p3, width = 12, height = 10, dpi = 600, bg = "white")
ggsave(plot3_bottom_file, plot = p3_bottom, width = 12, height = 11, device = "pdf")
ggsave(plot3_bottom_png, plot = p3_bottom, width = 12, height = 11, dpi = 600, bg = "white")

# Save p4
ggsave(plot4_file, plot = p4, width = 10, height = 10, device = "pdf")
ggsave(plot4_png, plot = p4, width = 10, height = 10, dpi = 600, bg = "white")
ggsave(plot4_bottom_file, plot = p4_bottom, width = 10, height = 11, device = "pdf")
ggsave(plot4_bottom_png, plot = p4_bottom, width = 10, height = 11, dpi = 600, bg = "white")

# Save p5
ggsave(plot5_file, plot = p5, width = 12, height = 12, device = "pdf")
ggsave(plot5_png, plot = p5, width = 12, height = 12, dpi = 600, bg = "white")
ggsave(plot5_bottom_file, plot = p5_bottom, width = 12, height = 13, device = "pdf")
ggsave(plot5_bottom_png, plot = p5_bottom, width = 12, height = 13, dpi = 600, bg = "white")

log_message(paste0("Plots saved to: ", OUTPUT_ANALYSIS))

# ---------------------------------------------------------------------------------------------------------------------
# 6. Wrap up 
# ---------------------------------------------------------------------------------------------------------------------

# Generate README file
# README
readme_file <- file.path(OUTPUT_README, "README_PPI_Visualization.md")
readme_content <- c(
  "# PPI Network Visualization Analysis Report / PPI",
  "",
  paste("- **Date / :**", Sys.Date),
  "- **Author / :** Yaoxin",
  paste("- **Input File / :**", INPUT_FILE),
  paste("- **Output Directory / :**", OUTPUT_BASE),
  "",
  "## Description / ",
  "This directory contains the results of Protein-Protein Interaction (PPI) network visualization.",
  "-（PPI）. ",
  "",
  "## Directory Structure / ",
  "- `analysis/`: Contains network plots (PDF/PNG). （PDF/PNG）. ",
  "- `logs/`: Contains execution logs and session info. . ",
  "- `readme/`: Contains network statistics, column mappings, and this README file. 、. ",
  "- `scripts/`: Contains a copy of the script used for this analysis. . ",
  "",
  "## Analysis Results / ",
  "Network plots successfully generated in multiple layouts (Cytoscape Style, Concentric, Edge Bundling).",
  "（Cytoscape、、）. ",
  "Network statistics including Degree, Betweenness, and Closeness are saved in `network_statistics.csv`.",
  "、 `network_statistics.csv` . "
)
writeLines(readme_content, readme_file)
log_message("README file generated. README. ")

# Generate Image Description MD
# MD
image_desc_file <- file.path(OUTPUT_README, "Image_Description.md")
image_desc_content <- c(
  "# PPI Network Visualizations Description / PPI",
  "",
  "-（PPI）. ",
  "This document describes the three generated visualizations of the Protein-Protein Interaction (PPI) network and their biological and graph-theoretical significance.",
  "",
  "## 1. PPI_Network_Cytoscape_Style",
  "- **Layout :** Fruchterman-Reingold (, ). ",
  "- **Nodes :** /. (Degree), Hub. ",
  "- **Colors :** (Cluster). Louvain, . (Modularity), . ",
  "- **Edges :** . STRING(combined_score), . ",
  "",
  "## 2. PPI_Network_Concentric",
  "- **Layout :** Concentric Circles . ",
  "- **Nodes :** ** (Betweenness Centrality)**. , “”. ",
  "- **Colors :** , . ",
  "- **Size :** (Degree). ",
  "",
  "## 3. PPI_Network_EdgeBundling",
  "- **Layout :** Hierarchical Edge Bundling , . ",
  "- **Grouping :** LouvainCluster, . ",
  "- **Edges :** , . Cluster. "
)
writeLines(image_desc_content, image_desc_file)
log_message("Image description MD file generated. MD. ")

# Calculate run time
end_time <- Sys.time()
run_time <- round(difftime(end_time, start_time, units = "mins"), 2)
log_message(paste0("Total Run Time: ", run_time, " mins. : ", run_time, " . "))

log_message("Analysis completed successfully. . ")
log_message(paste0("Outputs are in: ", OUTPUT_BASE))

# Print session info for reproducibility
session_info_file <- file.path(OUTPUT_LOGS, "session_info.txt")
writeLines(capture.output(sessionInfo()), session_info_file)
log_message(paste0("Session info saved to: ", session_info_file))
