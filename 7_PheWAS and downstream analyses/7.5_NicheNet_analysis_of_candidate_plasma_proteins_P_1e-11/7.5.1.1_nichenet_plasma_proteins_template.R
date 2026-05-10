#!/usr/bin/env Rscript
# ==============================================================================
# [Script]: 7.5.1.1_nichenet_plasma_proteins_template.R
# [Method]: NicheNet Analysis
# [Step]: NicheNet Network Prediction and Visualization
# 
# [Function]:
# Predict intercellular communication (ligand-receptor-target) using NicheNet.
# 
# [Input]:
#   --gene_list    : CSV file containing the list of genes/proteins
#   --nichenet_dir : Path to directory containing NicheNet prior models (e.g. lr_network.rds)
#   --out_dir      : Output directory path
# 
# [Output]:
# NicheNet predictions (CSV), Chord diagrams, Network graphs, Heatmaps.
# ==============================================================================

rm(list = ls())
gc()

# 1. Load Dependencies
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(nichenetr)
  library(tidyverse)
  library(scales)
  library(igraph)
  library(ggraph)
  library(circlize)
  library(ComplexHeatmap)
})

# 2. Parse Command Line Arguments
option_list <- list(
  make_option(c("--gene_list"), type="character", default=NULL, help="Input CSV with gene/protein list"),
  make_option(c("--nichenet_dir"), type="character", default=NULL, help="Directory containing NicheNet RDS models"),
  make_option(c("--out_dir"), type="character", default="./NicheNet_Results", help="Output directory path")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$gene_list) || is.null(opt$nichenet_dir)) {
  print_help(opt_parser)
  stop("Missing required input arguments.")
}

INPUT_GENES <- opt$gene_list
NICHENET_DIR <- opt$nichenet_dir
OUTPUT_DIR <- opt$out_dir

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

#  (Resolve global variable binding warnings for lintr)
invisible(utils::globalVariables(c(
  "Ligand", "Receptor", "Target_Gene", "Regulatory_Potential",
  "name", "Role", "Edge_Weight", "Edge_Type", "Degree",
  "Max_Score", "Mean_Score", "from", "to", "n",
  "Chain_Count", "Total_Strength", "Peak_Strength", "Hit_Count",
  "x", "y", "x_from", "y_from", "x_to", "y_to", "label",
  "Ligand_Receptor", "Anchor_Score", "Order_Index", "Target_Total_Strength",
  "Sector_ID", "Label", "label_x", "label_y", "label_angle", "hjust",
  "Role_Level", "radius", "is_leaf", "index", "node_count"
)))

# ---  (Configuration Paths) ---
INPUT_GENE_FILE <- "./"

# NicheNet
NICHENET_DB_DIR <- "./"

OUTPUT_BASE_DIR <- "./"

CURRENT_SCRIPT_PATH <- "./"

#  (Create output directory with timestamp)
TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTPUT_DIR <- file.path(OUTPUT_BASE_DIR, paste0("NicheNet_Plasma_Proteins_Result_", TIMESTAMP))
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)

#  (Create logs directory)
LOG_DIR <- file.path(OUTPUT_DIR, "logs")
if (!dir.exists(LOG_DIR)) dir.create(LOG_DIR, recursive = TRUE)
LOG_FILE <- file.path(LOG_DIR, paste0("run_nichenet_", TIMESTAMP, ".log"))

#  (Create figures directory)
FIGURE_DIR <- file.path(OUTPUT_DIR, "figures")
if (!dir.exists(FIGURE_DIR)) dir.create(FIGURE_DIR, recursive = TRUE)

#  (Define logging function)
log_msg <- function(msg) {
  timestamp_log <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  message(paste(timestamp_log, msg))
  cat(paste(timestamp_log, msg, "\n"), file = LOG_FILE, append = TRUE)
}

#  (Backup the script)
SCRIPT_BACKUP_DIR <- file.path(OUTPUT_DIR, "scripts")
if (!dir.exists(SCRIPT_BACKUP_DIR)) dir.create(SCRIPT_BACKUP_DIR, recursive = TRUE)
if (file.exists(CURRENT_SCRIPT_PATH)) {
  file.copy(CURRENT_SCRIPT_PATH, file.path(SCRIPT_BACKUP_DIR, basename(CURRENT_SCRIPT_PATH)))
  log_msg(" scripts  (Script backed up to output directory).")
}

log_msg("=== NicheNet (Starting NicheNet Internal Network Analysis) ===")
log_msg(paste(" (Output Directory):", OUTPUT_DIR))

# 2.  (Load Input Data) --------------------------------------------
log_msg("... (Loading input gene list...)")
if (!file.exists(INPUT_GENE_FILE)) {
  stop(paste(" (Input file not found):", INPUT_GENE_FILE))
}

gene_data <- read_csv(INPUT_GENE_FILE, show_col_types = FALSE)
plasma_proteins <- unique(gene_data$gene)
log_msg(paste("", length(plasma_proteins), " (Loaded", length(plasma_proteins), "plasma protein genes)."))

# 3. NicheNet (Load NicheNet Databases) ------------------------------
log_msg("NicheNet... (Loading NicheNet prior knowledge networks...)")
lr_network_file <- file.path(NICHENET_DB_DIR, "lr_network_human_21122021.rds")
ligand_target_file <- file.path(NICHENET_DB_DIR, "ligand_target_matrix_nsga2r_final.rds")

if (!file.exists(lr_network_file) || !file.exists(ligand_target_file)) {
  stop("NicheNet, .  (NicheNet database files not found.)")
}

lr_network <- readRDS(lr_network_file)
ligand_target_matrix <- readRDS(ligand_target_file)
log_msg("NicheNet.  (NicheNet databases loaded successfully.)")


# 4.  (Internal Network Check for Plasma Proteins) ------------------
log_msg("--- :  (Part: Internal Network Extraction) ---")

count_input_genes <- length(plasma_proteins)
count_internal_lr <- 0
count_internal_lrt <- 0
count_internal_lr_any_t <- 0

# 4.1 -  (Ligand and Receptor both in plasma proteins list - )
log_msg("-... (Extracting internal Ligand-Receptor interactions...)")

# NicheNet
internal_LR_raw <- lr_network %>%
  filter(from %in% plasma_proteins & to %in% plasma_proteins) %>%
  rename(Ligand = from, Receptor = to) %>%
  distinct(Ligand, Receptor)

log_msg(paste("NicheNet", nrow(internal_LR_raw), "-. "))

log_msg("OmniPath API, /... (Downloading OmniPath data...)")
omnipath_url <- "https://omnipathdb.org/interactions?fields=genesymbols"
omnipath_data <- tryCatch({
  read_tsv(omnipath_url, show_col_types = FALSE)
}, error = function(e) {
  log_msg(paste(": OmniPath (Error downloading OmniPath):", e$message))
  return(NULL)
})

if (!is.null(omnipath_data)) {
  # OmniPath
  omnipath_lr <- omnipath_data %>%
    filter(!is.na(source_genesymbol) & !is.na(target_genesymbol)) %>%
    select(Ligand = source_genesymbol, Receptor = target_genesymbol, 
           is_stimulation, is_inhibition) %>%
    group_by(Ligand, Receptor) %>%
    summarise(
      is_stimulation = any(is_stimulation == TRUE, na.rm = TRUE),
      is_inhibition = any(is_inhibition == TRUE, na.rm = TRUE),
      .groups = "drop"
    )
  
  # OmniPath, stimulationinhibitionTRUE
  internal_LR <- internal_LR_raw %>%
    inner_join(omnipath_lr, by = c("Ligand", "Receptor")) %>%
    filter(is_stimulation == TRUE | is_inhibition == TRUE) %>%
    arrange(Ligand, Receptor)
  
  log_msg("OmniPath: is_stimulationis_inhibitionTRUE. ")
} else {
  log_msg(": OmniPath, NicheNet. ")
  internal_LR <- internal_LR_raw %>% arrange(Ligand, Receptor)
}

if (nrow(internal_LR) > 0) {
  count_internal_lr <- nrow(internal_LR)
  write_csv(internal_LR, file.path(OUTPUT_DIR, "1_Internal_Ligand_Receptor.csv"))
  log_msg(paste("", count_internal_lr, "- (Found", count_internal_lr, "internal L-R interactions)."))
} else {
  log_msg("- (No internal L-R interactions found).")
}

# 4.2 --  (Ligand, Receptor, and Target ALL in plasma proteins list)
# 4.3 -  (Internal L-R predicting ANY target genes)
log_msg("--... (Extracting internal L-R-T and L-R-Any-T full chains...)")
internal_LRT_list <- list()
internal_LR_Any_T_list <- list()

valid_ligands_matrix <- intersect(plasma_proteins, colnames(ligand_target_matrix))

for (lig in valid_ligands_matrix) {
  # （, OmniPathinternal_LR）
  recs_in_P <- internal_LR %>% filter(Ligand == lig) %>% pull(Receptor) %>% unique()
  if (length(recs_in_P) == 0) next
  
  scores <- ligand_target_matrix[, lig]
  top_targets <- scores[scores > 0.01]
  
  if (length(top_targets) > 0) {
    # : 、, **>0.01 (L-R -> Any Target)
    grid_any_t <- expand.grid(Ligand = lig, Receptor = recs_in_P, Target_Gene = names(top_targets), stringsAsFactors = FALSE)
    grid_any_t$Regulatory_Potential <- as.numeric(top_targets[grid_any_t$Target_Gene])
    internal_LR_Any_T_list[[length(internal_LR_Any_T_list) + 1]] <- grid_any_t
    
    #  (L-R -> Target in list)
    targets_in_P <- intersect(names(top_targets), plasma_proteins)
    
    if (length(targets_in_P) > 0) {
      grid_lrt <- expand.grid(Ligand = lig, Receptor = recs_in_P, Target_Gene = targets_in_P, stringsAsFactors = FALSE)
      grid_lrt$Regulatory_Potential <- as.numeric(top_targets[grid_lrt$Target_Gene])
      internal_LRT_list[[length(internal_LRT_list) + 1]] <- grid_lrt
    }
  }
}

#  4.2  (L-R-T in list)
internal_LRT <- bind_rows(internal_LRT_list)
if (nrow(internal_LRT) > 0) {
  internal_LRT <- internal_LRT %>% 
    distinct() %>%
    arrange(desc(Regulatory_Potential))
  
  count_internal_lrt <- nrow(internal_LRT)
  write_csv(internal_LRT, file.path(OUTPUT_DIR, "2_Internal_Ligand_Receptor_Target.csv"))
  log_msg(paste("", count_internal_lrt, "-- (Found", count_internal_lrt, "internal L-R-T full chains)."))
} else {
  log_msg("-- (No internal L-R-T full chains found).")
}

#  4.3  (L-R in list -> Any Target)
internal_LR_Any_T <- bind_rows(internal_LR_Any_T_list)
if (nrow(internal_LR_Any_T) > 0) {
  internal_LR_Any_T <- internal_LR_Any_T %>% 
    distinct() %>%
    arrange(desc(Regulatory_Potential))
  
  count_internal_lr_any_t <- nrow(internal_LR_Any_T)
  write_csv(internal_LR_Any_T, file.path(OUTPUT_DIR, "3_Internal_LR_Any_Target.csv"))
  log_msg(paste("", count_internal_lr_any_t, "- (Found", count_internal_lr_any_t, "internal L-R predicting Any Target chains)."))
} else {
  log_msg("- (No internal L-R predicting Any Target chains found).")
}

# 4.4  (Visualization Output) ----------------------------------------
log_msg("... (Generating concentric network plot and checking network coverage...)")

prepare_visualization_components <- function(network_df) {
  network_df <- network_df %>%
    distinct(Ligand, Receptor, Target_Gene, .keep_all = TRUE) %>%
    arrange(desc(Regulatory_Potential), Ligand, Receptor, Target_Gene)

  backbone_df <- network_df

  if (nrow(backbone_df) == 0) {
    return(NULL)
  }

  ligand_nodes <- backbone_df %>%
    group_by(Ligand) %>%
    summarise(
      Chain_Count = n(),
      Total_Strength = sum(Regulatory_Potential),
      Peak_Strength = max(Regulatory_Potential),
      .groups = "drop"
    ) %>%
    arrange(desc(Total_Strength), desc(Chain_Count), desc(Peak_Strength), Ligand) %>%
    mutate(
      name = Ligand,
      Role = "Ligand",
      Label = Ligand,
      Sector_ID = paste0("Ligand::", Ligand),
      Order_Index = row_number(),
      Role_Level = 1
    )
  
  receptor_nodes <- backbone_df %>%
    group_by(Receptor) %>%
    summarise(
      Chain_Count = n(),
      Total_Strength = sum(Regulatory_Potential),
      Peak_Strength = max(Regulatory_Potential),
      .groups = "drop"
    ) %>%
    arrange(desc(Total_Strength), desc(Chain_Count), desc(Peak_Strength), Receptor) %>%
    mutate(
      name = Receptor,
      Role = "Receptor",
      Label = Receptor,
      Sector_ID = paste0("Receptor::", Receptor),
      Order_Index = row_number(),
      Role_Level = 2
    ) %>%
    select(-Receptor)

  target_nodes <- backbone_df %>%
    group_by(Target_Gene) %>%
    summarise(
      Chain_Count = n(),
      Total_Strength = sum(Regulatory_Potential),
      Peak_Strength = max(Regulatory_Potential),
      .groups = "drop"
    ) %>%
    arrange(desc(Total_Strength), desc(Chain_Count), desc(Peak_Strength), Target_Gene) %>%
    mutate(
      name = Target_Gene,
      Role = "Target",
      Label = Target_Gene,
      Sector_ID = paste0("Target::", Target_Gene),
      Order_Index = row_number(),
      Role_Level = 3
    )

  node_df <- bind_rows(ligand_nodes, receptor_nodes, target_nodes) %>%
    mutate(Role = factor(Role, levels = c("Ligand", "Receptor", "Target"))) %>%
    arrange(Role_Level, Order_Index)

  edge_lr <- backbone_df %>%
    group_by(Ligand, Receptor) %>%
    summarise(
      Edge_Weight = sum(Regulatory_Potential),
      Chain_Count = n(),
      .groups = "drop"
    ) %>%
    mutate(
      Edge_Type = "Ligand-Receptor",
      from = paste0("Ligand::", Ligand),
      to = paste0("Receptor::", Receptor)
    )
  
  edge_rt <- backbone_df %>%
    group_by(Receptor, Target_Gene) %>%
    summarise(
      Edge_Weight = sum(Regulatory_Potential),
      Chain_Count = n(),
      .groups = "drop"
    ) %>%
    mutate(
      Edge_Type = "Receptor-Target",
      from = paste0("Receptor::", Receptor),
      to = paste0("Target::", Target_Gene)
    )

  edge_df <- bind_rows(edge_lr, edge_rt) %>%
    arrange(desc(Edge_Weight), Edge_Type, from, to)

  list(
    backbone_df = backbone_df,
    node_df = node_df,
    edge_df = edge_df
  )
}

build_concentric_layout <- function(node_df) {
  radius_map <- c("Ligand" = 1.25, "Receptor" = 2.35, "Target" = 3.45)
  
  node_df %>%
    arrange(Role_Level, desc(Total_Strength), desc(Chain_Count), name) %>%
    group_by(Role) %>%
    mutate(
      node_count = n(),
      Relative_Strength = if(max(Total_Strength) > min(Total_Strength)) {
        (Total_Strength - min(Total_Strength)) / (max(Total_Strength) - min(Total_Strength))
      } else {
        0.5
      },
      angle = {
        role_n <- dplyr::n()
        if (role_n == 1) {
          rep(pi / 2, role_n)
        } else {
          seq(pi / 2, pi / 2 - 2 * pi, length.out = role_n + 1)[seq_len(role_n)]
        }
      },
      radius = radius_map[as.character(Role)],
      x = radius * cos(angle),
      y = radius * sin(angle),
      label_x = (radius + 0.34) * cos(angle),
      label_y = (radius + 0.34) * sin(angle),
      normalized_angle = (angle * 180 / pi) %% 360,
      label_angle = ifelse(normalized_angle > 90 & normalized_angle < 270, normalized_angle + 180, normalized_angle),
      hjust = ifelse(cos(angle) >= -1e-5, 0, 1)
    ) %>%
    ungroup()
}

join_edge_coordinates <- function(edge_df, node_layout) {
  edge_df %>%
    left_join(
      node_layout %>%
        select(Sector_ID, x_from = x, y_from = y),
      by = c("from" = "Sector_ID")
    ) %>%
    left_join(
      node_layout %>%
        select(Sector_ID, x_to = x, y_to = y),
      by = c("to" = "Sector_ID")
    )
}

validate_concentric_coverage <- function(network_df, components) {
  if (is.null(components) || nrow(network_df) == 0) {
    return(NULL)
  }
  
  network_df_unique <- network_df %>%
    distinct(Ligand, Receptor, Target_Gene, .keep_all = TRUE)
  
  expected_ligands <- sort(unique(network_df_unique$Ligand))
  expected_receptors <- sort(unique(network_df_unique$Receptor))
  expected_targets <- sort(unique(network_df_unique$Target_Gene))
  
  actual_ligands <- components$node_df %>%
    filter(Role == "Ligand") %>%
    pull(name) %>%
    unique() %>%
    sort()
  
  actual_receptors <- components$node_df %>%
    filter(Role == "Receptor") %>%
    pull(name) %>%
    unique() %>%
    sort()
  
  actual_targets <- components$node_df %>%
    filter(Role == "Target") %>%
    pull(name) %>%
    unique() %>%
    sort()
  
  expected_lr_edges <- network_df_unique %>%
    distinct(Ligand, Receptor) %>%
    transmute(edge_key = paste(Ligand, Receptor, sep = "||")) %>%
    pull(edge_key) %>%
    sort()
  
  actual_lr_edges <- components$edge_df %>%
    filter(Edge_Type == "Ligand-Receptor") %>%
    transmute(
      edge_key = paste(
        sub("^Ligand::", "", from),
        sub("^Receptor::", "", to),
        sep = "||"
      )
    ) %>%
    pull(edge_key) %>%
    unique() %>%
    sort()
  
  expected_rt_edges <- network_df_unique %>%
    distinct(Receptor, Target_Gene) %>%
    transmute(edge_key = paste(Receptor, Target_Gene, sep = "||")) %>%
    pull(edge_key) %>%
    sort()
  
  actual_rt_edges <- components$edge_df %>%
    filter(Edge_Type == "Receptor-Target") %>%
    transmute(
      edge_key = paste(
        sub("^Receptor::", "", from),
        sub("^Target::", "", to),
        sep = "||"
      )
    ) %>%
    pull(edge_key) %>%
    unique() %>%
    sort()
  
  missing_ligands <- setdiff(expected_ligands, actual_ligands)
  missing_receptors <- setdiff(expected_receptors, actual_receptors)
  missing_targets <- setdiff(expected_targets, actual_targets)
  missing_lr_edges <- setdiff(expected_lr_edges, actual_lr_edges)
  missing_rt_edges <- setdiff(expected_rt_edges, actual_rt_edges)
  
  list(
    expected_chain_count = nrow(network_df_unique),
    expected_ligand_count = length(expected_ligands),
    expected_receptor_count = length(expected_receptors),
    expected_target_count = length(expected_targets),
    expected_lr_edge_count = length(expected_lr_edges),
    expected_rt_edge_count = length(expected_rt_edges),
    displayed_ligand_count = length(actual_ligands),
    displayed_receptor_count = length(actual_receptors),
    displayed_target_count = length(actual_targets),
    displayed_lr_edge_count = length(actual_lr_edges),
    displayed_rt_edge_count = length(actual_rt_edges),
    missing_ligands = missing_ligands,
    missing_receptors = missing_receptors,
    missing_targets = missing_targets,
    missing_lr_edges = missing_lr_edges,
    missing_rt_edges = missing_rt_edges,
    node_edge_coverage_complete = all(
      length(missing_ligands) == 0,
      length(missing_receptors) == 0,
      length(missing_targets) == 0,
      length(missing_lr_edges) == 0,
      length(missing_rt_edges) == 0
    ),
    row_level_triplets_drawn_individually = FALSE
  )
}

plot_concentric_network <- function(components, output_dir) {
  if (is.null(components) || nrow(components$edge_df) == 0) {
    return(invisible(NULL))
  }
  
  role_palette <- c(
    "Ligand" = "#8E2C2C",
    "Receptor" = "#2F5D8A",
    "Target" = "#2C7A5B"
  )
  edge_palette <- c(
    "Ligand-Receptor" = "#C7CCD6",
    "Receptor-Target" = "#355C7D"
  )
  
  node_layout <- build_concentric_layout(components$node_df)
  edge_layout <- join_edge_coordinates(components$edge_df, node_layout)
  
  ring_path <- function(radius_value, role_name) {
    tibble(
      angle = seq(0, 2 * pi, length.out = 360),
      x = radius_value * cos(angle),
      y = radius_value * sin(angle),
      Role = role_name
    )
  }
  
  ring_df <- bind_rows(
    ring_path(1.25, "Ligand"),
    ring_path(2.35, "Receptor"),
    ring_path(3.45, "Target")
  )
  
  concentric_plot <- ggplot() +
    geom_path(
      data = ring_df,
      aes(x = x, y = y, group = Role),
      colour = "#E5E7EB",
      linewidth = 0.4
    ) +
    geom_curve(
      data = edge_layout %>% filter(Edge_Type == "Ligand-Receptor"),
      aes(
        x = x_from,
        y = y_from,
        xend = x_to,
        yend = y_to,
        linewidth = Edge_Weight,
        alpha = Edge_Weight
      ),
      curvature = 0.22,
      lineend = "round",
      colour = edge_palette["Ligand-Receptor"]
    ) +
    geom_curve(
      data = edge_layout %>% filter(Edge_Type == "Receptor-Target"),
      aes(
        x = x_from,
        y = y_from,
        xend = x_to,
        yend = y_to,
        linewidth = Edge_Weight,
        alpha = Edge_Weight
      ),
      curvature = 0.14,
      lineend = "round",
      colour = edge_palette["Receptor-Target"]
    ) +
    geom_point(
      data = node_layout,
      aes(x = x, y = y, size = Relative_Strength, fill = Role),
      shape = 21,
      colour = "#FFFFFF",
      stroke = 0.9
    ) +
    geom_text(
      data = node_layout,
      aes(
        x = label_x,
        y = label_y,
        label = Label,
        angle = label_angle,
        hjust = hjust,
        colour = Role
      ),
      family = "sans",
      size = 3.05,
      show.legend = FALSE
    ) +
    scale_fill_manual(values = role_palette) +
    scale_colour_manual(values = role_palette, guide = guide_legend(override.aes = list(size = 4))) +
    scale_size_continuous(range = c(2.5, 12.5), guide = guide_legend(override.aes = list(colour = "#666666"))) +
    scale_linewidth_continuous(range = c(0.38, 2.2)) +
    scale_alpha_continuous(range = c(0.18, 0.82)) +
    coord_equal(clip = "off") +
    labs(
      fill = "Node role",
      colour = "Node role",
      size = "Relative\nstrength",
      linewidth = "Edge\nstrength",
      alpha = "Edge\nstrength"
    ) +
    theme_void(base_family = "sans") +
    theme(
      legend.position = "right",
      legend.title = element_text(size = 11, face = "bold", colour = "#111111"),
      legend.text = element_text(size = 10, colour = "#1A1A1A"),
      legend.spacing.y = unit(0.3, "cm"),
      plot.margin = margin(10, 10, 10, 10),
      text = element_text(family = "sans")
    )
  
  concentric_plot_bottom <- concentric_plot + theme(legend.position = "bottom", legend.box = "horizontal")
  
  ggplot2::ggsave(
    filename = file.path(output_dir, "5_Internal_LRT_Network_Concentric.pdf"),
    plot = concentric_plot,
    width = 12.8,
    height = 11.6,
    dpi = 600,
    bg = "white"
  )
  ggplot2::ggsave(
    filename = file.path(output_dir, "5_Internal_LRT_Network_Concentric.png"),
    plot = concentric_plot,
    width = 12.8,
    height = 11.6,
    dpi = 600,
    bg = "white"
  )
  ggplot2::ggsave(
    filename = file.path(output_dir, "5_Internal_LRT_Network_Concentric_bottom_legend.pdf"),
    plot = concentric_plot_bottom,
    width = 12.8,
    height = 12.6,
    dpi = 600,
    bg = "white"
  )
  ggplot2::ggsave(
    filename = file.path(output_dir, "5_Internal_LRT_Network_Concentric_bottom_legend.png"),
    plot = concentric_plot_bottom,
    width = 12.8,
    height = 12.6,
    dpi = 600,
    bg = "white"
  )
}

scale_rows <- function(input_matrix) {
  scaled_matrix <- t(scale(t(input_matrix)))
  scaled_matrix[is.na(scaled_matrix)] <- 0
  scaled_matrix
}


visualization_source <- internal_LRT
if (nrow(visualization_source) == 0 && nrow(internal_LR_Any_T) > 0) {
  visualization_source <- internal_LR_Any_T
}

network_coverage_check <- NULL

if (nrow(visualization_source) > 0) {
  visualization_components <- prepare_visualization_components(visualization_source)
  plot_concentric_network(visualization_components, FIGURE_DIR)
  if (nrow(internal_LRT) > 0) {
    network_coverage_check <- validate_concentric_coverage(internal_LRT, visualization_components)
    if (isTRUE(network_coverage_check$node_edge_coverage_complete)) {
      log_msg(" 2_Internal_Ligand_Receptor_Target.csv  (Concentric plot covers all nodes and edges from 2_Internal_Ligand_Receptor_Target.csv).")
    } else {
      log_msg(":  2_Internal_Ligand_Receptor_Target.csv  (Warning: coverage is incomplete).")
    }
  }
  log_msg(" (Concentric network plot generated).")
} else {
  log_msg(",  (No figures generated due to unavailable chains).")
}


# 5. CSV (Generate Summary CSV) ------------------
log_msg("... (Generating detailed summary statistics CSV...)")

summary_stats <- data.frame(
  Metric_En = c(
    "Total input plasma proteins",
    "Internal Ligand-Receptor pairs found (Prior knowledge)",
    "Internal Ligand-Receptor-Target chains found (All in list)",
    "Internal L-R predicting ANY Target chains found"
  ),
  Metric_Cn = c(
    "",
    "- ",
    "(--) ",
    "- "
  ),
  Value = c(
    count_input_genes,
    count_internal_lr,
    count_internal_lrt,
    count_internal_lr_any_t
  ),
  stringsAsFactors = FALSE
)

# CSV
write_csv(summary_stats, file.path(OUTPUT_DIR, "4_Analysis_Summary_Statistics.csv"))
log_msg("CSV (Summary statistics CSV generated).")

log_msg("===  (Analysis Completed) ===")

#  (Record execution time)
log_msg(paste(" (Completion Time):", Sys.time()))
