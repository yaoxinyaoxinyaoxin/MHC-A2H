# ==============================================================================
# [Script]: 5.3.1.1_run_nichenet_analysis_template.R
# [Method]: Intercellular Communication Prediction
# [Step]: 5.3.1.1_run_nichenet_analysis
#
# [Function]:
# Performs a comprehensive NicheNet analysis based on MR/Colocalization results (Upstream & Downstream).
#
# [Usage]: 
#   Rscript 5.3.1.1_run_nichenet_analysis_template.R \
#     --input_mr <path> \
#     --cell_map <path> \
#     --bg_genes <path> \
#     --nichenet_db <path> \
#     --out_dir <path>
# ==============================================================================

# 1. Setup Environment & Configuration --------------------------------------
suppressPackageStartupMessages({
  library(optparse)
  library(nichenetr)
  library(tidyverse)
  library(dplyr)
  library(ggplot2)
  library(circlize)
  library(RColorBrewer)
  library(cowplot)
  library(ggpubr)
  library(readxl) # For reading the mapping file
})

# Clear environment variables
rm(list = ls())
gc()

# Define command line arguments
option_list <- list(
  make_option(c("--input_mr"), type="character", default=NULL, help="MR Results CSV file"),
  make_option(c("--cell_map"), type="character", default=NULL, help="Cell Type Mapping Excel File"),
  make_option(c("--bg_genes"), type="character", default=NULL, help="Background Genes Expression CSV"),
  make_option(c("--nichenet_db"), type="character", default=NULL, help="NicheNet DB Directory"),
  make_option(c("--out_dir"), type="character", default="./NicheNet_Analysis", help="Output directory path")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_mr) || is.null(opt$cell_map) || is.null(opt$bg_genes) || is.null(opt$nichenet_db)) {
  print_help(opt_parser)
  stop("Missing required arguments. / . ", call.=FALSE)
}

# --- Configuration Paths ---
INPUT_FILE <- opt$input_mr
CELL_MAPPING_FILE <- opt$cell_map
BACKGROUND_GENES_FILE <- opt$bg_genes
NICHENET_DB_DIR <- opt$nichenet_db
OUTPUT_BASE_DIR <- opt$out_dir

# Create Output Directory with Timestamp
TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTPUT_DIR <- file.path(OUTPUT_BASE_DIR, paste0("NicheNet_Analysis_", TIMESTAMP))
if (!dir.exists(OUTPUT_DIR)) dir.create(OUTPUT_DIR, recursive = TRUE)


# Logging
LOG_FILE <- file.path(OUTPUT_DIR, paste0("nichenet_analysis_", TIMESTAMP, ".log"))

log_msg <- function(msg) {
  timestamp_log <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  message(paste(timestamp_log, msg))
  cat(paste(timestamp_log, msg, "\n"), file = LOG_FILE, append = TRUE)
}

log_msg("Starting NicheNet Analysis Pipeline (Upstream & Downstream)")
log_msg(paste("Input File:", INPUT_FILE))
log_msg(paste("Cell Mapping File:", CELL_MAPPING_FILE))
log_msg(paste("Output Directory:", OUTPUT_DIR))

# 2. Load NicheNet Networks -------------------------------------------------
log_msg("Loading NicheNet Prior Knowledge Networks...")

lr_network <- readRDS(file.path(NICHENET_DB_DIR, "lr_network_human_21122021.rds"))
ligand_target_matrix <- readRDS(file.path(NICHENET_DB_DIR, "ligand_target_matrix_nsga2r_final.rds"))
weighted_networks <- readRDS(file.path(NICHENET_DB_DIR, "weighted_networks_nsga2r_final.rds"))

log_msg("NicheNet Networks loaded successfully.")

# 3. Load and Process Input Data & Cell Mapping -----------------------------

# Load Cell Mapping
log_msg("Loading Cell Type Mapping Table...")
if (!file.exists(CELL_MAPPING_FILE)) {
  stop(paste("Cell mapping file not found:", CELL_MAPPING_FILE))
}

cell_map_df <- read_excel(CELL_MAPPING_FILE, skip = 1) %>%
  dplyr::select(, ) %>%
  rename(Abbr = , FullName = ) %>%
  mutate(Abbr = gsub("[^[:alnum:]]", "", Abbr)) %>% # Remove invisible characters/spaces
  mutate(Abbr_Lower = tolower(Abbr)) %>% # Convert to lowercase for matching
  na.omit()

# Create mapping vector: Lowercase Abbr -> Full Name
cell_map_vec <- setNames(cell_map_df$FullName, cell_map_df$Abbr_Lower)

# Manual fixes for mismatches between MR data abbreviations and Excel abbreviations
# MR uses: nk, nkr, plasma, dc
# Excel (after processing) has: nkcells, nkrecruitingcells, plasmacells, dcs
cell_map_vec["nk"] <- cell_map_vec["nkcells"]
cell_map_vec["nkr"] <- cell_map_vec["nkrecruitingcells"]
cell_map_vec["plasma"] <- cell_map_vec["plasmacells"]
cell_map_vec["dc"] <- cell_map_vec["dcs"]

log_msg(paste("Loaded cell mapping for", length(cell_map_vec), "cell types."))
log_msg("Mapping Preview:")
print(head(cell_map_vec))

# Helper to apply mapping
apply_cell_mapping <- function(cell_names) {
  # Convert input to lowercase to ensure matching
  cell_names_lower <- tolower(cell_names)
  mapped <- cell_map_vec[cell_names_lower]
  
  # If no match found, keep original name (but warn if needed)
  final_names <- ifelse(is.na(mapped), cell_names, mapped)
  return(final_names)
}

# Load MR Results
log_msg("Loading MR/Colocalization Results...")
if (!file.exists(INPUT_FILE)) stop(paste("Input file not found:", INPUT_FILE))
mr_data <- read.csv(INPUT_FILE)
log_msg(paste("Loaded", nrow(mr_data), "rows from MR results."))

# Apply Mapping to MR Data (Receiver Cells)
mr_data$cell <- apply_cell_mapping(mr_data$cell)

# Check if CD4ET (now mapped name) is present
cd4et_name <- cell_map_vec["cd4et"]
if (!is.na(cd4et_name) && cd4et_name %in% mr_data$cell) {
  log_msg(paste("Confirmed: Effector CD4 T cells (", cd4et_name, ") present in MR data."))
} else {
  log_msg("Warning: Effector CD4 T cells not found in mapped MR data. Checking raw data...")
  if ("cd4et" %in% tolower(mr_data$cell) || "CD4ET" %in% mr_data$cell) {
     log_msg("  Found 'cd4et' in raw data but mapping failed. Check mapping table.")
  }
}


# Load Background Data
log_msg("Loading Background Genes (Expression Data)...")
if (!file.exists(BACKGROUND_GENES_FILE)) stop(paste("Background file not found:", BACKGROUND_GENES_FILE))
bg_data <- read.csv(BACKGROUND_GENES_FILE)

# Apply Mapping to Background Data (Sender Cells)
bg_data$cell <- apply_cell_mapping(bg_data$cell)

# Map Ensembl IDs to Symbols if needed
if (!"gene_symbol" %in% colnames(bg_data)) {
  log_msg("Mapping Ensembl IDs to Gene Symbols for background genes...")
  map_genes <- function(ens_ids) {
    if (requireNamespace("EnsDb.Hsapiens.v75", quietly = TRUE) && requireNamespace("AnnotationDbi", quietly = TRUE)) {
      edb <- EnsDb.Hsapiens.v75::EnsDb.Hsapiens.v75
      AnnotationDbi::mapIds(edb, keys = ens_ids, column = "SYMBOL", keytype = "GENEID")
    } else if (requireNamespace("org.Hs.eg.db", quietly = TRUE) && requireNamespace("AnnotationDbi", quietly = TRUE)) {
      AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = ens_ids, column = "SYMBOL", keytype = "ENSEMBL")
    } else {
      nichenetr::convert_human_gene_ids(ens_ids)
    }
  }
  
  unique_ensg <- unique(bg_data$gene_id)
  gene_map <- map_genes(unique_ensg)
  gene_map_df <- data.frame(gene_id = names(gene_map), gene_symbol = as.character(gene_map), stringsAsFactors = FALSE)
  
  bg_data <- bg_data %>%
    left_join(gene_map_df, by = "gene_id") %>%
    filter(!is.na(gene_symbol) & gene_symbol != "")
}

log_msg(paste("Background genes ready. Total unique expressed symbols:", length(unique(bg_data$gene_symbol))))

# Define Universe of Expressed Genes per Cell
expressed_genes_per_cell <- split(bg_data$gene_symbol, bg_data$cell)
all_expressed_genes <- unique(bg_data$gene_symbol)

# Define Potential Ligands (Global)
ligands <- lr_network %>% pull(from) %>% unique()
expressed_ligands_global <- intersect(ligands, all_expressed_genes)
log_msg(paste("Number of potential ligands (expressed in dataset):", length(expressed_ligands_global)))


# =============================================================================
# PART A: Upstream Analysis (Input Gene = Target)
# Question: Which ligands regulate my MR genes?
# =============================================================================
log_msg("--- Starting PART A: Upstream Analysis (Input as Target) ---")

# Extract Targets: Gene Name (Symbol) and Receiver Cell
mr_targets <- mr_data %>%
  dplyr::select(gene_name, cell) %>%
  distinct() %>%
  filter(!is.na(gene_name) & gene_name != "") %>%
  rename(target_gene = gene_name, receiver_cell = cell)

upstream_results_list <- list()
upstream_ligand_activities <- list()

receiver_cells <- unique(mr_targets$receiver_cell)

for (receiver_cell in receiver_cells) {
  log_msg(paste("  [Upstream] Analyzing Receiver Cell:", receiver_cell))
  
  # 1. Define Genes of Interest (Targets)
  genes_of_interest <- mr_targets %>%
    filter(receiver_cell == !!receiver_cell) %>%
    pull(target_gene) %>%
    unique()
  
  # Filter targets present in Ligand-Target Matrix
  genes_of_interest <- genes_of_interest[genes_of_interest %in% rownames(ligand_target_matrix)]
  
  if (length(genes_of_interest) == 0) {
    log_msg("    Skipping: No valid genes of interest found in NicheNet DB.")
    
    # DIAGNOSTIC for specific cells (like CD4ET)
    # Check if this cell is related to CD4ET (using mapped name)
    cd4et_mapped <- cell_map_vec["cd4et"]
    if (!is.na(cd4et_mapped) && receiver_cell == cd4et_mapped) {
       raw_genes <- mr_targets %>% 
         filter(receiver_cell == !!receiver_cell) %>% 
         pull(target_gene) %>% 
         unique()
       log_msg(paste("    [DEBUG] Raw genes for", receiver_cell, ":", paste(raw_genes, collapse=", ")))
      log_msg(paste("    [DEBUG] None of these genes are in NicheNet ligand_target_matrix rownames (targets)."))
    }
    next
  }
  
  # 2. Define Background Genes for Receiver
  if (!receiver_cell %in% names(expressed_genes_per_cell)) {
    log_msg(paste("    Warning: Receiver cell", receiver_cell, "not found in background data. Using all expressed genes as background."))
    expressed_genes_receiver <- all_expressed_genes
  } else {
    expressed_genes_receiver <- expressed_genes_per_cell[[receiver_cell]]
  }
  
  # Ensure background covers targets
  background_genes <- unique(c(expressed_genes_receiver, genes_of_interest))
  background_genes <- background_genes[background_genes %in% rownames(ligand_target_matrix)]
  
  # 3. Predict Ligand Activities
  ligand_activities <- predict_ligand_activities(
    geneset = genes_of_interest,
    background_expressed_genes = background_genes,
    ligand_target_matrix = ligand_target_matrix,
    potential_ligands = expressed_ligands_global
  )
  
  if (is.null(ligand_activities) || nrow(ligand_activities) == 0) {
    log_msg("    No active ligands found.")
    next
  }
  
  # Store Activities
  ligand_activities <- as.data.frame(ligand_activities)
  ligand_activities$receiver_cell <- receiver_cell
  upstream_ligand_activities[[receiver_cell]] <- ligand_activities
  
  # 4. Infer Full Chain: Sender -> Ligand -> Receptor -> Target
  
  # Select Top Ligands (e.g., top 30 or Pearson > 0.05)
  top_ligands <- ligand_activities %>% 
    arrange(desc(pearson)) %>% 
    top_n(30, pearson) %>% 
    pull(test_ligand)
  
  if (length(top_ligands) == 0) next
  
  # A. Ligand -> Target Links
  active_ligand_target_links <- genes_of_interest %>%
    lapply(function(target) {
      scores <- ligand_target_matrix[target, top_ligands]
      data.frame(ligand = names(scores), target_gene = target, regulatory_potential = as.numeric(scores))
    }) %>%
    bind_rows() %>%
    filter(regulatory_potential > 0.05)
  
  if (nrow(active_ligand_target_links) == 0) next
  
  # B. Ligand -> Receptor Links
  receptors <- lr_network %>% pull(to) %>% unique()
  expressed_receptors <- intersect(receptors, expressed_genes_receiver)
  
  lr_links <- lr_network %>%
    filter(from %in% top_ligands & to %in% expressed_receptors) %>%
    dplyr::select(from, to) %>%
    rename(ligand = from, receptor = to) %>%
    distinct()
  
  if (nrow(lr_links) == 0) next
  
  # C. Sender -> Ligand Links
  sender_links_list <- list()
  for (ligand in top_ligands) {
    senders <- bg_data %>%
      filter(gene_symbol == ligand) %>%
      pull(cell) %>%
      unique()
    
    if (length(senders) > 0) {
      sender_links_list[[ligand]] <- data.frame(sender_cell = senders, ligand = ligand)
    }
  }
  sender_links <- bind_rows(sender_links_list)
  
  if (nrow(sender_links) == 0) next
  
  # D. Merge All
  full_chain <- sender_links %>%
    inner_join(lr_links, by = "ligand") %>%
    inner_join(active_ligand_target_links, by = "ligand")
  
  full_chain$receiver_cell <- receiver_cell
  
  full_chain <- full_chain %>%
    left_join(ligand_activities %>% dplyr::select(test_ligand, pearson), by = c("ligand" = "test_ligand")) %>%
    rename(ligand_activity_pearson = pearson)
  
  upstream_results_list[[receiver_cell]] <- full_chain
}

# Save Upstream Results
if (length(upstream_results_list) > 0) {
  final_upstream_results <- bind_rows(upstream_results_list)
  final_upstream_activities <- bind_rows(upstream_ligand_activities)
  
  write.csv(final_upstream_results, file.path(OUTPUT_DIR, "Upstream_Analysis_InputAsTarget_Chain.csv"), row.names = FALSE)
  write.csv(final_upstream_activities, file.path(OUTPUT_DIR, "Upstream_Analysis_InputAsTarget_Activities.csv"), row.names = FALSE)
  log_msg("Upstream Analysis results saved.")
} else {
  log_msg("No significant Upstream results found.")
}


# =============================================================================
# PART B: Downstream Analysis (Input Gene = Ligand)
# Question: If my MR gene is a ligand, what are its receptors and downstream targets?
# =============================================================================
log_msg("--- Starting PART B: Downstream Analysis (Input as Ligand) ---")

# Identify Input Genes that are known Ligands
mr_ligands <- mr_data %>%
  dplyr::select(gene_name, cell) %>%
  distinct() %>%
  rename(ligand = gene_name, sender_cell = cell) %>%
  filter(ligand %in% ligands) # Must be in NicheNet ligand database

log_msg(paste("  Found", nrow(mr_ligands), "potential ligand-sender pairs from input."))

downstream_results_list <- list()

if (nrow(mr_ligands) > 0) {
  
  # Get all possible receiver cells
  all_receiver_cells <- names(expressed_genes_per_cell)
  
  for (i in 1:nrow(mr_ligands)) {
    curr_ligand <- mr_ligands$ligand[i]
    curr_sender <- mr_ligands$sender_cell[i]
    
    # Check if ligand is actually expressed in the claimed sender cell in our background data
    # (Optional validation, but good to have. If not expressed, we still proceed as MR evidence is strong)
    is_expressed <- curr_ligand %in% expressed_genes_per_cell[[curr_sender]]
    if (!is_expressed) {
      # log_msg(paste("    Note:", curr_ligand, "not detected in expression data of", curr_sender, "- proceeding based on MR evidence."))
    }
    
    # Find Receptors for this Ligand
    potential_receptors <- lr_network %>%
      filter(from == curr_ligand) %>%
      pull(to) %>%
      unique()
    
    if (length(potential_receptors) == 0) next
    
    # Iterate over all possible Receiver Cells
    for (receiver_cell in all_receiver_cells) {
      
      # Find Receptors expressed in this Receiver
      expressed_genes_recv <- expressed_genes_per_cell[[receiver_cell]]
      active_receptors <- intersect(potential_receptors, expressed_genes_recv)
      
      if (length(active_receptors) == 0) next
      
      # Find Top Targets in this Receiver
      # We look for targets in the ligand_target_matrix for this ligand
      if (curr_ligand %in% colnames(ligand_target_matrix)) {
        scores <- ligand_target_matrix[, curr_ligand]
        
        # Filter targets expressed in receiver
        valid_targets <- intersect(names(scores), expressed_genes_recv)
        scores <- scores[valid_targets]
        
        # Get top 20 targets
        top_targets <- sort(scores, decreasing = TRUE)[1:20]
        top_targets <- top_targets[top_targets > 0.01] # Threshold
        
        if (length(top_targets) > 0) {
          
          # Create Result Rows
          # We expand grid of Receptor x Target
          
          res_df <- expand.grid(
            receptor = active_receptors,
            target_gene = names(top_targets),
            stringsAsFactors = FALSE
          )
          
          res_df$sender_cell <- curr_sender
          res_df$ligand <- curr_ligand
          res_df$receiver_cell <- receiver_cell
          res_df$target_regulatory_potential <- top_targets[res_df$target_gene]
          
          downstream_results_list[[paste(curr_sender, curr_ligand, receiver_cell, sep="_")]] <- res_df
        }
      } else {
        # Ligand not in matrix (no target predictions available), but Receptors exist
        # Output just Ligand-Receptor info
         res_df <- data.frame(
            sender_cell = curr_sender,
            ligand = curr_ligand,
            receptor = active_receptors,
            receiver_cell = receiver_cell,
            target_gene = NA,
            target_regulatory_potential = NA
         )
         downstream_results_list[[paste(curr_sender, curr_ligand, receiver_cell, sep="_")]] <- res_df
      }
    }
  }
}

# Save Downstream Results
if (length(downstream_results_list) > 0) {
  final_downstream_results <- bind_rows(downstream_results_list) %>%
    dplyr::select(sender_cell, ligand, receptor, receiver_cell, target_gene, target_regulatory_potential)
  
  write.csv(final_downstream_results, file.path(OUTPUT_DIR, "Downstream_Analysis_InputAsLigand_Chain.csv"), row.names = FALSE)
  log_msg("Downstream Analysis results saved.")
} else {
  log_msg("No significant Downstream results found (Input genes may not be known ligands).")
}

# 5. Save and Visualize (Consolidated) ---------------------------------------
log_msg("Generating Visualizations...")

# Define plot metadata if not present (was in original script but missing definition in snippet, 
# assuming CELL_META was global or needs definition for Circos. 
# Re-implementing simplified Circos plotting without external CELL_META dependency for safety)

# --- Visualizations for UPSTREAM Analysis ---
if (exists("final_upstream_results") && nrow(final_upstream_results) > 0) {
  
  # Plot 1: Sender - Receiver Connectivity (Chord Diagram)
  circos_df <- final_upstream_results %>%
    dplyr::select(sender_cell, receiver_cell, ligand, ligand_activity_pearson) %>%
    distinct() %>%
    group_by(sender_cell, receiver_cell) %>%
    summarise(weight = sum(ligand_activity_pearson), .groups = 'drop')
  
  if (nrow(circos_df) > 0) {
    circos_mat <- circos_df %>%
      pivot_wider(names_from = receiver_cell, values_from = weight, values_fill = 0) %>%
      column_to_rownames("sender_cell") %>%
      as.matrix()
    
    pdf(file.path(OUTPUT_DIR, "Upstream_Sender_Receiver_Circos.pdf"), width = 10, height = 10)
    circos.clear()
    chordDiagram(circos_mat, transparency = 0.5, annotationTrack = "grid", preAllocateTracks = 1)
    
    # Simple labels
    circos.track(track.index = 1, panel.fun = function(x, y) {
      xlim = get.cell.meta.data("xlim")
      ylim = get.cell.meta.data("ylim")
      sector.name = get.cell.meta.data("sector.index")
      
      # Use smaller font for long names, or wrap text if possible
      # Here we just adjust cex and use the full name (sector.name)
      circos.text(mean(xlim), ylim[1] + 0.1, sector.name, 
                  facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 0.6)
    }, bg.border = NA)
    
    title("Upstream: Predicted Cell-Cell Communication")
    dev.off()
  }
  
  # Plot 2: Ligand - Target Heatmap
  top_ligands_global <- final_upstream_activities %>% 
    group_by(receiver_cell) %>% 
    top_n(5, pearson) %>% 
    ungroup() %>% 
    pull(test_ligand) %>% 
    unique()
  
  plot_data_lt <- final_upstream_results %>%
    filter(ligand %in% top_ligands_global) %>%
    dplyr::select(ligand, target_gene, regulatory_potential, receiver_cell) %>%
    distinct()
  
  if (nrow(plot_data_lt) > 0) {
    p_lt <- ggplot(plot_data_lt, aes(x = target_gene, y = ligand, fill = regulatory_potential)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "red") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      facet_wrap(~receiver_cell, scales = "free") +
      labs(title = "Upstream: Ligand-Target Regulatory Potential", x = "Target Gene (MR Hit)", y = "Ligand")
    
    ggsave(file.path(OUTPUT_DIR, "Upstream_Ligand_Target_Heatmap.pdf"), p_lt, width = 15, height = 12)
  }
}

# --- Visualizations for DOWNSTREAM Analysis ---
if (exists("final_downstream_results") && nrow(final_downstream_results) > 0) {
  
  # Plot 1: Dotplot of Ligand (Sender) -> Receptor (Receiver)
  # Aggregating to count unique targets or just presence
  
  ds_plot_data <- final_downstream_results %>%
    dplyr::select(sender_cell, ligand, receptor, receiver_cell) %>%
    distinct()
  
  if (nrow(ds_plot_data) > 0) {
    ds_plot_data <- ds_plot_data %>%
      mutate(interaction = paste0(ligand, " -> ", receptor))
    
    p_ds <- ggplot(ds_plot_data, aes(x = receiver_cell, y = interaction, color = sender_cell)) +
      geom_point(size = 3) +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      labs(title = "Downstream: Ligand-Receptor Interactions", 
           subtitle = "Ligand (MR Hit) -> Receptor (Receiver)",
           x = "Receiver Cell", y = "Ligand -> Receptor Pair", color = "Sender Cell")
    
    ggsave(file.path(OUTPUT_DIR, "Downstream_Ligand_Receptor_Dotplot.pdf"), p_ds, width = 12, height = 10)
  }
  
  # Plot 2: Top Targets of Downstream Ligands (Heatmap)
  ds_target_data <- final_downstream_results %>%
    filter(!is.na(target_regulatory_potential)) %>%
    group_by(ligand, receiver_cell) %>%
    top_n(10, target_regulatory_potential) %>%
    ungroup()
  
  if (nrow(ds_target_data) > 0) {
     p_ds_target <- ggplot(ds_target_data, aes(x = target_gene, y = ligand, fill = target_regulatory_potential)) +
      geom_tile() +
      scale_fill_gradient(low = "white", high = "blue") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5)) +
      facet_wrap(~receiver_cell, scales = "free") +
      labs(title = "Downstream: Predicted Top Targets", x = "Predicted Target", y = "Ligand (MR Hit)")
     
     ggsave(file.path(OUTPUT_DIR, "Downstream_Predicted_Targets_Heatmap.pdf"), p_ds_target, width = 15, height = 12)
  }
}

log_msg("Analysis Completed.")
