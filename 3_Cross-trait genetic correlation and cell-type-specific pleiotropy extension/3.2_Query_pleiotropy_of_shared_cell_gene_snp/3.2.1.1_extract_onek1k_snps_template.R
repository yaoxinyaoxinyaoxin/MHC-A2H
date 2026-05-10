# ==============================================================================
# [Script]: extract_onek1k_snps_M2_ultra.R
# [Method]: cell-gene-snp
# [Step]: OneK1K sc-eQTL
# 
# [Function]:
# Execute the 'OneK1K sc-eQTL' step within the 'cell-gene-snp' analytical framework.
# 
# [Parameters / ]:
# Standard predefined thresholds.
# 
# [Steps / ]:
#   1. Data loading and initialization / 
#   2. Core analytical execution / 
#   3. Results formatting and output / 
# ==============================================================================

# -------------------------------------------------------------------------
# Script Name: extract_onek1k_snps.R
# Description: 
#   Extracts SNP data from OneK1K single-cell eQTL summary statistics based on 
#   a provided list of cell types and lead SNPs.
#
# Usage: 
#   Rscript extract_onek1k_snps.R
#
# Dependencies: 
#   data.table, dplyr, stringr, fs, readr, glue, ggplot2, ggalluvial
#
# Steps:
#   1. Initialize environment and paths. . 
#   2. Load input CSV containing target cell types and SNPs. SNPCSV. 
#   3. Clean and filter input data (remove bulk, summary rows). （bulk）. 
#   4. Iterate through each unique cell type. . 
#   5. Read the corresponding .tsv.gz file from source directory. .tsv.gz. 
#   6. Filter for target SNPs. SNP. 
#   7. Aggregate results. . 
#   8. Filter results by FDR <= 0.05 and save. FDR <= 0.05. 
#   9. Create unique best P-value version from FDR filtered data. FDRP. 
#   10. Generate visualizations (Sankey Diagram & Dot Plot). . 
#   11. Save outputs and generate report. . 
#
# Inputs:
#   - Forest_Plot_Hierarchical_Cell_Gene_SNP_Shared_aging_vs_Herpes_Zoster.csv: Contains 'Display_Name' (Cell Type) and 'Lead_SNP'.
#   - OneK1K source data: .tsv.gz files in specified directory.
#
# Outputs:
#   - Extracted SNP data (CSV). SNP. 
#   - Extracted SNP data filtered by FDR <= 0.05 (CSV). FDR <= 0.05SNP. 
#   - Unique Best P-value data (CSV). P. 
#   - Plots (PDF/PNG) showing Cell-SNP to Gene relationships. Cell-SNPGene. 
#   - Statistics.
# -------------------------------------------------------------------------

# 1. Load Libraries and Set Options ---------------------------------------
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(fs)
  library(readr)
  library(glue)
  library(AnnotationDbi)
  library(EnsDb.Hsapiens.v75)
})

# ==============================================================================
# 0. Command Line Arguments / 
# ==============================================================================
option_list <- list(
  make_option(c("--input_csv"), type="character", default=NULL,
              help="Path to the input CSV file containing target cell types and SNPs / SNPCSV"),
  make_option(c("--source_dir"), type="character", default=NULL,
              help="Directory containing OneK1K source .tsv.gz files / OneK1K.tsv.gz"),
  make_option(c("--out_dir"), type="character", default="./Extracted_OneK1K",
              help="Output directory path / ")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_csv) || is.null(opt$source_dir)) {
  print_help(opt_parser)
  stop("Input CSV and Source Directory must be provided. / CSV. ")
}

# Check and load visualization libraries
if (!require("ggplot2")) {
  message("Installing ggplot2...")
  install.packages("ggplot2", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  library(ggplot2)
}
if (!require("ggalluvial")) {
  message("Installing ggalluvial...")
  install.packages("ggalluvial", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
  library(ggalluvial)
}

# Enable multithreading for data.table
# data.table
setDTthreads(threads = 0) # Use all available cores

# 2. Define Paths and Constants -------------------------------------------

# Input File Path
input_csv_path <- opt$input_csv

# Source Data Directory
source_data_dir <- opt$source_dir

# Base Output Directory
base_output_dir <- opt$out_dir

# Create Timestamped Output Directory
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- file.path(base_output_dir, paste0("extracted_snps_", timestamp))
dir_create(output_dir)

# Create Logs Directory
logs_dir <- file.path(output_dir, "logs")
dir_create(logs_dir)

# Log File Path
log_file <- file.path(logs_dir, "extraction_log.txt")

# 3. Helper Functions -----------------------------------------------------

# Logging function
log_message <- function(message) {
  timestamped_msg <- paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", message)
  cat(timestamped_msg, "\n")
  write(timestamped_msg, file = log_file, append = TRUE)
}

# 4. Main Execution -------------------------------------------------------

log_message("Starting OneK1K SNP extraction script...")
log_message(paste("Input CSV:", input_csv_path))
log_message(paste("Source Data Directory:", source_data_dir))
log_message(paste("Output Directory:", output_dir))

# 4.1 Load and Process Input CSV
# CSV
log_message("Loading input CSV...")

if (!file_exists(input_csv_path)) {
  stop(paste("Input file not found:", input_csv_path))
}

dt_input <- fread(input_csv_path)

# --- Parse Hierarchical Gene Information ---
log_message("Parsing hierarchical Gene information...")
dt_input[, Parsed_Gene := NA_character_]
current_gene <- NA_character_
gene_vector <- character(nrow(dt_input))

for (i in 1:nrow(dt_input)) {
  display_name <- dt_input$Display_Name[i]
  lead_snp <- dt_input$Lead_SNP[i]
  
  # Check if this is a gene header row
  # Usually gene header has text like "GeneName (location)" and empty SNP info
  if (!is.na(display_name) && !grepl("^\\s", display_name) && (is.na(lead_snp) | lead_snp == "")) {
    # Extract Gene Name (everything before first space or parenthesis)
    current_gene <- str_extract(display_name, "^[^\\s(]+")
  }
  gene_vector[i] <- current_gene
}
dt_input[, Parsed_Gene := gene_vector]
log_message("Gene parsing complete.")
# -------------------------------------------

# Filter out bulk data and empty SNPs
# bulkSNP
log_message(paste("Total rows in input:", nrow(dt_input)))

# Logic to extract cell type:
# In the provided file format, 'Display_Name' usually contains cell type for rows with valid SNP data.
# We filter rows where Lead_SNP is not empty and is_bulk is FALSE (or NA).
# Also filter out summary rows (where is_summary is TRUE).
# : , SNP, 'Display_Name'. 
# Lead_SNPis_bulkFALSE（NA）. 
# （is_summaryTRUE）. 

dt_filtered <- dt_input[
  !is.na(Lead_SNP) & Lead_SNP != "" & 
  (is.na(is_bulk) | is_bulk == FALSE) &
  (is.na(is_summary) | is_summary == FALSE)
]

log_message(paste("Rows after filtering bulk and summary:", nrow(dt_filtered)))

# Clean cell type names
dt_filtered[, clean_cell_type := str_trim(Display_Name)]

# Select relevant columns
targets <- dt_filtered[, .(cell_type = clean_cell_type, snp = Lead_SNP)] %>% unique()

log_message(paste("Unique cell type - SNP pairs to extract:", nrow(targets)))

# Get unique cell types to iterate over
unique_cells <- unique(targets$cell_type)
log_message(paste("Unique cell types found:", paste(unique_cells, collapse = ", ")))

# 4.2 Extract Data
results_list <- list()
stats_list <- list()

start_time <- Sys.time()

for (cell in unique_cells) {
  log_message(paste("Processing cell type:", cell))
  
  # Construct file path
  file_name <- paste0(cell, ".tsv.gz")
  file_path <- file.path(source_data_dir, file_name)
  
  if (!file_exists(file_path)) {
    log_message(paste("WARNING: Source file not found for cell type:", cell, "- Path:", file_path))
    stats_list[[cell]] <- data.table(cell_type = cell, status = "File Not Found", snps_found = 0, snps_requested = 0)
    next
  }
  
  # Get target SNPs for this cell
  # SNP
  target_snps <- targets[cell_type == cell, snp]
  log_message(paste("  Target SNPs count:", length(target_snps)))
  
  # Read source file
  tryCatch({
    # Use fread with select to optimize if possible.
    # fread. OneK1K. 
    
    dt_source <- fread(file_path)
    
    # Filter by RSID
    # RSID
    if ("RSID" %in% names(dt_source)) {
      dt_extracted <- dt_source[RSID %in% target_snps]
    } else {
      log_message(paste("ERROR: 'RSID' column not found in file:", file_name))
      stats_list[[cell]] <- data.table(cell_type = cell, status = "Column Missing", snps_found = 0, snps_requested = length(target_snps))
      next
    }
    
    # Add input cell type info
    dt_extracted[, Input_Cell_Type := cell]
    
    count_found <- nrow(dt_extracted)
    log_message(paste("  SNPs found:", count_found))
    
    results_list[[cell]] <- dt_extracted
    stats_list[[cell]] <- data.table(cell_type = cell, status = "Success", snps_found = count_found, snps_requested = length(target_snps))
    
  }, error = function(e) {
    log_message(paste("ERROR reading file for cell:", cell, "- Message:", e$message))
    stats_list[[cell]] <- data.table(cell_type = cell, status = "Read Error", snps_found = 0, snps_requested = length(target_snps))
  })
}

# 4.3 Combine and Save Results
log_message("Combining results...")

if (length(results_list) > 0) {
  final_result <- rbindlist(results_list, fill = TRUE)
  
  # Save original result
  output_file <- file.path(output_dir, "OneK1K_extracted_snps.csv")
  fwrite(final_result, output_file)
  log_message(paste("Saved extracted data to:", output_file))
  
  # 4.4 Filter by FDR <= 0.05
  #  FDR <= 0.05
  log_message("Filtering data by FDR <= 0.05...")
  
  if ("FDR" %in% names(final_result)) {
    # Ensure FDR is numeric
    final_result[, FDR := as.numeric(FDR)]
    
    # Save FDR filtered result
    fdr_filtered_result <- final_result[FDR <= 0.05]
    output_file_fdr <- file.path(output_dir, "OneK1K_extracted_snps_FDR_0.05.csv")
    fwrite(fdr_filtered_result, output_file_fdr)
    log_message(paste("Saved FDR filtered data to:", output_file_fdr))
    log_message(paste("Rows with FDR <= 0.05:", nrow(fdr_filtered_result)))
    
    # 4.5 Create Unique Best P-value Datasets
    # P
    log_message("Filtering unique best P-value per Cell-SNP-Gene...")
    
    if ("P_VALUE" %in% names(final_result) & "GENE" %in% names(final_result)) {
      # Ensure P_VALUE is numeric
      final_result[, P_VALUE := as.numeric(P_VALUE)]
      
      # --- Version 1: Filtered by FDR <= 0.05 (Original) ---
      unique_best_result <- fdr_filtered_result[order(P_VALUE), .SD[1], by = .(Input_Cell_Type, RSID, GENE)]
      output_file_unique <- file.path(output_dir, "OneK1K_extracted_snps_FDR_0.05_Unique_Best.csv")
      fwrite(unique_best_result, output_file_unique)
      log_message(paste("Saved unique best P-value data (FDR<=0.05) to:", output_file_unique))
      log_message(paste("Rows after deduplication (FDR<=0.05):", nrow(unique_best_result)))

      # --- Version 2: All Genes (FDR > 0.05 included) ---
      # Deduplicate ALL data first
      unique_best_all_result <- final_result[order(P_VALUE), .SD[1], by = .(Input_Cell_Type, RSID, GENE)]
      output_file_unique_all <- file.path(output_dir, "OneK1K_extracted_snps_All_FDR_Unique_Best.csv")
      fwrite(unique_best_all_result, output_file_unique_all)
      log_message(paste("Saved unique best P-value data (All FDR) to:", output_file_unique_all))
      log_message(paste("Rows after deduplication (All FDR):", nrow(unique_best_all_result)))
      
      # 5. Visualization (Dot Plot) ----------------------------
      
      log_message("Generating visualizations...")
      
      if (nrow(unique_best_all_result) > 0) {
        # Create Plots Directory
        plot_dir <- file.path(output_dir, "plots")
        dir_create(plot_dir)
        
        # Helper function to generate and save plot
        create_dot_plot <- function(plot_data, filename_suffix, is_all_fdr = FALSE) {
            
            # Create local copy
            viz_data <- copy(plot_data)
            
            # Determine Significance if not present
            if (!"Is_Significant" %in% names(viz_data)) {
                if (is_all_fdr) {
                    viz_data[, Is_Significant := ifelse(FDR <= 0.05, "Significant", "Non-Significant")]
                } else {
                    viz_data[, Is_Significant := "Significant"]
                }
            }
        
            # Create Label: Cell Type + SNP
            viz_data[, Cell_SNP := paste(Input_Cell_Type, RSID, sep = " - ")]
            
            # Determine Direction and Abs Beta
            if ("BETA" %in% names(viz_data)) {
              viz_data[, Direction := ifelse(BETA > 0, "Positive", "Negative")]
              viz_data[, Abs_Beta := abs(BETA)]
            } else {
              log_message("WARNING: BETA column not found for visualization. Using dummy direction.")
              viz_data[, Direction := "Unknown"]
              viz_data[, Abs_Beta := 1]
              viz_data[, BETA := 1]
            }
            
            # --- Apply Reference Ordering & Intersection Mapping ---
            # Use dt_filtered which already has Parsed_Gene and is cleaned
            if (exists("dt_filtered") && nrow(dt_filtered) > 0) {
              ref_data_clean <- copy(dt_filtered)
              ref_data_clean[, Clean_Cell := str_trim(Display_Name)]
              ref_data_clean <- ref_data_clean[!is.na(Clean_Cell) & !is.na(Lead_SNP) & !is.na(Parsed_Gene) & Clean_Cell != "" & Lead_SNP != ""]
              
              # Aggregate Outcomes
              intersection_map <- ref_data_clean[, .(Outcomes = paste(sort(unique(Outcome)), collapse = " & ")), 
                                                 by = .(Input_Cell_Type = Clean_Cell, RSID = Lead_SNP, GENE = Parsed_Gene)]
              
              # Merge into viz_data
              viz_data <- merge(viz_data, intersection_map, by = c("Input_Cell_Type", "RSID", "GENE"), all.x = TRUE)
              
              # Fill NA for Novel associations
              viz_data[is.na(Outcomes), Outcomes := "OneK1K Novel"]
      
              # Create composite key for ordering Y-axis (Cell_SNP)
              ref_order_df <- unique(ref_data_clean[, .(Cell = Clean_Cell, Lead_SNP)])
              ref_order_vec <- paste(ref_order_df$Cell, ref_order_df$Lead_SNP, sep = " - ")
              
              # Apply ordering to viz_data
              current_items <- unique(viz_data$Cell_SNP)
              ordered_levels <- intersect(ref_order_vec, current_items)
              remaining_items <- setdiff(current_items, ordered_levels)
              final_levels <- c(ordered_levels, remaining_items)
              
              viz_data[, Cell_SNP := factor(Cell_SNP, levels = rev(final_levels))]
              
            } else {
              viz_data[, Cell_SNP := factor(Cell_SNP, levels = unique(viz_data[order(Input_Cell_Type, P_VALUE), Cell_SNP]))]
              viz_data[, Outcomes := "OneK1K Novel"]
            }
            
            # --- Order Genes by Genomic Position (EnsDb) ---
            genes_unique <- unique(viz_data$GENE)
            
            tryCatch({
              edb <- EnsDb.Hsapiens.v75
              
              if ("GENE_ID" %in% names(viz_data)) {
                gene_ids_unique <- unique(viz_data$GENE_ID)
                gene_info <- genes(edb, filter = GeneIdFilter(gene_ids_unique), columns = c("gene_id", "symbol", "seq_name", "gene_seq_start"))
                gene_df <- as.data.frame(gene_info)
                
                chrom_levels <- c(as.character(1:22), "X", "Y", "MT")
                gene_df$seq_name <- factor(gene_df$seq_name, levels = chrom_levels)
                gene_df <- gene_df[order(gene_df$seq_name, gene_df$gene_seq_start), ]
                
                order_lookup <- data.table(GENE_ID = gene_df$gene_id, Order_Index = 1:nrow(gene_df))
                temp_order <- viz_data[, .(GENE, GENE_ID)] %>% unique()
                temp_order <- merge(temp_order, order_lookup, by = "GENE_ID", all.x = TRUE)
                setorder(temp_order, Order_Index, na.last = TRUE)
                final_gene_levels <- temp_order$GENE
                
              } else {
                 gene_info <- genes(edb, filter = SymbolFilter(genes_unique), columns = c("symbol", "seq_name", "gene_seq_start"))
                 gene_df <- as.data.frame(gene_info)
                 
                 chrom_levels <- c(as.character(1:22), "X", "Y", "MT")
                 gene_df$seq_name <- factor(gene_df$seq_name, levels = chrom_levels)
                 gene_df <- gene_df[order(gene_df$seq_name, gene_df$gene_seq_start), ]
                 
                 ordered_genes <- unique(gene_df$symbol)
                 missing_genes <- setdiff(genes_unique, ordered_genes)
                 final_gene_levels <- c(ordered_genes, sort(missing_genes))
              }
    
              viz_data[, GENE := factor(GENE, levels = final_gene_levels)]
            }, error = function(e) {
              viz_data[, GENE := factor(GENE, levels = sort(unique(GENE)))]
            })
            
            if (!is.factor(viz_data$Cell_SNP)) {
               viz_data[, Cell_SNP := factor(Cell_SNP, levels = unique(viz_data[order(Input_Cell_Type, P_VALUE), Cell_SNP]))]
            }
            
            # Create a clean 'Group' column for plotting shapes
            viz_data[, Group := Outcomes]
            viz_data[is.na(Group), Group := "OneK1K Novel"]
            
            # Handle Non-Significant Genes
            if (is_all_fdr) {
                viz_data[Is_Significant == "Non-Significant", Group := "Non-Significant"]
            }
    
            # Define shapes
            unique_groups <- unique(viz_data$Group)
            shape_values <- setNames(rep(16, length(unique_groups)), unique_groups) 
            
            for (grp in unique_groups) {
              grp_lower <- tolower(grp)
              if (grepl("non-significant", grp_lower)) {
                 shape_values[grp] <- 16 
              } else if (grepl("aging", grp_lower) && grepl("hz", grp_lower) && grepl("ra", grp_lower)) {
                shape_values[grp] <- 18 
              } else if (grepl("aging", grp_lower) && grepl("hz", grp_lower)) {
                shape_values[grp] <- 17 
              } else if (grepl("aging", grp_lower) && grepl("ra", grp_lower)) {
                 shape_values[grp] <- 15 
              } else if (grepl("hz", grp_lower) && grepl("ra", grp_lower)) {
                 shape_values[grp] <- 8  
              } else if (grp == "OneK1K Novel") {
                 shape_values[grp] <- 16 
              } else {
                 shape_values[grp] <- 16 
              }
            }
            
            viz_data_sig <- viz_data[Is_Significant == "Significant"]
             viz_data_nonsig <- viz_data[Is_Significant == "Non-Significant"]
             
             # Adjust plot dimensions for square cells
             # Assume 1 unit = 0.3 inches
             n_genes <- length(unique(viz_data$GENE))
             n_snps <- length(unique(viz_data$Cell_SNP))
             
             cell_size_inch <- 0.3
             margin_inch <- 2
             
             plot_height <- max(5, n_snps * cell_size_inch + margin_inch) 
             plot_width <- max(10, n_genes * cell_size_inch + margin_inch)
             
             p_dot <- ggplot()
             
             # Layer 1: Grid lines (Manually added to create cells between points)
             # Points are at 1, 2, 3...
             # Grid lines should be at 0.5, 1.5, 2.5...
             p_dot <- p_dot +
                 geom_vline(xintercept = seq(0.5, n_genes + 0.5, 1), color = "grey90", size = 0.5) +
                 geom_hline(yintercept = seq(0.5, n_snps + 0.5, 1), color = "grey90", size = 0.5)
             
             # Layer 2: Non-Significant (if any)
             if (nrow(viz_data_nonsig) > 0) {
               p_dot <- p_dot + 
                   geom_point(data = viz_data_nonsig, 
                          aes(x = GENE, y = Cell_SNP, shape = Group), 
                          color = "grey80", size = 3.5, alpha = 0.6)
             }
               
             # Layer 3: Significant
             if (nrow(viz_data_sig) > 0) {
                p_dot <- p_dot +
                   geom_point(data = viz_data_sig, 
                          aes(x = GENE, y = Cell_SNP, color = BETA, shape = Group), 
                          size = 3.5, alpha = 0.9)
             }
               
             p_dot <- p_dot +
               scale_color_gradient2(low = "#4575B4", mid = "white", high = "#D73027", midpoint = 0) +
               scale_shape_manual(values = shape_values) +
               coord_fixed(ratio = 1) + # Ensure square cells
               labs(
                    x = "Gene",
                    y = "Cell Type - Lead SNP",
                    color = "Beta",
                    shape = "Intersection Group") +
               theme_minimal(base_size = 14) +
               theme(
                 axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
                 axis.text.y = element_text(size = 10),
                 panel.grid = element_blank(), # Remove default grid
                 legend.position = "right",
                 plot.margin = margin(20, 20, 20, 20)
               )
            
            # FDR<=0.05
            if (filename_suffix == "DotPlot_Cell_SNP_Gene_3way_Intersection_FDR_0.05") {
               p_dot <- p_dot + theme(
                 legend.title = element_text(size = 16),
                 legend.text = element_text(size = 14),
                 axis.title = element_text(size = 18),
                 axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
                 axis.text.y = element_text(size = 14)
               )
            }
            
            ggsave(file.path(plot_dir, paste0(filename_suffix, ".pdf")), p_dot, width = plot_width, height = plot_height, limitsize = FALSE)
            ggsave(file.path(plot_dir, paste0(filename_suffix, ".png")), p_dot, width = plot_width, height = plot_height, dpi = 300, limitsize = FALSE)
            
            log_message(paste("Saved plot:", filename_suffix))
        }
        
        # --- Version 1: All FDR (Original Request) ---
        log_message("Creating Plot Version 1: All FDR...")
        create_dot_plot(unique_best_all_result, "DotPlot_Cell_SNP_Gene_All_FDR", is_all_fdr = TRUE)
        
        # --- Version 2: FDR <= 0.05 (Significant Only) ---
        log_message("Creating Plot Version 2: FDR <= 0.05 Only...")
        if (nrow(unique_best_result) > 0) {
            create_dot_plot(unique_best_result, "DotPlot_Cell_SNP_Gene_FDR_0.05", is_all_fdr = FALSE)
        } else {
            log_message("No significant results (FDR<=0.05) to plot.")
        }

        # --- Version 3: 3-way Intersection Only (Comprehensive Pipeline) ---
        log_message("Starting Comprehensive 3-way Intersection Analysis...")
        
        # 1. Identify 3-way intersection triplets from input
        # Note: dt_filtered is available in scope
        intersection_triplets <- dt_filtered[, .(Outcomes_Count = length(unique(Outcome)), 
                                                 Outcomes_List = list(unique(Outcome))), 
                                             by = .(Input_Cell_Type = clean_cell_type, RSID = Lead_SNP, GENE = Parsed_Gene)]
        
        target_outcomes <- c("aging", "HZ", "RA")
        # Check if all target outcomes are present in the list
        three_way_triplets <- intersection_triplets[sapply(Outcomes_List, function(x) all(target_outcomes %in% x))]
        
        log_message(paste("Number of 3-way intersection triplets found in input:", nrow(three_way_triplets)))
        
        if (nrow(three_way_triplets) > 0) {
            
            # 2. Extract Raw 3-way Data (from final_result)
            # Match on Cell and SNP ONLY to get ALL gene associations (including Novel) for these triplets
            # CellSNP, triplets（Novel）
            target_cell_snps <- unique(three_way_triplets[, .(Input_Cell_Type, RSID)])
            
            final_result_3way <- merge(final_result, target_cell_snps, 
                                       by = c("Input_Cell_Type", "RSID"))
            
            output_file_3way_raw <- file.path(output_dir, "OneK1K_extracted_snps_3way_intersection_Raw.csv")
            fwrite(final_result_3way, output_file_3way_raw)
            log_message(paste("Saved 3-way Raw data to:", output_file_3way_raw))
            
            # 3. Filter FDR <= 0.05 for 3-way
            if ("FDR" %in% names(final_result_3way)) {
                fdr_filtered_3way <- final_result_3way[FDR <= 0.05]
                output_file_3way_fdr <- file.path(output_dir, "OneK1K_extracted_snps_3way_intersection_FDR_0.05.csv")
                fwrite(fdr_filtered_3way, output_file_3way_fdr)
                log_message(paste("Saved 3-way FDR filtered data to:", output_file_3way_fdr))
                
                # 4. Unique Best P-value (FDR <= 0.05)
                unique_best_3way_sig <- fdr_filtered_3way[order(P_VALUE), .SD[1], by = .(Input_Cell_Type, RSID, GENE)]
                output_file_3way_unique_sig <- file.path(output_dir, "OneK1K_extracted_snps_3way_intersection_FDR_0.05_Unique_Best.csv")
                fwrite(unique_best_3way_sig, output_file_3way_unique_sig)
                log_message(paste("Saved 3-way Unique Best (FDR<=0.05) to:", output_file_3way_unique_sig))
            }
            
            # 5. Unique Best P-value (All FDR)
            unique_best_3way_all <- final_result_3way[order(P_VALUE), .SD[1], by = .(Input_Cell_Type, RSID, GENE)]
            output_file_3way_unique_all <- file.path(output_dir, "OneK1K_extracted_snps_3way_intersection_All_FDR_Unique_Best.csv")
            fwrite(unique_best_3way_all, output_file_3way_unique_all)
            log_message(paste("Saved 3-way Unique Best (All FDR) to:", output_file_3way_unique_all))
            
            # 6. Visualization
            if (nrow(unique_best_3way_all) > 0) {
                # Plot 1: All FDR
                log_message("Creating 3-way Intersection Plot (All FDR)...")
                create_dot_plot(unique_best_3way_all, "DotPlot_Cell_SNP_Gene_3way_Intersection_All_FDR", is_all_fdr = TRUE)
                
                # Plot 2: FDR <= 0.05
                log_message("Creating 3-way Intersection Plot (FDR <= 0.05)...")
                if (exists("unique_best_3way_sig") && nrow(unique_best_3way_sig) > 0) {
                   create_dot_plot(unique_best_3way_sig, "DotPlot_Cell_SNP_Gene_3way_Intersection_FDR_0.05", is_all_fdr = FALSE)
                } else {
                   log_message("No significant results (FDR<=0.05) for 3-way intersection to plot.")
                }
            }
        } else {
            log_message("No 3-way intersection triplets identified in input.")
        }

        log_message("Visualizations created.")
      } else {
        log_message("No unique results found to visualize.")
      }
    } else {
      log_message("WARNING: P_VALUE or GENE column missing. Skipping unique filtering.")
    }
  } else {
    log_message("WARNING: FDR column missing. Skipping FDR filtering.")
  }
} else {
  log_message("No results extracted.")
}

# 6. Generate Methodology Description (MD) ----------------------------------
log_message("Generating Methodology Description...")

method_md_content <- glue("
# Methodology: Investigation of SNP Pleiotropy in OneK1K Immune Cells

## Overview
This analysis investigates the pleiotropic effects of specific SNPs (identified from GWAS/eQTL studies) on gene expression across various immune cell types using the OneK1K single-cell eQTL dataset. The goal is to map the regulatory landscape of these SNPs and identify both known and novel target genes.

## Data Sources
1. **Input SNPs**: A hierarchical list of Cell Type - Lead SNP pairs derived from prior integration of GWAS and eQTL data (e.g., intersection of Aging, Herpes Zoster, and Rheumatoid Arthritis signals).
2. **OneK1K Dataset**: Single-cell eQTL summary statistics from the OneK1K cohort, providing cell-type-specific associations between SNPs and gene expression.

## Processing Steps

### 1. Data Extraction
- **Input Parsing**: The input list is parsed to identify target Cell Types and Lead SNPs.
- **Source Query**: For each target Cell Type, the corresponding OneK1K summary statistics file is queried.
- **SNP Filtering**: Rows matching the target Lead SNPs (RSID) are extracted.

### 2. Data Filtering and Deduplication
- **FDR Thresholding**: 
  - A 'Significant' set is defined by filtering for False Discovery Rate (FDR) ≤ 0.05.
  - A 'Complete' set includes all associations regardless of FDR (retaining non-significant results for background visualization).
- **Best P-value Selection**: 
  - To resolve potential duplicates (e.g., multiple probes or transcripts for the same gene), the data is grouped by `Cell Type`, `SNP`, and `Gene`.
  - Only the entry with the lowest P-value is retained for each unique triplet.

### 3. Annotation and Categorization
- **Intersection Mapping**: Extracted associations are cross-referenced with the input phenotypes (Aging, HZ, RA) to categorize them into intersection groups (e.g., 'Aging & HZ', '3-way Intersection').
- **Novelty Detection**: Associations found in OneK1K but not present in the initial input list are labeled as 'OneK1K Novel'.
- **Genomic Ordering**: Genes are ordered by their genomic coordinates (Chromosome + Start Position) using the `EnsDb.Hsapiens.v75` database (hg19 reference) to provide a biologically meaningful visualization layout.

### 4. Visualization
- **Dot Plot**:
  - **Axes**: X-axis represents Genes (ordered by position); Y-axis represents Cell-SNP pairs.
  - **Representation**: 
    - **Significant (FDR ≤ 0.05)**: Displayed as colored shapes. 
      - **Color**: Indicates Beta value (Effect Size), diverging from Blue (Negative) to Red (Positive).
      - **Shape**: Indicates the intersection group (e.g., Diamond for 3-way intersection, Triangle for 2-way).
    - **Non-Significant (FDR > 0.05)**: Displayed as **Gray Circles** to illustrate the testing background and specificity of effects.
  - **Aesthetics**: Titles are omitted for publication readiness.

## Output Files
- **CSV Data**: Detailed tables of extracted SNPs (All and FDR-filtered).
- **Plots**: High-resolution Dot Plots (PDF/PNG) showing the landscape of associations.
- **Statistics**: Summary counts of extracted and filtered associations.

Generated on: {Sys.time()}
")

write_file(method_md_content, file.path(output_dir, "Methodology_Description.md"))
log_message("Methodology Description saved.")

# 7. Final Statistics and Report --------------------------------------------
if (length(stats_list) > 0) {
  stats_df <- rbindlist(stats_list)
  stats_file <- file.path(output_dir, "extraction_statistics.csv")
  fwrite(stats_df, stats_file)
  log_message(paste("Saved statistics to:", stats_file))
}

end_time <- Sys.time()
run_time <- end_time - start_time

log_message("Script completed successfully.")
log_message(paste("Total Run Time:", round(run_time, 2), units(run_time)))
