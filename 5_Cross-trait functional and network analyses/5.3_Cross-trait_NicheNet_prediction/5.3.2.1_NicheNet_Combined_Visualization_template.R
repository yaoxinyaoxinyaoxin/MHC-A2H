#!/usr/bin/env Rscript
# ==============================================================================
# [Script]: 5.3.2.1_NicheNet_Combined_Visualization_template.R
# [Method]: NicheNet intercellular communication visualization
# [Step]: Combine and visualize NicheNet results across phenotypes
# 
# [Function]:
# Generates combined visualizations for NicheNet analysis across three phenotypes.
# 
# [Input]:
#   --aging_up   : Upstream NicheNet CSV for Aging
#   --aging_down : Downstream NicheNet CSV for Aging
#   --ra_up      : Upstream NicheNet CSV for RA
#   --ra_down    : Downstream NicheNet CSV for RA
#   --hz_up      : Upstream NicheNet CSV for HZ
#   --hz_down    : Downstream NicheNet CSV for HZ
#   --out_dir    : Output directory path
# 
# [Output]:
# Combined PDF/PNG plots and statistical tables.
# ==============================================================================

rm(list = ls())
gc()

# 1. Load Dependencies
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(ggalluvial)
  library(stringr)
  library(fs)
  library(RColorBrewer)
  library(patchwork)
  library(grid)
  library(cowplot)
  library(tidyr)
})

# 2. Parse Command Line Arguments
option_list <- list(
  make_option(c("--aging_up"), type="character", default=NULL, help="Upstream CSV for Aging"),
  make_option(c("--aging_down"), type="character", default=NULL, help="Downstream CSV for Aging"),
  make_option(c("--ra_up"), type="character", default=NULL, help="Upstream CSV for RA"),
  make_option(c("--ra_down"), type="character", default=NULL, help="Downstream CSV for RA"),
  make_option(c("--hz_up"), type="character", default=NULL, help="Upstream CSV for HZ"),
  make_option(c("--hz_down"), type="character", default=NULL, help="Downstream CSV for HZ"),
  make_option(c("--out_dir"), type="character", default="./NicheNet_Combined", help="Output directory path")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$aging_up) || is.null(opt$aging_down) || is.null(opt$ra_up) || is.null(opt$ra_down) || is.null(opt$hz_up) || is.null(opt$hz_down)) {
  print_help(opt_parser)
  stop("Missing input files. Please provide all 6 input CSV paths.")
}

# 3. Define Paths
AGING_UPSTREAM <- opt$aging_up
AGING_DOWNSTREAM <- opt$aging_down
RA_UPSTREAM <- opt$ra_up
RA_DOWNSTREAM <- opt$ra_down
HZ_UPSTREAM <- opt$hz_up
HZ_DOWNSTREAM <- opt$hz_down

TIMESTAMP <- format(Sys.time(), "%Y%m%d_%H%M%S")
OUTPUT_DIR <- file.path(opt$out_dir, paste0("Combined_Visualization_", TIMESTAMP))
dir_create(OUTPUT_DIR)
STATS_DIR <- file.path(OUTPUT_DIR, "stats")
dir_create(STATS_DIR)
SUPP_DIR <- file.path(OUTPUT_DIR, "supplement")
dir_create(SUPP_DIR)

# Initialize Log
LOG_FILE <- file.path(OUTPUT_DIR, "analysis.log")
log_message <- function(msg) {
  timestamp <- format(Sys.time(), "[%Y-%m-%d %H:%M:%S]")
  formatted_msg <- paste(timestamp, msg)
  cat(formatted_msg, "
")
  write(formatted_msg, file = LOG_FILE, append = TRUE)
}

log_message("Starting Combined NicheNet Visualization (mvAge / RA / HZ UKB) - Template")

# 3. Global Configuration
# ------------------------------------------------------------------------------
# Abbreviation Function (Copied from Source Script)
abbreviate_cell_types <- function(x) {
  x <- as.character(x)
  
  # Define mapping based on the provided Excel file
  mapping <- c(
    "CD4+ naive and central memory T cells" = "CD4NC",
    "CD8+ naive and central memory T cells" = "CD8NC",
    "CD4+ T cells with an effector memory or central memory phenotype" = "CD4ET",
    "CD4+ T cells with an effector memory phenotype" = "CD4ET", # Handling potential variation
    "CD8+ T cells with an effector memory phenotype" = "CD8ET",
    "CD8+ T cells expressing S100B" = "CD8S100B",
    "CD4+ T cells expressing SOX4" = "CD4SOX4",
    "Immature and naive B cells" = "BIN",
    "Memory B cells" = "BMem",
    "Plasma cells" = "Plasma cells",
    "Natural killer cells" = "NK cells",
    "NK recruiting cells" = "NK recruiting cells",
    "Classical monocytes" = "MonoC",
    "Nonclassical monocytes" = "MonoNC",
    "Dendritic cells" = "DCs"
  )
  
  # Use lookup, keep original if not found
  new_x <- mapping[x]
  ifelse(is.na(new_x), x, new_x)
}

# Helper to clean names for files
make_clean_names <- function(s) {
  s <- gsub(" ", "_", s)
  return(s)
}

# 4.1 Top-30 Lists & 3-way Intersection Stats / Top30
# ------------------------------------------------------------------------------
require_columns <- function(df, required_cols, label) {
  missing_cols <- setdiff(required_cols, colnames(df))
  if (length(missing_cols) > 0) {
    stop(paste0(label, " is missing required columns: ", paste(missing_cols, collapse = ", ")))
  }
}

compute_top30_upstream_ligands <- function(upstream_path, sort_by = c("Ligand_Activity", "Total_Ligand_Activity")) {
  sort_by <- match.arg(sort_by)
  if (!file_exists(upstream_path)) {
    stop(paste0("Upstream file not found: ", upstream_path))
  }
  upstream_data <- read_csv(upstream_path, show_col_types = FALSE)
  require_columns(
    upstream_data,
    c("ligand", "sender_cell", "receiver_cell", "regulatory_potential", "ligand_activity_pearson"),
    "Upstream CSV"
  )

  upstream_ligand_unique <- upstream_data %>%
    select(ligand, sender_cell, receiver_cell, regulatory_potential, ligand_activity_pearson) %>%
    distinct() %>%
    mutate(
      regulatory_potential = suppressWarnings(as.numeric(regulatory_potential)),
      ligand_activity_pearson = suppressWarnings(as.numeric(ligand_activity_pearson))
    )

  ligand_strength_summary <- upstream_ligand_unique %>%
    group_by(ligand) %>%
    summarise(
      Total_Count = n(),
      Total_Strength = sum(regulatory_potential, na.rm = TRUE),
      Mean_Strength = mean(regulatory_potential, na.rm = TRUE),
      Total_Ligand_Activity = sum(ligand_activity_pearson, na.rm = TRUE),
      .groups = "drop"
    )

  ligand_activity_summary <- upstream_data %>%
    select(ligand, ligand_activity_pearson) %>%
    distinct() %>%
    mutate(ligand_activity_pearson = suppressWarnings(as.numeric(ligand_activity_pearson))) %>%
    group_by(ligand) %>%
    summarise(
      Ligand_Activity = max(ligand_activity_pearson, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(Ligand_Activity = ifelse(is.infinite(Ligand_Activity), NA_real_, Ligand_Activity))

  ligand_summary <- ligand_strength_summary %>%
    left_join(ligand_activity_summary, by = "ligand")

  rank_table <- ligand_summary
  if (sort_by == "Ligand_Activity") {
    rank_table <- rank_table %>%
      arrange(desc(Ligand_Activity), desc(Total_Strength), desc(Total_Count), ligand)
  } else {
    rank_table <- rank_table %>%
      arrange(desc(Total_Ligand_Activity), desc(Ligand_Activity), desc(Total_Count), ligand)
  }

  top_ligands <- rank_table %>%
    slice_head(n = 30) %>%
    mutate(Rank = row_number()) %>%
    select(Rank, ligand, Ligand_Activity, Total_Ligand_Activity, Total_Strength, Total_Count, Mean_Strength)

  top_ligands
}

compute_top30_downstream_targets <- function(downstream_path, sort_by = c("Total_Strength", "Max_Single_Strength")) {
  sort_by <- match.arg(sort_by)
  if (!file_exists(downstream_path)) {
    stop(paste0("Downstream file not found: ", downstream_path))
  }
  downstream_data <- read_csv(downstream_path, show_col_types = FALSE)
  require_columns(
    downstream_data,
    c("target_gene", "receiver_cell", "target_regulatory_potential"),
    "Downstream CSV"
  )

  target_gene_long <- downstream_data %>%
    select(target_gene, receiver_cell, target_regulatory_potential) %>%
    mutate(target_regulatory_potential = suppressWarnings(as.numeric(target_regulatory_potential)))

  gene_summary <- target_gene_long %>%
    group_by(target_gene) %>%
    summarise(
      Total_Count = n(),
      Total_Strength = sum(target_regulatory_potential, na.rm = TRUE),
      Mean_Strength = mean(target_regulatory_potential, na.rm = TRUE),
      Max_Single_Strength = max(target_regulatory_potential, na.rm = TRUE),
      .groups = "drop"
    )

  rank_table <- gene_summary
  if (sort_by == "Total_Strength") {
    rank_table <- rank_table %>%
      arrange(desc(Total_Strength), desc(Total_Count), target_gene)
  } else {
    rank_table <- rank_table %>%
      arrange(desc(Max_Single_Strength), desc(Total_Strength), desc(Total_Count), target_gene)
  }

  top_genes <- rank_table %>%
    slice_head(n = 30) %>%
    mutate(Rank = row_number()) %>%
    select(Rank, target_gene, Max_Single_Strength, Total_Strength, Total_Count, Mean_Strength)

  top_genes
}

write_threeway_intersection_stats <- function(top_a, top_b, top_c, key_col, label_a, label_b, label_c, out_csv) {
  list_a <- top_a[[key_col]]
  list_b <- top_b[[key_col]]
  list_c <- top_c[[key_col]]

  inter <- Reduce(intersect, list(list_a, list_b, list_c))
  inter <- sort(unique(inter))

  if (length(inter) == 0) {
    out <- tibble::tibble(!!key_col := character(0))
    write_csv(out, out_csv)
    return(list(intersection = inter, table = out))
  }

  a_sub <- top_a %>% filter(.data[[key_col]] %in% inter) %>% mutate(Phenotype = label_a)
  b_sub <- top_b %>% filter(.data[[key_col]] %in% inter) %>% mutate(Phenotype = label_b)
  c_sub <- top_c %>% filter(.data[[key_col]] %in% inter) %>% mutate(Phenotype = label_c)

  merged <- bind_rows(a_sub, b_sub, c_sub) %>%
    arrange(.data[[key_col]], Phenotype)

  wide <- merged %>%
    pivot_wider(
      names_from = Phenotype,
      values_from = setdiff(colnames(merged), c(key_col, "Phenotype")),
      names_sep = "__"
    ) %>%
    arrange(.data[[key_col]])

  write_csv(wide, out_csv)
  list(intersection = inter, table = wide)
}

detect_categories <- function(symbols) {
  symbols <- unique(na.omit(as.character(symbols)))
  symbols_upper <- toupper(symbols)

  is_chemokine <- grepl("^(CCL|CXCL|XCL|CX3CL)", symbols_upper)
  is_type1_ifn <- grepl("^IFNA\\d+$", symbols_upper) | symbols_upper %in% c("IFNA", "IFNB1", "IFNB", "IFNE", "IFNK", "IFNW1", "IFNL1", "IFNL2", "IFNL3", "IFNL4")
  is_tnf_il1_il6_axis <- symbols_upper %in% c("TNF", "TNFSF10", "TNFSF11", "TNFSF13B", "IL1A", "IL1B", "IL1RN", "IL2", "IL6", "IL6R", "IL6ST", "IL16", "IL18", "IL18R1", "IL18RAP", "IL10", "IFNG", "TGFB1", "TGFB2", "TGFB3")
  is_ifn_stimulated_gene <- grepl("^(ISG|IFIT|MX|OAS|RSAD|IFI|IRF|STAT)", symbols_upper) | symbols_upper %in% c("BST2", "DDX58", "EIF2AK2", "TRIM22", "TRIM25", "USP18")
  is_antigen_presentation <- grepl("^HLA-", symbols_upper) | symbols_upper %in% c("B2M", "TAP1", "TAP2", "TAPBP", "CIITA")
  is_nfkb_pathway <- symbols_upper %in% c("NFKBIA", "NFKB1", "NFKB2", "RELA", "REL", "RELB", "IKBKB", "IKBKG", "RIPK2", "TNFAIP3", "TNIP1")
  is_cytotoxicity <- symbols_upper %in% c("GZMB", "GZMA", "PRF1", "GNLY", "NKG7")
  is_immune_receptors_adhesion <- grepl("^(CD|HLA-)", symbols_upper) | symbols_upper %in% c("SELPLG", "FCER2", "TNFRSF14")

  list(
    chemokines = symbols[is_chemokine],
    typeI_interferons = symbols[is_type1_ifn],
    tnf_il1_il6_axis = symbols[is_tnf_il1_il6_axis],
    ifn_stimulated_genes = symbols[is_ifn_stimulated_gene],
    antigen_presentation = symbols[is_antigen_presentation],
    nfkb_pathway = symbols[is_nfkb_pathway],
    cytotoxicity = symbols[is_cytotoxicity],
    immune_receptors_adhesion = symbols[is_immune_receptors_adhesion]
  )
}

write_intersection_interpretation_md <- function(
  md_path,
  upstream_intersection,
  downstream_intersection,
  upstream_csv,
  downstream_csv
) {
  upstream_intersection <- sort(unique(na.omit(as.character(upstream_intersection))))
  downstream_intersection <- sort(unique(na.omit(as.character(downstream_intersection))))

  up_cat <- detect_categories(upstream_intersection)
  down_cat <- detect_categories(downstream_intersection)

  collapse_or_none <- function(x) {
    if (length(x) == 0) return("None")
    paste(x, collapse = ", ")
  }

  paragraph_for_presence <- function(x, if_present, if_absent) {
    if (length(x) == 0) return(if_absent)
    if_present
  }

  md_lines <- c(
    "# Three-way Intersection of Top30 NicheNet Results (mvAge / RA / HZ UKB)",
    "",
    "## What is summarized",
    "",
    paste0("- **Upstream ligands**: three-way intersection among the **Top 30** upstream ligands of mvAge, RA, and HZ UKB."),
    paste0("- **Downstream target genes**: three-way intersection among the **Top 30** downstream targets of mvAge, RA, and HZ UKB."),
    paste0("- Detailed per-phenotype ranks and scores are saved as CSV:"),
    paste0("  - ", basename(upstream_csv)),
    paste0("  - ", basename(downstream_csv)),
    "",
    "## Results",
    "",
    paste0("### Upstream ligands (Top30 ∩ Top30 ∩ Top30), n = ", length(upstream_intersection)),
    "",
    if (length(upstream_intersection) == 0) "- (empty intersection)" else paste0("- ", collapse_or_none(upstream_intersection)),
    "",
    paste0("### Downstream targets (Top30 ∩ Top30 ∩ Top30), n = ", length(downstream_intersection)),
    "",
    if (length(downstream_intersection) == 0) "- (empty intersection)" else paste0("- ", collapse_or_none(downstream_intersection)),
    "",
    "## Interpretation (literature-grounded, pathway-level)",
    "",
    "This interpretation focuses on **shared biological themes** that can plausibly appear across mvAge (inflammaging), RA (chronic autoimmune inflammation), and HZ UKB (VZV reactivation and antiviral immunity). It is intentionally pathway-level to avoid over-interpreting individual genes without additional validation.",
    "",
    "### 1) Shared inflammatory cytokine axis",
    "",
    paste0("- Intersection cytokines (TNF/IL1/IL6/TGF/IL10-related): ", collapse_or_none(up_cat$tnf_il1_il6_axis)),
    paragraph_for_presence(
      up_cat$tnf_il1_il6_axis,
      "These shared cytokines support a common denominator of **immune activation and inflammatory signaling** that can plausibly span RA inflammation, age-associated immune dysregulation, and antiviral responses.",
      "No canonical TNF/IL-1/IL-6/TGF/IL-10 family cytokines are present in the three-way Top30 ligand intersection, suggesting the shared signal may be more **chemokine / receptor / innate-program** driven in this ranking."
    ),
    "",
    "### 2) Chemokine-driven immune cell recruitment",
    "",
    paste0("- Intersection chemokines: ", collapse_or_none(up_cat$chemokines)),
    paragraph_for_presence(
      up_cat$chemokines,
      "Chemokines indicate shared **leukocyte trafficking and tissue recruitment** programs. CXCL10 (IP-10) is well-described in RA synovium and can be induced by interferon-γ, aligning with inflammatory recruitment signatures and antiviral contexts.",
      "No chemokines are present in the three-way Top30 ligand intersection."
    ),
    "",
    "### 3) Antiviral/type I interferon programs",
    "",
    paste0("- Upstream type I IFNs: ", collapse_or_none(up_cat$typeI_interferons)),
    paste0("- Downstream IFN-stimulated genes: ", collapse_or_none(down_cat$ifn_stimulated_genes)),
    paragraph_for_presence(
      down_cat$ifn_stimulated_genes,
      "An IFN/ISG-enriched downstream signature supports a shared **viral-sensing / interferon-response** module. This is biologically expected for HZ UKB (VZV reactivation) and can overlap with immune activation states observed in RA and mvAge-related immune dysregulation.",
      "No obvious ISG/IFN-response targets are present in the three-way Top30 target intersection."
    ),
    "",
    "### 4) Antigen presentation / adaptive immunity",
    "",
    paste0("- Antigen presentation genes (e.g., HLA, B2M, TAP): ", collapse_or_none(down_cat$antigen_presentation)),
    paragraph_for_presence(
      down_cat$antigen_presentation,
      "These suggest shared changes in **antigen processing and presentation**, consistent with activation of adaptive immune responses that can appear in both autoimmunity (RA) and antiviral immunity (HZ UKB), and can be modulated with age.",
      "No clear antigen-processing/presentation targets are present in the three-way Top30 target intersection."
    ),
    "",
    "### 5) NF-κB and cytotoxicity signatures",
    "",
    paste0("- NF-κB pathway genes: ", collapse_or_none(down_cat$nfkb_pathway)),
    paste0("- Cytotoxicity genes (GZMB/PRF1/NKG7...): ", collapse_or_none(down_cat$cytotoxicity)),
    paragraph_for_presence(
      down_cat$nfkb_pathway,
      "Shared NF-κB regulators/mediators support a conserved **inflammatory transcriptional program**, which is central in RA pathobiology and can be engaged downstream of innate sensing and cytokine signaling in broader inflammatory settings.",
      "No NF-κB pathway genes are present in the three-way Top30 target intersection."
    ),
    "",
    "### 6) Immune receptor / adhesion / checkpoint-like ligands",
    "",
    paste0("- Immune receptors/adhesion (CD/HLA/SELPLG/FCER2/TNFRSF14): ", collapse_or_none(up_cat$immune_receptors_adhesion)),
    paragraph_for_presence(
      up_cat$immune_receptors_adhesion,
      "This category suggests shared **cell-cell recognition / adhesion / immune modulation** signals. For example, CD47 is a broadly expressed ligand for SIRPα that can modulate myeloid cell responses (often discussed as a myeloid 'don't-eat-me' checkpoint in other contexts), highlighting that core myeloid-interaction modules can emerge as shared nodes across phenotypes.",
      "No immune receptor/adhesion-like ligands are present in the three-way Top30 ligand intersection."
    ),
    "",
    "## References (key background)",
    "",
    "- Browaeys R, Saelens W, Saeys Y. NicheNet: modeling intercellular communication by linking ligands to target genes. Nat Methods. 2020. doi:10.1038/s41592-019-0667-5.",
    "- Franceschi C, Garagnani P, Parini P, Giuliani C, Santoro A. Inflammaging: a new immune–metabolic viewpoint for age-related diseases. Nat Rev Endocrinol. 2018. doi:10.1038/s41574-018-0059-4.",
    "- McInnes IB, Schett G. The pathogenesis of rheumatoid arthritis. N Engl J Med. 2011. doi:10.1056/NEJMra1004965.",
    "- McInnes IB, Schett G. Cytokines in the pathogenesis of rheumatoid arthritis. Nat Rev Immunol. 2007. doi:10.1038/nri2094.",
    "- Gershon AA, Breuer J, Cohen JI, et al. Varicella zoster virus infection. Nat Rev Dis Primers. 2015. doi:10.1038/nrdp.2015.16.",
    "- Laing KJ, et al. Immunobiology of Varicella-Zoster Virus Infection. J Infect Dis. 2018. (Free PMC review; see PubMed: PMID 30535082.)",
    "- Laragione T, et al. CXCL10 and its receptor CXCR3 regulate synovial fibroblast invasion in rheumatoid arthritis. Arthritis Rheum. 2011. doi:10.1002/art.30573.",
    "- Logtenberg MEW, Scheeren FA, Schumacher TN. The CD47-SIRPα immune checkpoint. Immunity. 2020. doi:10.1016/j.immuni.2020.04.011.",
    "- Croft M, Siegel RM. Beyond TNF: TNF superfamily cytokines as targets for the treatment of rheumatic diseases. Nat Rev Rheumatol. 2017. doi:10.1038/nrrheum.2017.22."
  )

  writeLines(md_lines, md_path, useBytes = TRUE)
}

# 4. Pre-scan for All Cell Types (To Ensure Unified Legend) / 
# ------------------------------------------------------------------------------
log_message("Pre-scanning input files for comprehensive cell type list...")

get_all_cells <- function(up_path, down_path) {
  if (!file_exists(up_path) || !file_exists(down_path)) return(character(0))
  
  up <- read_csv(up_path, show_col_types = FALSE)
  down <- read_csv(down_path, show_col_types = FALSE)
  
  cells <- unique(c(
    up$sender_cell, up$receiver_cell,
    down$sender_cell, down$receiver_cell
  ))
  
  return(cells)
}

all_cells_aging <- get_all_cells(AGING_UPSTREAM, AGING_DOWNSTREAM)
all_cells_ra <- get_all_cells(RA_UPSTREAM, RA_DOWNSTREAM)
all_cells_hz <- get_all_cells(HZ_UPSTREAM, HZ_DOWNSTREAM)

all_raw_cells <- unique(c(all_cells_aging, all_cells_ra, all_cells_hz))
all_abbr_cells <- unique(abbreviate_cell_types(all_raw_cells))
all_abbr_cells <- all_abbr_cells[!is.na(all_abbr_cells) & all_abbr_cells != "None"]

log_message(paste("Found", length(all_abbr_cells), "unique cell types across all projects."))

# Define Palette (Matching Source Scripts)
palette_base <- c(
  "CD4NC" = "#1f77b4",
  "CD8NC" = "#ff7f0e",
  "CD4ET" = "#2ca02c",
  "CD8ET" = "#d62728",
  "CD8S100B" = "#9467bd",
  "CD4SOX4" = "#8c564b",
  "BIN" = "#e377c2",
  "BMem" = "#7f7f7f",
  "Plasma cells" = "#bcbd22",
  "NK cells" = "#17becf",
  "NK recruiting cells" = "#aec7e8",
  "MonoC" = "#ffbb78",
  "MonoNC" = "#98df8a",
  "DCs" = "#ff9896"
)

# Handle Unknown Cells
unknown_cells <- setdiff(all_abbr_cells, names(palette_base))
if (length(unknown_cells) > 0) {
  log_message(paste("Generating colors for unknown cells:", paste(unknown_cells, collapse=", ")))
  extra_cols <- colorRampPalette(brewer.pal(12, "Paired"))(length(unknown_cells))
  palette_final <- c(palette_base, setNames(extra_cols, unknown_cells))
} else {
  palette_final <- palette_base
}

# Ensure all needed cells are in palette
cell_colors <- palette_final[all_abbr_cells]
cell_colors["None"] <- "transparent"

# 4.2 Compute 3-way Intersections (Top 30) / Top30
# ------------------------------------------------------------------------------
log_message("Computing Top30 intersections across phenotypes...")

top30_ligands_aging_activity <- compute_top30_upstream_ligands(AGING_UPSTREAM, sort_by = "Ligand_Activity")
top30_ligands_ra_activity <- compute_top30_upstream_ligands(RA_UPSTREAM, sort_by = "Ligand_Activity")
top30_ligands_hz_activity <- compute_top30_upstream_ligands(HZ_UPSTREAM, sort_by = "Ligand_Activity")

top30_ligands_aging_total_activity <- compute_top30_upstream_ligands(AGING_UPSTREAM, sort_by = "Total_Ligand_Activity")
top30_ligands_ra_total_activity <- compute_top30_upstream_ligands(RA_UPSTREAM, sort_by = "Total_Ligand_Activity")
top30_ligands_hz_total_activity <- compute_top30_upstream_ligands(HZ_UPSTREAM, sort_by = "Total_Ligand_Activity")

top30_ligands_aging <- top30_ligands_aging_activity
top30_ligands_ra <- top30_ligands_ra_activity
top30_ligands_hz <- top30_ligands_hz_activity

ligand_intersection_out <- file.path(STATS_DIR, "Top30_Upstream_Ligands_Intersection_mvAge_RA_HZ_UKB.csv")
ligand_intersection_stats <- write_threeway_intersection_stats(
  top_a = top30_ligands_aging,
  top_b = top30_ligands_ra,
  top_c = top30_ligands_hz,
  key_col = "ligand",
  label_a = "mvAge",
  label_b = "RA",
  label_c = "HZ UKB",
  out_csv = ligand_intersection_out
)
log_message(paste0("Three-way upstream ligand intersection size: ", length(ligand_intersection_stats$intersection)))
log_message(paste("Saved:", ligand_intersection_out))

top30_targets_aging_total_strength <- compute_top30_downstream_targets(AGING_DOWNSTREAM, sort_by = "Total_Strength")
top30_targets_ra_total_strength <- compute_top30_downstream_targets(RA_DOWNSTREAM, sort_by = "Total_Strength")
top30_targets_hz_total_strength <- compute_top30_downstream_targets(HZ_DOWNSTREAM, sort_by = "Total_Strength")

top30_targets_aging_single_strength <- compute_top30_downstream_targets(AGING_DOWNSTREAM, sort_by = "Max_Single_Strength")
top30_targets_ra_single_strength <- compute_top30_downstream_targets(RA_DOWNSTREAM, sort_by = "Max_Single_Strength")
top30_targets_hz_single_strength <- compute_top30_downstream_targets(HZ_DOWNSTREAM, sort_by = "Max_Single_Strength")

top30_targets_aging <- top30_targets_aging_total_strength
top30_targets_ra <- top30_targets_ra_total_strength
top30_targets_hz <- top30_targets_hz_total_strength

target_intersection_out <- file.path(STATS_DIR, "Top30_Downstream_Targets_Intersection_mvAge_RA_HZ_UKB.csv")
target_intersection_stats <- write_threeway_intersection_stats(
  top_a = top30_targets_aging,
  top_b = top30_targets_ra,
  top_c = top30_targets_hz,
  key_col = "target_gene",
  label_a = "mvAge",
  label_b = "RA",
  label_c = "HZ UKB",
  out_csv = target_intersection_out
)
log_message(paste0("Three-way downstream target intersection size: ", length(target_intersection_stats$intersection)))
log_message(paste("Saved:", target_intersection_out))

md_out <- file.path(OUTPUT_DIR, "Top30_Intersection_Interpretation_mvAge_RA_HZ_UKB.md")
write_intersection_interpretation_md(
  md_path = md_out,
  upstream_intersection = ligand_intersection_stats$intersection,
  downstream_intersection = target_intersection_stats$intersection,
  upstream_csv = ligand_intersection_out,
  downstream_csv = target_intersection_out
)
log_message(paste("Saved:", md_out))

top30_upstream_supp <- bind_rows(
  top30_ligands_aging_activity %>% mutate(Phenotype = "mvAge", SortBy = "Ligand_Activity"),
  top30_ligands_ra_activity %>% mutate(Phenotype = "RA", SortBy = "Ligand_Activity"),
  top30_ligands_hz_activity %>% mutate(Phenotype = "HZ UKB", SortBy = "Ligand_Activity"),
  top30_ligands_aging_total_activity %>% mutate(Phenotype = "mvAge", SortBy = "Total_Ligand_Activity"),
  top30_ligands_ra_total_activity %>% mutate(Phenotype = "RA", SortBy = "Total_Ligand_Activity"),
  top30_ligands_hz_total_activity %>% mutate(Phenotype = "HZ UKB", SortBy = "Total_Ligand_Activity")
) %>%
  select(Phenotype, SortBy, everything())

top30_downstream_supp <- bind_rows(
  top30_targets_aging_total_strength %>% mutate(Phenotype = "mvAge", SortBy = "Total_Strength"),
  top30_targets_ra_total_strength %>% mutate(Phenotype = "RA", SortBy = "Total_Strength"),
  top30_targets_hz_total_strength %>% mutate(Phenotype = "HZ UKB", SortBy = "Total_Strength"),
  top30_targets_aging_single_strength %>% mutate(Phenotype = "mvAge", SortBy = "Single_Strength"),
  top30_targets_ra_single_strength %>% mutate(Phenotype = "RA", SortBy = "Single_Strength"),
  top30_targets_hz_single_strength %>% mutate(Phenotype = "HZ UKB", SortBy = "Single_Strength")
) %>%
  select(Phenotype, SortBy, everything())

supp_upstream_csv <- file.path(SUPP_DIR, "Top30_Upstream_Ligands_Supplement.csv")
supp_downstream_csv <- file.path(SUPP_DIR, "Top30_Downstream_Targets_Supplement.csv")
write_csv(top30_upstream_supp, supp_upstream_csv)
write_csv(top30_downstream_supp, supp_downstream_csv)
log_message(paste("Saved:", supp_upstream_csv))
log_message(paste("Saved:", supp_downstream_csv))

calc_overlap <- function(df_a, df_b, key_col) {
  a <- unique(na.omit(as.character(df_a[[key_col]])))
  b <- unique(na.omit(as.character(df_b[[key_col]])))
  inter <- intersect(a, b)
  union_all <- union(a, b)
  tibble::tibble(
    Overlap_n = length(inter),
    Jaccard = ifelse(length(union_all) == 0, NA_real_, length(inter) / length(union_all)),
    OnlyA_n = length(setdiff(a, b)),
    OnlyB_n = length(setdiff(b, a))
  )
}

overlap_upstream <- bind_rows(
  calc_overlap(top30_ligands_aging_activity, top30_ligands_aging_total_activity, "ligand") %>% mutate(Phenotype = "mvAge"),
  calc_overlap(top30_ligands_ra_activity, top30_ligands_ra_total_activity, "ligand") %>% mutate(Phenotype = "RA"),
  calc_overlap(top30_ligands_hz_activity, top30_ligands_hz_total_activity, "ligand") %>% mutate(Phenotype = "HZ UKB")
) %>%
  mutate(Category = "Upstream_Ligand", SortA = "Ligand_Activity", SortB = "Total_Ligand_Activity") %>%
  select(Category, Phenotype, SortA, SortB, Overlap_n, Jaccard, OnlyA_n, OnlyB_n)

overlap_downstream <- bind_rows(
  calc_overlap(top30_targets_aging_total_strength, top30_targets_aging_single_strength, "target_gene") %>% mutate(Phenotype = "mvAge"),
  calc_overlap(top30_targets_ra_total_strength, top30_targets_ra_single_strength, "target_gene") %>% mutate(Phenotype = "RA"),
  calc_overlap(top30_targets_hz_total_strength, top30_targets_hz_single_strength, "target_gene") %>% mutate(Phenotype = "HZ UKB")
) %>%
  mutate(Category = "Downstream_Target", SortA = "Total_Strength", SortB = "Single_Strength") %>%
  select(Category, Phenotype, SortA, SortB, Overlap_n, Jaccard, OnlyA_n, OnlyB_n)

sorting_overlap_out <- file.path(SUPP_DIR, "Sorting_Scheme_Overlap_Summary.csv")
write_csv(bind_rows(overlap_upstream, overlap_downstream), sorting_overlap_out)
log_message(paste("Saved:", sorting_overlap_out))

manuscript_md <- file.path(OUTPUT_DIR, "Manuscript_Draft_Methods_Results_Discussion.md")
md_lines <- c(
  "# NicheNet cross-phenotype communication summary (mvAge / RA / HZ UKB)",
  "",
  "## Methods ",
  "",
  "### Study design ",
  "",
  "We integrated NicheNet-based intercellular communication results from three phenotypes (mvAge, Rheumatoid Arthritis [RA], and Herpes Zoster [HZ UKB]) and generated a harmonized 3-row figure. Each row contains: (i) upstream ligand ranking, (ii) a Sankey/alluvial diagram connecting upstream sender → bridge cell → downstream receiver, and (iii) downstream target gene ranking.",
  "",
  "（mvAge、RA、HZ UKB） NicheNet , . : 、Sankey/Alluvial 、. ",
  "",
  "### Input data and preprocessing ",
  "",
  "For each phenotype, we used two CSV files exported from the NicheNet workflow: an upstream table (ligand, sender cell, receiver cell, regulatory potential, ligand activity) and a downstream table (target gene, receiver cell, target regulatory potential). Cell types were mapped to a unified abbreviation set to ensure consistent legends across phenotypes.",
  "",
  " NicheNet  CSV: （、、、regulatory potential、ligand activity）（、、target regulatory potential）. , . ",
  "",
  "### Ranking metrics and plotting ",
  "",
  "**Upstream ligands (upward panel)**: we rank ligands by *Total Ligand Activity*, defined as the sum of ligand activity across unique (ligand, sender, receiver) triplets to avoid double counting through multiple target genes. We visualize per-sender contributions as stacked bars.",
  "",
  "**Upstream ligands (downward panel)**: we rank ligands by *Ligand Activity* and display, for each ligand, a single bar representing its ligand activity value. Bar color corresponds to the dominant sender cell (defined as the sender with the maximum regulatory potential for that ligand).",
  "",
  "**Downstream targets (upward panel)**: we rank target genes by *Total Strength* (the sum of target regulatory potential across receiver cells) and visualize per-receiver contributions as stacked bars.",
  "",
  "**Downstream targets (downward panel)**: we rank target genes by *Single Strength* (the maximum target regulatory potential across receiver cells) and display a single bar for each gene; the bar color corresponds to the receiver cell where this maximum occurs.",
  "",
  ": Ligand Activity ； regulatory potential . : Total Ligand Activity （ ligand-sender-receiver ）. : Single Strength（ target regulatory potential）；. : Total Strength. ",
  "",
  "### Cross-phenotype overlap and sensitivity ",
  "",
  paste0("We quantified the three-way intersection of Top30 upstream ligands and Top30 downstream targets across the three phenotypes, and additionally compared alternative ranking schemes. Overlap statistics (including Jaccard index) are saved in: `", basename(sorting_overlap_out), "`."),
  "",
  " Top30  Top30 , （ Jaccard ）, . ",
  "",
  "## Results ",
  "",
  paste0("### Three-way intersections "),
  paste0("- Upstream ligands intersection size: ", length(ligand_intersection_stats$intersection)),
  paste0("- Downstream targets intersection size: ", length(target_intersection_stats$intersection)),
  "",
  "The intersection suggests shared communication programs that recur across immune mvAge, chronic autoimmunity, and viral reactivation contexts. A pathway-level interpretation is provided in `Top30_Intersection_Interpretation_mvAge_RA_HZ_UKB.md`.",
  "",
  "、. . ",
  "",
  "### Supplementary Top30 tables Top30",
  "",
  paste0("- `", basename(supp_upstream_csv), "`: Top30 upstream ligands for each phenotype under two sorting schemes."),
  paste0("- `", basename(supp_downstream_csv), "`: Top30 downstream targets for each phenotype under two sorting schemes."),
  "",
  "（Supplementary Tables）, . ",
  "",
  "## Discussion ",
  "",
  "### Biological interpretation ",
  "",
  "Across phenotypes, shared upstream ligands and downstream targets likely reflect convergent inflammatory and immune-activation modules. In RA, persistent cytokine/chemokine signaling and myeloid–lymphoid crosstalk are expected; in HZ UKB, antiviral and interferon-linked programs can dominate; in mvAge, baseline inflammaging may shift both upstream ligand usage and downstream transcriptional responses. The observed overlaps provide a concise shortlist of candidates for mechanistic follow-up.",
  "",
  ": RA -；HZ UKB /；mvAge  inflammaging . . ",
  "",
  "### On the ranking strategy ",
  "",
  "Using *Single Strength* for the downstream upward panel prioritizes genes with a strong, cell-type-specific regulatory potential peak, which is useful for highlighting sharp, potentially actionable signals. In contrast, *Total Strength* emphasizes broad, multi-receiver contributions and captures signal breadth. Reporting both schemes (and their overlap) improves transparency and mitigates over-reliance on a single ranking definition.",
  "",
  " Single Strength “”, ； Total Strength , . . ",
  "",
  "### Limitations and next steps ",
  "",
  "These results depend on the upstream/downstream tables exported from the NicheNet workflow; thus, they inherit assumptions about ligand–target priors and the definition of ligand activity. Future work can validate key ligand–receiver–target triplets using independent datasets, perturbation evidence, and/or colocalization/MR-informed prioritization to connect genetic evidence with communication circuits.",
  "",
  " NicheNet . 、/MR  ligand–receiver–target , . "
)
writeLines(md_lines, manuscript_md, useBytes = TRUE)
log_message(paste("Saved:", manuscript_md))

# 5. Plot Generation Function / 
# ------------------------------------------------------------------------------
generate_row_plot <- function(upstream_path, downstream_path, row_title, show_legend_p0 = TRUE, show_legend_p1 = TRUE, show_legend_p2 = TRUE) {
  log_message(paste("----------------------------------------------------------------"))
  log_message(paste("Processing data for:", row_title))
  
  if (!file_exists(upstream_path) || !file_exists(downstream_path)) {
    log_message(paste("Error: Input files not found for", row_title))
    return(NULL)
  }
  
  # Read Data
  upstream_data <- read_csv(upstream_path, show_col_types = FALSE)
  downstream_data <- read_csv(downstream_path, show_col_types = FALSE)
  
  # ----------------------------------------------------------------------------
  # Data Processing (Copied from Source)
  # ----------------------------------------------------------------------------
  
  # Extract Upstream Columns
  upstream_clean <- upstream_data %>%
    select(sender_cell, target_gene, receiver_cell) %>%
    distinct() %>%
    rename(
      U_Sender = sender_cell,
      Bridge_Gene = target_gene,
      Bridge_Cell = receiver_cell
    ) %>%
    mutate(Source = "Upstream")

  # Extract Downstream Columns
  downstream_clean <- downstream_data %>%
    select(sender_cell, ligand, receiver_cell) %>%
    distinct() %>%
    rename(
      Bridge_Cell = sender_cell,
      Bridge_Gene = ligand,
      D_Receiver = receiver_cell
    ) %>%
    mutate(Source = "Downstream")

  # Extract Downstream Target Genes for Bar Plot
  target_gene_long <- downstream_data %>%
    select(target_gene, receiver_cell, target_regulatory_potential) %>%
    mutate(
      target_regulatory_potential = suppressWarnings(as.numeric(target_regulatory_potential)),
      D_Receiver_Abbr = abbreviate_cell_types(receiver_cell)
    )

  upstream_ligand_long <- upstream_data %>%
    select(ligand, sender_cell, regulatory_potential, ligand_activity_pearson) %>%
    mutate(
      regulatory_potential = suppressWarnings(as.numeric(regulatory_potential)),
      ligand_activity_pearson = suppressWarnings(as.numeric(ligand_activity_pearson)),
      U_Sender_Abbr = abbreviate_cell_types(sender_cell)
    )

  gene_summary <- target_gene_long %>%
    group_by(target_gene) %>%
    summarise(
      Total_Count = n(),
      Total_Strength = sum(target_regulatory_potential, na.rm = TRUE),
      Mean_Strength = mean(target_regulatory_potential, na.rm = TRUE),
      Max_Single_Strength = max(target_regulatory_potential, na.rm = TRUE),
      .groups = "drop"
    )

  # Fix: Ensure Ligand Activity is not over-counted by Target Genes
  # We first aggregate unique (Ligand, Sender, Receiver) triplets
  upstream_ligand_unique <- upstream_data %>%
    select(ligand, sender_cell, receiver_cell, regulatory_potential, ligand_activity_pearson) %>%
    distinct() %>%
    mutate(
      regulatory_potential = suppressWarnings(as.numeric(regulatory_potential)),
      ligand_activity_pearson = suppressWarnings(as.numeric(ligand_activity_pearson)),
      U_Sender_Abbr = abbreviate_cell_types(sender_cell)
    )

  ligand_strength_summary <- upstream_ligand_unique %>%
    group_by(ligand) %>%
    summarise(
      Total_Count = n(),
      Total_Strength = sum(regulatory_potential, na.rm = TRUE),
      Mean_Strength = mean(regulatory_potential, na.rm = TRUE),
      Total_Ligand_Activity = sum(ligand_activity_pearson, na.rm = TRUE),
      .groups = "drop"
    )

  ligand_activity_summary <- upstream_data %>%
    select(ligand, ligand_activity_pearson) %>%
    distinct() %>%
    mutate(ligand_activity_pearson = suppressWarnings(as.numeric(ligand_activity_pearson))) %>%
    group_by(ligand) %>%
    summarise(
      Ligand_Activity = max(ligand_activity_pearson, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(Ligand_Activity = ifelse(is.infinite(Ligand_Activity), NA_real_, Ligand_Activity))

  ligand_summary <- ligand_strength_summary %>%
    left_join(ligand_activity_summary, by = "ligand")

  # Merge Data (Full Join)
  chain_data <- full_join(
    upstream_clean,
    downstream_clean,
    by = c("Bridge_Cell", "Bridge_Gene"),
    relationship = "many-to-many"
  )

  # Apply Abbreviations and Handle NAs
  chain_data_viz <- chain_data %>%
    mutate(
      U_Sender_Abbr = abbreviate_cell_types(U_Sender),
      Bridge_Cell_Abbr = abbreviate_cell_types(Bridge_Cell),
      D_Receiver_Abbr = abbreviate_cell_types(D_Receiver),
      
      # Mark if connection exists
      Has_Upstream = !is.na(U_Sender),
      Has_Downstream = !is.na(D_Receiver)
    )

  # ----------------------------------------------------------------------------
  # Stats Calculation (Copied from Source)
  # ----------------------------------------------------------------------------
  
  # Upstream Stats
  upstream_link_stats <- upstream_data %>%
    group_by(sender_cell, receiver_cell) %>%
    summarise(Count = n(), .groups = "drop") %>%
    mutate(
      Total = sum(Count),
      Prop_U_Link = Count / Total
    ) %>%
    rename(U_Sender = sender_cell, Bridge_Cell = receiver_cell) %>%
    mutate(
      U_Sender_Abbr = abbreviate_cell_types(U_Sender),
      Bridge_Cell_Abbr = abbreviate_cell_types(Bridge_Cell)
    )

  # Downstream Stats
  downstream_sender_stats <- downstream_data %>%
    group_by(sender_cell) %>%
    summarise(Count = n(), .groups = "drop") %>%
    mutate(
      Total = sum(Count),
      Prop_D_Sender = Count / Total,
      Bridge_Cell_Abbr = abbreviate_cell_types(sender_cell)
    )

  downstream_link_stats <- downstream_data %>%
    group_by(sender_cell, receiver_cell) %>%
    summarise(Count = n(), .groups = "drop") %>%
    mutate(
      Total = sum(Count),
      Prop_D_Link = Count / Total
    ) %>%
    rename(Bridge_Cell = sender_cell, D_Receiver = receiver_cell) %>%
    mutate(
      Bridge_Cell_Abbr = abbreviate_cell_types(Bridge_Cell),
      D_Receiver_Abbr = abbreviate_cell_types(D_Receiver)
    )

  # ----------------------------------------------------------------------------
  # Visualization Construction (Copied from Source)
  # ----------------------------------------------------------------------------

  # Map Title to Abbreviation
  display_title <- if (row_title == "Rheumatoid Arthritis") "RA" else if (row_title == "HZ UKB") "HZ UKB" else row_title
  
  # --- Plot 1: Alluvial Plot (Sankey) ---
  
  # Prepare Data for Weight Assignment
  plot_data_prep <- chain_data_viz %>%
    filter(Has_Upstream | Has_Downstream) %>%
    mutate(id = row_number())

  # Join Proportions
  plot_data_prep <- plot_data_prep %>%
    left_join(upstream_link_stats %>% select(U_Sender_Abbr, Bridge_Cell_Abbr, Prop_U_Link), 
              by = c("U_Sender_Abbr", "Bridge_Cell_Abbr")) %>%
    left_join(downstream_link_stats %>% select(Bridge_Cell_Abbr, D_Receiver_Abbr, Prop_D_Link),
              by = c("Bridge_Cell_Abbr", "D_Receiver_Abbr"))

  # Calculate number of chains per link
  u_link_counts <- plot_data_prep %>%
    group_by(U_Sender_Abbr, Bridge_Cell_Abbr) %>%
    summarise(Chains_in_U_Link = n(), .groups = "drop")
  
  d_link_counts <- plot_data_prep %>%
    group_by(Bridge_Cell_Abbr, D_Receiver_Abbr) %>%
    summarise(Chains_in_D_Link = n(), .groups = "drop")
  
  plot_data_prep <- plot_data_prep %>%
    left_join(u_link_counts, by = c("U_Sender_Abbr", "Bridge_Cell_Abbr")) %>%
    left_join(d_link_counts, by = c("Bridge_Cell_Abbr", "D_Receiver_Abbr")) %>%
    mutate(
      # Apply fix for NA values to ensure robustness
      Prop_U_Link = replace_na(Prop_U_Link, 0),
      Prop_D_Link = replace_na(Prop_D_Link, 0),
      Chains_in_U_Link = replace_na(Chains_in_U_Link, 1),
      Chains_in_D_Link = replace_na(Chains_in_D_Link, 1),

      W1 = ifelse(Has_Upstream, Prop_U_Link / Chains_in_U_Link, 0),
      W2 = ifelse(Has_Downstream, Prop_D_Link / Chains_in_D_Link, 0),
      W3 = W2
    )

  # Create Lodes
  lodes_custom <- bind_rows(
    plot_data_prep %>% 
      select(id, U_Sender_Abbr, W1) %>% 
      rename(stratum = U_Sender_Abbr, y = W1) %>% 
      mutate(x = "Upstream"),
      
    plot_data_prep %>% 
      select(id, Bridge_Cell_Abbr, W2) %>% 
      rename(stratum = Bridge_Cell_Abbr, y = W2) %>% 
      mutate(x = "Bridge"),
      
    plot_data_prep %>% 
      select(id, D_Receiver_Abbr, W3) %>% 
      rename(stratum = D_Receiver_Abbr, y = W3) %>% 
      mutate(x = "Downstream")
  ) %>%
    mutate(
      x = factor(x, levels = c("Upstream", "Bridge", "Downstream"))
    ) %>%
    filter(!is.na(stratum), stratum != "None") %>%
    filter(y > 0)
    
  # Define Factor Levels
  bridge_levels <- downstream_sender_stats %>%
    arrange(desc(Prop_D_Sender)) %>%
    pull(Bridge_Cell_Abbr)
    
  lodes_custom <- lodes_custom %>%
    mutate(
      stratum = factor(stratum, levels = unique(c(levels(factor(lodes_custom$stratum)), bridge_levels)))
    )

  p1 <- ggplot(lodes_custom,
               aes(x = x, stratum = stratum, alluvium = id, y = y)) +
    geom_flow(aes(fill = stratum), width = 1/4, alpha = 0.6, curve_type = "cubic", na.rm = FALSE) +
    geom_stratum(aes(fill = stratum), width = 1/4, color = "grey30", alpha = 1, na.rm = FALSE) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 3, check_overlap = TRUE, na.rm = FALSE) +
    scale_x_discrete(limits = c("Upstream", "Bridge", "Downstream"), expand = c(.1, .1), labels = c("Upstream Sender", "Bridge Cell", "Downstream Receiver")) +
    scale_fill_manual(values = cell_colors, name = "Cell Type") +
    scale_y_continuous(labels = NULL, breaks = NULL) +
    theme_minimal(base_size = 14) +
    labs(
      title = NULL,
      subtitle = NULL,
      x = NULL,
      y = NULL
    ) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid = element_blank(),
      axis.text.x = element_text(face = "bold", size = 14),
      plot.title = element_text(face = "bold", size = 16)
    )
  
  # Save Sankey
  p1_filename <- paste0(make_clean_names(row_title), "_Sankey.pdf")
  ggsave(file.path(OUTPUT_DIR, p1_filename), p1, width = 10, height = 6)
  
  if (show_legend_p1) {
    p1 <- p1 + theme(legend.position = "bottom")
  } else {
    p1 <- p1 + theme(legend.position = "none")
  }

  # --- Plot 2: Downstream Target Bar Plot (Top 30 Bidirectional Independent) ---
  
  # 2.1. Top 30 UP (Single Strength - Single Bar)
  # Sort by Max Single Strength
  top_genes_up <- gene_summary %>%
    arrange(desc(Max_Single_Strength), desc(Total_Strength), desc(Total_Count), target_gene) %>%
    slice_head(n = 30)
    
  # Identify the Receiver Cell with the Max Strength for each Top Gene
  # If ties, take the first one (or based on count)
  gene_max_receiver <- target_gene_long %>%
    inner_join(top_genes_up %>% select(target_gene, Max_Single_Strength), by = "target_gene") %>%
    filter(target_regulatory_potential == Max_Single_Strength) %>%
    group_by(target_gene) %>%
    slice(1) %>% # Take one if ties
    select(target_gene, D_Receiver_Abbr)
    
  plot_data_bar_up <- top_genes_up %>%
    left_join(gene_max_receiver, by = "target_gene") %>%
    mutate(
      Target_Gene_Order = factor(target_gene, levels = top_genes_up$target_gene),
      y_value = Max_Single_Strength
    )

  label_bar_up <- top_genes_up %>%
    transmute(
      Target_Gene_Order = factor(target_gene, levels = top_genes_up$target_gene),
      label = target_gene,
      y_label = Max_Single_Strength
    )
  y_label_offset_up <- max(top_genes_up$Max_Single_Strength, na.rm = TRUE) * 0.05

  p2_up <- ggplot(plot_data_bar_up, aes(x = Target_Gene_Order, y = y_value, fill = D_Receiver_Abbr)) +
    geom_col(position = "identity", width = 0.7) +
    geom_text(
      data = label_bar_up,
      aes(x = Target_Gene_Order, y = y_label + y_label_offset_up, label = label),
      inherit.aes = FALSE,
      size = 2.5, fontface = "italic", angle = 90, hjust = 1, vjust = 0.5
    ) +
    scale_fill_manual(values = cell_colors) +
    scale_y_reverse(expand = expansion(mult = c(0, 0.3))) + 
    coord_cartesian(clip = "off") +
    labs(y = "Single Strength", x = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "none",
      plot.margin = margin(t = 5, b = 60),
      plot.title = element_blank(),
      plot.subtitle = element_blank()
    )

  # 2.2. Top 30 DOWN (Total Strength - Stacked)
  top_genes_down <- gene_summary %>%
    arrange(desc(Total_Strength), desc(Mean_Strength), desc(Total_Count), target_gene) %>%
    slice_head(n = 30)
    
  plot_data_bar_down <- target_gene_long %>%
    inner_join(top_genes_down, by = "target_gene") %>%
    filter(D_Receiver_Abbr %in% names(cell_colors)) %>%
    group_by(target_gene, D_Receiver_Abbr) %>%
    summarise(
      Cell_Count = n(),
      Cell_Total_Strength = sum(target_regulatory_potential, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Target_Gene_Order = factor(target_gene, levels = top_genes_down$target_gene),
      y_value = Cell_Total_Strength
    )

  label_bar_down <- top_genes_down %>%
    transmute(
      Target_Gene_Order = factor(target_gene, levels = top_genes_down$target_gene),
      label = target_gene,
      y_label = Total_Strength
    )
  y_label_offset_down <- max(top_genes_down$Total_Strength, na.rm = TRUE) * 0.1 # Increased offset

  p2_down <- ggplot(plot_data_bar_down, aes(x = Target_Gene_Order, y = y_value, fill = D_Receiver_Abbr)) +
    geom_col(position = "stack", width = 0.7) +
    geom_text(
      data = label_bar_down,
      aes(x = Target_Gene_Order, y = y_label + y_label_offset_down, label = label),
      inherit.aes = FALSE,
      size = 2.5, fontface = "italic", angle = 90, hjust = 0, vjust = 0.5
    ) +
    scale_fill_manual(values = cell_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.3))) + 
    coord_cartesian(clip = "off") + # Prevent clipping
    labs(y = "Total Strength", x = NULL, title = paste0("Top 30 Downstream Targets (", display_title, ")")) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "none",
      plot.margin = margin(t = 40, b = 2),
      plot.title = element_text(face = "bold", size = 14)
    )

  # Combine P2
  p2 <- p2_down / p2_up + plot_layout(heights = c(1, 1))

  # Save P2
  p2_filename <- paste0(make_clean_names(row_title), "_Downstream_Bar.pdf")
  ggsave(file.path(OUTPUT_DIR, p2_filename), p2, width = 8, height = 8)

  if (show_legend_p2) {
    p2 <- p2 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  } else {
    p2 <- p2 & theme(legend.position = "none")
  }

  # --- Plot 0: Upstream Ligand Bar Plot (Top 30 Bidirectional Independent) ---
  
  # 0.1. Top 30 UP (Single Ligand Activity - Single Bar)
  # Sort by Ligand Activity
  top_ligands_up <- ligand_summary %>%
    arrange(desc(Ligand_Activity), desc(Total_Strength), desc(Total_Count), ligand) %>%
    slice_head(n = 30)
  
  # Identify Dominant Sender for Color (Max Count or Max Strength)
  # We use Max Strength from a sender to determine color
  ligand_dominant_sender <- upstream_ligand_unique %>%
    group_by(ligand) %>%
    arrange(desc(regulatory_potential)) %>% # Use potential to pick dominant sender
    slice(1) %>%
    select(ligand, U_Sender_Abbr)
    
  plot_data_ligand_up <- top_ligands_up %>%
    left_join(ligand_dominant_sender, by = "ligand") %>%
    mutate(
      Ligand_Order = factor(ligand, levels = top_ligands_up$ligand),
      y_value = Ligand_Activity
    )

  label_ligand_up <- top_ligands_up %>%
    transmute(
      Ligand_Order = factor(ligand, levels = top_ligands_up$ligand),
      label = ligand,
      y_label = Ligand_Activity
    )
  y_label_offset_ligand_up <- max(top_ligands_up$Ligand_Activity, na.rm = TRUE) * 0.05

  p0_up <- ggplot(plot_data_ligand_up, aes(x = Ligand_Order, y = y_value, fill = U_Sender_Abbr)) +
    geom_col(position = "identity", width = 0.7) +
    geom_text(
      data = label_ligand_up,
      aes(x = Ligand_Order, y = y_label + y_label_offset_ligand_up, label = label),
      inherit.aes = FALSE,
      size = 2.5, fontface = "italic", angle = 90, hjust = 1, vjust = 0.5
    ) +
    scale_fill_manual(values = cell_colors) +
    scale_y_reverse(expand = expansion(mult = c(0, 0.3))) + 
    coord_cartesian(clip = "off") +
    labs(y = "Ligand Activity", x = NULL) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "none",
      plot.margin = margin(t = 5, b = 60),
      plot.title = element_blank(),
      plot.subtitle = element_blank()
    )

  # 0.2. Top 30 DOWN (Total Ligand Activity - Stacked)
  # Sort by Total Ligand Activity
  top_ligands_down <- ligand_summary %>%
    arrange(desc(Total_Ligand_Activity), desc(Ligand_Activity), desc(Total_Count), ligand) %>%
    slice_head(n = 30)
    
  plot_data_ligand_down <- upstream_ligand_unique %>%
    inner_join(top_ligands_down, by = "ligand") %>%
    filter(U_Sender_Abbr %in% names(cell_colors)) %>%
    group_by(ligand, U_Sender_Abbr) %>%
    summarise(
      Cell_Count = n(),
      Cell_Total_Activity = sum(ligand_activity_pearson, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(
      Ligand_Order = factor(ligand, levels = top_ligands_down$ligand),
      y_value = Cell_Total_Activity
    )

  label_ligand_down <- top_ligands_down %>%
    transmute(
      Ligand_Order = factor(ligand, levels = top_ligands_down$ligand),
      label = ligand,
      y_label = Total_Ligand_Activity
    )
  y_label_offset_ligand_down <- max(top_ligands_down$Total_Ligand_Activity, na.rm = TRUE) * 0.1

  p0_down <- ggplot(plot_data_ligand_down, aes(x = Ligand_Order, y = y_value, fill = U_Sender_Abbr)) +
    geom_col(position = "stack", width = 0.7) +
    geom_text(
      data = label_ligand_down,
      aes(x = Ligand_Order, y = y_label + y_label_offset_ligand_down, label = label),
      inherit.aes = FALSE,
      size = 2.5, fontface = "italic", angle = 90, hjust = 0, vjust = 0.5
    ) +
    scale_fill_manual(values = cell_colors) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.3))) + 
    coord_cartesian(clip = "off") +
    labs(y = "Total Ligand Activity", x = NULL, title = paste0("Top 30 Upstream Ligands (", display_title, ")")) +
    theme_minimal(base_size = 12) +
    theme(
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position = "none",
      plot.margin = margin(t = 40, b = 2),
      plot.title = element_text(face = "bold", size = 14)
    )

  # Combine P0
  p0 <- p0_down / p0_up + plot_layout(heights = c(1, 1))

  # Save P0
  p0_filename <- paste0(make_clean_names(row_title), "_Upstream_Bar.pdf")
  ggsave(file.path(OUTPUT_DIR, p0_filename), p0, width = 8, height = 8)

  if (show_legend_p0) {
    p0 <- p0 + plot_layout(guides = "collect") & theme(legend.position = "bottom")
  } else {
    p0 <- p0 & theme(legend.position = "none")
  }

  return(list(p0, p1, p2))
}

# 6. Generate Plots / 
# ------------------------------------------------------------------------------
# Aging (Hide Legend All)
plots_aging <- generate_row_plot(AGING_UPSTREAM, AGING_DOWNSTREAM, "mvAge", 
                                 show_legend_p0 = FALSE, show_legend_p1 = FALSE, show_legend_p2 = FALSE)

# RA (Hide Legend All)
plots_ra <- generate_row_plot(RA_UPSTREAM, RA_DOWNSTREAM, "Rheumatoid Arthritis", 
                              show_legend_p0 = FALSE, show_legend_p1 = FALSE, show_legend_p2 = FALSE)

# HZ (Show Legend Middle Only)
plots_hz <- generate_row_plot(HZ_UPSTREAM, HZ_DOWNSTREAM, "HZ UKB", 
                              show_legend_p0 = FALSE, show_legend_p1 = TRUE, show_legend_p2 = FALSE)

if (is.null(plots_aging) || is.null(plots_ra) || is.null(plots_hz)) {
  stop("One or more groups failed to generate plots. Check input paths.")
}

# 7. Compose Final Figure / 
# ------------------------------------------------------------------------------
log_message("Composing final multi-panel figure...")

# Layout: 3 Rows. 
# Row 1 (Aging): P0 + P1 + P2
# Row 2 (RA):    P0 + P1 + P2
# Row 3 (HZ):    P0 + P1 + P2
# Note: Legends for Row 1 and 2 are hidden in generate_row_plot. 
# Row 3 has legend.
# We do NOT use guides="collect" because we want the legend ONLY at the bottom (from Row 3) 
# and we've already manually controlled visibility.

# Create rows
row1 <- wrap_elements(plots_aging[[1]]) + plots_aging[[2]] + wrap_elements(plots_aging[[3]]) + plot_layout(widths = c(1, 2.5, 1))
row2 <- wrap_elements(plots_ra[[1]]) + plots_ra[[2]] + wrap_elements(plots_ra[[3]]) + plot_layout(widths = c(1, 2.5, 1))
row3 <- wrap_elements(plots_hz[[1]]) + plots_hz[[2]] + wrap_elements(plots_hz[[3]]) + plot_layout(widths = c(1, 2.5, 1))

# Combine Rows with wrap_elements to group them into single figures
combined_plot <- (wrap_elements(row1) / wrap_elements(row2) / wrap_elements(row3))

# 8. Save Output / 
# ------------------------------------------------------------------------------
pdf_file <- file.path(OUTPUT_DIR, "Combined_Phenotypes_Communication.pdf")
png_file <- file.path(OUTPUT_DIR, "Combined_Phenotypes_Communication.png")

plot_height <- 24 
plot_width <- 24

log_message(paste("Saving PDF to:", pdf_file))
ggsave(pdf_file, plot = combined_plot, width = plot_width, height = plot_height, limitsize = FALSE)

log_message(paste("Saving PNG to:", png_file))
ggsave(png_file, plot = combined_plot, width = plot_width, height = plot_height, dpi = 300, limitsize = FALSE)

log_message("Visualization complete.")
