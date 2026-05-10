#!/usr/bin/env Rscript
# ==============================================================================
# [Script]: 6.1.3.1_MR_Visualization_Compare_template.R
# [Method]: cis-MR Analysis (Two-Sample MR)
# [Step]: Visualization and intersection comparison
# 
# [Function]:
# Read and parse MR results for HZ and RA, extract significant results, find intersection, and plot.
# 
# [Input]:
#   --file_hz  : Path to HZ MR results CSV
#   --file_ra  : Path to RA MR results CSV
#   --out_dir  : Output directory path
# 
# [Output]:
# Volcano plots, forest plots, heatmaps, and summary statistics.
# ==============================================================================

rm(list = ls())
gc()

# 1. Load Dependencies
suppressPackageStartupMessages({
  library(optparse)
  library(tidyverse)
  library(ggplot2)
  library(reshape2)
  library(circlize)
  library(ComplexHeatmap)
  library(ggsci)
  library(showtext)
  library(ggrepel)
  library(magick)
  library(cowplot)
})

# 2. Parse Command Line Arguments
option_list <- list(
  make_option(c("--file_hz"), type="character", default=NULL, help="Path to HZ MR results CSV"),
  make_option(c("--file_ra"), type="character", default=NULL, help="Path to RA MR results CSV"),
  make_option(c("--out_dir"), type="character", default="./MR_Visualization", help="Output directory path")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$file_hz) || is.null(opt$file_ra)) {
  print_help(opt_parser)
  stop("Missing required input files.")
}

# Setup Paths
file_hz <- opt$file_hz
file_ra <- opt$file_ra
work_dir <- opt$out_dir

font_add("STSong", "/System/Library/Fonts/Supplemental/Songti.ttc")
showtext_auto()

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- file.path(work_dir, paste0("Visualization_Result_", timestamp))
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)


log_dir <- file.path(output_dir, "logs")
plots_dir <- file.path(output_dir, "plots")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(plots_dir, recursive = TRUE, showWarnings = FALSE)

log_file <- file.path(log_dir, paste0("run_log_", timestamp, ".txt"))
sink(log_file, split = TRUE)
start_time <- Sys.time()
cat("===================================================\n")
cat("MR Visualization Script Started at:", format(start_time, "%Y-%m-%d %H:%M:%S"), "\n")
cat("===================================================\n\n")

# 2. MR
file_hz <- "./"
file_ra <- "./"

cat("Reading input files...\n")
df_hz <- read.csv(file_hz, stringsAsFactors = FALSE)
df_ra <- read.csv(file_ra, stringsAsFactors = FALSE)

df_hz <- df_hz %>% mutate(exposure = str_replace(exposure, " protein levels$", ""))
df_ra <- df_ra %>% mutate(exposure = str_replace(exposure, " protein levels$", ""))

# 3. 
get_primary <- function(df) {
  df %>% filter(method %in% c("Inverse variance weighted", "Wald ratio")) %>%
    mutate(method_priority = ifelse(method == "Inverse variance weighted", 1, 2)) %>%
    group_by(exposure) %>% arrange(method_priority, pval) %>% slice(1) %>% ungroup() %>%
    select(-method_priority)
}

df_hz_primary <- get_primary(df_hz)
df_ra_primary <- get_primary(df_ra)

sig_hz <- df_hz_primary %>% filter(fdr <= 0.05) %>% pull(exposure)
sig_ra <- df_ra_primary %>% filter(fdr <= 0.05) %>% pull(exposure)
intersect_exp <- intersect(sig_hz, sig_ra)

cat("Significant exposures in HZ:", length(sig_hz), "\n")
cat("Significant exposures in RA:", length(sig_ra), "\n")
cat("Intersection exposures:", length(intersect_exp), "\n")

write.csv(data.frame(exposure = intersect_exp), file.path(output_dir, paste0("Intersection_Exposures_", timestamp, ".csv")), row.names = FALSE)

# 4.  (NPG : = #E64B35, = #4DBBD5)

plot_volcano <- function(df_primary, intersect_exp, title_str, out_prefix) {
  cat("Generating Volcano plot for", title_str, "...\n")
  df_primary <- df_primary %>%
    mutate(
      is_intersect = exposure %in% intersect_exp,
      sig = case_when(
        fdr <= 0.05 & or > 1 ~ "Significant Risk",
        fdr <= 0.05 & or < 1 ~ "Significant Protective",
        TRUE ~ "Not Significant"
      )
    )
  
  df_primary$sig <- factor(df_primary$sig, levels = c("Significant Risk", "Significant Protective", "Not Significant"))
  
  palette_sig <- c(
    "Significant Risk" = "#E64B35",
    "Significant Protective" = "#4DBBD5",
    "Not Significant" = "#D9D9D9"
  )
  
  p <- ggplot(df_primary, aes(x = log2(or), y = -log10(pval))) +
    geom_point(aes(fill = sig), shape = 21, color = "black", stroke = 0.5, size = 4.5, alpha = 0.85) +
    scale_fill_manual(values = palette_sig) +
    geom_hline(yintercept = -log10(max(df_primary$pval[df_primary$fdr <= 0.05], na.rm=TRUE)), linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "black", linewidth = 0.8) +
    geom_text_repel(
      data = df_primary %>% filter(is_intersect),
      aes(label = exposure),
      size = 4.5, max.overlaps = 100, box.padding = 1.0, fontface = "bold.italic", color = "black"
    ) +
    theme_classic(base_size = 14) +
    labs(
      x = "Log2(Odds Ratio)",
      y = "-Log10(P-value)",
      fill = "Significance"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 18),
      axis.title = element_text(size = 16, face = "bold"),
      axis.text = element_text(size = 14, color = "black"),
      legend.position = c(0.02, 0.98),
      legend.justification = c(0, 1),
      legend.title = element_text(face = "bold", size = 14),
      legend.text = element_text(size = 12),
      legend.background = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 1.0)
    )
  
  showtext_opts(dpi = 96)
  ggsave(file.path(plots_dir, paste0(out_prefix, ".pdf")), plot = p, width = 11, height = 10)
  showtext_opts(dpi = 300)
  ggsave(file.path(plots_dir, paste0(out_prefix, ".png")), plot = p, width = 11, height = 10, dpi = 300, bg="white")
  showtext_opts(dpi = 96)
  return(p)
}

plot_intersection_forest <- function(df_hz_primary, df_ra_primary, intersect_exp, out_prefix) {
  cat("Generating Dual-Outcome Intersection Forest plot...\n")
  
  df_hz_int <- df_hz_primary %>% filter(exposure %in% intersect_exp) %>% mutate(Outcome = "HZ (FinnGen)")
  df_ra_int <- df_ra_primary %>% filter(exposure %in% intersect_exp) %>% mutate(Outcome = "RA")
  
  df_int <- bind_rows(df_ra_int, df_hz_int)
  
  exp_order <- df_ra_int %>% arrange(desc(or)) %>% pull(exposure)
  df_int$exposure <- factor(df_int$exposure, levels = exp_order)
  
  # Set Outcome levels so RA is displayed above HZ
  # ggplot2 places earlier levels at the top in position_dodge when y is continuous-like, 
  # or rather: the order in factor levels determines dodge order.
  # For position_dodge with y axis, reversed factor levels often puts RA above HZ.
  df_int$Outcome <- factor(df_int$Outcome, levels = c("HZ (FinnGen)", "RA"))
  
  #  CSV 
  df_int_save <- df_int %>% 
    select(Outcome, exposure, method, nsnp, or, or_lci95, or_uci95, pval, fdr) %>%
    arrange(exposure, desc(Outcome))
  write_csv(df_int_save, file.path(output_dir, paste0("Intersection_Proteins_Data_", timestamp, ".csv")))
  
  # : xmax  xmin , 1
  max_x_limit <- 2.0
  min_x_limit <- 0.0
  
  df_int_plot <- df_int %>%
    mutate(
      or_plot = ifelse(or > max_x_limit, max_x_limit, ifelse(or < min_x_limit, min_x_limit, or)),
      or_lci95_plot = ifelse(or_lci95 < min_x_limit, min_x_limit, or_lci95),
      or_uci95_plot = ifelse(or_uci95 > max_x_limit, max_x_limit, or_uci95),
      is_truncated = or_uci95 > max_x_limit | or_lci95 < min_x_limit
    )
  
  max_dist <- max(abs(df_int_plot$or_uci95_plot - 1), abs(1 - df_int_plot$or_lci95_plot), na.rm = TRUE)
  axis_limit <- max_dist * 1.05
  
  num_exp <- length(levels(df_int_plot$exposure))
  bg_data <- data.frame(
    y = seq(1, num_exp, by = 2)
  )
  
  p <- ggplot(df_int_plot, aes(x = or_plot, y = exposure, color = Outcome)) +
    geom_rect(data = bg_data, aes(xmin = -Inf, xmax = Inf, ymin = y - 0.5, ymax = y + 0.5),
              fill = "grey92", alpha = 0.6, inherit.aes = FALSE) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray30", linewidth = 1.2) +
    geom_errorbar(aes(xmin = or_lci95_plot, xmax = or_uci95_plot), width = 0.4, linewidth = 1.2, position = position_dodge(width = 0.75)) +
    geom_segment(data = df_int_plot %>% filter(or_uci95 > max_x_limit), 
                 aes(x = or_plot, xend = max_x_limit, y = exposure, yend = exposure, color = Outcome),
                 arrow = arrow(length = unit(0.2, "cm"), type = "closed"), 
                 position = position_dodge(width = 0.75), linewidth = 1.2) +
    geom_point(aes(shape = Outcome, fill = Outcome), size = 6.0, position = position_dodge(width = 0.75), stroke = 1.0, color = "black") +
    scale_fill_manual(values = c("RA" = "#A65628", "HZ (FinnGen)" = "#4DAF4A"), breaks = c("RA", "HZ (FinnGen)")) +
    scale_color_manual(values = c("RA" = "#A65628", "HZ (FinnGen)" = "#4DAF4A"), breaks = c("RA", "HZ (FinnGen)")) +
    scale_shape_manual(values = c("RA" = 21, "HZ (FinnGen)" = 24), breaks = c("RA", "HZ (FinnGen)")) +
    coord_cartesian(xlim = c(1 - axis_limit, 1 + axis_limit)) +
    theme_classic(base_size = 24) +
    labs(
      x = "Odds Ratio (95% CI)",
      y = "Exposure (Protein)",
      fill = "Outcome",
      shape = "Outcome"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 28),
      axis.title = element_text(size = 26, face = "bold"),
      axis.text.y = element_text(size = 22, face = "bold.italic", color = "black"),
      axis.text.x = element_text(size = 22, color = "black"),
      legend.position = "right",
      legend.direction = "vertical",
      legend.title = element_text(face = "bold", size = 24, angle = -90, hjust = 0.5),
      legend.text = element_text(size = 22, angle = -90, hjust = 0.5),
      panel.border = element_blank(),
      panel.background = element_blank(),
      axis.line = element_line(colour = "black", linewidth = 1.0),
      legend.key.size = unit(1.8, "cm")
    ) +
    guides(color = "none") 
  
  fp_height <- max(12, length(intersect_exp) * 0.8)
  showtext_opts(dpi = 96)
  ggsave(file.path(plots_dir, paste0(out_prefix, ".pdf")), plot = p, width = 15, height = fp_height, limitsize = FALSE)
  showtext_opts(dpi = 300)
  ggsave(file.path(plots_dir, paste0(out_prefix, ".png")), plot = p, width = 15, height = fp_height, dpi = 300, bg="white", limitsize = FALSE)
  showtext_opts(dpi = 96)
}

plot_heatmap <- function(df_all, sig_exps, title_str, out_prefix) {
  cat("Generating Circular Heatmap for", title_str, "...\n")
  df_heat <- df_all %>% 
    filter(exposure %in% sig_exps) %>%
    filter(!method %in% c("Wald ratio")) %>%
    filter(nsnp > 2)
  
  if(nrow(df_heat) == 0) return(NULL)
  
  data_processed <- df_heat %>%
    mutate(log2_or = log2(or), distance_from_1 = abs(log2_or), TRAIT = exposure)
  
  ivw_data <- data_processed %>% filter(method == "Inverse variance weighted") %>%
    mutate(sort_value = case_when(or < 1 ~ distance_from_1, or > 1 ~ -distance_from_1, TRUE ~ 0)) %>%
    arrange(or < 1, sort_value)
  
  ordered_traits <- ivw_data$TRAIT
  missing_traits <- setdiff(unique(data_processed$TRAIT), ordered_traits)
  ordered_traits <- c(ordered_traits, missing_traits)
  
  data_processed <- data_processed %>% mutate(TRAIT = factor(TRAIT, levels = ordered_traits))
  start_degree <- 60
  
  or_mat <- acast(data_processed, TRAIT ~ method, value.var = "or")
  pval_mat <- acast(data_processed, TRAIT ~ method, value.var = "pval")
  log_or_mat <- log2(or_mat)
  
  pval_breaks <- c(0, 0.001, 0.005, 0.01, 0.02, 0.03, 0.04, 0.05, 0.051, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1)
  pval_colors <- c("#FF6666", "#FF8080", "#FF9999", "#FFB0B0", "#FFC2C2", "#FFD4D4", "#FFD4D4", "#FFE5E5", "#FFFFFF", "#E6F3FF", "#E3F0FD", "#E0EEFB", "#DDEBF9", "#DAE9F7", "#D7E6F5", "#D4E4F3", "#D1E1F1", "#CEDFEF", "#CBDCED", "#C8D9EB", "#C5D7E9", "#C2D4E7", "#BFD2E5", "#BCCFE3", "#B9CDE1", "#B6CADF", "#B3C8DD", "#B0C5DB")
  pval_col_fun <- colorRamp2(pval_breaks, pval_colors)
  
  or_breaks <- seq(0.5, 1.55, by = 0.01)
  n_breaks <- length(or_breaks)
  blue_colors <- colorRampPalette(c("#3C5488", "#4DBBD5", "#A5D6E6", "#D9EAF0"))(floor(n_breaks/2))
  red_colors <- colorRampPalette(c("#F7D6D0", "#F0A696", "#E64B35", "#B32E1C"))(ceiling(n_breaks/2))
  or_colors <- c(blue_colors, red_colors)
  or_col_fun <- colorRamp2(log2(or_breaks), or_colors)
  
  n_traits <- nrow(or_mat)
  label_cex <- dplyr::case_when(n_traits >= 220 ~ 1.0, n_traits >= 180 ~ 1.1, n_traits >= 140 ~ 1.3, n_traits >= 100 ~ 1.5, TRUE ~ 1.7)
  method_cex <- dplyr::case_when(n_traits >= 180 ~ 1.6, n_traits >= 120 ~ 1.8, TRUE ~ 2.0)
  plot_size <- max(20, min(34, 18 + n_traits / 18))
  
  draw_heatmap <- function(format) {
    if(format == "pdf") {
      pdf(file.path(plots_dir, paste0(out_prefix, ".pdf")), width = plot_size * 1.2, height = plot_size)
    } else {
      png(file.path(plots_dir, paste0(out_prefix, ".png")), width = plot_size * 1.2, height = plot_size, res = 300, units = "in", bg="white")
    }
    
    circos.clear()
    circos.par(start.degree = start_degree, gap.degree = 30, track.height = 0.18, track.margin = c(0.005, 0.005), cell.padding = c(0.01, 0, 0.01, 0), points.overflow.warning = FALSE)
    
    circos.heatmap(log_or_mat, col = or_col_fun, track.height = 0.34, cell.border = "gray90", cluster = FALSE, rownames.side = "outside", rownames.cex = label_cex)
    
    methods <- colnames(or_mat)
    n_methods <- length(methods)
    circos.track(track.index = get.current.track.index(), 
                 panel.fun = function(x, y) {
                   if(CELL_META$sector.numeric.index == 1) {
                     for(i in 1:n_methods) {
                       circos.text(CELL_META$cell.xlim[2] + mm_x(4), 
                                 n_methods - i + 0.5,
                                 methods[i],
                                 cex = method_cex,
                                 adj = c(0, 0.5),
                                 font = 2,
                                 facing = "inside")
                     }
                   }
                 })
    
    circos.heatmap(pval_mat, col = pval_col_fun, track.height = 0.22, cell.border = "gray90", cluster = FALSE, rownames.side = "none")
    
    # text(0, 1.10, paste("Odds Ratio Heatmap -", title_str), cex = 3.0, adj = c(0.5, 0), font = 2)
    # text(0, -1.10, "P-value Heatmap", cex = 3.0, adj = c(0.5, 1), font = 2)
    
    lgd_pval = Legend(title = "P-value", col_fun = pval_col_fun, at = c(1, 0.8, 0.6, 0.4, 0.2, 0.05, 0), labels = c("1.0", "0.8", "0.6", "0.4", "0.2", "0.05", ""), title_gp = gpar(fontsize = 26, fontface = "bold", just = "center"), labels_gp = gpar(fontsize = 24, just = "left"), grid_height = unit(8.0, "cm"), grid_width = unit(1.2, "cm"), legend_height = unit(8.0, "cm"), legend_width = unit(1.2, "cm"), title_gap = unit(1.0, "cm"))
    or_breaks_display <- c(0.50, 0.75, 1.00, 1.25, 1.5)
    lgd_or = Legend(title = "Odds Ratio", col_fun = or_col_fun, at = log2(or_breaks_display), labels = c("0.50", "0.75", "1.00", "1.25", "1.5"), title_gp = gpar(fontsize = 26, fontface = "bold", just = "center"), labels_gp = gpar(fontsize = 24, just = "left"), grid_height = unit(8.0, "cm"), grid_width = unit(1.2, "cm"), legend_height = unit(8.0, "cm"), legend_width = unit(1.2, "cm"), title_gap = unit(1.0, "cm"))
    
    pushViewport(viewport(x = 1.004, y = 0.99))
    draw(lgd_or, x = unit(-0.03, "npc"), y = unit(0.01, "npc"))
    draw(lgd_pval, x = unit(0.04, "npc"), y = unit(0.01, "npc"))
    upViewport()
    dev.off()
  }
  
  draw_heatmap("pdf")
  draw_heatmap("png")
}

# 5. 
p_v_hz <- plot_volcano(df_hz_primary, intersect_exp, "HZ (FinnGen)", "Volcano_Plot_HZ")
p_v_ra <- plot_volcano(df_ra_primary, intersect_exp, "RA", "Volcano_Plot_RA")

plot_intersection_forest(df_hz_primary, df_ra_primary, intersect_exp, "Forest_Plot_Intersection_Dual")

plot_heatmap(df_hz, sig_hz, "HZ (FinnGen)", "Circular_Heatmap_HZ_vertical")
plot_heatmap(df_ra, sig_ra, "RA", "Circular_Heatmap_RA_vertical")

# 6. 4: RA, HZ；, 
cat("Combining Heatmaps and Volcano plots (Left: RA, Right: HZ)...\n")
tryCatch({
  p_h_ra <- ggdraw() + draw_image(file.path(plots_dir, "Circular_Heatmap_RA_vertical.png"))
  p_h_hz <- ggdraw() + draw_image(file.path(plots_dir, "Circular_Heatmap_HZ_vertical.png"))
  
  # plot_grid by row:
  # row 1: RA Volcano, HZ Volcano
  # row 2: RA Heatmap, HZ Heatmap
  # This makes Left column = RA, Right column = HZ. Top row = Volcano, Bottom row = Heatmap.
  combined <- plot_grid(
    p_v_ra, p_v_hz, 
    p_h_ra, p_h_hz, 
    ncol = 2, nrow = 2, align = "v", labels = NULL
  )
  
  showtext_opts(dpi = 96)
  ggsave(file.path(plots_dir, "Combined_Heatmap_Volcano.pdf"), combined, width = 22, height = 24, bg = "white")
  showtext_opts(dpi = 300)
  ggsave(file.path(plots_dir, "Combined_Heatmap_Volcano.png"), combined, width = 22, height = 24, bg = "white", dpi = 300)
  showtext_opts(dpi = 96)
  cat("Combined plots saved successfully.\n")
}, error = function(e) {
  cat("Error in combining plots:", e$message, "\n")
})

# 7. 
cat("Saving statistics...\n")
stats_info <- paste(
  "=========================================\n",
  "MR Analysis Visualization Compare Statistics\n",
  "=========================================\n",
  "Project: UKB Plasma Protein cis-MR (250kb) HZ vs RA\n",
  "Time:", Sys.time(), "\n",
  "-----------------------------------------\n",
  "Inputs:\n",
  "HZ Results:", file_hz, "\n",
  "RA Results:", file_ra, "\n",
  "-----------------------------------------\n",
  "Significant exposures in HZ (FDR<0.05):", length(sig_hz), "\n",
  "Significant exposures in RA (FDR<0.05):", length(sig_ra), "\n",
  "Intersection exposures:", length(intersect_exp), "\n",
  "-----------------------------------------\n",
  "Outputs generated in:", output_dir, "\n",
  "=========================================\n"
)
writeLines(stats_info, file.path(output_dir, paste0("UKB_MR_Vis_Compare_Stats_", timestamp, ".txt")))

current_script <- file.path(work_dir, "MR_Visualization_Compare.R")
if(file.exists(current_script)) {
  file.copy(current_script, file.path(output_dir, paste0("MR_Visualization_Compare_Copy_", timestamp, ".R")))
}

end_time <- Sys.time()
run_time <- difftime(end_time, start_time, units = "mins")
cat("\nScript execution completed successfully.\n")
cat("Total run time:", round(run_time, 2), "minutes.\n")
sink()

cat("MR Visualization Compare script has been completed successfully.\n")
