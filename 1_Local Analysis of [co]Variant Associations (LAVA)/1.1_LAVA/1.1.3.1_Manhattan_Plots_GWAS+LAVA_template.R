#!/usr/bin/env Rscript
# ==============================================================================
# [Script]: 1.1.3.1_Manhattan_Plots_GWAS+LAVA_template.R
# [Method]: LAVA & GWAS Visualization
# [Step]: + (Visualization + Manhattan Plots)
# 
# [Function]:
# Generates high-quality, publication-ready Manhattan plots for multiple GWAS 
#       datasets (e.g., Aging, RA, HZ) and combines them with LAVA genetic correlation 
#       results. Creates a multi-panel figure highlighting the MHC region.
#       . MHC. 
# 
# [Data Availability / ]:
# Requires summary statistics for the analyzed traits and LAVA bivariate results.
# ==============================================================================

# ---------------------------------------------------------------------------------------------------------------------
# 1. Environment Setup & Arguments 
# ---------------------------------------------------------------------------------------------------------------------
set.seed(123)
rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(gridExtra)
  library(scales)
  library(patchwork)
  library(ragg)
})

option_list <- list(
  make_option(c("--gwas_mvAge"), type="character", default=NULL, help="Path to mvAge GWAS summary stats"),
  make_option(c("--gwas_RA"), type="character", default=NULL, help="Path to RA GWAS summary stats"),
  make_option(c("--gwas_HZ_UKB"), type="character", default=NULL, help="Path to HZ (UKB) GWAS summary stats"),
  make_option(c("--gwas_HZ_FinnGen"), type="character", default=NULL, help="Path to HZ (FinnGen) GWAS summary stats"),
  make_option(c("--lava_results"), type="character", default=NULL, help="Path to LAVA bivar results (mhc_bivar_results.txt)"),
  make_option(c("--out_dir"), type="character", default="./plots_output", help="Output directory path")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (any(sapply(list(opt$gwas_mvAge, opt$gwas_RA, opt$gwas_HZ_UKB, opt$gwas_HZ_FinnGen, opt$lava_results), is.null))) {
  print_help(opt_parser)
  stop("Missing required input paths for GWAS or LAVA results.")
}

out_dir <- opt$out_dir
dir_plots <- file.path(out_dir, "Plots")
if (!dir.exists(dir_plots)) {
  dir.create(dir_plots, recursive = TRUE)
}

log_file <- file.path(out_dir, "analysis_log.txt")
log_conn <- file(log_file, open = "wt")

write_log <- function(message) {
  cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", message, "\n"))
  cat(paste0("[", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "] ", message, "\n"), file = log_conn)
}

write_log("Starting visualization script...")

# ---------------------------------------------------------------------------------------------------------------------
# 2. Define Datasets 
# ---------------------------------------------------------------------------------------------------------------------
datasets <- list(
  list(name = "mvAge", path = opt$gwas_mvAge, cols = c(SNP="SNP", CHR="CHR", BP="BP", P="Pvalue", EAF="eaf"), is_gz=grepl("\\.gz$", opt$gwas_mvAge)),
  list(name = "RA", path = opt$gwas_RA, cols = c(SNP="rsid", CHR="chromosome", BP="base_pair_location", P="p_value", EAF="effect_allele_frequency"), is_gz=grepl("\\.gz$", opt$gwas_RA)),
  list(name = "HZ (UKB)", path = opt$gwas_HZ_UKB, cols = c(SNP="rsid", CHR="chromosome", BP="base_pair_location", P="p_value", EAF="effect_allele_frequency"), is_gz=grepl("\\.gz$", opt$gwas_HZ_UKB)),
  list(name = "HZ (FinnGen)", path = opt$gwas_HZ_FinnGen, cols = c(SNP="rsids", CHR="#chrom", BP="pos", P="pval", EAF="af_alt"), is_gz=grepl("\\.gz$", opt$gwas_HZ_FinnGen))
)

# ---------------------------------------------------------------------------------------------------------------------
# 3. Processing Function 
# ---------------------------------------------------------------------------------------------------------------------
generate_plots <- function(dataset_info, show_x_title = TRUE) {
  write_log(paste0("Processing dataset: ", dataset_info$name))
  
  if (!file.exists(dataset_info$path)) {
    write_log(paste0("ERROR: File not found: ", dataset_info$path))
    return(NULL)
  }
  
  if (dataset_info$is_gz) {
    header_line <- system(paste("gzcat", shQuote(dataset_info$path), "| head -n 1"), intern = TRUE)
  } else {
    header_line <- system(paste("head -n 1", shQuote(dataset_info$path)), intern = TRUE)
  }
  
  col_names <- strsplit(header_line, "\t")[[1]]
  col_idx <- match(dataset_info$cols, col_names)
  names(col_idx) <- names(dataset_info$cols)
  
  if (any(is.na(col_idx))) {
    missing_cols <- names(col_idx)[is.na(col_idx)]
    stop(paste("Columns not found in", dataset_info$name, ":", paste(missing_cols, collapse=", ")))
  }
  
  awk_script <- sprintf(
    "BEGIN {FS=\"\\t\"; OFS=\"\\t\"} NR==1 {print \"SNP\\tCHR\\tBP\\tP\\tEAF\"; next} { eaf=$%d; maf=(eaf<0.5?eaf:1-eaf); if(maf>0.01) print $%d, $%d, $%d, $%d, $%d }",
    col_idx["EAF"], col_idx["SNP"], col_idx["CHR"], col_idx["BP"], col_idx["P"], col_idx["EAF"]
  )
  
  cmd <- paste(ifelse(dataset_info$is_gz, "gzcat", "cat"), shQuote(dataset_info$path), "| awk", shQuote(awk_script))
  df_plot <- fread(cmd = cmd, header = TRUE, data.table = FALSE)
  
  suppressWarnings({
    df_plot$CHR <- as.numeric(as.character(df_plot$CHR))
    df_plot$BP <- as.numeric(as.character(df_plot$BP))
    df_plot$P <- as.numeric(as.character(df_plot$P))
  })
  
  df_plot <- df_plot %>% filter(!is.na(CHR) & !is.na(BP) & !is.na(P) & CHR %in% 1:22)
  write_log(paste0("Data loaded and filtered. Rows: ", nrow(df_plot)))
  
  shared_y_limit <- max(-log10(df_plot$P), na.rm = TRUE) * 1.05
  high_signal_cutoff <- 0.01
  max_background_points <- 25000 
  
  df_high <- df_plot %>% filter(P < high_signal_cutoff)
  df_low <- df_plot %>% filter(P >= high_signal_cutoff)
  
  df_low_sampled <- if (nrow(df_low) > max_background_points) sample_n(df_low, max_background_points) else df_low
  df_wg <- bind_rows(df_high, df_low_sampled)
  
  data_cum <- df_wg %>% group_by(CHR) %>% summarise(max_bp = max(BP)) %>% mutate(bp_add = lag(cumsum(max_bp), default = 0)) %>% select(CHR, bp_add)
  df_wg <- df_wg %>% inner_join(data_cum, by = "CHR") %>% mutate(bp_cum = BP + bp_add)
  axis_set <- df_wg %>% group_by(CHR) %>% summarize(center = mean(bp_cum))
  
  df_wg <- df_wg %>%
    mutate(
      is_mhc = (CHR == 6 & BP >= 25000000 & BP <= 34000000),
      color_group = case_when(is_mhc ~ "MHC", CHR %% 2 == 1 ~ "Odd", TRUE ~ "Even")
    )
  df_wg$color_group <- factor(df_wg$color_group, levels = c("Odd", "Even", "MHC"))
  df_wg <- df_wg %>% arrange(color_group) 
  sig_threshold <- 5e-8
  
  p_wg <- ggplot(df_wg, aes(x = bp_cum, y = -log10(P))) +
    geom_point(aes(color = color_group), alpha = 0.8, size = 1.0) +
    scale_color_manual(values = c("Odd" = "grey60", "Even" = "black", "MHC" = "#8E44AD")) +
    scale_x_continuous(label = axis_set$CHR, breaks = axis_set$center) +
    scale_y_continuous(limits = c(0, shared_y_limit)) +
    geom_hline(yintercept = -log10(sig_threshold), color = "red", linetype = "dashed") +
    labs(x = ifelse(show_x_title, "Chromosome", ""), y = "-log10(P)") +
    theme_minimal(base_size = 14) +
    theme(legend.position = "none", panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(angle = 0, size = 10), axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 14, face = "bold"), plot.title = element_blank())
  
  df_mhc <- df_plot %>% filter(CHR == 6 & BP >= 25000000 & BP <= 34000000)
  df_mhc$Phenotype <- dataset_info$name
    
  if (nrow(df_mhc) > 0) {
    p_mhc <- ggplot(df_mhc, aes(x = BP/1e6, y = -log10(P))) +
      geom_point(color = "#8E44AD", alpha = 0.8, size = 1.5) +
      scale_y_continuous(limits = c(0, shared_y_limit)) +
      geom_hline(yintercept = -log10(sig_threshold), color = "red", linetype = "dashed") +
      labs(x = ifelse(show_x_title, "Chromosome 6", ""), y = "-log10(P)") +
      facet_grid(Phenotype ~ .) +
      theme_minimal(base_size = 14) +
      theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(),
            legend.position = "none", panel.grid.minor.x = element_blank(), axis.text.x = element_text(size = 10),
            axis.title = element_text(size = 14, face = "bold"), plot.title = element_blank(),
            strip.background = element_rect(fill = "grey90", color = NA), strip.text = element_text(face = "bold", size = 12, color = "black"))
  } else {
    p_mhc <- ggplot(data.frame(Phenotype = dataset_info$name)) + geom_blank() + facet_grid(Phenotype ~ .) + 
      labs(title = "No SNPs in MHC") + theme_void() +
      theme(strip.background = element_rect(fill = "grey90", color = NA), strip.text = element_text(face = "bold", size = 12, color = "black"))
  }
  
  return(list(wg = p_wg, mhc = p_mhc))
}

# ---------------------------------------------------------------------------------------------------------------------
# 3.1 LAVA Manhattan Plot ( LAVA )
# ---------------------------------------------------------------------------------------------------------------------
generate_lava_manhattan_plot <- function(annotate_mhc = FALSE) {
  if (!file.exists(opt$lava_results)) stop(paste0("LAVA result file not found: ", opt$lava_results))

  res_lava <- fread(opt$lava_results)
  res_lava[, pos_mb := (START + STOP) / 2 / 1e6]
  res_lava[, logp := -log10(p)]

  pheno_map <- c(herpes_zoster = "HZ (UKB)", aging = "mvAge", RA = "RA", FinnGen_HZ = "HZ (FinnGen)")
  res_lava[, phen1_short := pheno_map[phen1]]
  res_lava[, phen2_short := pheno_map[phen2]]
  res_lava[is.na(phen1_short), phen1_short := as.character(phen1)]
  res_lava[is.na(phen2_short), phen2_short := as.character(phen2)]

  pheno_priority <- c("mvAge" = 1, "RA" = 2, "HZ (UKB)" = 3, "HZ (FinnGen)" = 4)
  res_lava[, rank1 := pheno_priority[phen1_short]]
  res_lava[, rank2 := pheno_priority[phen2_short]]
  res_lava[is.na(rank1), rank1 := 99]
  res_lava[is.na(rank2), rank2 := 99]
  
  res_lava[, pair := ifelse(rank1 < rank2, paste0(phen1_short, " vs ", phen2_short), paste0(phen2_short, " vs ", phen1_short))]
  
  desired_order <- c("mvAge vs RA", "mvAge vs HZ (UKB)", "mvAge vs HZ (FinnGen)", "RA vs HZ (UKB)", "RA vs HZ (FinnGen)", "HZ (UKB) vs HZ (FinnGen)")
  final_levels <- c(intersect(desired_order, unique(res_lava$pair)), setdiff(unique(res_lava$pair), desired_order))
  res_lava[, pair := factor(pair, levels = final_levels)]

  theme_pub_lava <- theme_bw(base_size = 14) +
    theme(text = element_text(family = "sans"), panel.grid.minor = element_blank(), panel.grid.major = element_line(color = "grey92"),
          panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), axis.ticks = element_line(colour = "black"),
          legend.position = "right", legend.background = element_rect(fill = "white", color = NA), plot.title = element_blank(),
          plot.subtitle = element_blank(), strip.background = element_rect(fill = "grey90", color = NA), strip.text = element_text(face = "bold", size = 12, color = "black"),
          axis.title = element_text(size = 14, face = "bold"))

  p <- ggplot(res_lava, aes(x = pos_mb, y = logp))
  
  if (annotate_mhc) {
    mhc_regions <- data.frame(xmin = c(29.672373, 31.511125, 32.224068), xmax = c(31.511124, 32.224067, 33.148800), Region = c("Class I", "Class III", "Class II"))
    mhc_regions$Region <- factor(mhc_regions$Region, levels = c("Class I", "Class III", "Class II"))
    p <- p + geom_rect(data = mhc_regions, aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf, fill = Region), inherit.aes = FALSE, alpha = 0.2) +
             scale_fill_manual(values = c("Class I" = "#D9D9D9", "Class III" = "#BDBDBD", "Class II" = "#969696"), name = "MHC Region")
  }

  p <- p + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red", alpha = 0.7) +
    geom_point(aes(color = rho, size = abs(rho)), alpha = 0.9) +
    scale_color_gradient2(low = "#2166AC", mid = "#F7F7F7", high = "#B2182B", midpoint = 0, limit = c(-1, 1), oob = scales::squish, name = "Genetic Correlation (Rho)") +
    scale_size_continuous(range = c(2, 6), name = "|Rho|") +
    facet_grid(pair ~ ., scales = "free_y") + labs(x = "Chromosome 6", y = "-log10(P-value)") +
    guides(color = guide_colorbar(order = 1), size = guide_legend(order = 2), fill = guide_legend(order = 3)) +
    theme_pub_lava
  return(p)
}

# ---------------------------------------------------------------------------------------------------------------------
# 4. Main Execution 
# ---------------------------------------------------------------------------------------------------------------------
all_plots <- list()
for (i in seq_along(datasets)) {
  info <- datasets[[i]]
  res <- generate_plots(info, show_x_title = (i == 4))
  if (!is.null(res)) all_plots[[info$name]] <- res
}

if (length(all_plots) == 4) {
  write_log("Combining plots...")
  row1 <- all_plots[["mvAge"]]$wg + all_plots[["mvAge"]]$mhc + plot_layout(ncol = 2, widths = c(2, 1))
  row2 <- all_plots[["RA"]]$wg + all_plots[["RA"]]$mhc + plot_layout(ncol = 2, widths = c(2, 1))
  row3 <- all_plots[["HZ (UKB)"]]$wg + all_plots[["HZ (UKB)"]]$mhc + plot_layout(ncol = 2, widths = c(2, 1))
  row4 <- all_plots[["HZ (FinnGen)"]]$wg + all_plots[["HZ (FinnGen)"]]$mhc + plot_layout(ncol = 2, widths = c(2, 1))
  gwas_plot <- (row1 / row2 / row3 / row4)
  
  lava_plot_annotated <- generate_lava_manhattan_plot(annotate_mhc = TRUE)
  final_plot_annotated <- (gwas_plot | lava_plot_annotated) + plot_layout(ncol = 2, widths = c(3, 2))
  final_plot_annotated_labeled <- final_plot_annotated + plot_annotation(tag_levels = 'a') & theme(plot.tag = element_text(size = 32, face = "bold"))
  
  outfile_base <- file.path(dir_plots, "Combined_Manhattan_Plots_Aging_RA_HZ_FinnGen_GWAS_LAVA")
  ggsave(paste0(outfile_base, "_GWAS_Only.pdf"), plot = gwas_plot, width = 14.4, height = 20, bg = "white")
  ggsave(paste0(outfile_base, "_LAVA_Annotated_Only.pdf"), plot = lava_plot_annotated, width = 9.6, height = 20, bg = "white")
  ggsave(paste0(outfile_base, "_Annotated_Combined.pdf"), plot = final_plot_annotated, width = 24, height = 20, bg = "white")
  write_log("Saved all combined PDF figures.")
} else {
  write_log("ERROR: Not all datasets were processed successfully. Skipping combination.")
}

close(log_conn)
print(paste0("Analysis finished. Output in: ", out_dir))
