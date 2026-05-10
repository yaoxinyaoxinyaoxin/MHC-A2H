#!/usr/bin/env Rscript
# ==============================================================================
# [Script]: 7.9.3.1_plot_manhattan_template.R
# [Method]: GWAS/PheWAS Manhattan Plot
# [Step]: Visual Alignment of Manhattan Plots
# 
# [Function]:
# Generate Manhattan plots for visually aligning signals across traits/cohorts.
# 
# [Input]:
#   --input_file   : Path to GWAS/PheWAS summary statistics
#   --out_dir      : Output directory path
#   --title        : Plot title
# 
# [Output]:
# High-quality Manhattan plot (PDF/PNG).
# ==============================================================================

rm(list = ls())
gc()

# 1. Load Dependencies
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(ggplot2)
  library(CMplot)
  library(gridExtra)
})

# 2. Parse Command Line Arguments
option_list <- list(
  make_option(c("--input_file"), type="character", default=NULL, help="Input summary stats file"),
  make_option(c("--out_dir"), type="character", default="./Manhattan_Plots", help="Output directory path"),
  make_option(c("--title"), type="character", default="Manhattan Plot", help="Title of the plot")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_file)) {
  print_help(opt_parser)
  stop("Missing required input file.")
}

INPUT_FILE <- opt$input_file
OUTPUT_DIR <- opt$out_dir
PLOT_TITLE <- opt$title

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

packages <- c("data.table", "ggplot2", "dplyr", "stringr", "parallel", "patchwork", "ragg")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
    library(pkg, character.only = TRUE)
  }
}

# , 
# : , 
setDTthreads(0) # 

# ------------------------------------------------------------------------------
# Step 1: Define Paths and Setup Directories
# : 1: 
# ------------------------------------------------------------------------------
#  (Input path)
input_dir <- "./"
# manifest (FinnGen manifest file)
manifest_file <- "./"
#  (Base output path)
base_out_dir <- OUTPUT_DIR
# PLINK (PLINK path)
plink_bin <- "./"
# LD  (LD reference panel path - hg19/GRCh37)
ld_panel <- "./"

# ==========================================
# 1. /SNP (GRCh37)
# ==========================================
target_snp <- "rs1800628"

# : chr6:31546850 (GRCh37)
target_chr <- 6
target_pos <- 31546850
cat(sprintf("Using hardcoded GRCh37 coordinates for %s: chr%d:%d\n", target_snp, target_chr, target_pos))

#  ( 1Mb)
region_start <- target_pos - 1000000
region_end <- target_pos + 1000000
region_title <- sprintf("Region %s±1Mb (chr%d:%d-%d)", target_snp, target_chr, region_start, region_end)

#  (Create timestamped output folder)
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path(base_out_dir, paste0("1_Manhattan_FinnGen_Plots_", timestamp))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# log
log_dir <- file.path(out_dir, "logs")
scripts_dir <- file.path(out_dir, "scripts")
ld_out_dir <- file.path(out_dir, "LD_results")
dir.create(log_dir, showWarnings = FALSE)
dir.create(scripts_dir, showWarnings = FALSE)
dir.create(ld_out_dir, showWarnings = FALSE)

#  (Initialize log file)
log_file <- file.path(log_dir, paste0("run_log_", timestamp, ".txt"))
log_msg <- function(msg) {
  time_stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", time_stamp, msg), file = log_file, append = TRUE)
  cat(sprintf("[%s] %s\n", time_stamp, msg))
}
start_time <- Sys.time()
log_msg("Script started.")

args <- commandArgs(trailingOnly = FALSE)
script_name <- sub("--file=", "", args[grep("--file=", args)])
if(length(script_name) > 0 && file.exists(script_name)) {
  file.copy(script_name, file.path(scripts_dir, basename(script_name)))
  log_msg("Script backed up to scripts directory.")
} else if(file.exists(file.path(base_out_dir, "1_plot_manhattan.R"))) {
  file.copy(file.path(base_out_dir, "1_plot_manhattan.R"), file.path(scripts_dir, "1_plot_manhattan.R"))
  log_msg("Script backed up to scripts directory.")
}

# ------------------------------------------------------------------------------
# Step 2: Read Mapping Files and Identify Target FinnGen Files
# : 2: FinnGen
# ------------------------------------------------------------------------------
log_msg("Step 2: Identifying target FinnGen files and phenotypes...")

manifest <- fread(manifest_file)
#  .gz 
gwas_files <- list.files(input_dir, pattern = "^finngen_R12_.*\\.gz$", full.names = TRUE)

target_targets <- data.table(File = character(), Phenocode = character(), Phenotype = character())

for (f in gwas_files) {
  filename <- basename(f)
  #  phenocode ( finngen_R12_AUTOIMMUNE.gz -> AUTOIMMUNE)
  curr_phenocode <- str_replace(filename, "^finngen_R12_", "")
  curr_phenocode <- str_replace(curr_phenocode, "\\.gz$", "")
  
  #  manifest 
  matched_pheno <- manifest[phenocode == curr_phenocode, phenotype]
  if (length(matched_pheno) > 0) {
    pheno_name <- matched_pheno[1]
  } else {
    pheno_name <- curr_phenocode #  curr_phenocode
  }
  
  target_targets <- rbind(target_targets, data.table(File = f, Phenocode = curr_phenocode, Phenotype = pheno_name))
}

log_msg(sprintf("Identified %d target FinnGen files to process.", nrow(target_targets)))
target_targets[, Disease_Name := Phenotype]
target_targets[, Disease_Name := str_replace(Disease_Name, "\\(FG\\)", "")]
target_targets[, Disease_Name := str_trim(Disease_Name)]
target_targets[, Disease_Name := paste0(Disease_Name, " (FinnGen)")]
target_targets[Phenocode == "SLE_FG", Disease_Name := "Systemic lupus erythematosus (FinnGen)"]
target_targets[Phenocode == "M13_SJOGREN", Disease_Name := "Sicca syndrome (FinnGen)"]
target_targets[Phenocode == "DERMATOPOLY_FG", Disease_Name := "Dermatopolymyositis (FinnGen)"]

# ------------------------------------------------------------------------------
# Step 3: Define Plotting Functions and Colors
# : 3: 
# ------------------------------------------------------------------------------
ld_colors <- c(
  "R2 > 0.8" = "#D62828", # 
  "0.6 < R2 <= 0.8" = "#F77F00", # 
  "0.4 < R2 <= 0.6" = "#FCBF49", # 
  "0.2 < R2 <= 0.4" = "#88D49E", # 
  "R2 <= 0.2" = "#118AB2", # 
  "Unknown" = "#E0E0E0" # 
)

plot_manhattan <- function(df, title = NULL, region = NULL, use_ld = TRUE, xlab_text = NULL, custom_y_limit = NULL) {
  #  (-log10(P))
  df_plot <- df %>%
    dplyr::mutate(logP = -log10(P)) %>%
    dplyr::mutate(logP = ifelse(is.infinite(logP), max(c(10, logP[is.finite(logP)]), na.rm=TRUE) + 5, logP))

  if (is.null(region)) {
    df_plot <- df_plot %>%
      dplyr::arrange(CHR, BP) %>%
      dplyr::mutate(CHR = as.integer(CHR)) %>%
      dplyr::filter(!is.na(CHR))
      
    don <- df_plot %>% 
      dplyr::group_by(CHR) %>% 
      dplyr::summarise(chr_len = max(BP, na.rm=TRUE)) %>% 
      dplyr::mutate(tot = cumsum(as.numeric(chr_len)) - chr_len) %>%
      dplyr::select(-chr_len)
    
    df_plot <- dplyr::left_join(df_plot, don, by="CHR") %>%
      dplyr::arrange(CHR, BP) %>%
      dplyr::mutate(BPcum = BP + tot)
      
    axisdf <- df_plot %>% dplyr::group_by(CHR) %>% dplyr::summarize(center = (max(BPcum) + min(BPcum)) / 2)
    
    if (is.null(custom_y_limit)) {
      max_y <- max(df_plot$logP, na.rm = TRUE)
      y_limit <- max_y * 1.1
    } else {
      y_limit <- custom_y_limit
    }
    
    if(use_ld && "R2_group" %in% colnames(df_plot)) {
      df_plot <- df_plot %>%
        dplyr::mutate(plot_color_group = dplyr::case_when(
          R2_group != "Unknown" ~ as.character(R2_group),
          CHR %% 2 == 0 ~ "Chr_Light", 
          TRUE ~ "Chr_Dark"
        ))
        
      df_plot$plot_color_group <- factor(df_plot$plot_color_group, 
        levels = c("Chr_Dark", "Chr_Light", "R2 <= 0.2", "0.2 < R2 <= 0.4", "0.4 < R2 <= 0.6", "0.6 < R2 <= 0.8", "R2 > 0.8"))
      
      df_plot <- df_plot %>% dplyr::arrange(plot_color_group)
      
      wg_colors <- c(
        "Chr_Dark" = "#A9A9A9",
        "Chr_Light" = "#E0E0E0",
        "R2 <= 0.2" = "#118AB2",
        "0.2 < R2 <= 0.4" = "#88D49E",
        "0.4 < R2 <= 0.6" = "#FCBF49",
        "0.6 < R2 <= 0.8" = "#F77F00",
        "R2 > 0.8" = "#D62828"
      )
      
      p <- ggplot(df_plot, aes(x = BPcum, y = logP, color = plot_color_group)) +
        geom_point(alpha = 0.8, size = 1.5) +
        scale_color_manual(values = wg_colors, breaks = c("R2 > 0.8", "0.6 < R2 <= 0.8", "0.4 < R2 <= 0.6", "0.2 < R2 <= 0.4", "R2 <= 0.2")) +
        scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, max(c(10, y_limit)))) +
        geom_hline(yintercept = -log10(5e-8), color = "#e63946", linetype = "dashed", linewidth = 0.8) +
        geom_hline(yintercept = -log10(1e-5), color = "#F77F00", linetype = "dashed", linewidth = 0.8) +
        labs(title = NULL, x = ifelse(is.null(xlab_text), "Chromosome", xlab_text), 
             y = expression(-log[10](italic(P))), color = bquote(R^2 ~ "with " ~ .(target_snp))) +
        theme_classic(base_size = 14) +
        theme(
          legend.position = c(0.01, 0.99),
          legend.justification = c("left", "top"),
          legend.background = element_blank(),
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(size = 10, vjust = 0.5, angle = 0),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 14, face = "bold"),
          plot.title = element_blank(),
          plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10)
        )
    } else {
      p <- ggplot(df_plot, aes(x = BPcum, y = logP)) +
        geom_point(aes(color = as.factor(CHR)), alpha = 0.8, size = 1.5) +
        scale_color_manual(values = rep(c("#457b9d", "#1d3557"), 22)) +
        scale_x_continuous(label = axisdf$CHR, breaks = axisdf$center) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, max(c(10, y_limit)))) +
        geom_hline(yintercept = -log10(5e-8), color = "#e63946", linetype = "dashed", linewidth = 0.8) +
        geom_hline(yintercept = -log10(1e-5), color = "#F77F00", linetype = "dashed", linewidth = 0.8) +
        labs(title = NULL, x = ifelse(is.null(xlab_text), "Chromosome", xlab_text), y = expression(-log[10](italic(P)))) +
        theme_classic(base_size = 14) +
        theme(
          legend.position = "none",
          panel.border = element_blank(),
          panel.grid.major.x = element_blank(),
          panel.grid.minor.x = element_blank(),
          axis.text.x = element_text(size = 10, vjust = 0.5, angle = 0),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 14, face = "bold"),
          plot.title = element_blank(),
          plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
        )
    }
  } else {
    if (is.null(custom_y_limit)) {
      max_y <- max(df_plot$logP, na.rm = TRUE)
      y_limit <- max_y * 1.1
    } else {
      y_limit <- custom_y_limit
    }
    
    if(use_ld && "R2_group" %in% colnames(df_plot)) {
      df_plot$R2_group <- factor(df_plot$R2_group, levels = c("Unknown", "R2 <= 0.2", "0.2 < R2 <= 0.4", "0.4 < R2 <= 0.6", "0.6 < R2 <= 0.8", "R2 > 0.8"))
      df_plot <- df_plot %>% dplyr::arrange(R2_group)
      
      p <- ggplot(df_plot, aes(x = BP / 1e6, y = logP, color = R2_group)) + 
        geom_point(alpha = 0.8, size = 2.5) +
        scale_color_manual(values = ld_colors, breaks = c("R2 > 0.8", "0.6 < R2 <= 0.8", "0.4 < R2 <= 0.6", "0.2 < R2 <= 0.4", "R2 <= 0.2")) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, y_limit)) +
        geom_hline(yintercept = -log10(5e-8), color = "#D62828", linetype = "dashed", linewidth = 0.8) +
        geom_hline(yintercept = -log10(1e-5), color = "#F77F00", linetype = "dashed", linewidth = 0.8) +
        labs(title = NULL, 
             x = ifelse(is.null(xlab_text), paste0("Chromosome ", region$chr, " Position (Mb)"), xlab_text), 
             y = expression(-log[10](italic(P))),
             color = bquote(R^2 ~ "with " ~ .(target_snp))) +
        theme_classic(base_size = 14) +
        theme(
          legend.position = c(0.01, 0.99),
          legend.justification = c("left", "top"),
          legend.background = element_blank(),
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 14, face = "bold"),
          plot.title = element_blank(),
          plot.margin = margin(t = 20, r = 20, b = 20, l = 20),
          legend.title = element_text(size = 12, face = "bold"),
          legend.text = element_text(size = 10)
        )
    } else {
      p <- ggplot(df_plot, aes(x = BP / 1e6, y = logP)) + 
        geom_point(color = "#457b9d", alpha = 0.8, size = 2.0) +
        scale_y_continuous(expand = c(0, 0), limits = c(0, y_limit)) +
        geom_hline(yintercept = -log10(5e-8), color = "#e63946", linetype = "dashed", linewidth = 0.8) +
        geom_hline(yintercept = -log10(1e-5), color = "#F77F00", linetype = "dashed", linewidth = 0.8) +
        labs(title = NULL, 
             x = ifelse(is.null(xlab_text), paste0("Chromosome ", region$chr, " Position (Mb)"), xlab_text), 
             y = expression(-log[10](italic(P)))) +
        theme_classic(base_size = 14) +
        theme(
          legend.position = "none",
          axis.text.x = element_text(size = 12),
          axis.text.y = element_text(size = 12),
          axis.title = element_text(size = 14, face = "bold"),
          plot.title = element_blank(),
          plot.margin = margin(t = 20, r = 20, b = 20, l = 20)
        )
    }
  }
  
  return(p)
}

# ------------------------------------------------------------------------------
# Step 4: Extract SNPs using awk concurrently
# : 4: awkSNP
# ------------------------------------------------------------------------------
log_msg(sprintf("Step 4: Extracting SNPs in chr%d:%d-%d using awk concurrently...", target_chr, region_start, region_end))

target_files <- target_targets$File

num_awk_cores <- min(16, length(target_files))
log_msg(sprintf("Using %d cores for parallel awk extraction...", num_awk_cores))

extract_snps_awk <- function(file_path, chr, start, end) {
  awk_script <- sprintf("
  BEGIN { FS=\"\\t\"; OFS=\"\\t\" }
  NR==1 { 
    for(i=1;i<=NF;i++) { 
      if($i==\"#chrom\" || $i==\"chromosome\") chr_col=i; 
      if($i==\"pos\" || $i==\"base_pair_location\") bp_col=i; 
      if($i==\"rsids\" || $i==\"rsid\") rsid_col=i;
    } 
  }
  NR>1 { 
    if($chr_col == %d && $bp_col >= %d && $bp_col <= %d) { 
      snp_name = (rsid_col != \"\" && $rsid_col != \"\" && $rsid_col != \"NA\") ? $rsid_col : $chr_col\":\"$bp_col;
      print snp_name
    } 
  }
  ", chr, start, end)
  
  cmd <- sprintf("gzcat '%s' | awk '%s'", file_path, awk_script)
  snps <- system(cmd, intern = TRUE)
  return(snps)
}

cl_awk <- makeCluster(num_awk_cores)
clusterExport(cl_awk, c("extract_snps_awk", "target_chr", "region_start", "region_end"))
region_snps_list <- parLapply(cl_awk, target_files, function(f) {
  extract_snps_awk(f, target_chr, region_start, region_end)
})
stopCluster(cl_awk)

# ------------------------------------------------------------------------------
# Step 5: Calculate LD using extracted SNPs
# : 5: SNPLD
# ------------------------------------------------------------------------------
log_msg(sprintf("Step 5: Calculating LD for union SNPs in chr%d:%d-%d (±1Mb of %s)...", target_chr, region_start, region_end, target_snp))

ld_map <- data.table(SNP = character(), R2 = numeric())

if(length(region_snps_list) > 0) {
  union_snps <- unique(unlist(region_snps_list))
  union_snps <- union_snps[union_snps != "" & union_snps != "rsid" & union_snps != "rsids"]
  log_msg(sprintf("Found %d union SNPs in the target region.", length(union_snps)))
  
  if(length(union_snps) > 0) {
    snp_list_file <- file.path(ld_out_dir, "union_snps.txt")
    write.table(union_snps, file = snp_list_file, quote = FALSE, row.names = FALSE, col.names = FALSE)
    
    if(!(target_snp %in% union_snps)) {
      write.table(target_snp, file = snp_list_file, quote = FALSE, row.names = FALSE, col.names = FALSE, append = TRUE)
      log_msg(sprintf("Added target SNP %s to the list.", target_snp))
    }
    
    log_msg(sprintf("Calculating LD with %s using PLINK...", target_snp))
    plink_out_prefix <- file.path(ld_out_dir, "ld_results")
    
    plink_cmd <- sprintf("'%s' --bfile '%s' --r2 --ld-snp %s --ld-window 99999 --ld-window-kb 2000 --ld-window-r2 0 --extract '%s' --out '%s'", 
                         plink_bin, ld_panel, target_snp, snp_list_file, plink_out_prefix)
    
    system(plink_cmd)
    
    ld_res_file <- paste0(plink_out_prefix, ".ld")
    if(file.exists(ld_res_file)) {
      log_msg("LD calculation successful. Processing results...")
      ld_res <- fread(ld_res_file)
      
      fwrite(ld_res, file.path(ld_out_dir, "ld_results.csv"))
      log_msg("Saved LD results as CSV format.")
      
      ld_map <- ld_res[, .(SNP = SNP_B, R2)]
      ld_map <- rbind(ld_map, data.table(SNP = target_snp, R2 = 1))
      ld_map <- unique(ld_map, by = "SNP")
      
      ld_map[, R2_group := cut(R2, 
                               breaks = c(-Inf, 0.2, 0.4, 0.6, 0.8, Inf), 
                               labels = c("R2 <= 0.2", "0.2 < R2 <= 0.4", "0.4 < R2 <= 0.6", "0.6 < R2 <= 0.8", "R2 > 0.8"))]
      
      strong_ld_snps <- ld_map[R2 > 0.8, SNP]
      write.table(strong_ld_snps, file = file.path(ld_out_dir, "strong_ld_snps_R2_gt_0.8.txt"), 
                  quote = FALSE, row.names = FALSE, col.names = FALSE)
      log_msg(sprintf("Found %d SNPs in strong LD (R2 > 0.8) with %s.", length(strong_ld_snps), target_snp))
      
    } else {
      log_msg("Error: PLINK LD calculation failed or output not found.")
    }
  } else {
    log_msg("No union SNPs found. Skipping LD calculation.")
  }
} else {
  log_msg("No data available to find union SNPs.")
}

# ------------------------------------------------------------------------------
# Step 6: Reading GWAS datasets, calculating LD, plotting, and combining images using parallel processing
# : 6: 、LD
# ------------------------------------------------------------------------------
log_msg("Step 6: Reading GWAS datasets, calculating LD, plotting, and combining images using parallel processing...")

process_and_plot_worker <- function(i, target_targets, out_dir, target_chr, region_start, region_end, target_snp, ld_map, plot_manhattan, ld_colors) {
  library(data.table)
  library(ggplot2)
  library(dplyr)
  library(stringr)
  library(patchwork)
  setDTthreads(1)
  
  file_path <- target_targets$File[i]
  phenocode <- target_targets$Phenocode[i]
  phenotype_name <- target_targets$Phenotype[i]
  disease_name <- target_targets$Disease_Name[i]
  
  safe_pheno_name <- str_replace_all(phenotype_name, "[^A-Za-z0-9_]", "_")
  safe_pheno_name <- paste0(safe_pheno_name, "_", phenocode)
  pheno <- safe_pheno_name
  
  if(!file.exists(file_path)) {
    return(list(pheno = pheno, status = "Failed", error = "File not found", stats = NULL))
  }
  
  tryCatch({
    header <- fread(cmd = paste("gzcat", shQuote(file_path), "| head -n 1"), nrows = 0)
    cols_to_read <- intersect(c("#chrom", "chromosome", "pos", "base_pair_location", "pval", "p_value", "rsids", "rsid", "af_alt"), colnames(header))
    
    dat <- fread(cmd = paste("gzcat", shQuote(file_path)), select = cols_to_read)
    
    # MAF 
    if ("af_alt" %in% colnames(dat)) {
      dat[, MAF := pmin(af_alt, 1 - af_alt)]
      dat <- dat[MAF > 0.01]
      dat[, c("af_alt", "MAF") := NULL]
    }
    gc()
    
    if ("#chrom" %in% colnames(dat)) setnames(dat, "#chrom", "CHR")
    if ("chromosome" %in% colnames(dat)) setnames(dat, "chromosome", "CHR")
    if ("pos" %in% colnames(dat)) setnames(dat, "pos", "BP")
    if ("base_pair_location" %in% colnames(dat)) setnames(dat, "base_pair_location", "BP")
    if ("pval" %in% colnames(dat)) setnames(dat, "pval", "P")
    if ("p_value" %in% colnames(dat)) setnames(dat, "p_value", "P")
    
    if ("rsids" %in% colnames(dat) && !all(is.na(dat$rsids))) {
      dat[, SNP := rsids]
    } else if ("rsid" %in% colnames(dat) && !all(is.na(dat$rsid))) {
      dat[, SNP := rsid]
    } else {
      dat[, SNP := paste0(CHR, ":", BP)]
    }
    
    dat <- dat[!is.na(P) & !is.na(CHR) & !is.na(BP)]
    dat <- dat[CHR %in% 1:22]
    
    stats_row <- data.frame(
      Phenotype = phenotype_name,
      Phenocode = phenocode,
      Total_SNPs = nrow(dat),
      Min_P = min(dat$P, na.rm = TRUE),
      Sig_SNPs_5e8 = nrow(dat[P < 5e-8])
    )
    
    dat_sig <- dat[P < 0.05 | (CHR == 6 & BP >= 25000000 & BP <= 35000000)]
    rm(dat)
    gc()
    
    if(nrow(ld_map) > 0) {
      dat_ld <- merge(dat_sig, ld_map, by = "SNP", all.x = TRUE)
      dat_ld[is.na(R2_group), R2_group := "Unknown"]
    } else {
      dat_ld <- dat_sig
      dat_ld[, R2_group := "Unknown"]
    }
    
    rm(dat_sig)
    gc()
    
    df_plot_tmp <- dat_ld %>% dplyr::mutate(logP = -log10(P)) %>% 
      dplyr::mutate(logP = ifelse(is.infinite(logP), max(c(10, logP[is.finite(logP)]), na.rm=TRUE) + 5, logP))
    global_max_y <- max(df_plot_tmp$logP, na.rm = TRUE)
    shared_y_limit <- global_max_y * 1.1
    rm(df_plot_tmp)
    gc()
    
    pheno_out_dir <- file.path(out_dir, pheno)
    dir.create(pheno_out_dir, showWarnings = FALSE)
    
    plot_width <- 12
    plot_height <- 6
    plot_dpi <- 1500
    
    p_whole <- NULL
    p_mhc <- NULL
    p_spec <- NULL
    
    dat_plot_whole <- dat_ld[P < 0.05]
    if(nrow(dat_plot_whole) > 0) {
      p_whole <- plot_manhattan(dat_plot_whole, title = paste0("Whole Genome - ", pheno), use_ld = TRUE, custom_y_limit = shared_y_limit)
      ggsave(file.path(pheno_out_dir, paste0(pheno, "_whole_genome_LD_colored.png")), plot = p_whole, width = plot_width, height = plot_height, dpi = plot_dpi, device = ragg::agg_png)
      ggsave(file.path(pheno_out_dir, paste0(pheno, "_whole_genome_LD_colored.pdf")), plot = p_whole, width = plot_width, height = plot_height)
    }
    
    dat_mhc <- dat_ld[CHR == 6 & BP >= 25000000 & BP <= 35000000]
    if(nrow(dat_mhc) > 0) {
      p_mhc <- plot_manhattan(dat_mhc, title = paste0("MHC Region (chr6:25M-35M) - ", pheno), 
                              region = list(chr = 6, start = 25000000, end = 35000000), use_ld = TRUE, custom_y_limit = NULL)
      ggsave(file.path(pheno_out_dir, paste0(pheno, "_MHC_chr6_25M_35M_LD_colored.png")), plot = p_mhc, width = 9, height = plot_height, dpi = plot_dpi, device = ragg::agg_png)
      ggsave(file.path(pheno_out_dir, paste0(pheno, "_MHC_chr6_25M_35M_LD_colored.pdf")), plot = p_mhc, width = 9, height = plot_height)
    }
    
    dat_spec <- dat_ld[CHR == target_chr & BP >= region_start & BP <= region_end]
    if(nrow(dat_spec) > 0) {
      region_title <- sprintf("Specific Region (%s ± 1Mb) - %s", target_snp, pheno)
      p_spec <- plot_manhattan(dat_spec, title = region_title, 
                               region = list(chr = target_chr, start = region_start, end = region_end), use_ld = TRUE, custom_y_limit = NULL)
      ggsave(file.path(pheno_out_dir, paste0(pheno, "_specific_region_1Mb_LD_colored.png")), plot = p_spec, width = 9, height = plot_height, dpi = plot_dpi, device = ragg::agg_png)
      ggsave(file.path(pheno_out_dir, paste0(pheno, "_specific_region_1Mb_LD_colored.pdf")), plot = p_spec, width = 9, height = plot_height)
    }
    
    if (!is.null(p_whole) && !is.null(p_mhc) && !is.null(p_spec)) {
      p_label <- ggplot() + 
        annotate("text", x = 0, y = 0, label = disease_name, angle = 270, color = "black", size = 6, fontface = "bold") +
        coord_cartesian(xlim = c(-1, 1), ylim = c(-1, 1)) +
        theme_void() +
        theme(
          plot.background = element_blank(),
          panel.background = element_rect(fill = "#D3D3D3", color = NA),
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0)
        )
      
      p_mhc_shared <- plot_manhattan(dat_mhc, title = paste0("MHC Region (chr6:25M-35M) - ", pheno), 
                                     region = list(chr = 6, start = 25000000, end = 35000000), use_ld = TRUE, custom_y_limit = shared_y_limit)
      p_spec_shared <- plot_manhattan(dat_spec, title = region_title, 
                                      region = list(chr = target_chr, start = region_start, end = region_end), use_ld = TRUE, custom_y_limit = shared_y_limit)
      
      p_mhc_v1 <- p_mhc_shared + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none", plot.margin = margin(t = 20, r = 10, b = 20, l = 5)) +
        labs(title = "MHC Region")
      p_spec_v1 <- p_spec_shared + theme(axis.title.y = element_blank(), axis.text.y = element_blank(), axis.ticks.y = element_blank(), legend.position = "none", plot.margin = margin(t = 20, r = 2, b = 20, l = 5)) +
        labs(title = "Target Region")
      p_whole_v1 <- p_whole + theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 20))
      p_label_v1 <- p_label + theme(plot.margin = margin(t = 20, r = 20, b = 20, l = 0))
      
      p_combined_v1 <- p_whole_v1 + p_mhc_v1 + p_spec_v1 + p_label_v1 + plot_layout(ncol = 4, widths = c(2, 1, 1, 0.08))
      ggsave(file.path(pheno_out_dir, paste0(pheno, "_combined_v1_horizontal.png")), plot = p_combined_v1, width = 20.5, height = 6, dpi = plot_dpi, device = ragg::agg_png)
      ggsave(file.path(pheno_out_dir, paste0(pheno, "_combined_v1_horizontal.pdf")), plot = p_combined_v1, width = 20.5, height = 6)
      
      p_mhc_v2 <- p_mhc + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none", plot.margin = margin(t = 5, r = 2, b = 5, l = 5))
      p_spec_v2 <- p_spec + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none", plot.margin = margin(t = 5, r = 5, b = 5, l = 2))
      p_inset <- p_mhc_v2 + p_spec_v2 + plot_layout(ncol = 2) + plot_annotation(theme = theme(plot.background = element_rect(fill = "white", color = "black", linewidth = 1)))
      
      p_whole_v2 <- p_whole + theme(plot.margin = margin(t = 20, r = 2, b = 20, l = 20))
      p_label_v2 <- p_label + theme(plot.margin = margin(t = 20, r = 20, b = 20, l = 0))
      
      p_combined_v2_base <- p_whole_v2 + inset_element(p_inset, left = 0.55, bottom = 0.55, right = 0.99, top = 0.99, align_to = "panel")
      p_combined_v2 <- p_combined_v2_base + p_label_v2 + plot_layout(ncol = 2, widths = c(1, 0.03))
      ggsave(file.path(pheno_out_dir, paste0(pheno, "_combined_v2_inset.png")), plot = p_combined_v2, width = 15.5, height = 8, dpi = plot_dpi, device = ragg::agg_png)
      ggsave(file.path(pheno_out_dir, paste0(pheno, "_combined_v2_inset.pdf")), plot = p_combined_v2, width = 15.5, height = 8)
    }
    
    rm(dat_ld, dat_plot_whole, dat_mhc, dat_spec)
    if(exists("p_whole")) rm(p_whole)
    if(exists("p_mhc")) rm(p_mhc)
    if(exists("p_spec")) rm(p_spec)
    if(exists("p_label")) rm(p_label)
    if(exists("p_combined_v1")) rm(p_combined_v1)
    if(exists("p_combined_v2")) rm(p_combined_v2)
    if(exists("p_combined_v2_base")) rm(p_combined_v2_base)
    gc()
    
    return(list(pheno = pheno, status = "Success", error = NA, stats = stats_row))
  }, error = function(e) {
    return(list(pheno = pheno, status = "Failed", error = e$message, stats = NULL))
  })
}

num_plot_workers <- min(16, nrow(target_targets))
log_msg(sprintf("Starting PSOCK cluster with %d workers for plotting...", num_plot_workers))
cl_plot <- makeCluster(num_plot_workers)

clusterExport(cl_plot, c("process_and_plot_worker", "target_targets", "out_dir", "target_chr", "region_start", "region_end", "target_snp", "ld_map", "plot_manhattan", "ld_colors"))

parallel_results <- parLapply(cl_plot, 1:nrow(target_targets), function(i) {
  process_and_plot_worker(i, target_targets, out_dir, target_chr, region_start, region_end, target_snp, ld_map, plot_manhattan, ld_colors)
})

stopCluster(cl_plot)
log_msg("Parallel plotting tasks completed.")

stats_list <- list()
for (res in parallel_results) {
  if (!is.null(res$stats)) {
    stats_list[[length(stats_list) + 1]] <- res$stats
  }
  if (res$status == "Failed") {
    log_msg(sprintf("Task Failed for %s: %s", res$pheno, res$error))
  } else {
    log_msg(sprintf("Task Success for %s", res$pheno))
  }
}

# ------------------------------------------------------------------------------
# Step 7: Save Statistics and Finalize
# : 7: 
# ------------------------------------------------------------------------------
if(length(stats_list) > 0) {
  stats_df <- bind_rows(stats_list)
  write.csv(stats_df, file.path(out_dir, paste0("FinnGen_Manhattan_Stats_", timestamp, ".csv")), row.names = FALSE)
  log_msg("Statistics saved.")
}

end_time <- Sys.time()
run_time <- round(as.numeric(difftime(end_time, start_time, units = "mins")), 2)
log_msg(sprintf("Script completed successfully. Total run time: %s minutes.", run_time))
