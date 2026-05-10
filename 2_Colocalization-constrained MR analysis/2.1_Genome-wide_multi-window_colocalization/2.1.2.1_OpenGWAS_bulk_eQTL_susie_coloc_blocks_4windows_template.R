#!/usr/bin/env Rscript
# ==============================================================================
# [Script]: 2.1.2.1_OpenGWAS_bulk_eQTL_susie_coloc_blocks_4windows_template.R
# [Method]: SuSiE-based Colocalization 
# [Step]:  ( OpenGWAS Bulk eQTL)
# 
# [Function]:
# Conduct SuSiE-based colocalization analysis to identify shared causal variants
#       between Bulk eQTL and GWAS traits. Execute analysis in 4 windows per LD block:
#       1. Block-based, 2. Gene +/- 100kb, 3. Gene +/- 50kb, 4. Gene +/- 10kb.
#        LD block  4 : 
#       1.  Block, 2.  100kb, 3.  50kb, 4.  10kb. 
# 
# [Data Availability / ]:
# Requires preprocessed eQTL, preprocessed GWAS, Gene info with LD blocks, and UKB LD reference.
# ==============================================================================

suppressWarnings(suppressMessages({
  rm(list = ls()); gc()
}))

start_time <- Sys.time()

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(future)
  library(future.apply)
  library(RhpcBLASctl)
  library(stringr)
  library(tools)
  library(susieR)
  library(coloc)
  library(Matrix)
  library(ieugwasr)
  library(ggplot2)
  library(cowplot)
  library(patchwork)
  library(bigsnpr)
  library(bigreadr)
})

if (!requireNamespace("locuscomparer", quietly = TRUE)) {
  if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools", repos="https://cloud.r-project.org")
  devtools::install_github("caleblare/locuscomparer")
}
library(locuscomparer)

# Data is hg19, so prefer EnsDb.Hsapiens.v75
if (requireNamespace("EnsDb.Hsapiens.v75", quietly = TRUE)) {
  library(EnsDb.Hsapiens.v75)
  ens_db <- EnsDb.Hsapiens.v75
} else if (requireNamespace("EnsDb.Hsapiens.v86", quietly = TRUE)) {
  library(EnsDb.Hsapiens.v86)
  ens_db <- EnsDb.Hsapiens.v86
}

option_list <- list(
  make_option(c("--work_dir"), type="character", default=getwd(), help="Working directory"),
  make_option(c("--eqtl_dir"), type="character", default=NULL, help="Path to preprocessed eQTL directory"),
  make_option(c("--gwas_file"), type="character", default=NULL, help="Path to preprocessed GWAS file"),
  make_option(c("--gene_info"), type="character", default=NULL, help="Path to Gene info with LD blocks"),
  make_option(c("--ld_dir"), type="character", default=NULL, help="Path to merged UKB LD Reference Panel"),
  make_option(c("--block_ld_dir"), type="character", default=NULL, help="Path to block-based UKB LD Reference Panel"),
  make_option(c("--out_dir"), type="character", default=NULL, help="Output directory"),
  make_option(c("--tasks"), type="integer", default=11, help="Number of parallel workers"),
  make_option(c("--threads_per_task"), type="integer", default=2, help="Number of threads per task"),
  make_option(c("--resume"), type="logical", default=TRUE, help="Enable in-place resume mode"),
  make_option(c("--test_n"), type="integer", default=0, help="Test limit for number of genes (0 for all)"),
  make_option(c("--use_block_ld"), type="logical", default=TRUE, help="Use block-based LD panel"),
  make_option(c("--chunk_size"), type="integer", default=10, help="Genes per parallel chunk")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (any(sapply(list(opt$eqtl_dir, opt$gwas_file, opt$gene_info, opt$block_ld_dir), is.null))) {
  print_help(opt_parser)
  stop("Missing required arguments.")
}

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
BASE_OUT <- if (!is.null(opt$out_dir)) opt$out_dir else file.path(opt$work_dir, paste0("susie_coloc_4windows_", timestamp))

DIR_RESULTS <- file.path(BASE_OUT, "results")
DIR_LOGS <- file.path(BASE_OUT, "logs")
DIR_SCRIPTS <- file.path(BASE_OUT, "scripts")
DIR_TEMP <- file.path(BASE_OUT, "temp")

for (d in c(BASE_OUT, DIR_RESULTS, DIR_LOGS, DIR_SCRIPTS, DIR_TEMP)) {
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

LOG_FILE <- file.path(DIR_LOGS, paste0("run_log_", timestamp, ".txt"))
ERR_FILE <- file.path(DIR_LOGS, paste0("error_log_", timestamp, ".txt"))
STRONG_SIGNAL_LOG <- file.path(DIR_LOGS, paste0("strong_signals_realtime_", timestamp, ".txt"))
CHECKPOINT_FILE <- file.path(DIR_LOGS, "checkpoint_completed_tasks.txt")
SUMMARY_FILE <- file.path(BASE_OUT, "summary_coloc_results.csv")

log_msg <- function(...) {
  msg <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", paste(..., collapse = " "))
  cat(msg, "\n"); write(msg, file = LOG_FILE, append = TRUE)
}
log_err <- function(...) {
  msg <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", paste(..., collapse = " "))
  cat(msg, "\n"); write(msg, file = ERR_FILE, append = TRUE)
}

if (!file.exists(STRONG_SIGNAL_LOG)) {
  header <- paste("Timestamp", "Cell", "EnsID", "Gene", "Block", "Window", "SuSiE_H4", "ABF_H4", "LeadSNP", sep="\t")
  write(header, file = STRONG_SIGNAL_LOG)
}

get_completed_tasks <- function() {
  if (file.exists(CHECKPOINT_FILE)) {
    tryCatch({
      dt <- fread(CHECKPOINT_FILE, header = FALSE, sep = "\t", fill = TRUE, quote = "")
      if (ncol(dt) >= 2) return(unique(trimws(as.character(dt[[2]]))))
      return(character(0))
    }, error = function(e) character(0))
  } else character(0)
}

# ------------------------------------------------------------------------------
# Data Loading
# ------------------------------------------------------------------------------
log_msg("Reading gene info:", opt$gene_info)
gene_info <- fread(opt$gene_info)
if (!"cell" %in% colnames(gene_info)) gene_info[, cell := "Bulk"]
setnames(gene_info, c("gene_id", "chromosome"), c("ens_gene_id", "chr"), skip_absent=TRUE)

if (!"gene_name" %in% colnames(gene_info) && exists("ens_db")) {
  all_ens <- unique(gene_info$ens_gene_id)
  mapped_syms <- tryCatch(mapIds(ens_db, keys=all_ens, column="SYMBOL", keytype="GENEID", multiVals="first"), error=function(e) NULL)
  if (!is.null(mapped_syms)) {
    map_dt <- data.table(ens_gene_id = names(mapped_syms), gene_name = as.character(mapped_syms))
    gene_info <- merge(gene_info, map_dt, by="ens_gene_id", all.x=TRUE)
    gene_info[is.na(gene_name) | gene_name=="", gene_name := ens_gene_id]
  } else {
    gene_info[, gene_name := ens_gene_id]
  }
} else if (!"gene_name" %in% colnames(gene_info)) {
  gene_info[, gene_name := ens_gene_id]
}

filter_palindromic <- function(dt, a1_col = "A1", a2_col = "A2", verbose = FALSE) {
  if (is.null(dt) || nrow(dt) == 0) return(dt)
  a1 <- toupper(trimws(dt[[a1_col]])); a2 <- toupper(trimws(dt[[a2_col]]))
  is_palindromic <- (a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="G" & a2=="C") | (a1=="C" & a2=="G")
  is_palindromic[is.na(is_palindromic)] <- FALSE
  if (sum(is_palindromic) > 0) return(dt[!is_palindromic])
  return(dt)
}

log_msg("Reading GWAS file:", opt$gwas_file)
gwas_dt <- fread(opt$gwas_file)
setkey(gwas_dt, chr, pos)
gwas_dt <- filter_palindromic(gwas_dt, verbose = TRUE)

find_block_ld_file <- function(chr, start, end, block_ld_dir) {
  chr_dir <- file.path(block_ld_dir, paste0("chr", chr))
  if (!dir.exists(chr_dir)) return(NULL)
  block_files <- list.files(chr_dir, pattern = "\\.bed$", full.names = FALSE)
  block_files <- block_files[grepl(paste0("chr", chr, "_"), block_files)]
  if (length(block_files) == 0) return(NULL)
  
  block_info <- data.table(file = block_files)
  block_info[, start_pos := as.integer(sub(paste0(".*chr", chr, "_([0-9]+)_.*"), "\\1", file))]
  block_info[, end_pos := as.integer(sub(paste0(".*_", chr, "_([0-9]+)\\.bed"), "\\1", file))]
  overlapping <- block_info[(start_pos <= end) & (end_pos >= start)]
  if (nrow(overlapping) == 0) return(NULL)
  overlapping[, overlap := pmin(end_pos, end) - pmax(start_pos, start)]
  best_block <- overlapping[which.max(overlap)]
  return(file.path(chr_dir, sub("\\.bed$", "", best_block$file)))
}

compute_ld_mem <- function(snps, ld_res, chr, start, end, block_ld_dir) {
  if (ld_res$type != "split") return(NULL)
  block_path <- find_block_ld_file(chr, start, end, block_ld_dir)
  if (is.null(block_path) || !block_path %in% names(ld_res$cache)) return(NULL)
  
  bigsnp <- ld_res$cache[[block_path]]
  map <- bigsnp$map
  idx <- match(snps, map$marker.ID)
  valid <- which(!is.na(idx))
  if (length(valid) < 2) return(NULL)
  
  matched_idx <- idx[valid]
  matched_snps <- snps[valid]
  corr_sparse <- bigsnpr::snp_cor(bigsnp$genotypes, ind.col = matched_idx, size = length(matched_idx))
  R <- as.matrix(corr_sparse)
  rownames(R) <- colnames(R) <- matched_snps
  return(R)
}

run_susie <- function(dataset, R) {
  z <- dataset$beta / sqrt(dataset$varbeta)
  z[!is.finite(z)] <- 0
  n <- suppressWarnings(as.integer(stats::median(dataset$N, na.rm = TRUE)))
  if (is.na(n) || n <= 0) n <- length(z)
  
  snp_order <- intersect(dataset$snp, rownames(R))
  if (length(snp_order) == 0) return(NULL)
  R <- R[snp_order, snp_order, drop = FALSE]
  z <- z[match(snp_order, dataset$snp)]
  
  sus <- tryCatch(susieR::susie_rss(z = z, R = R, n = n, L = 10, coverage = 0.95, max_iter = 200, check_input = FALSE),
                  error = function(e) tryCatch(susieR::susie_rss(z = z, LD = R, n = n, L = 10, coverage = 0.95, max_iter = 200), error=function(e2) NULL))
  
  if (is.null(sus)) return(NULL)
  sus <- coloc::annotate_susie(sus, snp_order, LD = R)
  sus$snp_names <- snp_order
  return(sus)
}

create_locuscompare_plot <- function(eqtl, gwas, gene_id, gene_name, chr, start, end, out_dir, file_prefix, lead_snp=NULL) {
  e_df <- data.frame(rsid=eqtl$snp, pval=eqtl$pval)
  g_df <- data.frame(rsid=gwas$snp, pval=gwas$pval)
  e_df <- e_df[!is.na(e_df$pval) & e_df$pval > 0, ]
  g_df <- g_df[!is.na(g_df$pval) & g_df$pval > 0, ]
  if (nrow(e_df) == 0 || nrow(g_df) == 0) return(NULL)
  
  temp_d <- file.path(tempdir(), paste0("lc_", gene_id, "_", sample(1:10000,1)))
  dir.create(temp_d, showWarnings=FALSE)
  fn1 <- file.path(temp_d, "eqtl.tsv"); fn2 <- file.path(temp_d, "gwas.tsv")
  write.table(e_df, fn1, sep="\t", quote=FALSE, row.names=FALSE)
  write.table(g_df, fn2, sep="\t", quote=FALSE, row.names=FALSE)
  
  plot_file <- file.path(out_dir, paste0(file_prefix, "_locuscompare.png"))
  tryCatch({
    p <- locuscomparer::locuscompare(in_fn1=fn2, in_fn2=fn1, title1="GWAS", title2=paste0(gene_name, " eQTL"), snp=lead_snp, genome="hg19")
    if (inherits(p, "ggplot")) ggsave(plot_file, p, width=12, height=10, bg="white")
  }, error=function(e) NULL, finally=unlink(temp_d, recursive=TRUE))
  return(plot_file)
}

execute_window_analysis <- function(tsk, eqtl, ld_res, gwas_data, super_data = NULL) {
  RhpcBLASctl::blas_set_num_threads(opt$threads_per_task)
  ens <- tsk$ens; cell <- tsk$cell; chr <- tsk$chr; s <- tsk$start; e <- tsk$end
  bid <- tsk$block_id; win <- tsk$window_tag
  task_id <- paste(cell, ens, bid, win, sep="_")
  
  gene_name <- if ("gene_name" %in% colnames(gene_info)) gene_info[ens_gene_id == ens, gene_name][1] else ens
  if (is.na(gene_name)) gene_name <- ens
  
  out_gene_dir <- file.path(DIR_RESULTS, cell, ens)
  if (!dir.exists(out_gene_dir)) dir.create(out_gene_dir, recursive = TRUE)
  out_prefix <- paste0(cell, "_", gene_name, "_block", bid, "_", win)
  full_out_prefix <- file.path(out_gene_dir, out_prefix)
  
  merged <- NULL; R <- NULL
  if (!is.null(super_data)) {
    merged <- super_data$merged[pos_eqtl >= s & pos_eqtl <= e]
    if (nrow(merged) < 10) return(list(task_id=task_id, status="insufficient_data"))
    common_snps <- intersect(merged$snp, rownames(super_data$R))
    if (length(common_snps) < 10) return(list(task_id=task_id, status="low_overlap"))
    merged <- merged[snp %in% common_snps]
    R <- super_data$R[common_snps, common_snps]
  } else {
    eqtl_sub <- filter_palindromic(eqtl[chr == tsk$chr & pos >= s & pos <= e])
    gwas_sub <- filter_palindromic(gwas_data[chr == tsk$chr & pos >= s & pos <= e])
    if ("pval" %in% names(gwas_sub)) gwas_sub <- gwas_sub[order(pval)]
    gwas_sub <- unique(gwas_sub, by = "snp")
    if (nrow(eqtl_sub) < 10 || nrow(gwas_sub) < 10) return(list(task_id=task_id, status="insufficient_data"))
    
    merged <- merge(eqtl_sub, gwas_sub, by="snp", suffixes=c("_eqtl", "_gwas"))
    flip <- (merged$A1_eqtl == merged$A2_gwas & merged$A2_eqtl == merged$A1_gwas)
    merged$beta_gwas[flip] <- -merged$beta_gwas[flip]
    merged <- merged[(merged$A1_eqtl == merged$A1_gwas & merged$A2_eqtl == merged$A2_gwas) | flip]
    if (nrow(merged) < 10) return(list(task_id=task_id, status="low_overlap"))
    
    R <- compute_ld_mem(unique(merged$snp), ld_res, chr, s, e, opt$block_ld_dir)
    if (is.null(R)) return(list(task_id=task_id, status="ld_fail"))
    common <- intersect(merged$snp, rownames(R))
    merged <- merged[snp %in% common]
    R <- R[common, common]
  }
  
  d_eqtl <- list(beta=merged$beta_eqtl, varbeta=merged$se_eqtl^2, MAF=pmin(merged$eaf_eqtl, 1-merged$eaf_eqtl), N=merged$N_eqtl, type="quant", snp=merged$snp)
  d_gwas <- list(beta=merged$beta_gwas, varbeta=merged$se_gwas^2, MAF=pmin(merged$eaf_gwas, 1-merged$eaf_gwas), N=merged$N_gwas, type="cc", snp=merged$snp)
  
  fwrite(merged[, .(snp, chr=chr_eqtl, pos=pos_eqtl, A1=A1_eqtl, A2=A2_eqtl, beta_eqtl, se_eqtl, pval_eqtl, eaf_eqtl, cell, ens_gene_id=ens, gene_name, block_id=bid, window=win)], paste0(full_out_prefix, "_harmonized_eqtl.csv"))
  fwrite(merged[, .(snp, chr=chr_gwas, pos=pos_gwas, A1=A1_gwas, A2=A2_gwas, beta_gwas, se_gwas, pval_gwas, eaf_gwas, cell, ens_gene_id=ens, gene_name, block_id=bid, window=win)], paste0(full_out_prefix, "_harmonized_gwas.csv"))
  
  coloc_abf <- tryCatch(coloc::coloc.abf(d_eqtl, d_gwas), error=function(e) NULL)
  pp_h4_abf <- if(!is.null(coloc_abf)) coloc_abf$summary["PP.H4.abf"] else NA
  
  sus_eqtl <- tryCatch(run_susie(d_eqtl, R), error=function(e) NULL)
  sus_gwas <- tryCatch(run_susie(d_gwas, R), error=function(e) NULL)
  res_susie <- if (!is.null(sus_eqtl) && !is.null(sus_gwas)) tryCatch(coloc::coloc.susie(sus_eqtl, sus_gwas), error=function(e) NULL) else NULL
  pp_h4_susie_max <- if (!is.null(res_susie) && !is.null(res_susie$summary)) max(res_susie$summary$PP.H4.abf, na.rm=TRUE) else NA
  
  lead_snp <- if (!is.na(pp_h4_susie_max)) as.character(res_susie$summary$hit1[which.max(res_susie$summary$PP.H4.abf)]) else
              if (!is.na(pp_h4_abf) && pp_h4_abf > 0.8) as.character(coloc_abf$results$snp[which.max(coloc_abf$results$SNP.PP.H4)]) else
              merged$snp[which.min(merged$pval_eqtl * merged$pval_gwas)]
              
  if (!is.null(sus_eqtl$sets$cs)) fwrite(as.data.table(sus_eqtl$sets$cs)[, `:=`(cell=cell, ens_gene_id=ens, gene_name=gene_name, block_id=bid, window=win)], paste0(full_out_prefix, "_susie_eqtl_cs.csv"))
  if (!is.null(sus_gwas$sets$cs)) fwrite(as.data.table(sus_gwas$sets$cs)[, `:=`(cell=cell, ens_gene_id=ens, gene_name=gene_name, block_id=bid, window=win)], paste0(full_out_prefix, "_susie_gwas_cs.csv"))
  if (!is.null(res_susie$summary)) fwrite(as.data.table(res_susie$summary)[, `:=`(cell=cell, ens_gene_id=ens, gene_name=gene_name, block_id=bid, window=win, chr=chr, start=s, end=e, n_snps=nrow(merged), lead_snp=lead_snp)], paste0(full_out_prefix, "_coloc_summary.csv"))
  if (!is.null(coloc_abf)) fwrite(as.data.table(t(coloc_abf$summary))[, `:=`(cell=cell, ens_gene_id=ens, gene_name=gene_name, block_id=bid, window=win, chr=chr, start=s, end=e, lead_snp=lead_snp)], paste0(full_out_prefix, "_coloc_abf.csv"))
  
  if ((!is.na(pp_h4_susie_max) && pp_h4_susie_max > 0.8) || (!is.na(pp_h4_abf) && pp_h4_abf > 0.8)) {
    cat(paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), cell, ens, gene_name, bid, win, round(pp_h4_susie_max,4), round(pp_h4_abf,4), lead_snp, sep="\t"), "\n", file=STRONG_SIGNAL_LOG, append=TRUE)
    create_locuscompare_plot(merged[,.(snp, pval=pval_eqtl)], merged[,.(snp, pval=pval_gwas)], ens, gene_name, chr, s, e, out_gene_dir, out_prefix, lead_snp)
  }
  
  return(list(task_id=task_id, status="ok", pp_h4=pp_h4_susie_max, pp_h4_abf=pp_h4_abf))
}

process_block_chunk <- function(chunk, gwas_chr_data) {
  RhpcBLASctl::omp_set_num_threads(opt$threads_per_task)
  if (!exists("Global_LD_Cache")) Global_LD_Cache <<- list(block_key = NULL, ld_res = NULL)
  
  block_key <- chunk$block_key
  if (is.null(Global_LD_Cache$block_key) || Global_LD_Cache$block_key != block_key) {
    parts <- strsplit(block_key, "_")[[1]]
    block_file <- find_block_ld_file(parts[1], as.integer(parts[2]), as.integer(parts[3]), opt$block_ld_dir)
    if (!is.null(block_file) && file.exists(paste0(block_file, ".rds"))) {
      Global_LD_Cache <<- list(block_key = block_key, ld_res = list(type="split", cache=setNames(list(bigsnpr::snp_attach(paste0(block_file, ".rds"))), block_file)))
    } else Global_LD_Cache <<- list(block_key = NULL, ld_res = NULL)
  }
  
  chunk_results <- list()
  for (tsk_info in chunk$tasks) {
    eqtl_file <- file.path(opt$eqtl_dir, paste0(tsk_info$ens, ".csv"))
    if (!file.exists(eqtl_file)) eqtl_file <- file.path(opt$eqtl_dir, tsk_info$cell, paste0(tsk_info$ens, ".csv"))
    eqtl_dt <- if (file.exists(eqtl_file)) tryCatch(fread(eqtl_file), error=function(e) NULL) else NULL
    if (is.null(eqtl_dt)) next
    
    s_min <- min(sapply(tsk_info$windows, function(w) w$start))
    e_max <- max(sapply(tsk_info$windows, function(w) w$end))
    super_data <- NULL
    tryCatch({
      eqtl_sup <- eqtl_dt[chr == tsk_info$chr & pos >= s_min & pos <= e_max]
      gwas_sup <- gwas_chr_data[chr == tsk_info$chr & pos >= s_min & pos <= e_max]
      if (nrow(eqtl_sup) >= 10 && nrow(gwas_sup) >= 10) {
        merged <- merge(eqtl_sup, unique(gwas_sup[order(pval)], by="snp"), by="snp", suffixes=c("_eqtl", "_gwas"))
        flip <- (merged$A1_eqtl == merged$A2_gwas & merged$A2_eqtl == merged$A1_gwas)
        merged$beta_gwas[flip] <- -merged$beta_gwas[flip]
        merged <- merged[(merged$A1_eqtl == merged$A1_gwas & merged$A2_eqtl == merged$A2_gwas) | flip]
        if (nrow(merged) >= 10 && !is.null(Global_LD_Cache$ld_res)) {
          R <- compute_ld_mem(unique(merged$snp), Global_LD_Cache$ld_res, tsk_info$chr, s_min, e_max, opt$block_ld_dir)
          if (!is.null(R)) {
            common <- intersect(merged$snp, rownames(R))
            super_data <- list(merged=merged[snp %in% common], R=R[common, common])
          }
        }
      }
    }, error=function(e) NULL)
    
    for (w in tsk_info$windows) {
      res <- tryCatch(execute_window_analysis(list(ens=tsk_info$ens, cell=tsk_info$cell, chr=tsk_info$chr, start=w$start, end=w$end, block_id=w$block_id, window_tag=w$window_tag), eqtl_dt, Global_LD_Cache$ld_res, gwas_chr_data, super_data), error=function(e) list(status="error"))
      chunk_results[[length(chunk_results)+1]] <- res
      if (res$status == "ok") cat(paste0("0\t", res$task_id, "\n"), file=CHECKPOINT_FILE, append=TRUE)
    }
  }
  return(chunk_results)
}

completed_tasks <- get_completed_tasks()
log_msg("Completed tasks:", length(completed_tasks))

plan(multisession, workers = opt$tasks)
unique_chrs <- sort(unique(gene_info$chr))
all_stats <- list()
task_count <- 0

for (chrom in unique_chrs) {
  chr_genes <- gene_info[chr == chrom]
  gwas_chr_subset <- gwas_dt[chr == chrom]
  tasks_flat <- list()
  
  for (i in seq_len(nrow(chr_genes))) {
    row <- chr_genes[i]
    for (b in 1:row$num_ld_blocks) {
      bs <- row[[paste0("block", b, "_start")]]; be <- row[[paste0("block", b, "_end")]]
      if (is.na(bs)) next
      pending_windows <- list()
      for (wd in list(list(tag="block", s=bs, e=be), list(tag="100kb", s=max(bs, row$start_hg19-1e5), e=min(be, row$end_hg19+1e5)), list(tag="50kb", s=max(bs, row$start_hg19-5e4), e=min(be, row$end_hg19+5e4)), list(tag="10kb", s=max(bs, row$start_hg19-1e4), e=min(be, row$end_hg19+1e4)))) {
        if (!paste(row$cell, row$ens_gene_id, b, wd$tag, sep="_") %in% completed_tasks) pending_windows[[length(pending_windows)+1]] <- list(block_id=b, window_tag=wd$tag, start=wd$s, end=wd$e)
      }
      if (length(pending_windows) > 0) tasks_flat[[length(tasks_flat)+1]] <- list(block_key=paste0(chrom,"_",bs,"_",be), ens=row$ens_gene_id, cell=row$cell, chr=chrom, windows=pending_windows)
    }
  }
  
  if (length(tasks_flat) == 0) next
  tasks_dt <- rbindlist(lapply(tasks_flat, function(x) data.table(block_key=x$block_key, idx=1L)), idcol="list_idx")
  tasks_dt[, count := .N, by=block_key]
  
  job_chunks <- list()
  for (bk in unique(tasks_dt[order(-count)]$block_key)) {
    block_tasks <- tasks_flat[tasks_dt[block_key == bk, list_idx]]
    for (k in 1:ceiling(length(block_tasks)/opt$chunk_size)) job_chunks[[length(job_chunks)+1]] <- list(chunk_id=paste0(bk,"_p",k), block_key=bk, tasks=block_tasks[((k-1)*opt$chunk_size+1):min(k*opt$chunk_size, length(block_tasks))])
  }
  
  if (opt$test_n > 0) job_chunks <- head(job_chunks, ceiling(opt$test_n / opt$chunk_size))
  
  res_chr <- future_lapply(job_chunks, function(chunk) {
    library(data.table); library(susieR); library(coloc); library(bigsnpr); if(requireNamespace("locuscomparer", quietly=TRUE)) library(locuscomparer)
    process_block_chunk(chunk, gwas_chr_subset)
  }, future.seed=TRUE)
  
  all_stats <- c(all_stats, unlist(res_chr, recursive=FALSE))
  task_count <- task_count + sum(sapply(job_chunks, function(x) length(x$tasks)))
  if (opt$test_n > 0 && task_count >= opt$test_n) break
}

fwrite(rbindlist(all_stats, fill=TRUE), SUMMARY_FILE)
log_msg("Pipeline completed in", round(as.numeric(difftime(Sys.time(), start_time, units="mins")), 2), "minutes.")
