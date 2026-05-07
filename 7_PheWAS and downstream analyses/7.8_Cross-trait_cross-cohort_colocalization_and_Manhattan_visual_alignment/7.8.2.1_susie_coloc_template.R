#!/usr/bin/env Rscript
# ==============================================================================
# [Script]: 7.8.2.1_susie_coloc_template.R
# [Method]: SuSiE Colocalization Analysis
# [Step]: Run pairwise SuSiE coloc
# 
# [Function]:
# Run SuSiE-based colocalization analysis to identify shared causal variants.
# 
# [Input]:
#   --dataset1_dir : Directory containing formatted dataset 1
#   --dataset2_dir : Directory containing formatted dataset 2
#   --out_dir      : Output directory path
#   --ld_matrix    : Path to local LD matrix or reference panel
# 
# [Output]:
# SuSiE colocalization results (H4 posteriors, shared variants).
# ==============================================================================

rm(list = ls())
gc()

# 1. Load Dependencies
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(coloc)
  library(susieR)
})

# 2. Parse Command Line Arguments
option_list <- list(
  make_option(c("--dataset1_dir"), type="character", default=NULL, help="Directory for Dataset 1"),
  make_option(c("--dataset2_dir"), type="character", default=NULL, help="Directory for Dataset 2"),
  make_option(c("--out_dir"), type="character", default="./Coloc_Results", help="Output directory path"),
  make_option(c("--ld_matrix"), type="character", default=NULL, help="Path to LD reference panel")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$dataset1_dir) || is.null(opt$dataset2_dir)) {
  print_help(opt_parser)
  stop("Missing required input arguments.")
}

DATASET1_DIR <- opt$dataset1_dir
DATASET2_DIR <- opt$dataset2_dir
OUTPUT_DIR <- opt$out_dir
LD_MATRIX <- opt$ld_matrix

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

}))

start_time <- Sys.time()
timestamp <- format(start_time, "%Y%m%d_%H%M%S")

## 1. Environment & Paths Configuration / 
WORK_DIR <- "./"
BASE_OUT <- file.path(WORK_DIR, paste0("1_Pairwise_Coloc_rs1800628_50kb_", timestamp))

# Tools & Ref Panel
PLINK19_BIN <- "./"
LD_DIR <- "./"

# Target Region /  (rs1800628 in GRCh37)
TARGET_CHR <- 6
TARGET_POS <- 31546850
TARGET_SNP <- "rs1800628"
WINDOW_SIZES <- c(50000, 10000)
MAX_WINDOW <- max(WINDOW_SIZES)

# Input List
LIST_PATH <- "./"

# Create Output Dirs
DIR_RESULTS <- file.path(BASE_OUT, "results")
DIR_LOGS <- file.path(BASE_OUT, "logs")
DIR_README <- file.path(BASE_OUT, "readme")
DIR_SCRIPTS <- file.path(BASE_OUT, "scripts")
DIR_TEMP <- file.path(BASE_OUT, "temp")
for (d in c(BASE_OUT, DIR_RESULTS, DIR_LOGS, DIR_README, DIR_SCRIPTS, DIR_TEMP)) {
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
}

# Logger
LOG_FILE <- file.path(DIR_LOGS, paste0("run_log_", timestamp, ".txt"))
log_msg <- function(...) {
  msg <- paste(format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "-", paste(..., collapse = " "))
  cat(msg, "\n")
  write(msg, file = LOG_FILE, append = TRUE)
}

log_msg("Started Pairwise Colocalization Analysis Script.")
log_msg(paste("Target Region:", TARGET_CHR, ":", TARGET_POS, "with windows", paste(WINDOW_SIZES/1000, "kb", collapse=", ")))

## 2. Load Libraries / 
safe_pkg <- function(pkg, github = NULL) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    if (is.null(github)) install.packages(pkg, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/") else {
      if (!requireNamespace("remotes")) install.packages("remotes")
      remotes::install_github(github)
    }
  }
  suppressPackageStartupMessages(library(pkg, character.only = TRUE))
}

for (p in c("data.table", "future", "future.apply", "stringr", "Matrix")) safe_pkg(p)
safe_pkg("susieR")
safe_pkg("coloc")
safe_pkg("ggplot2")
safe_pkg("patchwork")

if (!requireNamespace("locuscomparer", quietly = TRUE)) {
  log_msg("Installing locuscomparer...")
  if (!requireNamespace("devtools")) install.packages("devtools")
  devtools::install_github("caleblare/locuscomparer")
}
library(locuscomparer)

## 3. Merge Lists / 
log_msg("Reading dataset list...")
unified <- fread(LIST_PATH)

valid_idx <- sapply(unified$FilePath, file.exists)
unified <- unified[valid_idx]
# Clean phenotype names: Remove "(UKB data field ...)" suffixes
unified$Phenotype <- trimws(gsub("\\s*\\(UKB data field [^)]+\\)\\s*", "", unified$Phenotype, ignore.case = TRUE))
# Remove "protein levels" for plasma proteins
unified$Phenotype <- trimws(gsub("protein levels", "", unified$Phenotype, ignore.case = TRUE))
# Standardize specific phenotype names
unified$Phenotype <- trimws(gsub("(?i)\\s*\\(FinnGen\\)", "", unified$Phenotype))
unified$Phenotype <- trimws(gsub("(?i)\\s*\\(FG\\)", "", unified$Phenotype))
unified$Phenotype[unified$Phenotype == "Dermatopolymyositis"] <- "Dermatopolymyositis"
unified$Phenotype[unified$ID == "GCST90474052" | grepl("(?i)^systemic lupus erythematosus", unified$Phenotype) | grepl("(?i)^SLE", unified$Phenotype)] <- "SLE"
log_msg(paste("Total valid datasets after merge:", nrow(unified)))

fwrite(unified, file.path(DIR_README, paste0("Unified_Datasets_List_", timestamp, ".csv")))

## 4. Generate Pairwise Tasks / 
  # NOTE: Set TEST_MODE to TRUE for test analysis
  TEST_MODE <- FALSE
  if (TEST_MODE) {
    log_msg("TEST MODE ON: Selecting first 5 datasets.")
    unified <- unified[1:5]
  }

clean_name <- function(x) {
  x <- trimws(x)
  x <- gsub("[^A-Za-z0-9]+", "_", x)
  x <- sub("_$", "", x)
  return(x)
}

n_ds <- nrow(unified)
tasks <- list()
for (i in 1:(n_ds - 1)) {
  for (j in (i + 1):n_ds) {
    t_id <- paste0(unified$ID[i], "_", clean_name(unified$Phenotype[i]), "_", unified$Source[i], "_vs_", 
                   unified$ID[j], "_", clean_name(unified$Phenotype[j]), "_", unified$Source[j])
    tasks[[length(tasks) + 1]] <- list(
      task_id = t_id,
      ds1 = unified[i, ],
      ds2 = unified[j, ]
    )
  }
}
log_msg(paste("Generated", length(tasks), "pairwise colocalization tasks."))

task_df <- rbindlist(lapply(tasks, function(t) {
  data.table(Task_ID = t$task_id, 
             ID1 = t$ds1$ID, Phenotype1 = t$ds1$Phenotype, Source1 = t$ds1$Source,
             ID2 = t$ds2$ID, Phenotype2 = t$ds2$Phenotype, Source2 = t$ds2$Source)
}))
fwrite(task_df, file.path(DIR_README, paste0("Pairwise_Tasks_List_", timestamp, ".csv")))

## 5. Helper Functions / 
filter_palindromic <- function(dt, a1_col = "A1", a2_col = "A2") {
  if (is.null(dt) || nrow(dt) == 0) return(dt)
  a1 <- toupper(trimws(dt[[a1_col]]))
  a2 <- toupper(trimws(dt[[a2_col]]))
  is_pal <- (a1 == "A" & a2 == "T") | (a1 == "T" & a2 == "A") | (a1 == "G" & a2 == "C") | (a1 == "C" & a2 == "G")
  is_pal[is.na(is_pal)] <- FALSE
  return(dt[!is_pal])
}

compute_ld_plink <- function(snps, chr, ld_dir, temp_dir) {
  bfile_prefix <- file.path(ld_dir, "UKB_LD_Reference_Panel_EUR_Merged")
  if (!file.exists(paste0(bfile_prefix, ".bed"))) {
     bfile_prefix <- file.path(ld_dir, paste0("chr", chr))
  }
  if (!file.exists(paste0(bfile_prefix, ".bed"))) return(NULL)
  
  tmp_d <- file.path(temp_dir, paste0("ld_", Sys.getpid(), "_", sample(1:100000, 1)))
  dir.create(tmp_d, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(tmp_d, recursive = TRUE))
  
  snps_file <- file.path(tmp_d, "snps.txt")
  writeLines(snps, snps_file)
  out_prefix <- file.path(tmp_d, "ld_out")
  
  cmd <- paste(
    shQuote(PLINK19_BIN),
    "--bfile", shQuote(bfile_prefix),
    "--extract", shQuote(snps_file),
    "--r", "square", "spaces",
    "--write-snplist",
    "--out", shQuote(out_prefix),
    "--memory", 2000,
    "--threads", 1
  )
  system(cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  ld_file <- paste0(out_prefix, ".ld")
  snp_list_file <- paste0(out_prefix, ".snplist")
  if (file.exists(ld_file) && file.exists(snp_list_file)) {
      actual_snps <- readLines(snp_list_file)
      R <- tryCatch(as.matrix(fread(ld_file, header = FALSE)), error=function(e) NULL)
      if (is.null(R) || nrow(R) == 0) return(NULL)
      rownames(R) <- actual_snps
      colnames(R) <- actual_snps
      R[is.na(R)] <- 0
      return(R)
  }
  return(NULL)
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
  
  sus <- tryCatch({
    susieR::susie_rss(z = z, R = R, n = n, L = 10, coverage = 0.95, max_iter = 200, check_input = FALSE)
  }, error = function(e) {
    tryCatch({
       susieR::susie_rss(z = z, LD = R, n = n, L = 10, coverage = 0.95, max_iter = 200)
    }, error = function(e2) NULL)
  })
  
  if (is.null(sus)) return(NULL)
  sus <- coloc::annotate_susie(sus, snp_order, LD = R)
  sus$snp_names <- snp_order
  return(sus)
}

create_locuscompare_plot <- function(d1, d2, name1, name2, lead_snp, out_file) {
  suppressPackageStartupMessages(library(locuscomparer))
  suppressPackageStartupMessages(library(ggplot2))
  suppressPackageStartupMessages(library(cowplot))
  
  # Override locuscomparer colors
  my_assign_color <- function (rsid, snp, ld) {
      ld = ld[ld$SNP_A == snp, ]
      ld$color = as.character(cut(ld$R2, breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), 
          labels = c("#118AB2", "#88D49E", "#FCBF49", "#F77F00", "#D62828"), include.lowest = TRUE))
      color = data.frame(rsid, stringsAsFactors = FALSE)
      color = merge(color, ld[, c("SNP_B", "color")], by.x = "rsid", 
          by.y = "SNP_B", all.x = TRUE)
      color[is.na(color$color), "color"] = "#118AB2"
      if (snp %in% color$rsid) {
          color[rsid == snp, "color"] = "purple"
      } else {
          color = rbind(color, data.frame(rsid = snp, color = "purple"))
      }
      res = color$color
      names(res) = color$rsid
      return(res)
  }
  
  my_make_scatterplot <- function (merged, title1, title2, color, shape, size, legend = TRUE, 
      legend_position = c("bottomright", "topright", "topleft")) {
      p = ggplot(merged, aes(logp1, logp2)) + geom_point(aes(fill = rsid, 
          size = rsid, shape = rsid), alpha = 0.8) + geom_point(data = merged[merged$label != 
          "", ], aes(logp1, logp2, fill = rsid, size = rsid, shape = rsid)) + 
          xlab(bquote(.(title1) ~ -log[10] * "(P)")) + ylab(bquote(.(title2) ~ 
          -log[10] * "(P)")) + scale_fill_manual(values = color, 
          guide = "none") + scale_shape_manual(values = shape, 
          guide = "none") + scale_size_manual(values = size, guide = "none") + 
          ggrepel::geom_text_repel(aes(label = label)) + theme_classic()
      if (legend == TRUE) {
          legend_position = match.arg(legend_position)
          if (legend_position == "bottomright") {
              legend_box = data.frame(x = 0.8, y = seq(0.4, 0.2, -0.05))
          } else if (legend_position == "topright") {
              legend_box = data.frame(x = 0.8, y = seq(0.8, 0.6, -0.05))
          } else {
              legend_box = data.frame(x = 0.2, y = seq(0.8, 0.6, -0.05))
          }
          fill_colors = rev(c("#118AB2", "#88D49E", "#FCBF49", "#F77F00", "#D62828"))
          p = ggdraw(p) + geom_rect(data = legend_box, aes(xmin = x, 
              xmax = x + 0.05, ymin = y, ymax = y + 0.05), color = "black", 
              fill = fill_colors) + 
              draw_label("0.8", x = legend_box$x[1] + 0.05, y = legend_box$y[1], hjust = -0.3, size = 10) + 
              draw_label("0.6", x = legend_box$x[2] + 0.05, y = legend_box$y[2], hjust = -0.3, size = 10) + 
              draw_label("0.4", x = legend_box$x[3] + 0.05, y = legend_box$y[3], hjust = -0.3, size = 10) + 
              draw_label("0.2", x = legend_box$x[4] + 0.05, y = legend_box$y[4], hjust = -0.3, size = 10) + 
              draw_label(parse(text = "r^2"), x = legend_box$x[1] + 0.05, y = legend_box$y[1], vjust = -2, size = 10)
      }
      return(p)
  }
  
  assignInNamespace("assign_color", my_assign_color, ns="locuscomparer")
  assignInNamespace("make_scatterplot", my_make_scatterplot, ns="locuscomparer")

  df1 <- data.frame(rsid=d1$snp, pval=d1$pval)
  df2 <- data.frame(rsid=d2$snp, pval=d2$pval)
  df1 <- df1[!is.na(df1$pval) & df1$pval > 0, ]
  df2 <- df2[!is.na(df2$pval) & df2$pval > 0, ]
  if (nrow(df1) == 0 || nrow(df2) == 0) return(NULL)
  
  tmp_d <- file.path(tempdir(), paste0("lc_", sample(1:10000,1)))
  dir.create(tmp_d, showWarnings=FALSE)
  on.exit(unlink(tmp_d, recursive=TRUE))
  fn1 <- file.path(tmp_d, "d1.tsv"); fn2 <- file.path(tmp_d, "d2.tsv")
  write.table(df1, fn1, sep="\t", quote=FALSE, row.names=FALSE)
  write.table(df2, fn2, sep="\t", quote=FALSE, row.names=FALSE)
  
  tryCatch({
    p <- locuscomparer::locuscompare(
      in_fn1 = fn1, in_fn2 = fn2, 
      title1 = name1, title2 = name2,
      snp = lead_snp, genome = "hg19"
    )
    if (inherits(p, "ggplot")) {
       #  LZW  tiff 
       ggsave(out_file, p, width = 7.2, height = 6.5, dpi = 1500, bg = "white", compression = "lzw")
       return(out_file)
    } else {
       write(paste("Plot is not ggplot. Class:", paste(class(p), collapse=", ")), file=paste0(out_file, "_not_ggplot.log"))
    }
  }, error=function(e) {
       write(paste("LocusCompare error:", e$message), file=paste0(out_file, "_error.log"))
       return(NULL)
  })
  return(NULL)
}

## 6. Execution Engine / 
process_task <- function(tsk) {
  tid <- tsk$task_id
  ds1 <- tsk$ds1
  ds2 <- tsk$ds2
  
  out_task_dir <- file.path(DIR_RESULTS, tid)
  dir.create(out_task_dir, recursive = TRUE, showWarnings = FALSE)
  
  # Read Data (Maximum window)
  dt1_full <- fread(ds1$FilePath)[chr == TARGET_CHR & pos >= (TARGET_POS - MAX_WINDOW) & pos <= (TARGET_POS + MAX_WINDOW)]
  dt2_full <- fread(ds2$FilePath)[chr == TARGET_CHR & pos >= (TARGET_POS - MAX_WINDOW) & pos <= (TARGET_POS + MAX_WINDOW)]
  
  dt1_full <- filter_palindromic(dt1_full)
  dt2_full <- filter_palindromic(dt2_full)
  
  summary_list <- list()
  details_list <- list()
  
  for (win_size in WINDOW_SIZES) {
      win_kb <- win_size / 1000
      out_prefix <- file.path(out_task_dir, paste0(tid, "_", win_kb, "kb"))
      
      dt1 <- dt1_full[pos >= (TARGET_POS - win_size) & pos <= (TARGET_POS + win_size)]
      dt2 <- dt2_full[pos >= (TARGET_POS - win_size) & pos <= (TARGET_POS + win_size)]
      
      if (nrow(dt1) < 10 || nrow(dt2) < 10) {
          summary_list[[paste0(win_kb)]] <- data.table(Task_ID = tid, Window = paste0(win_kb, "kb"), status = "insufficient_data")
          next
      }
      
      dt1 <- unique(dt1[order(pval)], by="snp")
      dt2 <- unique(dt2[order(pval)], by="snp")
      
      merged <- merge(dt1, dt2, by="snp", suffixes=c("_1", "_2"))
      if (nrow(merged) < 10) {
          summary_list[[paste0(win_kb)]] <- data.table(Task_ID = tid, Window = paste0(win_kb, "kb"), status = "low_overlap")
          next
      }
      
      # Harmonize
      flip <- (merged$A1_1 == merged$A2_2 & merged$A2_1 == merged$A1_2)
      merged$beta_2[flip] <- -merged$beta_2[flip]
      merged$eaf_2[flip] <- 1 - merged$eaf_2[flip]
      match_ok <- (merged$A1_1 == merged$A1_2 & merged$A2_1 == merged$A2_2) | flip
      merged <- merged[match_ok]
      if (nrow(merged) < 10) {
          summary_list[[paste0(win_kb)]] <- data.table(Task_ID = tid, Window = paste0(win_kb, "kb"), status = "low_overlap_after_harm")
          next
      }
      
      # LD
      snps <- unique(merged$snp)
      R <- compute_ld_plink(snps, TARGET_CHR, LD_DIR, DIR_TEMP)
      if (is.null(R)) {
          summary_list[[paste0(win_kb)]] <- data.table(Task_ID = tid, Window = paste0(win_kb, "kb"), status = "ld_fail")
          next
      }
      
      common <- intersect(merged$snp, rownames(R))
      merged <- merged[snp %in% common]
      R <- R[common, common]
      
      # Format for coloc
      get_coloc_list <- function(ds, beta, se, eaf, n_col, snp) {
        type <- "quant"
        s <- NULL
        if (!is.na(ds$N_case) && ds$N_case > 0 && !is.na(ds$N_control) && ds$N_control > 0) {
          type <- "cc"
          s <- ds$N_case / (ds$N_case + ds$N_control)
        }
        res <- list(beta=beta, varbeta=se^2, MAF=pmin(eaf, 1-eaf), N=n_col, type=type, snp=snp)
        if (type == "cc") res$s <- s
        return(res)
      }
      
      d1_list <- get_coloc_list(ds1, merged$beta_1, merged$se_1, merged$eaf_1, merged$N_1, merged$snp)
      d2_list <- get_coloc_list(ds2, merged$beta_2, merged$se_2, merged$eaf_2, merged$N_2, merged$snp)
      
      # Run SuSiE
      sus1 <- tryCatch(run_susie(d1_list, R), error=function(e) NULL)
      sus2 <- tryCatch(run_susie(d2_list, R), error=function(e) NULL)
      
      pp_h4_susie <- NA
      res_susie <- NULL
      if (!is.null(sus1) && !is.null(sus2)) {
        res_susie <- tryCatch(coloc::coloc.susie(sus1, sus2), error=function(e) NULL)
        if (!is.null(res_susie) && !is.null(res_susie$summary)) {
            pp_h4_susie <- max(res_susie$summary$PP.H4.abf, na.rm=TRUE)
        }
      }
      
      # Run ABF as backup
      coloc_abf <- tryCatch(coloc::coloc.abf(d1_list, d2_list), error=function(e) NULL)
      pp_h4_abf <- if(!is.null(coloc_abf)) coloc_abf$summary["PP.H4.abf"] else NA
      
      # Process SuSiE summary and evaluate Signal Sets
      if (!is.null(res_susie) && !is.null(res_susie$summary) && nrow(res_susie$summary) > 0) {
          susie_summ <- as.data.table(res_susie$summary)
          # Add LD with rs1800628
          susie_summ$r2_with_rs1800628 <- NA_real_
          if ("rs1800628" %in% rownames(R)) {
              for (i in seq_len(nrow(susie_summ))) {
                  h1 <- as.character(susie_summ$hit1[i])
                  if (h1 %in% rownames(R)) {
                      susie_summ$r2_with_rs1800628[i] <- R["rs1800628", h1]^2
                  }
              }
          }
          susie_summ$Is_rs1800628_Linked <- ifelse(!is.na(susie_summ$r2_with_rs1800628) & susie_summ$r2_with_rs1800628 > 0.6, "Yes", "No")
        
        # Assign Signal Sets (group by r2 > 0.6)
        susie_summ$Signal_Set <- NA_integer_
        set_id <- 1
        for (i in seq_len(nrow(susie_summ))) {
            if (is.na(susie_summ$Signal_Set[i])) {
                susie_summ$Signal_Set[i] <- set_id
                h1_i <- as.character(susie_summ$hit1[i])
                for (j in seq_len(nrow(susie_summ))) {
                    if (i != j && is.na(susie_summ$Signal_Set[j])) {
                        h1_j <- as.character(susie_summ$hit1[j])
                        if (h1_i %in% rownames(R) && h1_j %in% rownames(R)) {
                            r2_ij <- R[h1_i, h1_j]^2
                            if (!is.na(r2_ij) && r2_ij > 0.6) {
                                susie_summ$Signal_Set[j] <- set_id
                            }
                          } else if (h1_i == h1_j) {
                              susie_summ$Signal_Set[j] <- set_id
                          }
                      }
                  }
                  set_id <- set_id + 1
              }
          }
          res_susie$summary <- as.data.frame(susie_summ)
      }
      
      # Lead SNP selection (prioritize rs1800628 linked strong signals)
    lead_snp <- NULL
    if (!is.null(res_susie) && !is.null(res_susie$summary) && nrow(res_susie$summary) > 0) {
        susie_summ <- as.data.table(res_susie$summary)
        strong_linked <- susie_summ[PP.H4.abf > 0.8 & Is_rs1800628_Linked == "Yes"]
        if (nrow(strong_linked) > 0) {
            best_linked_idx <- which.max(strong_linked$PP.H4.abf)
            lead_snp <- as.character(strong_linked$hit1[best_linked_idx])
        } else {
            best_idx <- which.max(susie_summ$PP.H4.abf)
            if (length(best_idx) > 0) lead_snp <- as.character(susie_summ$hit1[best_idx])
        }
    }
    if (is.null(lead_snp) && !is.null(coloc_abf) && "SNP.PP.H4" %in% colnames(coloc_abf$results)) {
         max_abf_idx <- which.max(coloc_abf$results$SNP.PP.H4)
         if (length(max_abf_idx) > 0 && coloc_abf$results$SNP.PP.H4[max_abf_idx] > 0.8) lead_snp <- as.character(coloc_abf$results$snp[max_abf_idx])
    }
      if (is.null(lead_snp)) lead_snp <- merged$snp[which.min(merged$pval_1 * merged$pval_2)]
      
      # Visualize ( tiff )
      viz_file_lead <- file.path(out_task_dir, paste0(tid, "_", win_kb, "kb_locuscompare_lead.tiff"))
      create_locuscompare_plot(
          d1 = merged[, .(snp, pval=pval_1)], 
          d2 = merged[, .(snp, pval=pval_2)], 
          name1 = paste0(ds1$Phenotype, " (", ds1$Source, ")"), 
          name2 = paste0(ds2$Phenotype, " (", ds2$Source, ")"), 
          lead_snp = lead_snp, out_file = viz_file_lead
      )
      
      viz_file_rs1800628 <- file.path(out_task_dir, paste0(tid, "_", win_kb, "kb_locuscompare_rs1800628.tiff"))
      create_locuscompare_plot(
          d1 = merged[, .(snp, pval=pval_1)], 
          d2 = merged[, .(snp, pval=pval_2)], 
          name1 = paste0(ds1$Phenotype, " (", ds1$Source, ")"), 
          name2 = paste0(ds2$Phenotype, " (", ds2$Source, ")"), 
          lead_snp = "rs1800628", out_file = viz_file_rs1800628
      )
      
      # Save SuSiE results summary
    if (!is.null(res_susie) && !is.null(res_susie$summary) && nrow(res_susie$summary) > 0) {
        susie_summ <- as.data.table(res_susie$summary)
        # Get strong signals (> 0.8) or just the best one
        strong_signals <- susie_summ[PP.H4.abf > 0.8]
        if (nrow(strong_signals) == 0) {
            strong_signals <- susie_summ[which.max(susie_summ$PP.H4.abf)]
        }
          # Sort: rs1800628 linked first, then by PP.H4.abf
          strong_signals <- strong_signals[order(Is_rs1800628_Linked == "Yes", PP.H4.abf, decreasing = TRUE)]
          
          res_rows <- lapply(seq_len(nrow(strong_signals)), function(i) {
              data.table(
                  Task_ID = tid, Window = paste0(win_kb, "kb"), 
                  ID1 = ds1$ID, Phenotype1 = ds1$Phenotype, Source1 = ds1$Source,
                  ID2 = ds2$ID, Phenotype2 = ds2$Phenotype, Source2 = ds2$Source,
                  n_snps = nrow(merged), 
                  lead_snp = as.character(strong_signals$hit1[i]),
                  hit2 = as.character(strong_signals$hit2[i]),
                  Method = "SuSiE-coloc",
                  PP.H4.SuSiE = strong_signals$PP.H4.abf[i],
                  Signal_Set = strong_signals$Signal_Set[i],
                  Is_rs1800628_Linked = strong_signals$Is_rs1800628_Linked[i],
                  r2_with_rs1800628 = strong_signals$r2_with_rs1800628[i],
                  status = "ok"
              )
          })
          res_row <- rbindlist(res_rows)
      } else {
          res_row <- data.table(
              Task_ID = tid, Window = paste0(win_kb, "kb"), 
              ID1 = ds1$ID, Phenotype1 = ds1$Phenotype, Source1 = ds1$Source,
              ID2 = ds2$ID, Phenotype2 = ds2$Phenotype, Source2 = ds2$Source,
              n_snps = nrow(merged), lead_snp = lead_snp, hit2 = NA_character_,
              Method = "coloc-abf",
              PP.H4.SuSiE = pp_h4_abf,
              Signal_Set = NA_integer_,
              Is_rs1800628_Linked = NA_character_,
              r2_with_rs1800628 = NA_real_,
              status = "ok"
          )
      }
      
      fwrite(res_row, paste0(out_prefix, "_summary.csv"))
      summary_list[[paste0(win_kb)]] <- res_row
      
      if (!is.null(res_susie) && !is.null(res_susie$summary)) {
          susie_details <- as.data.table(res_susie$summary)
          susie_details$Method <- "SuSiE-coloc"
          susie_details$Task_ID <- tid
          susie_details$Window <- paste0(win_kb, "kb")
          susie_details$ID1 <- ds1$ID
          susie_details$Phenotype1 <- ds1$Phenotype
          susie_details$Source1 <- ds1$Source
          susie_details$ID2 <- ds2$ID
          susie_details$Phenotype2 <- ds2$Phenotype
          susie_details$Source2 <- ds2$Source
          fwrite(susie_details, paste0(out_prefix, "_coloc_susie_details.csv"))
          details_list[[paste0(win_kb)]] <- susie_details
      }
  }
  
  return(list(summary = rbindlist(summary_list, fill=TRUE), details = rbindlist(details_list, fill=TRUE)))
}

## 7. Parallel Execution / 
log_msg(paste("Starting parallel execution with 16 workers..."))
plan(multisession, workers = 16)

results_list <- future_lapply(tasks, function(tsk) {
  tryCatch({
    process_task(tsk)
  }, error = function(e) {
    list(summary = data.table(Task_ID = tsk$task_id, status = paste("error:", e$message)), details = NULL)
  })
}, future.seed = TRUE)

## 8. Summary & Cleanup / 
summary_list <- lapply(results_list, function(x) if("summary" %in% names(x)) x$summary else NULL)
details_list <- lapply(results_list, function(x) if("details" %in% names(x)) x$details else NULL)

final_res <- rbindlist(summary_list, fill = TRUE)
fwrite(final_res, file.path(DIR_README, paste0("All_Pairwise_Coloc_Summary_", timestamp, ".csv")))

final_details <- rbindlist(details_list, fill = TRUE)
if (nrow(final_details) > 0) {
  # reorder columns to put context first
  setcolorder(final_details, c("Task_ID", "ID1", "Phenotype1", "Source1", "ID2", "Phenotype2", "Source2", "Method", "Window"))
  fwrite(final_details, file.path(DIR_README, paste0("All_Pairwise_Coloc_Details_", timestamp, ".csv")))
}

end_time <- Sys.time()
duration <- round(difftime(end_time, start_time, units = "mins"), 2)
log_msg(paste("Analysis completed in", duration, "minutes."))

# Write README
readme_content <- c(
  "# Pairwise SuSiE Colocalization Analysis Report",
  paste("Date:", Sys.time()),
  paste("Target Region:", TARGET_CHR, ":", TARGET_POS, "(", TARGET_SNP, "with windows", paste(WINDOW_SIZES/1000, "kb", collapse=", "), ")"),
  paste("Total Input Datasets:", n_ds),
  paste("Total Pairwise Tasks:", length(tasks)),
  paste("Total Completed Successfully:", sum(final_res$status == "ok", na.rm=TRUE)),
  paste("Runtime:", duration, "minutes"),
  "",
  "## Directories",
  "- `results/`: Individual coloc details and locuscompare plots for each pair.",
  "- `logs/`: Execution logs.",
  "- `readme/`: Unified list and summary CSV.",
  "- `Detailed_Stats_*/`: Contains statistical summaries based on signal type and phenotype."
)
writeLines(readme_content, file.path(DIR_README, "README.md"))
log_msg("Saved README.md. Done.")

# ==========================================
# 9.  (Auto-generate Detailed Stats)
# ==========================================
log_msg("Invoking generate_detailed_stats.R to generate statistical summaries...")
stats_script <- file.path(WORK_DIR, "generate_detailed_stats.R")
if (file.exists(stats_script)) {
    #  summary （, ）
    tryCatch({
        system(paste("Rscript", shQuote(stats_script), shQuote(BASE_OUT)))
        log_msg("Successfully generated detailed statistics.")
    }, error = function(e) {
        log_msg(paste("Failed to run generate_detailed_stats.R:", e$message))
    })
} else {
    log_msg("Warning: generate_detailed_stats.R not found in working directory.")
}
