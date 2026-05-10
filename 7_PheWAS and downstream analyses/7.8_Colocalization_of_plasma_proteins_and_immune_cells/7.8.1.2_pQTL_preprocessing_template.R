#!/usr/bin/env Rscript
# ==============================================================================
# [Script]: 7.8.1.2_pQTL_preprocessing_template.R
# [Method]: SuSiE Colocalization Preprocessing
# [Step]: Preprocess GWAS/pQTL data for SuSiE coloc
# 
# [Function]:
# Preprocess summary statistics to the required format for SuSiE colocalization.
# 
# [Input]:
#   --input_file   : Raw GWAS/pQTL summary statistics file
#   --out_dir      : Output directory path
# 
# [Output]:
# Formatted summary statistics ready for colocalization analysis.
# ==============================================================================

rm(list = ls())
gc()

# 1. Load Dependencies
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(readr)
})

# 2. Parse Command Line Arguments
option_list <- list(
  make_option(c("--input_file"), type="character", default=NULL, help="Input summary stats file"),
  make_option(c("--out_dir"), type="character", default="./Prep_Results", help="Output directory path")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_file)) {
  print_help(opt_parser)
  stop("Missing required input arguments.")
}

INPUT_FILE <- opt$input_file
OUTPUT_DIR <- opt$out_dir

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

options(stringsAsFactors = FALSE)
Sys.setenv(TZ = "UTC")

# 1)  / Configuration
pQTL_src_dir  <- "./"
ld_snp_info   <- "./"
output_dir <- OUTPUT_DIR
csv_input     <- "./"

#  / Fixed arguments
cores_arg     <- 10       # 10, CPUIO
test_mode     <- FALSE    # 5 -> 
resume_arg    <- TRUE

# 2)  / Load libraries
suppressPackageStartupMessages({
  cran_mirror <- "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"
  if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table", repos = cran_mirror)
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr", repos = cran_mirror)
  if (!requireNamespace("stringr", quietly = TRUE)) install.packages("stringr", repos = cran_mirror)
  if (!requireNamespace("readr", quietly = TRUE)) install.packages("readr", repos = cran_mirror)
  if (!requireNamespace("tibble", quietly = TRUE)) install.packages("tibble", repos = cran_mirror)
  if (!requireNamespace("parallel", quietly = TRUE)) install.packages("parallel", repos = cran_mirror)
  
  library(data.table)
  library(dplyr)
  library(stringr)
  library(readr)
  library(tibble)
  library(parallel)
})

# 3)  / Common config
available_cores <- parallel::detectCores()
cores <- max(1, min(cores_arg, available_cores - 2)) # Leave cores for system
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

out_root <- file.path(output_dir, paste0("UKB_pQTL_preprocessing_", timestamp))
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

#  / Subfolders
dir.create(file.path(out_root, "results"), showWarnings = FALSE)
dir.create(file.path(out_root, "stats"), showWarnings = FALSE)
dir.create(file.path(out_root, "logs"), showWarnings = FALSE)
dir.create(file.path(out_root, "checkpoints"), showWarnings = FALSE)
dir.create(file.path(out_root, "tmp"), showWarnings = FALSE)
dir.create(file.path(out_root, "scripts"), showWarnings = FALSE)

#  / Backup script
this_script_path <- file.path(output_dir, "UKB_pQTL_preprocessing_for_SuSiE_finetuned.R")
try({ file.copy(from = this_script_path, to = file.path(out_root, "scripts", "UKB_pQTL_preprocessing_for_SuSiE_finetuned.R"), overwrite = TRUE) }, silent = TRUE)

# 4)  / Initialize logging
log_file <- file.path(out_root, "logs", paste0("run_log_", timestamp, ".txt"))
sink(log_file, append = TRUE, split = TRUE)
cat("========================================\n")
cat("UKB pQTL Preprocessing for SuSiE/coloc (Genome-wide)\n")
cat("UKBpQTL SuSiE/coloc \n")
cat("========================================\n")
cat("Start Time:", as.character(Sys.time()), "\n")
cat("Cores:", cores, "(available:", available_cores, ")\n")
cat("Output Root:", out_root, "\n")
cat("pQTL Source Dir:", pQTL_src_dir, "\n")
cat("LD SNP Info:", ld_snp_info, "\n")
cat("CSV Input:", csv_input, "\n")
cat("Test Mode:", test_mode, "\n")
cat("Resume:", resume_arg, "\n\n")

# 5)  / Utility functions
safe_upper <- function(x) { toupper(gsub("\n|\r", "", as.character(x))) }

# pQTL  / pQTL specific column mapping
detect_sceqtl_cols <- function(dt) {
  cols <- colnames(dt)
  map <- list()
  map$rsid <- intersect(cols, c("rs_id","SNP","rsid","RSID","variant","MarkerName"))[1]
  map$chr  <- intersect(cols, c("chromosome","CHR","chr","Chromosome","chr.exposure"))[1]
  map$pos  <- intersect(cols, c("base_pair_location","POS","pos","BP","position","Position","Position_BP","pos.exposure"))[1]
  map$ea   <- intersect(cols, c("effect_allele","EA","A1","Allele1","ALT","effect_allele.exposure"))[1]
  map$oa   <- intersect(cols, c("other_allele","NEA","A2","Allele2","REF","other_allele.exposure"))[1]
  map$beta <- intersect(cols, c("beta","BETA","effect_size","ES","beta.exposure"))[1]
  map$se   <- intersect(cols, c("standard_error","SE","se","se.exposure"))[1]
  map$pval <- intersect(cols, c("p_value","P_VALUE","pval","p","P","pvalue","pval.exposure"))[1]
  map$eaf  <- intersect(cols, c("effect_allele_frequency","eaf","EAF","freq","AF","eaf.exposure"))[1]
  map$n    <- intersect(cols, c("n","N","samplesize","samplesize.exposure","Sample size"))[1]
  map
}

validate_mapping <- function(map) {
  required <- c("rsid","chr","pos","ea","oa","beta","se","pval","n")
  all(!is.na(map[required]))
}

# LD  / Read LD panel
read_ld_panel <- function(ld_path) {
  ld <- data.table::fread(ld_path, sep = "\t")
  setnames(ld, old = c("SNP_ID","Chromosome","Position_BP","Allele1","Allele2","Allele_Direction"), 
               new = c("rsid","chr","pos","ld_a1","ld_a2","ld_dir"))
  ld$rsid <- as.character(ld$rsid)
  ld$ld_a1 <- safe_upper(ld$ld_a1)
  ld$ld_a2 <- safe_upper(ld$ld_a2)
  # ,  / Genome-wide, removed position filtering
  ld
}

# LD / Allele alignment to LD panel
# , 、beta、eaf
align_alleles_to_ld <- function(dt, ld_dt, map) {
  dt <- dt[!is.na(dt[[map$rsid]]),]
  
  dt[[map$ea]] <- safe_upper(dt[[map$ea]])
  dt[[map$oa]] <- safe_upper(dt[[map$oa]])
  dt <- dt[dt[[map$ea]] %in% c("A","C","G","T") & dt[[map$oa]] %in% c("A","C","G","T"), ]

  merged <- dt %>% dplyr::inner_join(ld_dt, by = c(setNames("rsid", map$rsid)))

  ea <- merged[[map$ea]]; oa <- merged[[map$oa]]
  ld_a1 <- merged$ld_a1; ld_a2 <- merged$ld_a2
  
  status <- ifelse(is.na(ld_a1) | is.na(ld_a2), "missing_ld",
             ifelse(is.na(ea) | is.na(oa), "missing_data",
               ifelse(ea == ld_a1 & oa == ld_a2, "consistent",
                 ifelse(ea == ld_a2 & oa == ld_a1, "flipped", "inconsistent"))))

  merged$alignment_status <- status

  beta <- as.numeric(merged[[map$beta]])
  eaf  <- suppressWarnings(as.numeric(merged[[map$eaf]]))
  
  #  flipped, beta, eaf1-eaf
  merged$final_beta <- ifelse(status == "flipped", -beta, beta)
  merged$final_eaf  <- ifelse(status == "flipped" & !is.na(eaf), 1 - eaf, eaf)

  # LD
  out <- data.table::data.table(
    snp = merged[[map$rsid]],
    chr = suppressWarnings(as.integer(merged[[map$chr]])),
    pos = suppressWarnings(as.integer(merged[[map$pos]])),
    A1  = ifelse(status %in% c("consistent","flipped"), merged$ld_a1, NA_character_),
    A2  = ifelse(status %in% c("consistent","flipped"), merged$ld_a2, NA_character_),
    beta = merged$final_beta,
    se   = suppressWarnings(as.numeric(merged[[map$se]])),
    pval = suppressWarnings(as.numeric(merged[[map$pval]])),
    eaf  = merged$final_eaf,
    N    = suppressWarnings(as.numeric(merged[[map$n]]))
  )

  out <- out[status %in% c("consistent","flipped")]
  out <- out[!is.na(snp) & !is.na(chr) & !is.na(pos) & !is.na(beta) & !is.na(se) & !is.na(pval)]
  
  list(out = out, stats = list(
    total = nrow(merged),
    consistent = sum(status == "consistent", na.rm = TRUE),
    flipped = sum(status == "flipped", na.rm = TRUE),
    excluded = sum(status %in% c("missing_ld","missing_data","inconsistent"), na.rm = TRUE),
    aligned = nrow(out)
  ))
}

#  / Process single protein file
process_protein <- function(task_id, ld_dt, out_root, ld_rsids_file, awk_script_file) {
  data.table::setDTthreads(1)
  
  cell_dir <- file.path(out_root, "results")
  out_file <- file.path(cell_dir, paste0(task_id, ".tsv.gz"))
  
  if (file.exists(out_file)) {
    return(list(id = task_id, processed = TRUE, skipped = TRUE, n_aligned = NA))
  }
  
  src_file <- file.path(pQTL_src_dir, paste0(task_id, ".tsv.gz"))
  if (!file.exists(src_file)) {
    return(list(id = task_id, processed = FALSE, skipped = TRUE, error = "Source file not found", n_aligned = 0))
  }

  cat(sprintf("[%s] Processing %s\n", format(Sys.time(), "%H:%M:%S"), task_id))

  #  awk ,  (Extract SNPs present in LD panel)
  tmp_file <- file.path(out_root, "tmp", paste0(task_id, "_filtered.tsv"))
  
  awk_cmd <- sprintf(
    "bash -c \"gunzip -c '%s' | awk -f '%s' '%s' - > '%s'\"",
    src_file, awk_script_file, ld_rsids_file, tmp_file
  )
  
  system(awk_cmd, ignore.stdout = TRUE, ignore.stderr = TRUE)
  
  if (!file.exists(tmp_file) || file.info(tmp_file)$size == 0) {
    return(list(id = task_id, processed = FALSE, skipped = TRUE, error = "Awk extraction failed or empty", n_aligned = 0))
  }

  dt <- tryCatch({ data.table::fread(tmp_file) }, error = function(e) NULL)
  unlink(tmp_file)
  
  if (is.null(dt) || nrow(dt) == 0) {
    return(list(id = task_id, processed = FALSE, skipped = TRUE, error = "Empty, unreadable file", n_aligned = 0))
  }
  
  source_total <- nrow(dt)
  map <- detect_sceqtl_cols(dt)
  
  mapping_rec <- data.frame(
    task_id = task_id,
    rsid = map$rsid, chr = map$chr, pos = map$pos,
    ea = map$ea, oa = map$oa, beta = map$beta, se = map$se, pval = map$pval, eaf = map$eaf, N = map$n,
    stringsAsFactors = FALSE
  )
  mapping_file <- file.path(out_root, "stats", paste0(task_id, "_mapping.csv"))
  data.table::fwrite(mapping_rec, mapping_file)

  if (!validate_mapping(map)) {
    return(list(id = task_id, processed = FALSE, skipped = TRUE, error = "Missing key columns", n_aligned = 0))
  }

  align_res <- align_alleles_to_ld(dt, ld_dt, map)
  out <- align_res$out
  stats <- align_res$stats
  
  if (!is.null(out$eaf)) {
    out$maf <- ifelse(out$eaf > 0.5, 1 - out$eaf, out$eaf)
    out <- subset(out, maf > 0.01)
    out$maf <- NULL
  }
  
  dup_snps <- out$snp[duplicated(out$snp)]
  if (length(dup_snps) > 0) {
    dup_info <- data.frame(
      task_id = task_id,
      n_duplicates = length(dup_snps),
      example_dups = paste(head(unique(dup_snps), 5), collapse = ";"),
      stringsAsFactors = FALSE
    )
    dup_file <- file.path(out_root, "stats", paste0("dup_", task_id, ".csv"))
    data.table::fwrite(dup_info, dup_file)
  }
  
  data.table::fwrite(out, out_file, sep = "\t", compress = "gzip")
  
  cp_file <- file.path(out_root, "checkpoints", "processed_genes.txt")
  write(x = task_id, file = cp_file, append = TRUE)
  
  rm(dt, align_res, out)
  gc()
  
  list(id = task_id, processed = TRUE, skipped = FALSE, n_aligned = stats$aligned, source_total = source_total, intersection_total = stats$total, stats = stats)
}

# 6) LDRSIDawk
cat("LD: ", ld_snp_info, "\n")
ld_dt <- read_ld_panel(ld_snp_info)
cat("LD SNP: ", nrow(ld_dt), "\n")

ld_rsids_file <- file.path(out_root, "tmp", "ld_rsids.txt")
cat("LD rsidsawk...\n")
data.table::fwrite(ld_dt[, .(rsid)], ld_rsids_file, col.names = FALSE, quote = FALSE)

awk_script_file <- file.path(out_root, "tmp", "filter.awk")
awk_code <- "
BEGIN { FS=\"\\t\"; OFS=\"\\t\" }
NR==FNR { a[$1]; next }
FNR==1 {
  col=0;
  for(i=1; i<=NF; i++) {
    if($i==\"rs_id\" || $i==\"SNP\" || $i==\"rsid\" || $i==\"RSID\" || $i==\"variant\" || $i==\"MarkerName\") {
      col=i; break;
    }
  }
  print $0; next
}
col > 0 && ($col in a) { print $0 }
"
writeLines(awk_code, awk_script_file)

# 7) CSV
cat(": ", csv_input, "\n")
csv_data <- data.table::fread(csv_input)
# GCST
if (!"GCST" %in% colnames(csv_data)) stop("CSVGCST")

# GCST
unique_gcst <- unique(csv_data$GCST)
unique_gcst <- unique_gcst[!is.na(unique_gcst) & unique_gcst != ""]

task_list <- data.frame(
  task_id = unique_gcst,
  stringsAsFactors = FALSE
)

cat("CSV: ", nrow(task_list), "\n")

processed <- character(0)
if (resume_arg) {
  cp_file <- file.path(out_root, "checkpoints", "processed_genes.txt")
  if (file.exists(cp_file)) processed <- readLines(cp_file)
}

to_process <- task_list[!task_list$task_id %in% processed, , drop = FALSE]

if (test_mode) {
  to_process <- utils::head(to_process, 5)
  cat(": 5. \n")
}

cat(": ", nrow(to_process), "\n")

# 8)  
start_time <- Sys.time()
res_list <- list()

if (nrow(to_process) > 0) {
  chunk_size <- 50
  n_tasks <- nrow(to_process)
  n_chunks <- ceiling(n_tasks / chunk_size)

  cat(sprintf("\n %d  %d ,  %d . \n", n_tasks, n_chunks, chunk_size))

  data.table::setDTthreads(1)

  for (chunk_idx in seq_len(n_chunks)) {
    start_idx <- (chunk_idx - 1) * chunk_size + 1
    end_idx <- min(chunk_idx * chunk_size, n_tasks)
    
    cat(sprintf("\n---  %d/%d  ( %d  %d) ---\n", chunk_idx, n_chunks, start_idx, end_idx))
    gc()
    
    chunk_indices <- start_idx:end_idx
    
    #  mclapply (mc.preschedule = TRUE) macOSFD
    chunk_res <- parallel::mclapply(chunk_indices, function(i) {
      rec <- to_process[i, , drop = FALSE]
      tryCatch({
        process_protein(rec$task_id, ld_dt, out_root, ld_rsids_file, awk_script_file)
      }, error = function(e) {
        list(id = rec$task_id, processed = FALSE, skipped = TRUE, n_aligned = 0, error = as.character(e))
      })
    }, mc.cores = cores, mc.preschedule = TRUE)
    
    res_list <- c(res_list, chunk_res)
    rm(chunk_res)
    gc()
  }
} else {
  cat(". \n")
}

# 10) 
summary_dt <- data.table::rbindlist(lapply(res_list, function(x) {
  if (is.null(x)) return(NULL)
  data.table::data.table(
    task_id = x$id,
    processed = isTRUE(x$processed),
    skipped = isTRUE(x$skipped),
    error = ifelse(is.null(x$error), NA_character_, x$error),
    source_snp_total = ifelse(is.null(x$source_total), NA_integer_, as.integer(x$source_total)),
    intersection_snp_total = ifelse(is.null(x$intersection_total), NA_integer_, as.integer(x$intersection_total)),
    n_aligned = ifelse(is.null(x$n_aligned), NA_integer_, as.integer(x$n_aligned)),
    consistent = ifelse(is.null(x$stats), NA_integer_, as.integer(x$stats$consistent)),
    flipped = ifelse(is.null(x$stats), NA_integer_, as.integer(x$stats$flipped))
  )
}), fill = TRUE)

data.table::fwrite(summary_dt, file.path(out_root, "stats", "processing_summary.csv"))

dup_files <- list.files(file.path(out_root, "stats"), pattern = "^dup_.*\\.csv$", full.names = TRUE)
if (length(dup_files) > 0) {
  all_dups <- lapply(dup_files, data.table::fread)
  data.table::fwrite(data.table::rbindlist(all_dups, fill=TRUE), file.path(out_root, "stats", "all_duplicate_snps.csv"))
}

# 12) 
end_time <- Sys.time()
cat("\nEnd Time:", as.character(end_time), "\n")
cat("Duration:", as.character(difftime(end_time, start_time, units = "mins")), "mins\n")
sink(NULL)

# Cleanup
try({ unlink(file.path(out_root, "tmp"), recursive = TRUE) }, silent = TRUE)

if (test_mode) message("Test mode finished. . ")
