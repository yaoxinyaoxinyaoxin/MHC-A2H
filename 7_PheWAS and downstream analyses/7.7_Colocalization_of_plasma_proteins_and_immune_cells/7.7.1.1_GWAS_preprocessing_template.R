#!/usr/bin/env Rscript
# ==============================================================================
# [Script]: 7.7.1.1_GWAS_preprocessing_template.R
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

# 1)  / Parse CLI args
args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(flag, default = NULL) {
  idx <- which(args == flag)
  if (length(idx) == 1 && idx < length(args)) return(args[idx + 1])
  default
}

src_dir <- get_arg("--src_dir", 
  default = NULL)
ld_snp_info <- get_arg("--ld_snp_info", 
  default = NULL)
sample_size_file <- get_arg("--sample_size_file", 
  default = NULL)
output_dir <- get_arg("--output_dir", 
  default = NULL)
cores_arg <- as.integer(get_arg("--cores", default = "16"))
test_mode <- tolower(get_arg("--test_mode", default = "FALSE")) %in% c("true","1","t","yes")
ignore_history <- tolower(get_arg("--ignore_history", default = "FALSE")) %in% c("true","1","t","yes")

# 2)  / Load libraries
suppressPackageStartupMessages({
  cran_mirror <- "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"
  if (!requireNamespace("data.table", quietly = TRUE)) install.packages("data.table", repos = cran_mirror)
  if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr", repos = cran_mirror)
  if (!requireNamespace("future", quietly = TRUE)) install.packages("future", repos = cran_mirror)
  if (!requireNamespace("future.apply", quietly = TRUE)) install.packages("future.apply", repos = cran_mirror)
  if (!requireNamespace("parallel", quietly = TRUE)) install.packages("parallel", repos = cran_mirror)
  
  library(data.table)
  library(dplyr)
  library(future)
  library(future.apply)
  library(parallel)
})

# 3)  / Common config
available_cores <- parallel::detectCores()
cores <- max(1, min(cores_arg, available_cores))
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

base_work_dir <- if (is.null(output_dir) || output_dir == "") getwd() else output_dir
out_root <- file.path(base_work_dir, paste0("731_Immune_Cell_Preprocessed_", timestamp))
dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

#  / Subfolders
dir.create(file.path(out_root, "results"), showWarnings = FALSE)
dir.create(file.path(out_root, "stats"), showWarnings = FALSE)
dir.create(file.path(out_root, "logs"), showWarnings = FALSE)
dir.create(file.path(out_root, "readme"), showWarnings = FALSE)
dir.create(file.path(out_root, "checkpoints"), showWarnings = FALSE)
dir.create(file.path(out_root, "scripts"), showWarnings = FALSE)
dir.create(file.path(out_root, "temp"), showWarnings = FALSE) # LDrsid

#  scripts/ 
this_script <- "731_GRCh37_preprocessing_for_SuSiE_coloc.R"
try({ file.copy(from = file.path(getwd(), this_script), to = file.path(out_root, "scripts", this_script), overwrite = TRUE) }, silent = TRUE)

# 4)  / Initialize logging
log_file <- file.path(out_root, "logs", paste0("run_log_", timestamp, ".txt"))
sink(log_file, append = TRUE, split = TRUE)
cat("========================================\n")
cat("731 Immune Cell Preprocessing for SuSiE/coloc\n")
cat("731 (LD, Awk)\n")
cat("========================================\n")
cat("Start Time / :", as.character(Sys.time), "\n")
cat("Cores / :", cores, "(available:", available_cores, ")\n")
cat("Output Root / :", out_root, "\n")
cat("Source Dir / :", src_dir, "\n")
cat("LD SNP Info / LD:", ld_snp_info, "\n")
cat("Sample Size File / :", sample_size_file, "\n")
cat("Test Mode / :", test_mode, "\n\n")

# 5)  / Utility functions
safe_upper <- function(x) { toupper(gsub("\n|\r", "", as.character(x))) }

#  / Read Sample Size Data
read_sample_size <- function(path) {
  if (file.exists(path)) {
    ss <- data.table::fread(path)
    # Ensure STUDY and samplesize exist
    if ("STUDY" %in% names(ss) && "samplesize" %in% names(ss)) {
      cat(", : ", nrow(ss), "\n")
      return(ss[, .(STUDY, samplesize)])
    }
  }
  cat(": (STUDY, samplesize), N. \n")
  return(NULL)
}

#  / Auto column mapping
detect_cols <- function(dt) {
  cols <- colnames(dt)
  map <- list()
  map$rsid <- intersect(cols, c("rs_id","SNP","rsid","RSID","variant","MarkerName"))[1]
  map$chr  <- intersect(cols, c("chromosome","chr","CHR","Chromosome","chr.exposure"))[1]
  map$pos  <- intersect(cols, c("base_pair_location","pos","BP","position","Position","Position_BP","pos.exposure"))[1]
  map$ea   <- intersect(cols, c("effect_allele","EA","A1","Allele1","ALT","effect_allele.exposure"))[1]
  map$oa   <- intersect(cols, c("other_allele","NEA","A2","Allele2","REF","other_allele.exposure"))[1]
  map$beta <- intersect(cols, c("beta","BETA","effect_size","ES","beta.exposure"))[1]
  map$se   <- intersect(cols, c("standard_error","se","SE","se.exposure"))[1]
  map$pval <- intersect(cols, c("p_value","pval","p","P","pvalue","pval.exposure"))[1]
  map$eaf  <- intersect(cols, c("effect_allele_frequency","eaf","EAF","freq","AF","eaf.exposure"))[1]
  map$n    <- intersect(cols, c("N","n","samplesize","samplesize.exposure","Sample size"))[1]
  map
}

#  / Validate required columns
validate_mapping <- function(map) {
  required <- c("rsid","chr","pos","ea","oa","beta","se","pval","eaf")
  all(!is.na(map[required]))
}

# LDrsidawk / Extract LD rsids to temp file
create_ld_rsid_file <- function(ld_dt, temp_dir) {
  rsid_file <- file.path(temp_dir, "ld_rsids.txt")
  data.table::fwrite(ld_dt[, .(rsid)], rsid_file, col.names = FALSE, quote = FALSE)
  rsid_file
}

# （LD）/ Allele alignment to LD panel
align_alleles_to_ld <- function(dt, ld_dt, map) {
  dt[[map$ea]] <- safe_upper(dt[[map$ea]])
  dt[[map$oa]] <- safe_upper(dt[[map$oa]])

  # A/C/G/T
  dt <- dt[dt[[map$ea]] %in% c("A","C","G","T") & dt[[map$oa]] %in% c("A","C","G","T"), ]

  # LD
  merged <- dt %>% dplyr::inner_join(ld_dt, by = c(setNames("rsid", map$rsid)))

  # , 
  ea <- merged[[map$ea]]
  oa <- merged[[map$oa]]
  ld_a1 <- merged$ld_a1
  ld_a2 <- merged$ld_a2
  
  #  (Anchor to LD Panel):
  # consistent: EA == LD_A1 & OA == LD_A2
  # flipped: EA == LD_A2 & OA == LD_A1
  # inconsistent: 
  status <- ifelse(is.na(ld_a1) | is.na(ld_a2), "missing_ld",
             ifelse(is.na(ea) | is.na(oa), "missing_data",
               ifelse(ea == ld_a1 & oa == ld_a2, "consistent",
                 ifelse(ea == ld_a2 & oa == ld_a1, "flipped", "inconsistent"))))

  merged$alignment_status <- status

  # : , Beta, EAF
  beta <- as.numeric(merged[[map$beta]])
  eaf  <- suppressWarnings(as.numeric(merged[[map$eaf]]))
  merged$final_beta <- ifelse(status == "flipped", -beta, beta)
  merged$final_eaf  <- ifelse(status == "flipped" & !is.na(eaf), 1 - eaf, eaf)

  #  A1  ld_a1, A2  ld_a2 
  out <- data.table::data.table(
    snp  = merged[[map$rsid]],
    chr  = suppressWarnings(as.integer(merged[[map$chr]])),
    pos  = suppressWarnings(as.integer(merged[[map$pos]])),
    A1   = ifelse(status %in% c("consistent","flipped"), merged$ld_a1, NA_character_),
    A2   = ifelse(status %in% c("consistent","flipped"), merged$ld_a2, NA_character_),
    beta = merged$final_beta,
    se   = suppressWarnings(as.numeric(merged[[map$se]])),
    pval = suppressWarnings(as.numeric(merged[[map$pval]])),
    eaf  = merged$final_eaf,
    N    = suppressWarnings(as.numeric(merged[[map$n]])) # N, 
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

# 6) 、LD
cat(": ", src_dir, "\n")
all_files <- list.files(src_dir, pattern = "\\.tsv\\.gz$", full.names = TRUE)
if (length(all_files) == 0) stop(".tsv.gz: ", src_dir)

file_list_dt <- data.table(
  file_path = all_files,
  trait_id = sub("\\.tsv\\.gz$", "", basename(all_files))
)
cat(": ", nrow(file_list_dt), "\n")

# LD
cat("LD: ", ld_snp_info, "\n")
ld_dt <- data.table::fread(ld_snp_info, sep = "\t")
setnames(ld_dt, old = c("SNP_ID","Chromosome","Position_BP","Allele1","Allele2","Allele_Direction"), 
             new = c("rsid","chr","pos","ld_a1","ld_a2","ld_dir"))
ld_dt$rsid <- as.character(ld_dt$rsid)
ld_dt$ld_a1 <- safe_upper(ld_dt$ld_a1)
ld_dt$ld_a2 <- safe_upper(ld_dt$ld_a2)
cat("LDSNP: ", nrow(ld_dt), "\n")

# LD rsid Awk
temp_dir <- file.path(out_root, "temp")
ld_rsid_file <- create_ld_rsid_file(ld_dt, temp_dir)

sample_size_dt <- read_sample_size(sample_size_file)

# 7) 
cp_file <- file.path(out_root, "checkpoints", "processed_traits.txt")
processed <- if (file.exists(cp_file)) readLines(cp_file) else character(0)

if (!ignore_history) {
  prev_runs <- list.dirs(base_work_dir, recursive = FALSE, full.names = TRUE)
  prev_runs <- prev_runs[grepl("731_Immune_Cell_Preprocessed_\\d{8}_\\d{6}$", prev_runs)]
  for (pr in prev_runs) {
    pr_cp <- file.path(pr, "checkpoints", "processed_traits.txt")
    if (file.exists(pr_cp)) {
      processed <- unique(c(processed, readLines(pr_cp)))
    }
  }
}
to_process <- file_list_dt[!file_list_dt$trait_id %in% processed, ]
if (test_mode) {
  to_process <- utils::head(to_process, 20)
  cat(": 20. \n")
}
cat(": ", nrow(to_process), "\n")

# 8)  (Awk)
process_trait <- function(file_path, trait_id, ld_dt, ld_rsid_file, sample_size_dt, out_root) {
  trait_dir <- file.path(out_root, "results")
  out_file <- file.path(trait_dir, paste0(trait_id, ".tsv.gz"))
  
  if (file.exists(out_file)) {
    return(list(trait_id = trait_id, processed = TRUE, skipped = TRUE, n_aligned = NA))
  }
  
  #  Awk : LD reference  rsid, 
  # : ld_rsids, gz($1)
  cmd <- sprintf("gzcat '%s' | awk 'NR==FNR{a[$1]; next} FNR==1 {print; next} $1 in a {print}' '%s' -", file_path, ld_rsid_file)
  
  dt <- tryCatch({ data.table::fread(cmd = cmd) }, error = function(e) NULL)
  
  if (is.null(dt) || nrow(dt) == 0) {
    return(list(trait_id = trait_id, processed = FALSE, skipped = TRUE, n_aligned = 0))
  }
  
  source_total <- nrow(dt) - 1 # （, fread, nrow(dt)）
  
  map <- detect_cols(dt)
  
  mapping_rec <- data.frame(
    trait_id = trait_id,
    rsid = map$rsid, chr = map$chr, pos = map$pos,
    ea = map$ea, oa = map$oa, beta = map$beta, se = map$se, pval = map$pval, eaf = map$eaf, N = map$n,
    stringsAsFactors = FALSE
  )
  data.table::fwrite(mapping_rec, file.path(out_root, "stats", paste0(trait_id, "_column_mapping.csv")))
  
  if (!validate_mapping(map)) {
    return(list(trait_id = trait_id, processed = FALSE, skipped = TRUE, n_aligned = 0, error = ""))
  }
  
  #  (, )
  if (!is.null(sample_size_dt)) {
    ss_val <- sample_size_dt[STUDY == trait_id, samplesize]
    if (length(ss_val) > 0 && !is.na(ss_val[1])) {
      if (is.na(map$n)) map$n <- "N" # n
      dt[[map$n]] <- ss_val[1]
    }
  }
  
  align_res <- align_alleles_to_ld(dt, ld_dt, map)
  out <- align_res$out
  stats <- align_res$stats
  
  final_kept <- 0
  if (!is.null(out) && nrow(out) > 0) {
    final_kept <- nrow(out)
    data.table::fwrite(out, out_file, sep = "\t", compress = "gzip")
  }
  
  cp_file <- file.path(out_root, "checkpoints", "processed_traits.txt")
  write(x = trait_id, file = cp_file, append = TRUE)
  
  list(trait_id = trait_id, processed = TRUE, skipped = FALSE, n_aligned = stats$aligned, 
       source_total = source_total, intersection_total = stats$total, stats = stats, final_kept = final_kept)
}

# 9)  / Batch process
plan(multisession, workers = cores)

start_time <- Sys.time()
res_list <- list()

if (nrow(to_process) > 0) {
  cat(", :", cores, "\n")
  res_list <- future.apply::future_lapply(seq_len(nrow(to_process)), function(i) {
    rec <- to_process[i, ]
    tryCatch({
      process_trait(rec$file_path, rec$trait_id, ld_dt, ld_rsid_file, sample_size_dt, out_root)
    }, error = function(e) {
      tryCatch({ # Retry once
        process_trait(rec$file_path, rec$trait_id, ld_dt, ld_rsid_file, sample_size_dt, out_root)
      }, error = function(e2) {
        list(trait_id = rec$trait_id, processed = FALSE, skipped = TRUE, n_aligned = 0, error = as.character(e2))
      })
    })
  }, future.seed = TRUE)
}

# 10)  / Summaries and Cleanup
if (length(res_list) > 0) {
  summary_dt <- data.table::rbindlist(lapply(res_list, function(x) {
    if (is.null(x)) return(NULL)
    data.table::data.table(
      trait_id = x$trait_id,
      processed = isTRUE(x$processed),
      skipped = isTRUE(x$skipped),
      intersection_snp_total = ifelse(is.null(x$intersection_total), NA_integer_, as.integer(x$intersection_total)),
      consistent = ifelse(is.null(x$stats), NA_integer_, as.integer(x$stats$consistent)),
      flipped = ifelse(is.null(x$stats), NA_integer_, as.integer(x$stats$flipped)),
      excluded = ifelse(is.null(x$stats), NA_integer_, as.integer(x$stats$excluded)),
      final_kept = ifelse(is.null(x$final_kept), 0L, as.integer(x$final_kept)),
      error = ifelse(is.null(x$error), "", as.character(x$error))
    )
  }), fill = TRUE)
  
  data.table::fwrite(summary_dt, file.path(out_root, "stats", "processing_summary.csv"))
  data.table::fwrite(summary_dt, file.path(out_root, "logs", "traits_log_summary.csv"))
}

unlink(temp_dir, recursive = TRUE)

# 11)  README 
readme_path <- file.path(out_root, "readme", "README.txt")
readme_lines <- c(
  "731 Immune Cell Preprocessing for SuSiE/coloc",
  "731 (SuSiE/coloc )",
  "==================================================",
  paste("Timestamp / :", timestamp),
  paste("Source Dir / :", src_dir),
  paste("LD SNP Info / LD:", ld_snp_info),
  paste("Sample Size File / :", sample_size_file),
  "",
  "Description / :",
  "This pipeline extracts GWAS summary statistics for immune cells and rigidly aligns them to the UKB LD Reference Panel.",
  "GWAS, UKB LD. ",
  "",
  "Alignment Logic (Global Alignment Rule) / :",
  "1) Consistent : Source_EA == LD_A1 & Source_OA == LD_A2 -> Beta = Beta, EAF = EAF",
  "2) Flipped : Source_EA == LD_A2 & Source_OA == LD_A1 -> Beta = -Beta, EAF = 1 - EAF",
  "3) Output : Output_A1 MUST be LD_A1, Output_A2 MUST be LD_A2.",
  "   (This is critical: since Beta represents the effect of Output_A1, A1 MUST be forced to match the LD reference!)",
  "   (: BetaLD_A1, A1LD_A1！)",
  "",
  "Features / :",
  "- AWK stream filtering: Filter SNPs by LD rsids dynamically during read to reduce memory usage.",
  "- Sample size completion: Reads 'STUDY' and 'samplesize' from external CSV and overrides N column.",
  "- Removed MAF/Position filtering: Retains full genome, intersection defined strictly by LD panel.",
  "  (AWK；N；MAF/, LD)",
  "",
  "Outputs / :",
  "- results/*.tsv.gz: Aligned output files ready for coloc /  ",
  "- stats/*: Alignment summaries and mapping traces / ",
  "- logs/*: Running logs / ",
  "- checkpoints/*: Checkpoints for resume capability / "
)
writeLines(readme_lines, readme_path)

end_time <- Sys.time()
cat("\nEnd Time / :", as.character(end_time), "\n")
cat("Total Elapsed / :", as.character(difftime(end_time, start_time, units = "mins")), "mins\n")
cat("Output Directory / :", out_root, "\n")
if (test_mode) {
  cat("\n*** : 20.  ***\n")
}
sink(NULL)
