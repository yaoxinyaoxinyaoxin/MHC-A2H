#!/usr/bin/env Rscript
# ==============================================================================
# [Script]: 7.9.1.2_FinnGen_GWAS_preprocessing_template.R
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

packages <- c("data.table", "dplyr", "stringr", "parallel", "readr")
for (pkg in packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
    library(pkg, character.only = TRUE)
  }
}

setDTthreads(1) # Worker  1 ,  parallel 

# ------------------------------------------------------------------------------
#  1:  (Define Paths and Setup Directories)
# ------------------------------------------------------------------------------
base_dir <- "./"

finngen_gwas_dir <- "./"
finngen_manifest_path <- "./"

# LD 
ld_snp_info <- "./"

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path(base_dir, paste0("2_Preprocessed_FinnGen_SuSiE_", timestamp))
results_dir <- file.path(out_dir, "results")
logs_dir <- file.path(out_dir, "logs")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, showWarnings = FALSE)
dir.create(logs_dir, showWarnings = FALSE)

log_file <- file.path(logs_dir, paste0("run_log_", timestamp, ".txt"))
log_msg <- function(msg) {
  time_stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", time_stamp, msg), file = log_file, append = TRUE)
  cat(sprintf("[%s] %s\n", time_stamp, msg))
}

log_msg("Script started.")

# ------------------------------------------------------------------------------
#  2:  Manifest  (Gather Targets & Parse Manifest)
# ------------------------------------------------------------------------------
log_msg("Step 2: Gathering FinnGen GWAS datasets and matching manifest...")

# 2.1  .gz 
gwas_files <- list.files(finngen_gwas_dir, pattern = "\\.gz$", full.names = TRUE)
log_msg(sprintf("Found %d .gz files in the target directory.", length(gwas_files)))

#  phenocode
# : finngen_R12_SLE_FG.gz
file_names <- basename(gwas_files)
phenocodes <- gsub("^finngen_R12_", "", file_names)
phenocodes <- gsub("\\.gz$", "", phenocodes)

all_targets <- data.table(
  FilePath = gwas_files,
  Phenocode = phenocodes
)

# 2.2  Manifest 
manifest_df <- fread(finngen_manifest_path)

#  all_targets 
all_targets <- merge(all_targets, manifest_df[, .(phenocode, phenotype, num_cases, num_controls, ``)],
                     by.x = "Phenocode", by.y = "phenocode", all.x = TRUE)

setnames(all_targets, 
         old = c("phenotype", "num_cases", "num_controls", ""), 
         new = c("Phenotype", "N_case", "N_control", "N_total"))

#  GCST  Phenocode 
all_targets[, GCST := Phenocode]

input_list_file <- file.path(out_dir, paste0("Input_Datasets_List_", timestamp, ".csv"))
fwrite(all_targets, input_list_file)
log_msg(paste("Saved input dataset list to:", input_list_file))

# ============================  ============================
# , 
# all_targets <- head(all_targets, 3)
# ==================================================================

# ------------------------------------------------------------------------------
#  3:  (Parallel Preprocessing & Alignment)
# ------------------------------------------------------------------------------
log_msg("Step 3: Loading UKB LD Reference Panel...")
ld_dt <- fread(ld_snp_info, sep = "\t")
setnames(ld_dt, old = c("SNP_ID","Chromosome","Position_BP","Allele1","Allele2","Allele_Direction"), 
             new = c("rsid","chr","pos","ld_a1","ld_a2","ld_dir"))
ld_dt[, rsid := as.character(rsid)]
ld_dt[, ld_a1 := toupper(ld_a1)]
ld_dt[, ld_a2 := toupper(ld_a2)]
log_msg(sprintf("Loaded LD Reference Panel with %d SNPs.", nrow(ld_dt)))

process_gwas_worker <- function(i, all_targets, ld_dt, results_dir) {
  library(data.table)
  library(dplyr)
  setDTthreads(1)
  
  gcst <- all_targets$GCST[i]
  file_path <- all_targets$FilePath[i]
  n_total <- all_targets$N_total[i]
  
  out_file <- file.path(results_dir, paste0(gcst, ".csv"))
  if (file.exists(out_file)) {
    return(list(GCST = gcst, Status = "Skipped (Exists)", Error = NA, Aligned_SNPs = NA))
  }
  
  tryCatch({
    #  awk : MAF > 0.05
    #  FinnGen 
    # FinnGen : #chrom, pos, ref, alt, rsids, nearest_genes, pval, mlogp, beta, sebeta, af_alt
    awk_script <- "
    BEGIN { FS=\"\\t\"; OFS=\"\\t\" }
    NR==1 {
      for(i=1; i<=NF; i++) {
        if($i==\"#chrom\") chr=i;
        if($i==\"pos\") pos=i;
        if($i==\"alt\") ea=i;
        if($i==\"ref\") oa=i;
        if($i==\"beta\") beta=i;
        if($i==\"sebeta\") se=i;
        if($i==\"pval\") pval=i;
        if($i==\"af_alt\") eaf=i;
        if($i==\"rsids\") rsid=i;
      }
      print \"rsid\", \"chr\", \"pos\", \"ea\", \"oa\", \"beta\", \"se\", \"pval\", \"eaf\"
    }
    NR>1 {
      chr_val = (chr != \"\") ? $chr : \"\";
      pos_val = (pos != \"\") ? $pos : \"\";
      ea_val = (ea != \"\") ? $ea : \"\";
      oa_val = (oa != \"\") ? $oa : \"\";
      beta_val = (beta != \"\") ? $beta : \"\";
      se_val = (se != \"\") ? $se : \"\";
      pval_val = (pval != \"\") ? $pval : \"\";
      eaf_val = (eaf != \"\") ? $eaf : \"\";
      rsid_val = (rsid != \"\") ? $rsid : \"\";
      
      my_eaf = eaf_val;
      if (my_eaf != \"\" && my_eaf != \"NA\") {
        my_maf = (my_eaf > 0.5) ? (1 - my_eaf) : my_eaf;
        if (my_maf > 0.05) {
          #  rsid , 
          split(rsid_val, arr, \",\")
          my_rsid = (arr[1] != \"\" && arr[1] != \"NA\") ? arr[1] : chr_val\":\"pos_val
          print my_rsid, chr_val, pos_val, ea_val, oa_val, beta_val, se_val, pval_val, eaf_val
        }
      }
    }
    "
    
    #  awk 
    cmd <- sprintf("gzcat '%s' | awk '%s'", file_path, awk_script)
    dat <- fread(cmd = cmd)
    
    if (nrow(dat) == 0) {
      stop("No SNPs remained after MAF > 0.05 filtering or reading failed.")
    }
    
    if (is.na(n_total)) {
        dat[, n_final := 487448] # , 
    } else {
        dat[, n_final := as.numeric(n_total)]
    }
    
    dat[, ea := toupper(ea)]
    dat[, oa := toupper(oa)]
    
    #  ATCG
    dat <- dat[ea %in% c("A","C","G","T") & oa %in% c("A","C","G","T")]
    
    #  ( LD )
    merged <- merge(dat, ld_dt, by = "rsid", all.x = FALSE) # Inner join
    
    merged[, status := fcase(
      ea == ld_a1 & oa == ld_a2, "consistent",
      ea == ld_a2 & oa == ld_a1, "flipped",
      default = "inconsistent"
    )]
    
    #  consistent  flipped
    merged <- merged[status %in% c("consistent", "flipped")]
    
    #  Beta  EAF
    merged[, beta := as.numeric(beta)]
    merged[, eaf := as.numeric(eaf)]
    
    merged[, final_beta := fifelse(status == "flipped", -beta, beta)]
    merged[, final_eaf := fifelse(status == "flipped", 1 - eaf, eaf)]
    
    # : A1A2LD
    out <- data.table(
      snp = merged$rsid,
      chr = merged$chr.x,
      pos = merged$pos.x,
      A1  = merged$ld_a1,  # 
      A2  = merged$ld_a2,  # 
      beta = merged$final_beta,
      se   = as.numeric(merged$se),
      pval = as.numeric(merged$pval),
      eaf  = merged$final_eaf,
      N    = merged$n_final
    )
    
    out <- out[!is.na(beta) & !is.na(se) & !is.na(pval) & !is.na(eaf)]
    
    fwrite(out, out_file)
    
    return(list(GCST = gcst, Status = "Success", Error = NA, Aligned_SNPs = nrow(out)))
    
  }, error = function(e) {
    return(list(GCST = gcst, Status = "Failed", Error = e$message, Aligned_SNPs = NA))
  })
}

#  (16 workers)
num_cores <- 16
log_msg(sprintf("Starting PSOCK cluster with %d workers...", num_cores))
cl <- makeCluster(num_cores)

clusterExport(cl, c("process_gwas_worker", "all_targets", "ld_dt", "results_dir"))

results_list <- parLapply(cl, 1:nrow(all_targets), function(i) {
  process_gwas_worker(i, all_targets, ld_dt, results_dir)
})

stopCluster(cl)

res_dt <- rbindlist(results_list)
fwrite(res_dt, file.path(out_dir, paste0("Processing_Summary_", timestamp, ".csv")))

log_msg("Parallel processing completed. Summary:")
print(res_dt)

log_msg("Script finished successfully.")
