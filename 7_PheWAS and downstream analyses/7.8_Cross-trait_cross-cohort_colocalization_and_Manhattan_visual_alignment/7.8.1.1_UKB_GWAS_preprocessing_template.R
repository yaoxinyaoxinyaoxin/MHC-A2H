#!/usr/bin/env Rscript
# ==============================================================================
# [Script]: 7.8.1.1_UKB_GWAS_preprocessing_template.R
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

prot_gwas_dir <- "./"
disease_gwas_dir <- "./"
csv_834_path <- "./"

hub_file <- "./"
gene_list_file <- "./"
sig_file <- "./"
disease_pheno_csv_path <- "./"

# LD 
ld_snp_info <- "./"

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
out_dir <- file.path(base_dir, paste0("Preprocessed_SuSiE_", timestamp))
results_dir <- file.path(out_dir, "results")
logs_dir <- file.path(out_dir, "logs")
readme_dir <- file.path(out_dir, "readme")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(results_dir, showWarnings = FALSE)
dir.create(logs_dir, showWarnings = FALSE)
dir.create(readme_dir, showWarnings = FALSE)

log_file <- file.path(logs_dir, paste0("run_log_", timestamp, ".txt"))
log_msg <- function(msg) {
  time_stamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", time_stamp, msg), file = log_file, append = TRUE)
  cat(sprintf("[%s] %s\n", time_stamp, msg))
}

log_msg("Script started.")

# ------------------------------------------------------------------------------
#  2:  (Gather Target Datasets)
# ------------------------------------------------------------------------------
log_msg("Step 2: Gathering target datasets (Proteins and Diseases)...")

# --- 2.1  ---
hub_df <- fread(hub_file)
hub_genes <- hub_df[Type == "Hub", Gene]
gene_list <- fread(gene_list_file)
hub_phenotypes <- gene_list[gene %in% hub_genes, Phenotype]
sig_df <- fread(sig_file)
hub_gcst_map <- unique(sig_df[Phenotype %in% hub_phenotypes, .(GCST, Phenotype)])

manual_gcst_prot <- data.table(
  GCST = c("GCST90468120", "GCST90468136", "GCST90468138", "GCST90468112", "GCST90468147", 
           "GCST90468082", "GCST90468090", "GCST90468092", "GCST90468068", "GCST90468100", 
           "GCST90468095", "GCST90468098", "GCST90468076"),
  Phenotype = c("Coeliac disease", "Hyperthyroidism", "Hypothyroidism", "Asthma", "Psoriasis",
                "Lymphocyte count", "Monocyte count", "Neutrophill count", "Eosinophill count",
                "Reticulocyte count", "Platelet count", "Red blood cell erythrocyte count",
                "High light scatter reticulocyte count")
)
target_proteins <- rbind(hub_gcst_map, manual_gcst_prot)
target_proteins <- unique(target_proteins, by = "GCST")
target_proteins[, Source := "Proteins/Traits"]
target_proteins[, FilePath := file.path(prot_gwas_dir, paste0(GCST, ".tsv.gz"))]
target_proteins[, Disease_Name := Phenotype]

# --- 2.2  ---
disease_gcst_list <- c(
  "GCST90473161", "GCST90473888", "GCST90473152", "GCST90473157", "GCST90473169", 
  "GCST90473148", "GCST90473145", "GCST90474048", "GCST90474052", "GCST90473937", 
  "GCST90474055", "GCST90473712", "GCST90473130", "GCST90474159", "GCST90473113", 
  "GCST90473962", "GCST90473137", "GCST90473935", "GCST90473050", "GCST90473036", 
  "GCST90473953", "GCST90473167", "GCST90473375", "GCST90473599", "GCST90473115", 
  "GCST90473180", "GCST90473806", "GCST90473182", "GCST90473804", "GCST90473984", 
  "GCST90473921", "GCST90473701", "GCST90473235", "GCST90473430", "GCST90474176", 
  "GCST90473961", "GCST90474170",
  "GCST90474461", "GCST90474466", "GCST90474501", "GCST90474521", "GCST90474526", 
  "GCST90474531", "GCST90474536", "GCST90474541", "GCST90474576", "GCST90474601"
)
target_diseases <- data.table(GCST = disease_gcst_list, Phenotype = NA_character_)

if (file.exists(disease_pheno_csv_path)) {
  pheno_df <- fread(disease_pheno_csv_path, select = c("GCST", "Disease_Name"))
  pheno_df <- unique(pheno_df)
  target_diseases <- merge(target_diseases, pheno_df, by = "GCST", all.x = TRUE)
} else {
  target_diseases[, Disease_Name := GCST]
}
target_diseases[, Phenotype := Disease_Name]
target_diseases[, Source := "UKB_Diseases"]
target_diseases[, FilePath := file.path(disease_gwas_dir, paste0(GCST, "_hg37_converted.tsv.gz"))] # 

# （ _hg37_converted.tsv.gz  .tsv.gz）
target_diseases$FilePath <- sapply(target_diseases$GCST, function(x) {
  p1 <- file.path(disease_gwas_dir, paste0(x, "_hg37_converted.tsv.gz"))
  p2 <- file.path(disease_gwas_dir, paste0(x, ".tsv.gz"))
  if (file.exists(p1)) return(p1)
  if (file.exists(p2)) return(p2)
  return(p2) #  .tsv.gz
})

# ---  ---
all_targets <- rbind(target_proteins, target_diseases, fill = TRUE)
all_targets <- all_targets[file.exists(FilePath)]
log_msg(sprintf("Total target datasets found on disk: %d", nrow(all_targets)))

# ------------------------------------------------------------------------------
#  3:  834_UKB_.csv  (Parse Sample Size)
# ------------------------------------------------------------------------------
log_msg("Step 3: Parsing sample sizes from 834 UKB CSV...")
csv_834 <- read_csv(csv_834_path, show_col_types = FALSE)
setDT(csv_834)
csv_834 <- csv_834[!is.na(`STUDY ACCESSION`)]

#  cases  controls
# : "91,830 Non-Finnish European ancestry cases, 366,610 Non-Finnish European ancestry controls"
parse_sample_size <- function(size_str) {
  if (is.na(size_str)) return(list(case = NA_integer_, control = NA_integer_, total = NA_integer_))
  
  case_match <- str_match(size_str, "([0-9,]+)\\s+[^c]*cases")
  control_match <- str_match(size_str, "([0-9,]+)\\s+[^c]*controls")
  
  case_val <- ifelse(!is.na(case_match[,2]), as.integer(gsub(",", "", case_match[,2])), NA_integer_)
  control_val <- ifelse(!is.na(control_match[,2]), as.integer(gsub(",", "", control_match[,2])), NA_integer_)
  
  #  cases/controls,  (e.g., "400,000 Non-Finnish European ancestry individuals")
  total_val <- case_val + control_val
  if (is.na(total_val)) {
    ind_match <- str_match(size_str, "([0-9,]+)\\s+[^i]*individuals")
    if (!is.na(ind_match[,2])) {
      total_val <- as.integer(gsub(",", "", ind_match[,2]))
    }
  }
  
  return(list(case = case_val, control = control_val, total = total_val))
}

size_info <- rbindlist(lapply(csv_834$`INITIAL SAMPLE SIZE`, parse_sample_size))
csv_834 <- cbind(csv_834, size_info)

#  all_targets
all_targets <- merge(all_targets, csv_834[, .(`STUDY ACCESSION`, case, control, total)], 
                     by.x = "GCST", by.y = "STUDY ACCESSION", all.x = TRUE)

setnames(all_targets, c("case", "control", "total"), c("N_case", "N_control", "N_total"))

input_list_file <- file.path(readme_dir, paste0("Input_Datasets_List_", timestamp, ".csv"))
fwrite(all_targets, input_list_file)
log_msg(paste("Saved input dataset list to:", input_list_file))

# ============================  ============================
# , 
# all_targets <- head(all_targets, 3)
# ==================================================================

# ------------------------------------------------------------------------------
#  4:  (Parallel Preprocessing & Alignment)
# ------------------------------------------------------------------------------
log_msg("Step 4: Loading UKB LD Reference Panel...")
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
    awk_script <- "
    BEGIN { FS=\"\\t\"; OFS=\"\\t\" }
    NR==1 {
      for(i=1; i<=NF; i++) {
        if($i==\"chromosome\") chr=i;
        if($i==\"base_pair_location\") pos=i;
        if($i==\"effect_allele\") ea=i;
        if($i==\"other_allele\") oa=i;
        if($i==\"beta\") beta=i;
        if($i==\"standard_error\") se=i;
        if($i==\"p_value\") pval=i;
        if($i==\"effect_allele_frequency\" || $i==\"eaf\") eaf=i;
        if($i==\"rsid\" || $i==\"rs_id\") rsid=i;
        if($i==\"variant_id\") varid=i;
        if($i==\"n\") n=i;
      }
      print \"rsid\", \"chr\", \"pos\", \"ea\", \"oa\", \"beta\", \"se\", \"pval\", \"eaf\", \"n\"
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
      varid_val = (varid != \"\") ? $varid : \"\";
      n_val = (n != \"\") ? $n : \"\";
      
      my_eaf = eaf_val;
      if (my_eaf != \"\" && my_eaf != \"NA\") {
        my_maf = (my_eaf > 0.5) ? (1 - my_eaf) : my_eaf;
        if (my_maf > 0.05) {
          my_rsid = (rsid_val != \"\" && rsid_val != \"NA\") ? rsid_val : ((varid_val != \"\" && varid_val != \"NA\") ? varid_val : chr_val\":\"pos_val)
          print my_rsid, chr_val, pos_val, ea_val, oa_val, beta_val, se_val, pval_val, eaf_val, n_val
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
    
    #  N_total,  n 
    if (is.na(n_total)) {
      if (!all(is.na(dat$n))) {
        dat[, n_final := as.numeric(n)]
      } else {
        # ,  (UKB 394642  54219)
        dat[, n_final := 394642] # 
      }
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
fwrite(res_dt, file.path(readme_dir, paste0("Processing_Summary_", timestamp, ".csv")))

log_msg("Parallel processing completed. Summary:")
print(res_dt)

# ==============================================================================
#  Readme
# ==============================================================================
readme_text <- sprintf("
# GWAS Preprocessing for SuSiE/coloc
# : %s

##  (Description)
 GWAS , （ SuSiE, coloc）. 
（ LD ）. 

## 
1.  awk ,  MAF > 0.05 , . 
2.  RSID  UKB EUR LD . 
3. : 
   - Consistent (Source EA == LD A1): Beta  EAF . 
   - Flipped (Source EA == LD A2): Beta , EAF  (1 - EAF). 
4.  A1  A2  LD  ld_a1  ld_a2. 

## 
- : %d
- : %d
- : %d
- : `Processing_Summary_%s.csv`  `Input_Datasets_List_%s.csv`

## 
- : %s
- : %s
", 
timestamp, 
nrow(res_dt), sum(res_dt$Status == "Success"), sum(res_dt$Status == "Failed"),
timestamp, timestamp,
input_list_file, results_dir)

writeLines(readme_text, file.path(readme_dir, paste0("README_", timestamp, ".md")))

log_msg("Script finished successfully.")
