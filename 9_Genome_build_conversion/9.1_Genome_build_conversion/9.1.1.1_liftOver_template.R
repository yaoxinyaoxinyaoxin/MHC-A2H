#!/usr/bin/env Rscript
# ==============================================================================
# [Script]: 9.1.1.1_liftOver_template.R
# [Method]: Genome Assembly Conversion (LiftOver)
# [Step]: Convert genomic coordinates (e.g. hg38 to hg19)
# 
# [Function]:
# Convert GWAS/eQTL summary statistics between genome assemblies using liftOver.
# 
# [Input]:
#   --input_file   : Path to input summary statistics file
#   --chain_file   : Path to UCSC chain file
#   --out_dir      : Output directory path
# 
# [Output]:
# Converted summary statistics and unmapped records.
# ==============================================================================

rm(list = ls())
gc()

# 1. Load Dependencies
suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(dplyr)
  library(rtracklayer)
  library(GenomicRanges)
})

# 2. Parse Command Line Arguments
option_list <- list(
  make_option(c("--input_file"), type="character", default=NULL, help="Input summary stats file"),
  make_option(c("--chain_file"), type="character", default=NULL, help="Path to UCSC chain file (e.g. hg38ToHg19.over.chain)"),
  make_option(c("--out_dir"), type="character", default="./LiftOver_Results", help="Output directory path")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_file) || is.null(opt$chain_file)) {
  print_help(opt_parser)
  stop("Missing required input file or chain file.")
}

INPUT_FILE <- opt$input_file
CHAIN_FILE <- opt$chain_file
OUTPUT_DIR <- opt$out_dir

dir.create(OUTPUT_DIR, recursive = TRUE, showWarnings = FALSE)

options(stringsAsFactors = FALSE)

# Load libraries / 
suppressPackageStartupMessages({
  library(data.table)
  library(rtracklayer)
  library(GenomicRanges)
  library(R.utils)
})

# Define paths / 
project_dir <- "./"
input_file <- file.path(project_dir, "1_data/_/1__hg38/finngen_R12_AB1_ZOSTER.gz")
output_dir <- file.path(project_dir, "1_data/_/2_hg19")
output_file <- file.path(output_dir, "finngen_R12_AB1_ZOSTER.gz")
chain_file <- CHAIN_FILE
log_dir <- file.path(project_dir, "logs")

# Initialize log file / 
log_file <- file.path(log_dir, paste0("liftOver_hg38_to_hg19_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
sink(log_file, split = TRUE)

cat(sprintf("[%s] Script started.\n", Sys.time()))
cat(sprintf("Input file: %s\n", input_file))
cat(sprintf("Output file: %s\n", output_file))
cat(sprintf("Chain file: %s\n", chain_file))

# ---------------------------------------------------------------------------------------------------------------------
# 1. Read Data / 
# ---------------------------------------------------------------------------------------------------------------------
cat(sprintf("[%s] Reading input data...\n", Sys.time()))
# Read data using fread for efficiency / fread
dt <- fread(input_file)

# Display data structure / 
cat("Original data head:\n")
print(head(dt))
cat(sprintf("Original row count: %d\n", nrow(dt)))

# ---------------------------------------------------------------------------------------------------------------------
# 2. Prepare for LiftOver / LiftOver
# ---------------------------------------------------------------------------------------------------------------------
cat(sprintf("[%s] Preparing data for liftOver...\n", Sys.time()))

# Ensure chromosome column is character / 
dt[, `#chrom` := as.character(`#chrom`)]

# Add 'chr' prefix if missing (required for chain file) / 'chr'（chain）
# Check if chromosomes already have 'chr' / 'chr'
if (!any(grepl("^chr", dt$`#chrom`))) {
  dt[, chrom_u := paste0("chr", `#chrom`)]
} else {
  dt[, chrom_u := `#chrom`]
}

# Handle 'chr23', 'chr24', 'chrMT' mappings if necessary / 
# Finngen usually uses 23 for X, 24 for Y, 25 for XY, 26 for MT. Need to check specific mapping.
# For now, assume standard 1-22, X, Y, MT or numerical.
# If numeric 23 -> chrX, 24 -> chrY, 26 -> chrM
dt[chrom_u == "chr23", chrom_u := "chrX"]
dt[chrom_u == "chr24", chrom_u := "chrY"]
dt[chrom_u == "chr26", chrom_u := "chrM"]
dt[chrom_u == "chrMT", chrom_u := "chrM"]

# Create GRanges object / GRanges
gr <- GRanges(
  seqnames = dt$chrom_u,
  ranges = IRanges(start = dt$pos, end = dt$pos),
  mcols = dt[, .SD, .SDcols = !c("chrom_u")] # Keep other columns / 
)

# ---------------------------------------------------------------------------------------------------------------------
# 3. Perform LiftOver / LiftOver
# ---------------------------------------------------------------------------------------------------------------------
cat(sprintf("[%s] Importing chain file and performing liftOver...\n", Sys.time()))
chain <- import.chain(chain_file)
gr_lifted <- liftOver(gr, chain)

# Convert back to data.table / data.table
# liftOver returns a GRangesList. Unlist it. / liftOverGRangesList, . 
gr_lifted_unlist <- unlist(gr_lifted)

# Create a mapping index to original data / 
# elementMetadata contains the original columns if we attached them? No, we put them in mcols but liftOver might drop them if not careful?
# liftOver preserves mcols.
# But 'unlist' might reorder? No, unlist keeps order of elements found.
# Better approach: 
# The unlisted object has names corresponding to original indices if we set them, or we can use the mapping.
# Actually, it's safer to reconstruct data.table from gr_lifted_unlist directly since it carries mcols.

dt_lifted <- as.data.table(gr_lifted_unlist)

# 'seqnames' is the new chromosome, 'start' is the new position.
# mcols columns are preserved with prefix 'mcols.' or just as columns depending on conversion.
# Let's check names
cat("Lifted data columns:\n")
print(colnames(dt_lifted))

# ---------------------------------------------------------------------------------------------------------------------
# 4. Format Output / 
# ---------------------------------------------------------------------------------------------------------------------
cat(sprintf("[%s] Formatting output...\n", Sys.time()))

# Rename columns back to original standard / 
# Expected: #chrom, pos, ...
# dt_lifted has: seqnames, start, end, width, strand, mcols...
# We need to extract original columns which are inside 'mcols' (if as.data.table flattened it correctly).
# as.data.table on GRanges usually flattens mcols.
# Let's rename 'seqnames' to '#chrom' and 'start' to 'pos'.

# Check if original columns are present (they should be prefixed with mcols. or just names)
# If mcols were passed in GRanges constructor, they are columns in dt_lifted.
# Note: '#chrom' was in input, but we used 'chrom_u' for seqnames. The original '#chrom' might be in mcols if we included it.
# In step 2, we did: mcols = dt[, .SD, .SDcols = !c("chrom_u")]
# So original '#chrom' is in mcols.
# But we want the NEW coordinates. So we take new seqnames and new start.

# Update #chrom and pos
# Remove 'chr' prefix from new seqnames to match Finngen style (usually numbers)
new_chrom <- as.character(dt_lifted$seqnames)
new_chrom <- gsub("^chr", "", new_chrom)
new_chrom[new_chrom == "X"] <- "23"
new_chrom[new_chrom == "Y"] <- "24"
new_chrom[new_chrom == "M"] <- "26"
new_chrom[new_chrom == "MT"] <- "26"

# Rename mcols columns to original names / mcols
current_cols <- colnames(dt_lifted)
mcols_cols <- grep("^mcols\\.", current_cols, value = TRUE)
new_names <- gsub("^mcols\\.", "", mcols_cols)
# Handle special case if #chrom becomes .chrom or similar due to R naming
# In the log: "mcols..chrom". R usually replaces # with . or check.make.names
# Original was "#chrom". make.names("#chrom") -> "X.chrom". data.table might handle differently.
# Log says: "mcols..chrom" (two dots? or dot then hash?)
# Let's check the log output: "mcols..chrom".
# It seems `mcols.` + `#chrom` -> `mcols.#chrom` but printed/converted might be tricky.
# Let's rely on setnames by matching patterns.

# Remove 'mcols.' prefix
setnames(dt_lifted, mcols_cols, new_names)

# If there is a column named ".chrom", rename it to "#chrom" if that was the original
if (".chrom" %in% colnames(dt_lifted) && !("#chrom" %in% colnames(dt_lifted))) {
    setnames(dt_lifted, ".chrom", "#chrom")
}

# Update #chrom and pos with lifted values
dt_lifted[, `#chrom` := new_chrom]
dt_lifted[, pos := start]

# Keep only original columns / 
original_cols <- colnames(dt)[colnames(dt) != "chrom_u"]
# Check if all original cols are present in dt_lifted (some might be renamed if duplicate)
# data.table as.data.table(GRanges) usually keeps mcols names unless conflict.
# 'pos' and '#chrom' were in mcols? 
# In step 2: mcols = dt[, .SD, .SDcols = !c("chrom_u")]
# dt had '#chrom' and 'pos'.
# GRanges has 'seqnames', 'ranges' (start/end). 
# So '#chrom' and 'pos' from original dt are in mcols.
# But 'pos' conflicts with GRanges start? No, GRanges doesn't have 'pos'.
# '#chrom' is fine.
# Wait, if we included '#chrom' and 'pos' in mcols, they are now in dt_lifted.
# We should update them with the NEW values.
# But wait, dt_lifted already has 'mcols.#chrom' and 'mcols.pos' (or just '#chrom' and 'pos').
# If they exist, we overwrite them.

# Let's clean up.
# Identify columns to keep.
cols_to_keep <- original_cols
# Ensure they exist
missing_cols <- setdiff(cols_to_keep, colnames(dt_lifted))
if (length(missing_cols) > 0) {
  cat("Warning: Some original columns missing or renamed in lifted object. Checking...\n")
  print(colnames(dt_lifted))
  # Try to find them. maybe prefixed.
  # If we can't find them, we might have an issue.
  # But usually as.data.table works fine.
}

# Select and order columns
dt_final <- dt_lifted[, ..cols_to_keep]

# Sort by chrom and pos / 
# Convert chrom to numeric for sorting if possible, or keep character
# Finngen chrom is usually 1, 2...
dt_final[, chrom_num := as.numeric(`#chrom`)]
setorder(dt_final, chrom_num, pos)
dt_final[, chrom_num := NULL]

cat("Lifted data head:\n")
print(head(dt_final))
cat(sprintf("Lifted row count: %d (Original: %d, Lost: %d)\n", nrow(dt_final), nrow(dt), nrow(dt) - nrow(dt_final)))

# ---------------------------------------------------------------------------------------------------------------------
# 5. Save Output / 
# ---------------------------------------------------------------------------------------------------------------------
cat(sprintf("[%s] Saving output to %s...\n", Sys.time(), output_file))
fwrite(dt_final, output_file, sep = "\t", quote = FALSE, na = "NA", compress = "gzip")

# ---------------------------------------------------------------------------------------------------------------------
# 6. Finalize / 
# ---------------------------------------------------------------------------------------------------------------------
end_time <- Sys.time()
cat(sprintf("[%s] Script completed successfully.\n", end_time))
# cat(sprintf("Total execution time: %s\n", end_time - start_time)) # start_time not defined but Sys.time() diff works visually
sink()
