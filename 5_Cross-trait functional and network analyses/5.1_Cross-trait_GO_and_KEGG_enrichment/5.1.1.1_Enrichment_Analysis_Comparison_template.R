# ==============================================================================
# [Script]: 5.1.1.1_Enrichment_Analysis_Comparison_template.R
# [Method]: Cross-trait GO/KEGG Enrichment
# [Step]: 5.1.1.1_Enrichment_Analysis_Comparison
#
# [Function]:
# Conduct functional enrichment analysis (GO/KEGG) and construct Protein-Protein Interaction (PPI) networks.
#
# [Usage]: 
#   Rscript 5.1.1.1_Enrichment_Analysis_Comparison_template.R \
#     --input_aging <path> \
#     --input_ra <path> \
#     --input_hz <path> \
#     --out_dir <path>
# ==============================================================================


# -------------------------------------------------------------------------
# 1. Initialization / 
# -------------------------------------------------------------------------
# Clear environment / 
rm(list = ls())
gc()

# Set specific CRAN mirror / CRAN
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

# Add local library path for aPEAR

# Function to check and install packages
check_and_install <- function(pkg) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message(paste0("Installing missing package: ", pkg))
    tryCatch({
      install.packages(pkg)
    }, error = function(e) {
      message(paste0("Failed to install ", pkg, ": ", e$message))
    })
  }
}

# Check for ggVennDiagram
check_and_install("ggVennDiagram")

# Load required libraries / 
suppressPackageStartupMessages({
  library(optparse)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(readr)
  library(stringr)
  library(openxlsx)
  library(ggpubr)
  library(scales) # Load scales for oob handling
  
  if(requireNamespace("ggVennDiagram", quietly = TRUE)) {
    library(ggVennDiagram)
  } else {
    warning("ggVennDiagram package not found. Venn diagrams might fail.")
  }
  
  # Load aPEAR and arules (Required for network visualization)
  if(requireNamespace("aPEAR", quietly = TRUE)) {
    library(aPEAR)
    # Ensure arules is loaded (Best practice for efficient similarity calculation)
    if(requireNamespace("arules", quietly = TRUE)) {
      library(arules)
    } else {
      warning("arules package not found. Performance might be affected.")
    }
  } else {
    warning("aPEAR package not found. Network plots will be skipped.")
  }
})

# Define command line arguments
option_list <- list(
  make_option(c("--input_aging"), type="character", default=NULL,
              help="Input text file containing aging gene list / "),
  make_option(c("--input_ra"), type="character", default=NULL,
              help="Input text file containing RA gene list / RA"),
  make_option(c("--input_hz"), type="character", default=NULL,
              help="Input text file containing HZ gene list / HZ"),
  make_option(c("--out_dir"), type="character", default="./Enrichment_Comparison_Result",
              help="Output directory path / ")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_aging) || is.null(opt$input_ra) || is.null(opt$input_hz)) {
  print_help(opt_parser)
  stop("Missing required input files. / . ", call.=FALSE)
}

# Define paths / 
work_dir <- opt$out_dir

# Create output directory with timestamp / 
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- file.path(work_dir, paste0("Analysis_Result_", timestamp))
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

# Initialize log file / 
log_file <- file.path(output_dir, "analysis_log.txt")
sink(log_file, split = TRUE)

cat(paste0("Analysis started at: ", Sys.time(), "
"))
cat(paste0("Output directory: ", output_dir, "
"))

# -------------------------------------------------------------------------
# 2. Load Data / 
# -------------------------------------------------------------------------
cat("
[Step 2] Loading gene lists...
")

# File paths
file_aging <- opt$input_aging
file_ra <- opt$input_ra
file_shingles <- opt$input_hz


# Read files (assuming no header based on preview)
genes_aging_raw <- read_lines(file_aging)
genes_ra_raw <- read_lines(file_ra)
genes_shingles_raw <- read_lines(file_shingles)

# Clean gene symbols (remove empty strings)
genes_aging <- genes_aging_raw[genes_aging_raw != ""]
genes_ra <- genes_ra_raw[genes_ra_raw != ""]
genes_shingles <- genes_shingles_raw[genes_shingles_raw != ""]

cat(paste0("Loaded ", length(genes_aging), " genes for mvAge.\n"))
cat(paste0("Loaded ", length(genes_ra), " genes for Rheumatoid Arthritis.\n"))
cat(paste0("Loaded ", length(genes_shingles), " genes for HZ (UKB).\n"))

# Preview data / 
print(head(genes_aging))
print(head(genes_ra))
print(head(genes_shingles))

# ID Conversion (Symbol to Entrez) / ID
cat("\nConverting Gene Symbols to Entrez IDs...\n")

convert_ids <- function(genes, name) {
  tryCatch({
    eg <- bitr(genes, fromType="SYMBOL", toType="ENTREZID", OrgDb="org.Hs.eg.db")
    cat(paste0(name, ": Successfully converted ", nrow(eg), " genes out of ", length(genes), "\n"))
    return(eg)
  }, error = function(e) {
    cat(paste0("Error converting IDs for ", name, ": ", e$message, "\n"))
    return(NULL)
  })
}

eg_aging <- convert_ids(genes_aging, "mvAge")
eg_ra <- convert_ids(genes_ra, "Rheumatoid Arthritis")
eg_shingles <- convert_ids(genes_shingles, "HZ (UKB)")

if (is.null(eg_aging) || is.null(eg_ra) || is.null(eg_shingles)) {
  stop("Gene ID conversion failed for one or more lists. Please check input gene symbols.")
}

# -------------------------------------------------------------------------
# 3. Enrichment Analysis / 
# -------------------------------------------------------------------------
cat("\n[Step 3] Performing Enrichment Analysis...\n")

# Function for Enrichment / 
run_enrichment <- function(gene_df, name) {
  res_list <- list()
  
  # Retry helper function
  try_with_retry <- function(expr, max_attempts = 3, wait_seconds = 5) {
    attempt <- 1
    while (attempt <= max_attempts) {
      result <- tryCatch({
        expr
      }, error = function(e) {
        message(paste0("Attempt ", attempt, " failed: ", e$message))
        return(NULL)
      })
      
      if (!is.null(result)) return(result)
      
      attempt <- attempt + 1
      if (attempt <= max_attempts) {
        message(paste0("Retrying in ", wait_seconds, " seconds..."))
        Sys.sleep(wait_seconds)
      }
    }
    stop("Max attempts reached. Operation failed.")
  }
  
  # GO BP
  cat(paste0("Running GO BP for ", name, "...\n"))
  ego_bp <- try_with_retry({
    enrichGO(gene = gene_df$ENTREZID,
             OrgDb = org.Hs.eg.db,
             ont = "BP",
             pAdjustMethod = "BH",
             pvalueCutoff = 1, # Relaxed to include all results as requested / 
             qvalueCutoff = 1,
             readable = TRUE)
  })
  res_list[["GO_BP"]] <- ego_bp
  
  # GO CC
  cat(paste0("Running GO CC for ", name, "...\n"))
  ego_cc <- try_with_retry({
    enrichGO(gene = gene_df$ENTREZID,
             OrgDb = org.Hs.eg.db,
             ont = "CC",
             pAdjustMethod = "BH",
             pvalueCutoff = 1,
             qvalueCutoff = 1,
             readable = TRUE)
  })
  res_list[["GO_CC"]] <- ego_cc
  
  # GO MF
  cat(paste0("Running GO MF for ", name, "...\n"))
  ego_mf <- try_with_retry({
    enrichGO(gene = gene_df$ENTREZID,
             OrgDb = org.Hs.eg.db,
             ont = "MF",
             pAdjustMethod = "BH",
             pvalueCutoff = 1,
             qvalueCutoff = 1,
             readable = TRUE)
  })
  res_list[["GO_MF"]] <- ego_mf
  
  # KEGG
  cat(paste0("Running KEGG for ", name, "...\n"))
  kk <- tryCatch({
    try_with_retry({
      enrichKEGG(gene = gene_df$ENTREZID,
                 organism = 'hsa',
                 pAdjustMethod = "BH",
                 pvalueCutoff = 1,
                 qvalueCutoff = 1)
    }, max_attempts = 5, wait_seconds = 10) # More retries for KEGG
  }, error = function(e) {
    message(paste0("KEGG enrichment failed for ", name, ": ", e$message))
    return(NULL)
  })
  
  if (!is.null(kk)) {
    kk <- setReadable(kk, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  }
  res_list[["KEGG"]] <- kk
  
  return(res_list)
}

res_aging <- run_enrichment(eg_aging, "mvAge")
res_ra <- run_enrichment(eg_ra, "Rheumatoid Arthritis")
res_shingles <- run_enrichment(eg_shingles, "HZ (UKB)")

# -------------------------------------------------------------------------
# 4. Visualization & Comparison / 
# -------------------------------------------------------------------------
cat("\n[Step 4] Generating Visualizations...\n")

# Helper function for bubble plot / 
plot_bubble <- function(enrich_res, title, p_limits, count_limits, p_breaks, count_breaks) {
  if (is.null(enrich_res) || nrow(enrich_res) == 0) {
    return(ggplot() + annotate("text", x=1, y=1, label="No results") + theme_void())
  }
  
  # Filter for significant results for bubble plot / 
  enrich_res_sig <- enrich_res
  enrich_res_sig@result <- enrich_res@result[enrich_res@result$p.adjust < 0.05, ]
  
  if (nrow(enrich_res_sig) == 0) {
    return(ggplot() + annotate("text", x=1, y=1, label="No significant results (p.adj < 0.05)") + theme_void())
  }
  
  # Create dotplot
  p <- dotplot(enrich_res_sig, showCategory=10)
  
  # Customize
  p <- p + 
    scale_color_gradient(low="red", high="blue", limits=p_limits, breaks=p_breaks, guide = guide_colorbar(reverse = TRUE)) + # Red is low p (more significant)
    scale_size_continuous(limits=count_limits, breaks=count_breaks) +
    scale_x_continuous(n.breaks = 4) + # Limit x-axis ticks to 4 to prevent overlap
    scale_y_discrete(labels = function(x) stringr::str_wrap(x, width = 40)) + # Wrap long labels
    theme(
      plot.title = element_blank(),
      axis.text.x = element_text(angle = 45, hjust = 1), # Rotate x-axis labels to prevent overlap
      axis.text.y = element_text(size=8),
      # Legend at bottom-right inside the plot / 
      legend.position = c(1, 0), 
      legend.justification = c(1, 0),
      legend.background = element_rect(fill = alpha("white", 0.6), color = NA), # Transparent background
      legend.text = element_text(size = 6),
      legend.title = element_text(size = 7),
      legend.key.size = unit(0.3, "cm"),
      legend.margin = margin(2, 2, 2, 2)
    ) +
    labs(color = "p.adjust", size = "Count")
    
  return(p)
}

# Helper function for Venn diagram /  (Updated for 3 sets)
plot_venn_3way <- function(res1, res2, res3, name1="mvAge", name2="RA", name3="HZ (UKB)", category) {
  
  # Filter for significant results (p.adjust < 0.05) / 
  res1_sig <- if (!is.null(res1)) res1@result[res1@result$p.adjust < 0.05, ] else NULL
  res2_sig <- if (!is.null(res2)) res2@result[res2@result$p.adjust < 0.05, ] else NULL
  res3_sig <- if (!is.null(res3)) res3@result[res3@result$p.adjust < 0.05, ] else NULL

  # Extract IDs
  ids1 <- if (!is.null(res1_sig) && nrow(res1_sig) > 0) res1_sig$ID else c()
  ids2 <- if (!is.null(res2_sig) && nrow(res2_sig) > 0) res2_sig$ID else c()
  ids3 <- if (!is.null(res3_sig) && nrow(res3_sig) > 0) res3_sig$ID else c()
  
  venn_list <- list(
    mvAge = ids1,
    RA = ids2,
    HZ = ids3
  )
  names(venn_list) <- c(name1, name2, name3)
  
  if(requireNamespace("ggVennDiagram", quietly = TRUE)) {
    p <- ggVennDiagram(venn_list, label_alpha = 0) + 
      scale_fill_gradient(low = "white", high = "red") +
      theme(legend.position = "none", plot.title = element_blank())
    return(p)
  } else {
    return(ggplot() + annotate("text", x=1, y=1, label="ggVennDiagram not installed") + theme_void())
  }
}

# Generate plots
# Calculate global limits for consistent legends (Based on significant results only)
all_enrich <- c(res_aging, res_ra, res_shingles)
# Extract p.adjust and Count for significant results (p < 0.05)
all_p <- unlist(lapply(all_enrich, function(x) {
  if(!is.null(x) && nrow(x)>0) {
    p_vals <- x@result$p.adjust
    return(p_vals[p_vals < 0.05])
  } else {
    return(NULL)
  }
}))
all_count <- unlist(lapply(all_enrich, function(x) {
  if(!is.null(x) && nrow(x)>0) {
    counts <- x@result$Count
    p_vals <- x@result$p.adjust
    return(counts[p_vals < 0.05])
  } else {
    return(NULL)
  }
}))

min_p <- if(length(all_p)>0) min(all_p, na.rm=TRUE) else 0.05
max_p <- if(length(all_p)>0) max(all_p, na.rm=TRUE) else 0.05
min_count <- if(length(all_count)>0) min(all_count, na.rm=TRUE) else 0
max_count <- if(length(all_count)>0) max(all_count, na.rm=TRUE) else 10

p_limits <- c(min_p, max_p)
count_limits <- c(min_count, max_count)
p_breaks <- pretty(p_limits, n=4)
count_breaks <- pretty(count_limits, n=4)

plot_list <- list()

# Row 1: mvAge
plot_list[[1]] <- plot_bubble(res_aging[["GO_BP"]], NULL, p_limits, count_limits, p_breaks, count_breaks)
plot_list[[2]] <- plot_bubble(res_aging[["GO_CC"]], NULL, p_limits, count_limits, p_breaks, count_breaks)
plot_list[[3]] <- plot_bubble(res_aging[["GO_MF"]], NULL, p_limits, count_limits, p_breaks, count_breaks)
plot_list[[4]] <- plot_bubble(res_aging[["KEGG"]],  NULL, p_limits, count_limits, p_breaks, count_breaks)

# Row 2: Rheumatoid Arthritis
plot_list[[5]] <- plot_bubble(res_ra[["GO_BP"]], NULL, p_limits, count_limits, p_breaks, count_breaks)
plot_list[[6]] <- plot_bubble(res_ra[["GO_CC"]], NULL, p_limits, count_limits, p_breaks, count_breaks)
plot_list[[7]] <- plot_bubble(res_ra[["GO_MF"]], NULL, p_limits, count_limits, p_breaks, count_breaks)
plot_list[[8]] <- plot_bubble(res_ra[["KEGG"]],  NULL, p_limits, count_limits, p_breaks, count_breaks)

# Row 3: HZ (UKB)
plot_list[[9]] <- plot_bubble(res_shingles[["GO_BP"]], NULL, p_limits, count_limits, p_breaks, count_breaks)
plot_list[[10]] <- plot_bubble(res_shingles[["GO_CC"]], NULL, p_limits, count_limits, p_breaks, count_breaks)
plot_list[[11]] <- plot_bubble(res_shingles[["GO_MF"]], NULL, p_limits, count_limits, p_breaks, count_breaks)
plot_list[[12]] <- plot_bubble(res_shingles[["KEGG"]],  NULL, p_limits, count_limits, p_breaks, count_breaks)

# Row 4: Venn Diagrams (Comparison)
plot_list[[13]] <- plot_venn_3way(res_aging[["GO_BP"]], res_ra[["GO_BP"]], res_shingles[["GO_BP"]], category="GO BP")
plot_list[[14]] <- plot_venn_3way(res_aging[["GO_CC"]], res_ra[["GO_CC"]], res_shingles[["GO_CC"]], category="GO CC")
plot_list[[15]] <- plot_venn_3way(res_aging[["GO_MF"]], res_ra[["GO_MF"]], res_shingles[["GO_MF"]], category="GO MF")
plot_list[[16]] <- plot_venn_3way(res_aging[["KEGG"]],  res_ra[["KEGG"]],  res_shingles[["KEGG"]],  category="KEGG")

# Create label plots
create_label_plot <- function(label, angle=0, size=20, fontface="bold") {
  ggplot() + 
    annotate("text", x=0.5, y=0.5, label=label, angle=angle, size=size, fontface=fontface, family="sans") + 
    theme_void()
}

# Row labels (Left side, rotated 90 degrees)
r1 <- create_label_plot("mvAge", angle=90, size=7)
r2 <- create_label_plot("RA", angle=90, size=7)
r3 <- create_label_plot("HZ (UKB)", angle=90, size=7)
r4 <- create_label_plot("Comparison", angle=90, size=7)

# Column labels (Bottom)
c1 <- create_label_plot("GO BP", size=6)
c2 <- create_label_plot("GO CC", size=6)
c3 <- create_label_plot("GO MF", size=6)
c4 <- create_label_plot("KEGG", size=6)

# Combine plots with layout
# Row 1: Label + 4 Plots
row1 <- r1 + plot_list[[1]] + plot_list[[2]] + plot_list[[3]] + plot_list[[4]] + plot_layout(ncol=5, widths=c(0.2, 1, 1, 1, 1))
# Row 2: Label + 4 Plots
row2 <- r2 + plot_list[[5]] + plot_list[[6]] + plot_list[[7]] + plot_list[[8]] + plot_layout(ncol=5, widths=c(0.2, 1, 1, 1, 1))
# Row 3: Label + 4 Plots
row3 <- r3 + plot_list[[9]] + plot_list[[10]] + plot_list[[11]] + plot_list[[12]] + plot_layout(ncol=5, widths=c(0.2, 1, 1, 1, 1))
# Row 4: Label + 4 Plots
row4 <- r4 + plot_list[[13]] + plot_list[[14]] + plot_list[[15]] + plot_list[[16]] + plot_layout(ncol=5, widths=c(0.2, 1, 1, 1, 1))
# Column labels: Spacer + 4 Labels
col_row <- plot_spacer() + c1 + c2 + c3 + c4 + plot_layout(ncol=5, widths=c(0.2, 1, 1, 1, 1))

# Combine all
combined_plot <- row1 / row2 / row3 / row4 / col_row + plot_layout(heights=c(1,1,1,1,0.15))

# Save plot
ggsave(file.path(output_dir, "Combined_Enrichment_Plot_Optimized.pdf"), combined_plot, width = 22, height = 22) # Increased height
ggsave(file.path(output_dir, "Combined_Enrichment_Plot_Optimized.png"), combined_plot, width = 22, height = 22, dpi = 300)

cat("Plots saved to output directory.\n")

# -------------------------------------------------------------------------
# Generate Swapped Layout (Columns: Traits, Rows: Categories)
# -------------------------------------------------------------------------
cat("Generating Swapped Layout Plots...\n")

# Create labels for swapped layout
# Top labels (Traits) - Horizontal
t1 <- create_label_plot("mvAge", size=7)
t2 <- create_label_plot("RA", size=7)
t3 <- create_label_plot("HZ (UKB)", size=7)
t4 <- create_label_plot("Comparison", size=7)

# Left labels (Categories) - Vertical (Rotated 90)
l1 <- create_label_plot("GO BP", angle=90, size=7)
l2 <- create_label_plot("GO CC", angle=90, size=7)
l3 <- create_label_plot("GO MF", angle=90, size=7)
l4 <- create_label_plot("KEGG", angle=90, size=7)

# Combine plots for swapped layout
# Top Header: Spacer + mvAge + RA + HZ + Comparison
top_header <- plot_spacer() + t1 + t2 + t3 + t4 + plot_layout(ncol=5, widths=c(0.2, 1, 1, 1, 1))

# Row 1 (BP): Label + mvAge + RA + HZ + Comp
swapped_row1 <- l1 + plot_list[[1]] + plot_list[[5]] + plot_list[[9]] + plot_list[[13]] + plot_layout(ncol=5, widths=c(0.2, 1, 1, 1, 1))

# Row 2 (CC): Label + mvAge + RA + HZ + Comp
swapped_row2 <- l2 + plot_list[[2]] + plot_list[[6]] + plot_list[[10]] + plot_list[[14]] + plot_layout(ncol=5, widths=c(0.2, 1, 1, 1, 1))

# Row 3 (MF): Label + mvAge + RA + HZ + Comp
swapped_row3 <- l3 + plot_list[[3]] + plot_list[[7]] + plot_list[[11]] + plot_list[[15]] + plot_layout(ncol=5, widths=c(0.2, 1, 1, 1, 1))

# Row 4 (KEGG): Label + mvAge + RA + HZ + Comp
swapped_row4 <- l4 + plot_list[[4]] + plot_list[[8]] + plot_list[[12]] + plot_list[[16]] + plot_layout(ncol=5, widths=c(0.2, 1, 1, 1, 1))

# Combine all
combined_plot_swapped <- top_header / swapped_row1 / swapped_row2 / swapped_row3 / swapped_row4 + plot_layout(heights=c(0.15, 1, 1, 1, 1))

# Save swapped plot
ggsave(file.path(output_dir, "Combined_Enrichment_Plot_Swapped.pdf"), combined_plot_swapped, width = 22, height = 22)
ggsave(file.path(output_dir, "Combined_Enrichment_Plot_Swapped.png"), combined_plot_swapped, width = 22, height = 22, dpi = 300)

cat("Swapped plots saved to output directory.\n")

# -------------------------------------------------------------------------
# Step 4b: Network Visualization (aPEAR) / 
# -------------------------------------------------------------------------
cat("\n[Step 4b] Generating Network Plots with aPEAR...\n")

plot_enrich_network <- function(enrich_res, title, filename, show_labels = TRUE, show_title = TRUE, size_limits = NULL, size_breaks = NULL, legend_style = "standard") {
  if (is.null(enrich_res) || nrow(enrich_res) == 0) {
    cat(paste0("Skipping network plot for ", title, ": No results.\n"))
    return(NULL)
  }
  
  # aPEAR requires at least a few terms to cluster
  if (nrow(enrich_res) < 2) {
    cat(paste0("Skipping network plot for ", title, ": Too few results (<2).\n"))
    return(NULL)
  }

  tryCatch({
    cat(paste0("Generating network plot for ", title, "...\n"))
    
    # Convert to data frame for aPEAR
    enrich_df <- as.data.frame(enrich_res)
    
    # Safety Check: Limit to top 500 terms to prevent performance issues while maximizing inclusion
    # : 500, 
    if (nrow(enrich_df) > 500) {
      cat(paste0("Note: Result contains ", nrow(enrich_df), " terms. Limiting to top 500 for network visualization performance.\n"))
      enrich_df <- enrich_df[order(enrich_df$p.adjust), ][1:500, ]
    }
    
    # Calculate -log10(p.adjust) for visualization
    #  -log10(p.adjust) 
    enrich_df$LogP <- -log10(enrich_df$p.adjust)
    
    # Determine max limit to ensure consistent scale proportions
    limit_max <- 5.204
    
    v_grey_start <- 0
    v_grey_end   <- 1.29 / limit_max
    v_blue       <- 1.301 / limit_max  # Start of significant at 25%
    v_green      <- 2.0 / limit_max
    v_yellow     <- 3.0 / limit_max
    v_red        <- 4.0 / limit_max
    v_max        <- 1.0
    
    # Determine font size for labels (0 to hide)
    current_font_size <- if(show_labels) 2.25 else 0

    p <- enrichmentNetwork(enrich_df, 
                           simMethod = "jaccard", # Explicitly use Jaccard
                           clustMethod = "hier", # Changed to Hierarchical for stability with large/sparse data
                           colorBy = "LogP", 
                           nodeSize = "Count",
                           repelLabels = show_labels, # Controlled by parameter
                           drawEllipses = TRUE, # Visualize clusters
                           fontSize = current_font_size, # Hide labels if show_labels is FALSE
                           verbose = FALSE) +
          # Custom Color Scheme based on Reference Image
          scale_color_gradientn(
             colours = c("grey80", "grey80", "blue", "green", "yellow", "red", "red"),
             values = c(v_grey_start, v_grey_end, v_blue, v_green, v_yellow, v_red, v_max),
             limits = c(0, limit_max),
             oob = scales::squish, # Ensure values > limit_max are plotted as Red instead of Grey/NA
             breaks = c(1.301, 2, 3, 4),
             labels = c("p<=0.05", "p<0.01", "p<0.001", "p<0.0001"),
             name = "Log(p.adjust)",
             guide = guide_colorbar(
               title.position = "top",
               title.hjust = 0.5,
               barheight = 15,  # Increased height to prevent label overlap
               barwidth = 1.5,
               frame.colour = "black", 
               ticks.colour = "black",
               reverse = TRUE   # Top: Small value (Grey, p>0.05), Bottom: Large value (Red, p<0.0001)
             )
          ) +
          # Increase Pathway Size (Node Size) for better visibility
          scale_size_continuous(
             range = c(1, 4), # Reduced range (half of previous 2-8) as requested
             limits = size_limits,
             breaks = size_breaks,
             name = "Pathway size"
          ) +
          # Optimize Legends:
          # 1. Pathway size: Black circles (override color)
          # 2. Log(p.adjust): Color bar (handled by scale_color_gradientn)
          guides(
            size = guide_legend(
              title = "Pathway size",
              override.aes = list(color = "black", fill = "black"),
              order = 1
            ),
            color = guide_colorbar(order = 2)
          ) +
          labs(title = if(show_title) title else NULL, 
               subtitle = if(show_title) "Network Scheme: Jaccard + Hierarchical + Custom Gradient (Grey>0.05, Blue<=0.05 -> Red)" else NULL,
               size = "Pathway size") +
          theme(plot.title = element_text(hjust = 0.5, face="bold"),
                plot.subtitle = element_text(hjust = 0.5, size=8, color="grey30"))
    
    # Apply Legend Style
    if (legend_style == "inside") {
      # Add padding to the right to create space for legends "outside" the node area
      p <- p + scale_x_continuous(expand = expansion(mult = c(0.05, 0.35)))
      
      # Ensure main plot has transparent background so it can be placed on top of legend
      p <- p + theme(
        panel.background = element_rect(fill = "transparent", color = NA),
        plot.background = element_rect(fill = "transparent", color = NA)
      )
      
      # For "inside" style, we split the legends:
      p_temp <- p + 
        theme(
          legend.position = "right",
          legend.background = element_rect(fill = alpha("white", 0.0), color = NA), # Fully transparent
          legend.text = element_text(size = 8),   # Increased size for visibility
          legend.title = element_text(size = 9, face = "bold"),
          legend.key.size = unit(0.5, "cm"),      # Larger keys
          legend.margin = margin(5, 5, 5, 5)
        ) +
        guides(
          color = guide_colorbar(
            title.position = "top",
            title.hjust = 0.5,
            barheight = 6,
            barwidth = 0.8,
            frame.colour = "black", 
            ticks.colour = "black",
            reverse = TRUE,
            order = 2
          ),
          size = guide_legend(
            title = "Pathway size",
            override.aes = list(color = "black", fill = "black"),
            order = 1,
            direction = "vertical"
          )
        )
      
      p_size <- p_temp + guides(color = "none")
      leg_size <- tryCatch(ggpubr::get_legend(p_size), error=function(e) NULL)
      
      p_color <- p_temp + guides(size = "none")
      leg_color <- tryCatch(ggpubr::get_legend(p_color), error=function(e) NULL)
      
      p <- p + theme(legend.position = "none")
      
      if (!is.null(leg_size)) {
        p <- p + inset_element(leg_size, left=0.7, bottom=0.6, right=1.0, top=1.0, align_to = "panel", on_top = FALSE)
      }
      
      if (!is.null(leg_color)) {
        p <- p + inset_element(leg_color, left=0.7, bottom=0.0, right=1.0, top=0.4, align_to = "panel", on_top = FALSE)
      }
      
    } else if (legend_style == "none") {
      p <- p + theme(legend.position = "none")
    }
    
    # Save as PDF and PNG if filename is provided
    if (!is.null(filename)) {
      ggsave(file.path(output_dir, filename), p, width = 12, height = 10, bg = "white")
      ggsave(file.path(output_dir, gsub(".pdf", ".png", filename)), p, width = 12, height = 10, dpi=300, bg = "white")
    }
    
    return(p)
  }, error = function(e) {
    cat(paste0("Error generating network plot for ", title, ": ", e$message, "\n"))
    return(NULL)
  })
}

if ("package:aPEAR" %in% search()) {
  # -------------------------------------------------------------------------
  # Create Composite Plots (Bubble + Network)
  # -------------------------------------------------------------------------
  cat("\nGenerating Network Plots and Composite (Bubble + Network) Plots...\n")
  
  # Function to create empty placeholder if plot is NULL
  get_plot_or_empty <- function(p) {
    if (is.null(p)) return(plot_spacer())
    return(p)
  }
  
  # Initialize lists to capture network plots
  net_aging <- list()
  net_ra <- list()
  net_shingles <- list()
  
  # Generate network plots for mvAge (Save and Assign)
  net_aging[["GO_BP"]] <- plot_enrich_network(res_aging[["GO_BP"]], "mvAge - GO BP Network", "Network_mvAge_GO_BP.pdf", size_limits = count_limits, size_breaks = count_breaks)
  net_aging[["GO_CC"]] <- plot_enrich_network(res_aging[["GO_CC"]], "mvAge - GO CC Network", "Network_mvAge_GO_CC.pdf", size_limits = count_limits, size_breaks = count_breaks)
  net_aging[["GO_MF"]] <- plot_enrich_network(res_aging[["GO_MF"]], "mvAge - GO MF Network", "Network_mvAge_GO_MF.pdf", size_limits = count_limits, size_breaks = count_breaks)
  net_aging[["KEGG"]]  <- plot_enrich_network(res_aging[["KEGG"]],  "mvAge - KEGG Network",  "Network_mvAge_KEGG.pdf", size_limits = count_limits, size_breaks = count_breaks)
  
  # Generate network plots for RA (Save and Assign)
  net_ra[["GO_BP"]] <- plot_enrich_network(res_ra[["GO_BP"]], "RA - GO BP Network", "Network_RA_GO_BP.pdf", size_limits = count_limits, size_breaks = count_breaks)
  net_ra[["GO_CC"]] <- plot_enrich_network(res_ra[["GO_CC"]], "RA - GO CC Network", "Network_RA_GO_CC.pdf", size_limits = count_limits, size_breaks = count_breaks)
  net_ra[["GO_MF"]] <- plot_enrich_network(res_ra[["GO_MF"]], "RA - GO MF Network", "Network_RA_GO_MF.pdf", size_limits = count_limits, size_breaks = count_breaks)
  net_ra[["KEGG"]]  <- plot_enrich_network(res_ra[["KEGG"]],  "RA - KEGG Network",  "Network_RA_KEGG.pdf", size_limits = count_limits, size_breaks = count_breaks)
  
  # Generate network plots for HZ (UKB) (Save and Assign)
  net_shingles[["GO_BP"]] <- plot_enrich_network(res_shingles[["GO_BP"]], "HZ (UKB) - GO BP Network", "Network_HZ_UKB_GO_BP.pdf", size_limits = count_limits, size_breaks = count_breaks)
  net_shingles[["GO_CC"]] <- plot_enrich_network(res_shingles[["GO_CC"]], "HZ (UKB) - GO CC Network", "Network_HZ_UKB_GO_CC.pdf", size_limits = count_limits, size_breaks = count_breaks)
  net_shingles[["GO_MF"]] <- plot_enrich_network(res_shingles[["GO_MF"]], "HZ (UKB) - GO MF Network", "Network_HZ_UKB_GO_MF.pdf", size_limits = count_limits, size_breaks = count_breaks)
  net_shingles[["KEGG"]]  <- plot_enrich_network(res_shingles[["KEGG"]],  "HZ (UKB) - KEGG Network",  "Network_HZ_UKB_KEGG.pdf", size_limits = count_limits, size_breaks = count_breaks)
  
  cat("Individual Network plots saved.\n")

  # -------------------------------------------------------------------------
  # Re-generate network plots for Composite (No labels, No title)
  # -------------------------------------------------------------------------
  cat("Generating optimized network plots for composite figure (No labels, No title)...\n")
  
  net_aging_comp <- list()
  net_aging_comp[["GO_BP"]] <- plot_enrich_network(res_aging[["GO_BP"]], "mvAge - GO BP Network", NULL, show_labels=FALSE, show_title=FALSE, size_limits = count_limits, size_breaks = count_breaks, legend_style = "inside")
  net_aging_comp[["GO_CC"]] <- plot_enrich_network(res_aging[["GO_CC"]], "mvAge - GO CC Network", NULL, show_labels=FALSE, show_title=FALSE, size_limits = count_limits, size_breaks = count_breaks, legend_style = "inside")
  net_aging_comp[["GO_MF"]] <- plot_enrich_network(res_aging[["GO_MF"]], "mvAge - GO MF Network", NULL, show_labels=FALSE, show_title=FALSE, size_limits = count_limits, size_breaks = count_breaks, legend_style = "inside")
  net_aging_comp[["KEGG"]]  <- plot_enrich_network(res_aging[["KEGG"]],  "mvAge - KEGG Network",  NULL, show_labels=FALSE, show_title=FALSE, size_limits = count_limits, size_breaks = count_breaks, legend_style = "inside")
  
  net_ra_comp <- list()
  net_ra_comp[["GO_BP"]] <- plot_enrich_network(res_ra[["GO_BP"]], "RA - GO BP Network", NULL, show_labels=FALSE, show_title=FALSE, size_limits = count_limits, size_breaks = count_breaks, legend_style = "inside")
  net_ra_comp[["GO_CC"]] <- plot_enrich_network(res_ra[["GO_CC"]], "RA - GO CC Network", NULL, show_labels=FALSE, show_title=FALSE, size_limits = count_limits, size_breaks = count_breaks, legend_style = "inside")
  net_ra_comp[["GO_MF"]] <- plot_enrich_network(res_ra[["GO_MF"]], "RA - GO MF Network", NULL, show_labels=FALSE, show_title=FALSE, size_limits = count_limits, size_breaks = count_breaks, legend_style = "inside")
  net_ra_comp[["KEGG"]]  <- plot_enrich_network(res_ra[["KEGG"]],  "RA - KEGG Network",  NULL, show_labels=FALSE, show_title=FALSE, size_limits = count_limits, size_breaks = count_breaks, legend_style = "inside")

  net_shingles_comp <- list()
  net_shingles_comp[["GO_BP"]] <- plot_enrich_network(res_shingles[["GO_BP"]], "HZ (UKB) - GO BP Network", NULL, show_labels=FALSE, show_title=FALSE, size_limits = count_limits, size_breaks = count_breaks, legend_style = "inside")
  net_shingles_comp[["GO_CC"]] <- plot_enrich_network(res_shingles[["GO_CC"]], "HZ (UKB) - GO CC Network", NULL, show_labels=FALSE, show_title=FALSE, size_limits = count_limits, size_breaks = count_breaks, legend_style = "inside")
  net_shingles_comp[["GO_MF"]] <- plot_enrich_network(res_shingles[["GO_MF"]], "HZ (UKB) - GO MF Network", NULL, show_labels=FALSE, show_title=FALSE, size_limits = count_limits, size_breaks = count_breaks, legend_style = "inside")
  net_shingles_comp[["KEGG"]]  <- plot_enrich_network(res_shingles[["KEGG"]],  "HZ (UKB) - KEGG Network",  NULL, show_labels=FALSE, show_title=FALSE, size_limits = count_limits, size_breaks = count_breaks, legend_style = "inside")

  # Define composite cells (Left: Bubble, Right: Network)
  # mvAge Row
  comp_aging <- list()
  comp_aging[[1]] <- get_plot_or_empty(plot_list[[1]]) + theme(plot.margin=margin(r=0)) + get_plot_or_empty(net_aging_comp[["GO_BP"]]) + theme(plot.margin=margin(l=0)) + plot_layout(ncol=2, widths=c(0.6, 2.4))
  comp_aging[[2]] <- get_plot_or_empty(plot_list[[2]]) + theme(plot.margin=margin(r=0)) + get_plot_or_empty(net_aging_comp[["GO_CC"]]) + theme(plot.margin=margin(l=0)) + plot_layout(ncol=2, widths=c(0.6, 2.4))
  comp_aging[[3]] <- get_plot_or_empty(plot_list[[3]]) + theme(plot.margin=margin(r=0)) + get_plot_or_empty(net_aging_comp[["GO_MF"]]) + theme(plot.margin=margin(l=0)) + plot_layout(ncol=2, widths=c(0.6, 2.4))
  comp_aging[[4]] <- get_plot_or_empty(plot_list[[4]]) + theme(plot.margin=margin(r=0)) + get_plot_or_empty(net_aging_comp[["KEGG"]])  + theme(plot.margin=margin(l=0)) + plot_layout(ncol=2, widths=c(0.6, 2.4))

  # RA Row
  comp_ra <- list()
  comp_ra[[1]] <- get_plot_or_empty(plot_list[[5]]) + theme(plot.margin=margin(r=0)) + get_plot_or_empty(net_ra_comp[["GO_BP"]]) + theme(plot.margin=margin(l=0)) + plot_layout(ncol=2, widths=c(0.6, 2.4))
  comp_ra[[2]] <- get_plot_or_empty(plot_list[[6]]) + theme(plot.margin=margin(r=0)) + get_plot_or_empty(net_ra_comp[["GO_CC"]]) + theme(plot.margin=margin(l=0)) + plot_layout(ncol=2, widths=c(0.6, 2.4))
  comp_ra[[3]] <- get_plot_or_empty(plot_list[[7]]) + theme(plot.margin=margin(r=0)) + get_plot_or_empty(net_ra_comp[["GO_MF"]]) + theme(plot.margin=margin(l=0)) + plot_layout(ncol=2, widths=c(0.6, 2.4))
  comp_ra[[4]] <- get_plot_or_empty(plot_list[[8]]) + theme(plot.margin=margin(r=0)) + get_plot_or_empty(net_ra_comp[["KEGG"]])  + theme(plot.margin=margin(l=0)) + plot_layout(ncol=2, widths=c(0.6, 2.4))

  # Shingles Row
  comp_shingles <- list()
  comp_shingles[[1]] <- get_plot_or_empty(plot_list[[9]]) + theme(plot.margin=margin(r=0)) + get_plot_or_empty(net_shingles_comp[["GO_BP"]]) + theme(plot.margin=margin(l=0)) + plot_layout(ncol=2, widths=c(0.6, 2.4))
  comp_shingles[[2]] <- get_plot_or_empty(plot_list[[10]]) + theme(plot.margin=margin(r=0)) + get_plot_or_empty(net_shingles_comp[["GO_CC"]]) + theme(plot.margin=margin(l=0)) + plot_layout(ncol=2, widths=c(0.6, 2.4))
  comp_shingles[[3]] <- get_plot_or_empty(plot_list[[11]]) + theme(plot.margin=margin(r=0)) + get_plot_or_empty(net_shingles_comp[["GO_MF"]]) + theme(plot.margin=margin(l=0)) + plot_layout(ncol=2, widths=c(0.6, 2.4))
  comp_shingles[[4]] <- get_plot_or_empty(plot_list[[12]]) + theme(plot.margin=margin(r=0)) + get_plot_or_empty(net_shingles_comp[["KEGG"]])  + theme(plot.margin=margin(l=0)) + plot_layout(ncol=2, widths=c(0.6, 2.4))

  # Assemble Composite Plot
  # Row 1: Label + 4 Composite Plots (mvAge)
  comp_row1 <- r1 + comp_aging[[1]] + comp_aging[[2]] + comp_aging[[3]] + comp_aging[[4]] + plot_layout(ncol=5, widths=c(0.2, 1, 1, 1, 1))
  
  # Row 2: Label + 4 Composite Plots (RA)
  comp_row2 <- r2 + comp_ra[[1]] + comp_ra[[2]] + comp_ra[[3]] + comp_ra[[4]] + plot_layout(ncol=5, widths=c(0.2, 1, 1, 1, 1))
  
  # Row 3: Label + 4 Composite Plots (Shingles)
  comp_row3 <- r3 + comp_shingles[[1]] + comp_shingles[[2]] + comp_shingles[[3]] + comp_shingles[[4]] + plot_layout(ncol=5, widths=c(0.2, 1, 1, 1, 1))
  
  # Row 4: Label + 4 Venn Plots (Comparison)
  comp_row4 <- r4 + plot_list[[13]] + plot_list[[14]] + plot_list[[15]] + plot_list[[16]] + plot_layout(ncol=5, widths=c(0.2, 1, 1, 1, 1))
  
  # Column Labels
  comp_col_row <- plot_spacer() + c1 + c2 + c3 + c4 + plot_layout(ncol=5, widths=c(0.2, 1, 1, 1, 1))
  
  # Combine All
  final_composite <- comp_row1 / comp_row2 / comp_row3 / comp_row4 / comp_col_row + plot_layout(heights=c(1, 1, 1, 0.6, 0.1))
  
  # Save Composite Plot
  ggsave(file.path(output_dir, "Combined_Enrichment_Network_Composite.pdf"), final_composite, width = 40, height = 28, limitsize = FALSE)
  ggsave(file.path(output_dir, "Combined_Enrichment_Network_Composite.png"), final_composite, width = 40, height = 28, dpi = 300, limitsize = FALSE)
  
  cat("Composite plots saved.\n")

  # -------------------------------------------------------------------------
  # Generate Swapped Composite Layout (Columns: Traits, Rows: Categories)
  # -------------------------------------------------------------------------
  cat("Generating Swapped Composite Layout Plots...\n")
  
  # Define swapped rows for composite plots
   # Structure: Label (Vertical) | mvAge (Bubble+Net) | RA | HZ (Bubble+Net) | Spacer(0.2) | Comparison (Venn) | Spacer(0.4)
   # Widths: 0.2 : 3.5 : 3.5 : 3.5 : 0.2 : 1 : 0.4
   
   # Use inset_element on plot_spacer() to place Venn diagram on the top layer and enlarge it
   # while keeping the original layout widths exactly the same.
   wrap_venn <- function(p) {
     # (0.5)
     # right = 1.3  spacer(0.4) , 
     plot_spacer() + inset_element(p, left = -0.3, bottom = -0.2, right = 1.3, top = 1.2, align_to = "full", on_top = TRUE)
   }
   
   # Row 1 (BP)
   swapped_comp_row1 <- l1 + comp_aging[[1]] + comp_ra[[1]] + comp_shingles[[1]] + plot_spacer() + wrap_venn(plot_list[[13]]) + plot_spacer() + plot_layout(ncol=7, widths=c(0.2, 3.5, 3.5, 3.5, 0.2, 1, 0.4))
   
   # Row 2 (CC)
   swapped_comp_row2 <- l2 + comp_aging[[2]] + comp_ra[[2]] + comp_shingles[[2]] + plot_spacer() + wrap_venn(plot_list[[14]]) + plot_spacer() + plot_layout(ncol=7, widths=c(0.2, 3.5, 3.5, 3.5, 0.2, 1, 0.4))
   
   # Row 3 (MF)
   swapped_comp_row3 <- l3 + comp_aging[[3]] + comp_ra[[3]] + comp_shingles[[3]] + plot_spacer() + wrap_venn(plot_list[[15]]) + plot_spacer() + plot_layout(ncol=7, widths=c(0.2, 3.5, 3.5, 3.5, 0.2, 1, 0.4))
   
   # Row 4 (KEGG)
   swapped_comp_row4 <- l4 + comp_aging[[4]] + comp_ra[[4]] + comp_shingles[[4]] + plot_spacer() + wrap_venn(plot_list[[16]]) + plot_spacer() + plot_layout(ncol=7, widths=c(0.2, 3.5, 3.5, 3.5, 0.2, 1, 0.4))
   
   # Top Header: Spacer + mvAge + RA + HZ + Spacer(0.2) + Comparison + Spacer(0.4)
   swapped_comp_header <- plot_spacer() + t1 + t2 + t3 + plot_spacer() + t4 + plot_spacer() + plot_layout(ncol=7, widths=c(0.2, 3.5, 3.5, 3.5, 0.2, 1, 0.4))
  
  # Combine All
  final_composite_swapped <- swapped_comp_header / swapped_comp_row1 / swapped_comp_row2 / swapped_comp_row3 / swapped_comp_row4 + plot_layout(heights=c(0.15, 1, 1, 1, 1))
  
  # Save Swapped Composite Plot
  ggsave(file.path(output_dir, "Combined_Enrichment_Network_Composite_Swapped.pdf"), final_composite_swapped, width = 38, height = 32, limitsize = FALSE)
  ggsave(file.path(output_dir, "Combined_Enrichment_Network_Composite_Swapped.png"), final_composite_swapped, width = 38, height = 32, dpi = 300, limitsize = FALSE)
  
  cat("Swapped Composite plots saved.\n")
  
} else {
  cat("aPEAR package not loaded. Skipping network plots.\n")
}

# -------------------------------------------------------------------------
# 5. Save Results / 
# -------------------------------------------------------------------------
cat("\n[Step 5] Saving result tables...\n")

save_results <- function(res_list, name) {
  wb <- createWorkbook()
  
  for (category in names(res_list)) {
    if (!is.null(res_list[[category]]) && nrow(res_list[[category]]) > 0) {
      addWorksheet(wb, category)
      writeData(wb, category, as.data.frame(res_list[[category]]))
    }
  }
  
  saveWorkbook(wb, file.path(output_dir, paste0("Enrichment_Results_", name, ".xlsx")), overwrite = TRUE)
}

save_results(res_aging, "mvAge")
save_results(res_ra, "RA")
save_results(res_shingles, "HZ_UKB")

# Save merged summary (User requirement: "+")
# Create a comprehensive summary file
cat("Creating comprehensive summary file...\n")
wb_summary <- createWorkbook()

# Add Gene Lists
addWorksheet(wb_summary, "Gene_Lists")
max_len <- max(length(genes_aging), length(genes_ra), length(genes_shingles))
genes_df <- data.frame(
  mvAge_Genes = c(genes_aging, rep(NA, max_len - length(genes_aging))),
  RA_Genes = c(genes_ra, rep(NA, max_len - length(genes_ra))),
  HZ_UKB_Genes = c(genes_shingles, rep(NA, max_len - length(genes_shingles)))
)
writeData(wb_summary, "Gene_Lists", genes_df)

# Add Top Enrichment Results
add_top_results <- function(wb, res_list, prefix) {
  for (category in names(res_list)) {
    if (!is.null(res_list[[category]]) && nrow(res_list[[category]]) > 0) {
      sheet_name <- substr(paste0(prefix, "_", category), 1, 31) # Sheet name limit
      addWorksheet(wb, sheet_name)
      writeData(wb, sheet_name, as.data.frame(res_list[[category]]))
    }
  }
}

add_top_results(wb_summary, res_aging, "mvAge")
add_top_results(wb_summary, res_ra, "RA")
add_top_results(wb_summary, res_shingles, "HZ_UKB")

saveWorkbook(wb_summary, file.path(output_dir, "Summary_GeneLists_and_Enrichment.xlsx"), overwrite = TRUE)

# -------------------------------------------------------------------------
# NEW: Create a supplementary Excel file summarizing all significant results
# -------------------------------------------------------------------------
cat("\nCreating Supplementary Table for significant enrichment results...\n")

# Helper function to extract and format significant results
extract_sig_results <- function(res_list, pheno_name, category) {
  res_obj <- res_list[[category]]
  if (is.null(res_obj)) return(NULL)
  
  df <- as.data.frame(res_obj)
  if (nrow(df) == 0) return(NULL)
  
  df_sig <- df[df$p.adjust < 0.05, ]
  if (nrow(df_sig) == 0) return(NULL)
  
  # Add Phenotype and Category columns
  df_sig$Phenotype <- pheno_name
  df_sig$Category <- category
  return(df_sig)
}

# 1. Collect all GO results
all_go_list <- list()
for (pheno in c("mvAge", "RA", "HZ_UKB")) {
  if (pheno == "mvAge") res_list <- res_aging
  else if (pheno == "RA") res_list <- res_ra
  else if (pheno == "HZ_UKB") res_list <- res_shingles
  
  for (cat in c("GO_BP", "GO_CC", "GO_MF")) {
    df_cat <- extract_sig_results(res_list, pheno, cat)
    if (!is.null(df_cat)) all_go_list[[paste(pheno, cat, sep="_")]] <- df_cat
  }
}
all_go <- do.call(rbind, all_go_list)

# 2. Collect all KEGG results
all_kegg_list <- list()
for (pheno in c("mvAge", "RA", "HZ_UKB")) {
  if (pheno == "mvAge") res_list <- res_aging
  else if (pheno == "RA") res_list <- res_ra
  else if (pheno == "HZ_UKB") res_list <- res_shingles
  
  df_cat <- extract_sig_results(res_list, pheno, "KEGG")
  if (!is.null(df_cat)) all_kegg_list[[paste(pheno, "KEGG", sep="_")]] <- df_cat
}
all_kegg <- do.call(rbind, all_kegg_list)

# 3. Add Overlapping_Phenotypes column
add_overlap_col <- function(df) {
  if (is.null(df) || nrow(df) == 0) return(df)
  
  # Map Description to all phenotypes that have this description
  overlap_map <- aggregate(Phenotype ~ Description, data = df, FUN = function(x) paste(unique(x), collapse = ", "))
  names(overlap_map)[2] <- "Overlapping_Phenotypes"
  
  # Merge back
  df_merged <- merge(df, overlap_map, by = "Description", all.x = TRUE)
  
  # Ensure Phenotype factor ordering: mvAge -> RA -> HZ_UKB
  df_merged$Phenotype <- factor(df_merged$Phenotype, levels = c("mvAge", "RA", "HZ_UKB"))
  
  # Order by Phenotype, then Category, then p.adjust
  df_sorted <- df_merged[order(df_merged$Phenotype, df_merged$Category, df_merged$p.adjust), ]
  
  # Reorder columns to make it readable
  cols_first <- c("Phenotype", "Category", "ID", "Description", "Overlapping_Phenotypes")
  cols_rest <- setdiff(names(df_sorted), cols_first)
  df_sorted <- df_sorted[, c(cols_first, cols_rest)]
  
  return(df_sorted)
}

all_go_final <- add_overlap_col(all_go)
all_kegg_final <- add_overlap_col(all_kegg)

# 4. Save to Excel
wb_supp <- createWorkbook()
addWorksheet(wb_supp, "Input_Gene_Lists")
writeData(wb_supp, "Input_Gene_Lists", genes_df)

if (!is.null(all_go_final) && nrow(all_go_final) > 0) {
  addWorksheet(wb_supp, "GO_Enrichment")
  writeData(wb_supp, "GO_Enrichment", all_go_final)
}

if (!is.null(all_kegg_final) && nrow(all_kegg_final) > 0) {
  addWorksheet(wb_supp, "KEGG_Enrichment")
  writeData(wb_supp, "KEGG_Enrichment", all_kegg_final)
}

supp_file_path <- file.path(output_dir, "Supplementary_Table_Enrichment_Results.xlsx")
saveWorkbook(wb_supp, supp_file_path, overwrite = TRUE)
cat(paste("Saved Supplementary Table to:", supp_file_path, "\n"))

# Create Readme
readme_content <- c(
  "Project: Multi-omics and HZ (UKB)",
  "Analysis: Cross-trait Enrichment Analysis Comparison (mvAge, Rheumatoid Arthritis, HZ (UKB))",
  paste0("Date: ", Sys.Date()),
  "Author: AI Assistant",
  "",
  "Description:",
  "This folder contains enrichment analysis results for mvAge, Rheumatoid Arthritis (RA), and HZ (UKB) gene lists.",
  "1. Enrichment_Results_mvAge.xlsx: Detailed GO/KEGG results for mvAge.",
  "2. Enrichment_Results_RA.xlsx: Detailed GO/KEGG results for RA.",
  "3. Enrichment_Results_HZ_UKB.xlsx: Detailed GO/KEGG results for HZ (UKB).",
  "4. Summary_GeneLists_and_Enrichment.xlsx: Combined file with gene lists and results.",
  "5. Combined_Enrichment_Plot_Optimized.pdf/png: 4x4 grid visualization (Bubble Plots Only - Standard Layout).",
  "6. Combined_Enrichment_Plot_Swapped.pdf/png: 4x4 grid visualization (Bubble Plots Only - Swapped Layout).",
  "7. Combined_Enrichment_Network_Composite.pdf/png: Standard Layout (Bubble + Network Plots).",
  "8. Combined_Enrichment_Network_Composite_Swapped.pdf/png: Swapped Layout (Bubble + Network Plots).",
  "9. Network_*.pdf/png: Individual Enrichment network plots using aPEAR.",
  "",
  "Visualization Layout (Standard):",
  "Row 1: mvAge (Labels on Left)",
  "Row 2: Rheumatoid Arthritis (Labels on Left)",
  "Row 3: HZ (UKB) (Labels on Left)",
  "Row 4: Comparison (Labels on Left)",
  "Columns: GO BP, GO CC, GO MF, KEGG (Labels at Bottom)",
  "",
  "Visualization Layout (Swapped):",
  "Columns: mvAge, Rheumatoid Arthritis, HZ (UKB), Comparison (Labels at Top)",
  "Rows: GO BP, GO CC, GO MF, KEGG (Labels on Left)",
  "",
  "Visualization Interpretation (Comparison - Venn Diagrams):",
  "The Venn diagrams in the last row illustrate the intersection of significant enrichment results between mvAge, RA, and HZ (UKB).",
  "- Overlapping areas indicate shared biological mechanisms.",
  "",
  "Methods:",
  "Enrichment analysis performed using clusterProfiler (R).",
  "Significance cutoff: p.adjust < 0.05 (BH correction).",
  "Comparison metric: Intersection of significant terms."
)

writeLines(readme_content, file.path(output_dir, "README.txt"))

# -------------------------------------------------------------------------
# 6. Finalize / 
# -------------------------------------------------------------------------
end_time <- Sys.time()
cat("\nAnalysis completed successfully!\n")
cat(paste0("Start Time: ", timestamp, "\n"))
cat(paste0("End Time: ", end_time, "\n"))
print(sessionInfo())
sink()
