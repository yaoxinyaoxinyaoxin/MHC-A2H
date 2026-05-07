#!/usr/bin/env Rscript
# ==============================================================================
# [Script]: 2.2.2.1_Integrated_MR_Plots_Analysis_template.R
# [Method]: MR Analysis Visualization (MR)
# [Step]: (+)
# 
# [Function]:
# Integration script to replicate and combine circular MR plots from different 
#       analysis pipelines into a single figure.
#       （ a/b/c/d/e/f）. 
# 
# [Data Availability / ]:
# Requires 6 source R scripts that generate the circular/chord plots.
# ==============================================================================

rm(list = ls())
gc()

suppressPackageStartupMessages({
  library(optparse)
  library(circlize)
  library(ComplexHeatmap)
  library(grid)
  library(gridBase)
  library(cowplot)
  library(ggplot2)
  library(data.table)
  library(dplyr)
  library(stringr)
  library(RColorBrewer)
  library(colorspace)
})

start_time <- Sys.time()
cat("Script Started at:", as.character(start_time), "\n")

option_list <- list(
  make_option(c("--script_a"), type="character", default=NULL, help="Path to Script A (e.g. OpenGWAS Aging)"),
  make_option(c("--script_b"), type="character", default=NULL, help="Path to Script B (e.g. OpenGWAS RA)"),
  make_option(c("--script_c"), type="character", default=NULL, help="Path to Script C (e.g. OpenGWAS HZ)"),
  make_option(c("--script_d"), type="character", default=NULL, help="Path to Script D (e.g. OneK1K Aging)"),
  make_option(c("--script_e"), type="character", default=NULL, help="Path to Script E (e.g. OneK1K RA)"),
  make_option(c("--script_f"), type="character", default=NULL, help="Path to Script F (e.g. OneK1K HZ)"),
  make_option(c("--out_dir"), type="character", default="./Integrated_Plots", help="Output directory")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

scripts <- c(opt$script_a, opt$script_b, opt$script_c, opt$script_d, opt$script_e, opt$script_f)

if (any(sapply(scripts, is.null))) {
  print_help(opt_parser)
  stop("Missing one or more required script paths (--script_a to --script_f)")
}

for (s in scripts) {
  if (!file.exists(s)) stop("Script not found: ", s)
}

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- file.path(opt$out_dir, paste0("Integrated_Plots_", timestamp))
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

cat("Output Directory:", output_dir, "\n")

labels <- letters[1:6]

# ==============================================================================
# Execution Logic
# ==============================================================================
run_script_logic <- function(script_path, plot_label, new_out_dir) {
  cat("\n----------------------------------------------------------------\n")
  cat("Processing:", plot_label, "from", basename(script_path), "\n")
  cat("----------------------------------------------------------------\n")
  
  lines <- readLines(script_path, warn = FALSE)
  lines <- lines[!grepl("sink\\(", lines)]
  lines <- gsub("^output_dir\\s*<-.*", paste0('output_dir <- "', new_out_dir, '"'), lines)
  lines <- gsub("^output_dir\\s*=.*", paste0('output_dir = "', new_out_dir, '"'), lines)
  lines <- gsub("setwd\\(", "# setwd\\(", lines)
  
  new_pdf <- file.path(new_out_dir, paste0(plot_label, ".pdf"))
  new_png <- file.path(new_out_dir, paste0(plot_label, ".png"))
  
  # Check which function is used
  if (any(grepl("draw_circular_plot_optimized", lines))) {
    target_func <- "draw_circular_plot_optimized"
    idx <- grep(paste0("^", target_func, "\\s*\\("), lines)
    if (length(idx) > 0) lines[idx] <- paste0("# [Auto-Commented] ", lines[idx])
    
    injection <- c(
      "# --- Injected Plotting Calls ---",
      paste0('cat("Generating PDF for ', plot_label, '...\\n")'),
      paste0(target_func, '(to_file = TRUE, file_path = "', new_pdf, '", file_type = "pdf")'),
      paste0('cat("Generating PNG for ', plot_label, '...\\n")'),
      paste0(target_func, '(to_file = TRUE, file_path = "', new_png, '", file_type = "png")'),
      'circos.clear()'
    )
    lines <- c(lines, injection)
  } else if (any(grepl("draw_circular_heatmap_mirror", lines))) {
    target_func <- "draw_circular_heatmap_mirror"
    idx <- grep(paste0("^", target_func, "\\s*\\("), lines)
    if (length(idx) > 0) lines[idx] <- paste0("# [Auto-Commented] ", lines[idx])
    
    injection <- c(
      "# --- Injected Plotting Calls ---",
      paste0('cat("Generating PDF for ', plot_label, '...\\n")'),
      paste0(target_func, '(to_file = TRUE, file_path = "', new_pdf, '", file_type = "pdf", mirror_total_degree = 90)'),
      paste0('cat("Generating PNG for ', plot_label, '...\\n")'),
      paste0(target_func, '(to_file = TRUE, file_path = "', new_png, '", file_type = "png", mirror_total_degree = 90)'),
      'circos.clear()'
    )
    lines <- c(lines, injection)
  } else {
    warning("Neither draw_circular_plot_optimized nor draw_circular_heatmap_mirror found in ", basename(script_path))
  }
  
  env <- new.env()
  assign("output_dir", new_out_dir, envir = env)
  
  tryCatch({
    eval(parse(text = lines), envir = env)
    cat("Successfully processed:", plot_label, "\n")
  }, error = function(e) {
    cat("Error processing", plot_label, ":", e$message, "\n")
  })
  
  return(new_png)
}

png_files <- c()
for (i in 1:6) {
  res_png <- run_script_logic(scripts[i], labels[i], output_dir)
  png_files[labels[i]] <- res_png
}

# ==============================================================================
# Combine Plots
# ==============================================================================
cat("\n----------------------------------------------------------------\n")
cat("Combining plots into single figure...\n")
cat("----------------------------------------------------------------\n")

get_img_plot <- function(img_path, label) {
  if (!file.exists(img_path)) {
    warning("Image not found:", img_path)
    return(NULL)
  }
  ggdraw() + 
    draw_image(img_path, scale = 0.95) +
    draw_label(label, x = 0.02, y = 0.98, hjust = 0, vjust = 1, size = 30, fontface = "bold")
}

pa <- get_img_plot(png_files["a"], "a")
pb <- get_img_plot(png_files["b"], "b")
pc <- get_img_plot(png_files["c"], "c")
pd <- get_img_plot(png_files["d"], "d")
pe <- get_img_plot(png_files["e"], "e")
pf <- get_img_plot(png_files["f"], "f")

# Arrange in 2x3 grid
combined_plot <- plot_grid(pa, pb, pc, pd, pe, pf, ncol = 3, nrow = 2) + 
  theme(plot.background = element_rect(fill = "white", color = NA))

# Arrange in 3x2 grid
combined_plot_3x2 <- plot_grid(pa, pd, pb, pe, pc, pf, ncol = 2, nrow = 3) + 
  theme(plot.background = element_rect(fill = "white", color = NA))

# Save Final Figures
final_pdf <- file.path(output_dir, "Integrated_Figure_abcdef.pdf")
final_png <- file.path(output_dir, "Integrated_Figure_abcdef.png")
ggsave(final_pdf, combined_plot, width = 30, height = 20, limitsize = FALSE, bg = "white")
ggsave(final_png, combined_plot, width = 30, height = 20, dpi = 300, limitsize = FALSE, bg = "white")

final_pdf_3x2 <- file.path(output_dir, "Integrated_Figure_abcdef_3x2.pdf")
final_png_3x2 <- file.path(output_dir, "Integrated_Figure_abcdef_3x2.png")
ggsave(final_pdf_3x2, combined_plot_3x2, width = 20, height = 30, limitsize = FALSE, bg = "white")
ggsave(final_png_3x2, combined_plot_3x2, width = 20, height = 30, dpi = 300, limitsize = FALSE, bg = "white")

# Legend
legend_file <- file.path(output_dir, "Figure_Legend.md")
legend_content <- c(
  "# Figure 1 Legend /  1 ", "", "## English",
  "**Figure 1. Integrated Mendelian Randomization Analysis Results.**", "",
  "- **a-c**: Circular plots based on OpenGWAS bulk eQTL data for Aging (a), Rheumatoid Arthritis (b), and Herpes Zoster (c).",
  "- **d-f**: Chord diagrams based on OneK1K sc-eQTL data for Aging (d), Rheumatoid Arthritis (e), and Herpes Zoster (f).",
  "", "## ", "** 1. . **", "",
  "- **a-c**:  OpenGWAS bulk eQTL ,  (a)、 (b)  (c). ",
  "- **d-f**:  OneK1K sc-eQTL ,  (d)、 (e)  (f)."
)
writeLines(legend_content, legend_file)

end_time <- Sys.time()
cat("\n================================================================\n")
cat("All tasks completed.\n")
cat("Total Runtime:", round(difftime(end_time, start_time, units = "secs"), 2), "seconds\n")
cat("Results saved to:", output_dir, "\n")
cat("================================================================\n")
