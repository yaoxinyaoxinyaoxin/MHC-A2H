# ==============================================================================
# [Script]: 4.1.6.1_plot_3d_scatter_with_PCA_template.R
# [Method]: 3D PCA Scatter Plot
# [Step]: 4.1.6.1_plot_3d_scatter_with_PCA
#
# [Function]:
# Reads intersection signal data and visualizes their relationships in a 3D scatter plot.
#
# [Usage]: 
#   Rscript 4.1.6.1_plot_3d_scatter_with_PCA_template.R \
#     --input_dir <path> \
#     --cluster_file <path> \
#     --out_dir <path>
# ==============================================================================


# ------------------------------------------------------------------------------
# 1. Environment Setup / 
# ------------------------------------------------------------------------------
rm(list = ls) # Clear environment / 
options(stringsAsFactors = FALSE)

# Load required libraries / 
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(scatterplot3d)
  library(plotly)
  library(htmlwidgets)
  library(tidyr)
  library(RColorBrewer)
})

# Define command line arguments
option_list <- list(
  make_option(c("--input_dir"), type="character", default=NULL,
              help="Directory containing input intersection 3D data CSV files / 3DCSV"),
  make_option(c("--cluster_file"), type="character", default=NULL,
              help="Cluster summary CSV file path / CSV"),
  make_option(c("--out_dir"), type="character", default="./3D_Scatter_Output",
              help="Output directory path / ")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_dir) || is.null(opt$cluster_file)) {
  print_help(opt_parser)
  stop("Missing required arguments. / . ", call.=FALSE)
}

# Define paths / 
input_dir <- opt$input_dir
cluster_file <- opt$cluster_file
output_base_dir <- opt$out_dir

# Create output directory with timestamp / 
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- file.path(output_base_dir, timestamp)
if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)

# Create logs directory / 
log_dir <- file.path(output_dir, "logs")
if (!dir.exists(log_dir)) dir.create(log_dir, recursive = TRUE)


# Start logging / 
log_file <- file.path(log_dir, paste0("3D_Scatter_Plot_", timestamp, ".log"))
sink(log_file, split = TRUE)

cat(sprintf("[%s] Script started.\n", Sys.time()))
cat(sprintf("Input Directory: %s\n", input_dir))
cat(sprintf("Output Directory: %s\n", output_dir))

# ------------------------------------------------------------------------------
# 2. Data Loading and Processing / 
# ------------------------------------------------------------------------------
cat(sprintf("[%s] Loading data...\n", Sys.time()))


# Load Cluster Data
if (file.exists(cluster_file)) {
  cat(sprintf("Loading Cluster Summary from: %s\n", cluster_file))
  cluster_df <- read_csv(cluster_file, show_col_types = FALSE)
  
  # Process Cluster Data to create (Cell, SNP) -> ClusterID mapping
  cluster_map <- cluster_df %>%
    select(Exposure_Cluster_ID, Exposure_Cell, Lead_SNPs) %>%
    separate_rows(Lead_SNPs, sep = ";\\s*") %>%
    mutate(SNP = trimws(Lead_SNPs)) %>%
    select(cell = Exposure_Cell, SNP, Cluster_ID = Exposure_Cluster_ID)
    
  cat(sprintf("Loaded %d cluster mappings.\n", nrow(cluster_map)))
} else {
  warning(sprintf("Cluster file not found: %s\nProceeding without cluster coloring.", cluster_file))
  cluster_map <- NULL
}

# List all CSV files / CSV
files <- list.files(input_dir, pattern = "*intersection_3d_data.csv", full.names = TRUE)

# Exclude specific file (bulk_Signal_41) as requested previously (keeping this logic)
files <- files[basename(files) != "bulk_Signal_41_intersection_3d_data.csv"]

# 2. Version_2_No_Bulk Requirement: Exclude bulk data and specific cd4nc file
files <- files[!grepl("bulk", basename(files), ignore.case = TRUE)]
files <- files[basename(files) != "cd4nc_Signal_20_intersection_3d_data.csv"]

if (length(files) == 0) {
  stop("No input files found matching pattern '*intersection_3d_data.csv' in input directory.")
}

# Read and bind data / 
data_list <- lapply(files, function(f) {
  # Explicitly specify column types to prevent parsing errors
  df <- read_csv(f, show_col_types = FALSE, col_types = cols(
    cell = col_character(),
    gene = col_character(),
    SNP = col_character(),
    b_aging = col_double(),
    b_RA = col_double(),
    b_HZ = col_double(),
    z_aging = col_double(),
    z_RA = col_double(),
    z_HZ = col_double(),
    effect_allele = col_character(),
    other_allele = col_character(),
    beta_exposure = col_double(),
    eaf_exposure = col_double()
    # is_association_signal might be logical or character depending on file content, let readr guess or handle later
  ))

  # Check if required columns exist
  req_cols <- c("cell", "gene", "SNP", "b_aging", "b_RA", "b_HZ", "z_aging", "z_RA", "z_HZ", "effect_allele", "other_allele", "beta_exposure", "eaf_exposure")
  if (!all(req_cols %in% names(df))) {
    warning(sprintf("File %s is missing required columns. Skipping.", basename(f)))
    return(NULL)
  }
  
  # Extract Cluster_ID from filename
  fname <- basename(f)
  cluster_id <- sub("_intersection_3d_data\\.csv$", "", fname)
  
  # Assign Cluster_ID to all rows
  df$Cluster_ID <- cluster_id
  
  return(df)
})

# Remove NULLs and bind / 
all_data <- bind_rows(data_list)

if (nrow(all_data) == 0) {
  stop("No valid data loaded.")
}

# Clean cell names: map Bulk_MHC to bulk to match cluster file if needed
if ("Bulk_MHC" %in% all_data$cell) {
  cat("Mapping 'Bulk_MHC' to 'bulk'...\n")
  all_data$cell[all_data$cell == "Bulk_MHC"] <- "bulk"
}

# Filter out Unclustered
all_data <- all_data %>% filter(Cluster_ID != "Unclustered")
  
# --- Apply Jitter to coordinates to resolve exact overlaps ---
set.seed(123) # For reproducibility
jitter_amount <- 1e-6
all_data$b_aging <- all_data$b_aging + runif(nrow(all_data), -jitter_amount, jitter_amount)
all_data$b_RA <- all_data$b_RA + runif(nrow(all_data), -jitter_amount, jitter_amount)
all_data$b_HZ <- all_data$b_HZ + runif(nrow(all_data), -jitter_amount, jitter_amount)
cat("Applied jitter (1e-6) to coordinates to resolve exact overlaps.\n")

# Identify Bulk vs SC
all_data$is_bulk <- grepl("bulk", all_data$cell, ignore.case = TRUE) | grepl("bulk", all_data$Cluster_ID, ignore.case = TRUE)
cat(sprintf("Identified %d bulk rows and %d single-cell rows.\n", 
            sum(all_data$is_bulk), sum(!all_data$is_bulk)))

# --- Construct Hover Text /  ---
cat("Constructing hover text for better tooltip display...\n")

# Normalize is_association_signal
if ("is_association_signal" %in% names(all_data)) {
  all_data$is_association_signal <- all_data$is_association_signal %in% c(TRUE, "TRUE")
} else {
  all_data$is_association_signal <- FALSE
}

# Helper to replace NA
replace_na_char <- function(x) ifelse(is.na(x) | x == "NA", "-", as.character(x))
replace_na_num <- function(x) ifelse(is.na(x), "-", as.character(round(x, 4)))

all_data$hover_text <- paste0(
  "mvAge(b): ", round(all_data$b_aging, 3), ", RA(b): ", round(all_data$b_RA, 3), ", HZ(b): ", round(all_data$b_HZ, 3), "<br>",
  "mvAge(z): ", round(all_data$z_aging, 3), ", RA(z): ", round(all_data$z_RA, 3), ", HZ(z): ", round(all_data$z_HZ, 3), "<br>",
  "Cluster: ", all_data$Cluster_ID, "<br>",
  "Cell: ", all_data$cell, "<br>",
  "Gene: ", all_data$gene, "<br>",
  "SNP: ", all_data$SNP, "<br>",
  "EA: ", replace_na_char(all_data$effect_allele), " OA: ", replace_na_char(all_data$other_allele), "<br>",
  "Beta_Exp: ", replace_na_num(all_data$beta_exposure), " EAF: ", replace_na_num(all_data$eaf_exposure),
  ifelse(all_data$is_association_signal, "<br>(Assoc)", ""),
  ifelse(all_data$is_bulk, "<br>(Bulk)", "<br>(Single-Cell)")
)

cat(sprintf("Total rows loaded: %d\n", nrow(all_data)))
cat(sprintf("Unique cells: %d, Genes: %d, SNPs: %d, Clusters: %d\n", 
            n_distinct(all_data$cell), n_distinct(all_data$gene), 
            n_distinct(all_data$SNP), n_distinct(all_data$Cluster_ID)))


# ------------------------------------------------------------------------------
# 3. Global Configuration (Colors & Helpers) / 
# ------------------------------------------------------------------------------

# Define Colors based on Cluster_ID (Signal Set)
unique_clusters <- unique(all_data$Cluster_ID)
n_clusters <- length(unique_clusters)

# Generate a palette with enough colors
if (n_clusters <= 8) {
  cluster_colors <- setNames(RColorBrewer::brewer.pal(max(3, n_clusters), "Set1")[1:n_clusters], unique_clusters)
} else {
  cluster_colors <- setNames(rainbow(n_clusters), unique_clusters)
}
if ("Unclustered" %in% names(cluster_colors)) {
  cluster_colors["Unclustered"] <- "grey50"
}

# Assign colors to data
all_data$color <- cluster_colors[all_data$Cluster_ID]

# Generate Gene Colors (Global for consistency)
unique_genes <- sort(unique(all_data$gene))
n_genes <- length(unique_genes)
if (n_genes <= 12) {
  gene_colors <- setNames(RColorBrewer::brewer.pal(max(3, n_genes), "Set3")[1:n_genes], unique_genes)
} else {
  base_cols <- c(RColorBrewer::brewer.pal(8, "Set2"), RColorBrewer::brewer.pal(12, "Set3"), RColorBrewer::brewer.pal(9, "Set1"))
  gene_colors <- setNames(colorRampPalette(base_cols)(n_genes), unique_genes)
}

# Helper: Calculate symmetric limits
get_symmetric_limit <- function(vals) {
  max_val <- max(abs(vals), na.rm = TRUE)
  if (max_val == 0) max_val <- 1 
  limit <- max_val * 1.1 
  return(c(-limit, limit))
}

# Helper: Create mesh3d arrow data
create_arrow_mesh <- function(axis_dir, tip_pos, ranges, arrow_len_ratio=0.05, arrow_rad_ratio=0.02, N=12) {
    rx <- ranges[1]; ry <- ranges[2]; rz <- ranges[3]
    len_val <- ranges[which(c("x","y","z") == axis_dir)] * arrow_len_ratio
    rad_y <- ry * arrow_rad_ratio
    rad_z <- rz * arrow_rad_ratio
    rad_x <- rx * arrow_rad_ratio
    
    if (axis_dir == "x") {
      tip <- c(tip_pos, 0, 0)
      base_center <- c(tip_pos - len_val, 0, 0)
      theta <- seq(0, 2*pi, length.out=N+1)
      y_base <- 0 + rad_y * cos(theta)
      z_base <- 0 + rad_z * sin(theta)
      x_base <- rep(base_center[1], N+1)
      x_verts <- c(tip[1], x_base, base_center[1])
      y_verts <- c(tip[2], y_base, base_center[2])
      z_verts <- c(tip[3], z_base, base_center[3])
    } else if (axis_dir == "y") {
      tip <- c(0, tip_pos, 0)
      base_center <- c(0, tip_pos - len_val, 0)
      theta <- seq(0, 2*pi, length.out=N+1)
      x_base <- 0 + rad_x * cos(theta)
      z_base <- 0 + rad_z * sin(theta)
      y_base <- rep(base_center[2], N+1)
      x_verts <- c(tip[1], x_base, base_center[1])
      y_verts <- c(tip[2], y_base, base_center[2])
      z_verts <- c(tip[3], z_base, base_center[3])
    } else { # z
      tip <- c(0, 0, tip_pos)
      base_center <- c(0, 0, tip_pos - len_val)
      theta <- seq(0, 2*pi, length.out=N+1)
      x_base <- 0 + rad_x * cos(theta)
      y_base <- 0 + rad_y * sin(theta)
      z_base <- rep(base_center[3], N+1)
      x_verts <- c(tip[1], x_base, base_center[1])
      y_verts <- c(tip[2], y_base, base_center[2])
      z_verts <- c(tip[3], z_base, base_center[3])
    }
    
    i <- c(); j <- c(); k <- c()
    for (m in 1:N) { i <- c(i, 0); j <- c(j, m); k <- c(k, m+1) }
    center_idx <- N + 1 
    for (m in 1:N) { i <- c(i, center_idx); j <- c(j, m+1); k <- c(k, m) }
    
    return(list(x=x_verts, y=y_verts, z=z_verts, i=i, j=j, k=k))
}
mesh_idx <- list(i = c(0, 0), j = c(1, 2), k = c(2, 3))


# ------------------------------------------------------------------------------
# 4. Main Visualization Function / 
# ------------------------------------------------------------------------------
run_visualization <- function(data_subset, version_name, show_bulk_buttons) {
  
  cat(sprintf("\n[%s] Running visualization for: %s\n", Sys.time(), version_name))
  
  if (nrow(data_subset) == 0) {
    warning("No data provided for this version. Skipping.")
    return()
  }
  
  # Calculate Limits specific to this subset
  x <- data_subset$b_aging
  y <- data_subset$b_RA
  z <- data_subset$b_HZ
  
  xlim <- get_symmetric_limit(x)
  ylim <- get_symmetric_limit(y)
  zlim <- get_symmetric_limit(z)
  
  cat(sprintf("Axis limits for %s:\nX: %s to %s\nY: %s to %s\nZ: %s to %s\n", 
              version_name, xlim[1], xlim[2], ylim[1], ylim[2], zlim[1], zlim[2]))

  # Create Output Subdirectories
  version_dir <- file.path(output_dir, version_name)
  plots_static_dir <- file.path(version_dir, "Plots", "Static")
  plots_interactive_dir <- file.path(version_dir, "Plots", "Interactive")
  data_dir <- file.path(version_dir, "Data")
  readme_dir <- file.path(version_dir, "Readme")
  
  dir.create(plots_static_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(plots_interactive_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(data_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(readme_dir, recursive = TRUE, showWarnings = FALSE)

  # --- Static Plotting ---
  plot_3d_data <- function(d_sub, suffix = "", title_suffix = "") {
    if (nrow(d_sub) == 0) return()
    
    # Ensure color exists (it should be in d_sub)
    if (!"color" %in% names(d_sub)) d_sub$color <- cluster_colors[d_sub$Cluster_ID]
    
    pdf_file <- file.path(plots_static_dir, paste0("3D_Scatter_Plot", suffix, ".pdf"))
    png_file <- file.path(plots_static_dir, paste0("3D_Scatter_Plot", suffix, ".png"))
    
    draw_curr_plot <- function() {
      s3d <- scatterplot3d(
        x = d_sub$b_aging, y = d_sub$b_RA, z = d_sub$b_HZ,
        color = d_sub$color, pch = 16, cex.symbols = 1.2,     
        xlim = xlim, ylim = ylim, zlim = zlim,
        xlab = "mvAge", ylab = "RA", zlab = "HZ",
        main = paste("3D Scatter Plot", title_suffix),
        grid = TRUE, box = TRUE, type = "h", scale.y = 0.7          
      )
      
      uniq_subset_clusters <- unique(d_sub$Cluster_ID)
      if(length(uniq_subset_clusters) > 0) {
        if (length(uniq_subset_clusters) > 20) {
            legend_clusters <- head(uniq_subset_clusters, 20)
            legend_title <- "Clusters (Top 20)"
        } else {
            legend_clusters <- uniq_subset_clusters
            legend_title <- "Clusters"
        }
        legend("topright", legend = legend_clusters,
               col = cluster_colors[legend_clusters], pch = 16, inset = -0.1, xpd = TRUE, cex = 0.6, bty = "n", title = legend_title)
      }
      s3d$points3d(0, 0, 0, col = "black", pch = 3, cex = 2, lwd = 2)
    }
    
    pdf(pdf_file, width = 12, height = 10)
    draw_curr_plot()
    dev.off()
    
    png(png_file, width = 2400, height = 2000, res = 300)
    draw_curr_plot()
    dev.off()
  }
  
  plot_3d_data(data_subset, suffix = "_All", title_suffix = "(All Signals)")
  
  if ("is_association_signal" %in% names(data_subset)) {
    assoc_data <- data_subset %>% filter(is_association_signal == TRUE | is_association_signal == "TRUE")
    plot_3d_data(assoc_data, suffix = "_Associations_Only", title_suffix = "(Associations Only)")
  }
  
  if (show_bulk_buttons) {
      bulk_d <- data_subset[data_subset$is_bulk, ]
      sc_d <- data_subset[!data_subset$is_bulk, ]
      plot_3d_data(bulk_d, suffix = "_Bulk", title_suffix = "(Bulk Only)")
      plot_3d_data(sc_d, suffix = "_SingleCell", title_suffix = "(Single Cell Only)")
  }

  # --- Interactive Plotting ---
  cat(sprintf("[%s] Generating interactive plot for %s...\n", Sys.time(), version_name))
  
  # Calculate Regression Equations for each cluster
  cluster_eqs <- list()
  sorted_clusters <- sort(unique(data_subset$Cluster_ID))
  for (clus in sorted_clusters) {
    clus_data <- data_subset %>% filter(Cluster_ID == clus)
    # Use tryCatch to handle potential errors
    tryCatch({
      if (nrow(clus_data) >= 3) {
        # HZ ~ Aging + RA
        model <- lm(b_HZ ~ b_aging + b_RA, data = clus_data)
        coefs <- coef(model)
        r_sq <- summary(model)$r.squared
        
        # Format equation
        eq_str <- sprintf("HZ = %.3f*mvAge %+.3f*RA %+.3f", 
                          coefs["b_aging"], coefs["b_RA"], coefs["(Intercept)"])
        eq_str <- paste0(eq_str, " (R2=", round(r_sq, 3), ")")
        cluster_eqs[[clus]] <- eq_str
      } else {
        cluster_eqs[[clus]] <- "Not enough data (<3 points)"
      }
    }, error = function(e) {
      cluster_eqs[[clus]] <- "Regression calculation failed"
    })
  }
  
  tryCatch({
    trace_registry <- list()
    fig <- plot_ly()
    
    # 1. Static Traces (Planes, Axes)
    # Planes
    fig <- fig %>% add_trace(
      x = c(xlim[1], xlim[2], xlim[2], xlim[1]), y = c(ylim[1], ylim[1], ylim[2], ylim[2]), z = c(0, 0, 0, 0),
      i = mesh_idx$i, j = mesh_idx$j, k = mesh_idx$k, type = "mesh3d", opacity = 0.1, color = I("blue"), 
      lighting = list(ambient = 0.6, diffuse = 0.5, specular = 0.1), name = "Plane Z=0 (XY)", showlegend = TRUE, hoverinfo = "text", text = "Plane Z=0 (XY)"
    )
    trace_registry[[length(trace_registry) + 1]] <- list(type = "static")
    
    fig <- fig %>% add_trace(
      x = c(xlim[1], xlim[2], xlim[2], xlim[1]), y = c(0, 0, 0, 0), z = c(zlim[1], zlim[1], zlim[2], zlim[2]), 
      i = mesh_idx$i, j = mesh_idx$j, k = mesh_idx$k, type = "mesh3d", opacity = 0.1, color = I("green"), 
      lighting = list(ambient = 0.6, diffuse = 0.5, specular = 0.1), name = "Plane Y=0 (XZ)", showlegend = TRUE, hoverinfo = "text", text = "Plane Y=0 (XZ)"
    )
    trace_registry[[length(trace_registry) + 1]] <- list(type = "static")

    fig <- fig %>% add_trace(
      x = c(0, 0, 0, 0), y = c(ylim[1], ylim[2], ylim[2], ylim[1]), z = c(zlim[1], zlim[1], zlim[2], zlim[2]), 
      i = mesh_idx$i, j = mesh_idx$j, k = mesh_idx$k, type = "mesh3d", opacity = 0.1, color = I("red"), 
      lighting = list(ambient = 0.6, diffuse = 0.5, specular = 0.1), name = "Plane X=0 (YZ)", showlegend = TRUE, hoverinfo = "text", text = "Plane X=0 (YZ)"
    )
    trace_registry[[length(trace_registry) + 1]] <- list(type = "static")

    # Axes
    x_range_val <- xlim[2] - xlim[1]; y_range_val <- ylim[2] - ylim[1]; z_range_val <- zlim[2] - zlim[1]
    ranges <- c(x_range_val, y_range_val, z_range_val)
    
    # X Axis
    fig <- fig %>% add_trace(x = c(xlim[1], xlim[2]), y = c(0, 0), z = c(0, 0), type = "scatter3d", mode = "lines", line = list(color = "black", width = 5), name = "mvAge", showlegend = FALSE, hoverinfo = "none")
    trace_registry[[length(trace_registry) + 1]] <- list(type = "static")
    x_arrow <- create_arrow_mesh("x", xlim[2], ranges)
    fig <- fig %>% add_trace(x = x_arrow$x, y = x_arrow$y, z = x_arrow$z, i = x_arrow$i, j = x_arrow$j, k = x_arrow$k, type = "mesh3d", color = I("black"), opacity = 1.0, name = "X-axis Arrow", showlegend = FALSE, hoverinfo = "none")
    trace_registry[[length(trace_registry) + 1]] <- list(type = "static")
    fig <- fig %>% add_trace(x = c(xlim[2]), y = c(y_range_val * 0.05), z = c(z_range_val * 0.05), type = "scatter3d", mode = "text", text = c("mvAge"), textposition = "middle right", textfont = list(color = "black", size = 12, weight = "bold"), name = "X-axis Label", showlegend = FALSE, hoverinfo = "none")
    trace_registry[[length(trace_registry) + 1]] <- list(type = "static")
    
    # Y Axis
    fig <- fig %>% add_trace(x = c(0, 0), y = c(ylim[1], ylim[2]), z = c(0, 0), type = "scatter3d", mode = "lines", line = list(color = "black", width = 5), name = "RA", showlegend = FALSE, hoverinfo = "none")
    trace_registry[[length(trace_registry) + 1]] <- list(type = "static")
    y_arrow <- create_arrow_mesh("y", ylim[2], ranges)
    fig <- fig %>% add_trace(x = y_arrow$x, y = y_arrow$y, z = y_arrow$z, i = y_arrow$i, j = y_arrow$j, k = y_arrow$k, type = "mesh3d", color = I("black"), opacity = 1.0, name = "Y-axis Arrow", showlegend = FALSE, hoverinfo = "none")
    trace_registry[[length(trace_registry) + 1]] <- list(type = "static")
    fig <- fig %>% add_trace(x = c(x_range_val * 0.05), y = c(ylim[2]), z = c(z_range_val * 0.05), type = "scatter3d", mode = "text", text = c("RA"), textposition = "middle right", textfont = list(color = "black", size = 12, weight = "bold"), name = "Y-axis Label", showlegend = FALSE, hoverinfo = "none")
    trace_registry[[length(trace_registry) + 1]] <- list(type = "static")
    
    # Z Axis
    fig <- fig %>% add_trace(x = c(0, 0), y = c(0, 0), z = c(zlim[1], zlim[2]), type = "scatter3d", mode = "lines", line = list(color = "black", width = 5), name = "HZ", showlegend = FALSE, hoverinfo = "none")
    trace_registry[[length(trace_registry) + 1]] <- list(type = "static")
    z_arrow <- create_arrow_mesh("z", zlim[2], ranges)
    fig <- fig %>% add_trace(x = z_arrow$x, y = z_arrow$y, z = z_arrow$z, i = z_arrow$i, j = z_arrow$j, k = z_arrow$k, type = "mesh3d", color = I("black"), opacity = 1.0, name = "Z-axis Arrow", showlegend = FALSE, hoverinfo = "none")
    trace_registry[[length(trace_registry) + 1]] <- list(type = "static")
    fig <- fig %>% add_trace(x = c(x_range_val * 0.05), y = c(y_range_val * 0.05), z = c(zlim[2]), type = "scatter3d", mode = "text", text = c("HZ"), textposition = "middle right", textfont = list(color = "black", size = 12, weight = "bold"), name = "Z-axis Label", showlegend = FALSE, hoverinfo = "none")
    trace_registry[[length(trace_registry) + 1]] <- list(type = "static")

    # 2. Overview Traces
    sorted_clusters <- sort(unique(data_subset$Cluster_ID))
    
    for (clus in sorted_clusters) {
      clus_data <- data_subset %>% filter(Cluster_ID == clus)
      
      is_assoc <- if ("is_association_signal" %in% names(clus_data)) {
          clus_data$is_association_signal %in% c(TRUE, "TRUE")
      } else rep(FALSE, nrow(clus_data))
      
      is_bulk_rows <- clus_data$is_bulk
      
      subsets <- list(
          list(data = clus_data[is_assoc & is_bulk_rows, ], type = "assoc", is_bulk = TRUE),
          list(data = clus_data[is_assoc & !is_bulk_rows, ], type = "assoc", is_bulk = FALSE),
          list(data = clus_data[!is_assoc & is_bulk_rows, ], type = "non_assoc", is_bulk = TRUE),
          list(data = clus_data[!is_assoc & !is_bulk_rows, ], type = "non_assoc", is_bulk = FALSE)
      )
      
      for (sub in subsets) {
          if (nrow(sub$data) > 0) {
              is_assoc_sub <- sub$type == "assoc"
              is_bulk_sub <- sub$is_bulk
              
              marker_sym <- if (is_assoc_sub) 'diamond' else 'circle'
              marker_size <- if (is_assoc_sub) 6 else 5
              marker_opacity <- if (is_assoc_sub) 1 else 0.8
              
              fig <- fig %>% add_trace(
                  data = sub$data, x = ~b_aging, y = ~b_RA, z = ~b_HZ,
                  type = 'scatter3d', mode = 'markers',
                  marker = list(size = marker_size, opacity = marker_opacity, symbol = marker_sym, color = cluster_colors[clus]), 
                  text = ~hover_text,
                  hoverinfo = "text",
                  name = clus, legendgroup = clus, 
                  showlegend = is_assoc_sub, 
                  visible = TRUE,
                  meta = list(type = "overview", cluster = clus, subtype = sub$type, is_bulk = sub$is_bulk)
              )
              
              trace_registry[[length(trace_registry) + 1]] <- list(
                  type = "overview", cluster = clus, subtype = sub$type, 
                  is_bulk = sub$is_bulk, is_assoc = is_assoc_sub
              )
          }
      }
    }

    # 3. Detail Traces
    for (clus in sorted_clusters) {
      clus_data <- data_subset %>% filter(Cluster_ID == clus)
      clus_genes <- sort(unique(clus_data$gene))
      
      for (g in clus_genes) {
        gene_sub <- clus_data %>% filter(gene == g)
        
        is_assoc_g <- if ("is_association_signal" %in% names(gene_sub)) {
            gene_sub$is_association_signal %in% c(TRUE, "TRUE")
        } else rep(FALSE, nrow(gene_sub))
        
        is_bulk_g <- gene_sub$is_bulk
        
        subsets <- list(
            list(data = gene_sub[is_assoc_g & is_bulk_g, ], type = "assoc", is_bulk = TRUE),
            list(data = gene_sub[is_assoc_g & !is_bulk_g, ], type = "assoc", is_bulk = FALSE),
            list(data = gene_sub[!is_assoc_g & is_bulk_g, ], type = "non_assoc", is_bulk = TRUE),
            list(data = gene_sub[!is_assoc_g & !is_bulk_g, ], type = "non_assoc", is_bulk = FALSE)
        )
        
        curr_cluster_col <- cluster_colors[clus]
        curr_gene_col <- gene_colors[g]
        has_assoc <- any(is_assoc_g)
        
        for (sub in subsets) {
          if (nrow(sub$data) > 0) {
              is_assoc_sub <- sub$type == "assoc"
              is_bulk_sub <- sub$is_bulk
              
              marker_sym <- if (is_assoc_sub) 'diamond' else 'circle'
              marker_size <- if (is_assoc_sub) 6 else 5
              marker_opacity <- if (is_assoc_sub) 1 else 0.8
              do_show_leg <- if (has_assoc) is_assoc_sub else TRUE
              
              fig <- fig %>% add_trace(
                  data = sub$data, x = ~b_aging, y = ~b_RA, z = ~b_HZ,
                  type = 'scatter3d', mode = 'markers',
                  marker = list(size = marker_size, opacity = marker_opacity, symbol = marker_sym, color = curr_cluster_col), 
                  text = ~hover_text,
                  hoverinfo = "text",
                  customdata = ~paste(Cluster_ID, gene, SNP, sep="|"),
                  name = g, legendgroup = g, 
                  showlegend = do_show_leg, 
                  visible = FALSE,
                  meta = list(type = "detail", cluster = clus, gene = g, subtype = sub$type, is_bulk = sub$is_bulk, cluster_col = curr_cluster_col, gene_col = curr_gene_col)
              )
              
              trace_registry[[length(trace_registry) + 1]] <- list(
                  type = "detail", cluster = clus, gene = g, subtype = sub$type,
                  is_bulk = sub$is_bulk, is_assoc = is_assoc_sub,
                  cluster_col = curr_cluster_col, gene_col = curr_gene_col
              )
          }
        }
      }
    }
    
    # 4. Regression Lines (New: PCA-based)
    # 4a. Cluster-specific Lines
    for (clus in sorted_clusters) {
      clus_data <- data_subset %>% filter(Cluster_ID == clus)
      
      if (nrow(clus_data) >= 2) { # Need at least 2 points for a line
        # Extract coordinates
        coords <- clus_data %>% select(b_aging, b_RA, b_HZ) %>% as.matrix()
        
        # Calculate PCA (Principal Component Analysis)
        # The first principal component (PC1) is the line that minimizes the sum of squared orthogonal distances
        pca_res <- prcomp(coords, center = TRUE, scale. = FALSE)
        
        # Center of mass
        center_pt <- pca_res$center
        
        # Direction vector (PC1)
        direction_vec <- pca_res$rotation[, 1]
        
        # Project data points onto the line to find the range
        # t = (P - center) . direction
        centered_coords <- scale(coords, center = center_pt, scale = FALSE)
        t_vals <- centered_coords %*% direction_vec
        
        # Extend the line slightly beyond the data range (e.g., 10% padding)
        t_range <- range(t_vals)
        padding <- diff(t_range) * 0.1
        if (padding == 0) padding <- 0.1 # Handle case where all points are identical
        t_min <- t_range[1] - padding
        t_max <- t_range[2] + padding
        
        # Calculate start and end points of the line segment
        p_start <- center_pt + t_min * direction_vec
        p_end <- center_pt + t_max * direction_vec
        
        # Add line trace
        line_color <- cluster_colors[clus]
        
        # Format equation text for hover (approximate direction)
        # Line: P = Center + t * Direction
        eq_text <- sprintf("Line Direction (PC1):<br>mvAge: %.3f<br>RA: %.3f<br>HZ: %.3f<br>Center: (%.3f, %.3f, %.3f)",
                           direction_vec[1], direction_vec[2], direction_vec[3],
                           center_pt[1], center_pt[2], center_pt[3])

        fig <- fig %>% add_trace(
            x = c(p_start[1], p_end[1]), 
            y = c(p_start[2], p_end[2]), 
            z = c(p_start[3], p_end[3]),
            type = "scatter3d",
            mode = "lines",
            line = list(color = line_color, width = 6), # Thicker line
            hoverinfo = "text",
            text = eq_text,
            name = paste("Line:", clus),
            meta = list(type = "pca_line", cluster = clus), # Meta for JS control
            showlegend = FALSE,
            visible = FALSE
        )
        
        trace_registry[[length(trace_registry) + 1]] <- list(
            type = "regression_line", cluster = clus
        )

        # --- NEW: Gene Bubbles on PCA Line (PCA Mode) ---
        # Strategy A: Use Representative SNP (Max |PC1|) instead of Centroid
        # This aligns with the unified projection strategy to avoid LD collinearity issues.
        
        # 1. Add PC1 scores to clus_data
        # Ensure row order matches (pca_res was calculated from clus_data)
        clus_data$PC1 <- pca_res$x[, 1]
        
        # 2. Select Representative SNP for each gene
        gene_centroids <- clus_data %>%
          add_count(gene, name = "total_snps") %>%
          group_by(gene) %>%
          slice_max(order_by = abs(PC1), n = 1, with_ties = FALSE) %>%
          rename(
            mean_aging = b_aging, # Keep variable names for compatibility
            mean_RA = b_RA,
            mean_HZ = b_HZ,
            n_snps = total_snps
          ) %>%
          ungroup()
        
        # Project centroids onto the line
        # t = (P - Center) . Direction
        # P_proj = Center + t * Direction
        if (nrow(gene_centroids) > 0) {
            # Prepare vectors
            genes_vec <- gene_centroids$gene
            n_snps_vec <- gene_centroids$n_snps
            
            # Pre-allocate result vectors
            proj_x <- numeric(nrow(gene_centroids))
            proj_y <- numeric(nrow(gene_centroids))
            proj_z <- numeric(nrow(gene_centroids))
            hover_txts <- character(nrow(gene_centroids))
            bubble_colors <- character(nrow(gene_centroids))
            bubble_sizes <- numeric(nrow(gene_centroids))
            
            for (k in 1:nrow(gene_centroids)) {
                P <- c(gene_centroids$mean_aging[k], gene_centroids$mean_RA[k], gene_centroids$mean_HZ[k])
                P_centered <- P - center_pt
                t_proj <- sum(P_centered * direction_vec)
                P_proj <- center_pt + t_proj * direction_vec
                
                proj_x[k] <- P_proj[1]
                proj_y[k] <- P_proj[2]
                proj_z[k] <- P_proj[3]
                
                # Bubble Size (Log scale for better visualization of variance)
                # Base size 10, scaled by log of SNPs
                bubble_sizes[k] <- 8 + 3 * log(n_snps_vec[k])
                
                # Bubble Color (Use Gene Color)
                g_name <- genes_vec[k]
                bubble_colors[k] <- if (g_name %in% names(gene_colors)) gene_colors[g_name] else "grey"
                
                hover_txts[k] <- paste0(
                    "<b>Gene: ", g_name, "</b><br>",
                    "Cluster: ", clus, "<br>",
                    "SNPs: ", n_snps_vec[k], "<br>",
                    "Projected Pos: (", round(P_proj[1], 3), ", ", round(P_proj[2], 3), ", ", round(P_proj[3], 3), ")"
                )
            }
            
            # Add Bubble Trace (One trace per cluster for all genes to save traces, or one per gene?)
            # One trace per cluster is more efficient for JS toggling.
            fig <- fig %>% add_trace(
                x = proj_x, y = proj_y, z = proj_z,
                type = "scatter3d",
                mode = "markers",
                marker = list(
                    size = bubble_sizes, 
                    color = bubble_colors, 
                    opacity = 0.9, 
                    line = list(color = 'black', width = 1)
                ),
                text = hover_txts,
                hoverinfo = "text",
                name = paste("PCA Bubbles:", clus),
                meta = list(type = "pca_gene_bubble", cluster = clus),
                showlegend = FALSE,
                visible = FALSE
            )
            
            trace_registry[[length(trace_registry) + 1]] <- list(
                type = "pca_gene_bubble", cluster = clus
            )
        }

      }
    }

    # 4b. Global Regression Line (Overview - All Signals)
    # Only if we have enough data overall
    if (nrow(data_subset) >= 2) {
      coords_global <- data_subset %>% select(b_aging, b_RA, b_HZ) %>% as.matrix()
      pca_res_global <- prcomp(coords_global, center = TRUE, scale. = FALSE)
      center_pt_global <- pca_res_global$center
      direction_vec_global <- pca_res_global$rotation[, 1]
      centered_coords_global <- scale(coords_global, center = center_pt_global, scale = FALSE)
      t_vals_global <- centered_coords_global %*% direction_vec_global
      t_range_global <- range(t_vals_global)
      padding_global <- diff(t_range_global) * 0.1
      if (padding_global == 0) padding_global <- 0.1
      t_min_global <- t_range_global[1] - padding_global
      t_max_global <- t_range_global[2] + padding_global
      p_start_global <- center_pt_global + t_min_global * direction_vec_global
      p_end_global <- center_pt_global + t_max_global * direction_vec_global
      
      eq_text_global <- sprintf("Global Trend (Overview):<br>mvAge: %.3f<br>RA: %.3f<br>HZ: %.3f<br>Center: (%.3f, %.3f, %.3f)",
                         direction_vec_global[1], direction_vec_global[2], direction_vec_global[3],
                         center_pt_global[1], center_pt_global[2], center_pt_global[3])

      fig <- fig %>% add_trace(
          x = c(p_start_global[1], p_end_global[1]), 
          y = c(p_start_global[2], p_end_global[2]), 
          z = c(p_start_global[3], p_end_global[3]),
          type = "scatter3d",
          mode = "lines",
          line = list(color = "black", width = 8, dash = "solid"), # Distinct style
          hoverinfo = "text",
          text = eq_text_global,
          name = "Global Trend (Overview)",
          meta = list(type = "pca_line", cluster = "Overview"), # Meta for JS control
          showlegend = TRUE,
          visible = FALSE # Default Hidden
      )
      
      trace_registry[[length(trace_registry) + 1]] <- list(
          type = "global_regression_line"
      )

      # --- NEW: Global Gene Bubbles on Global PCA Line ---
      # Strategy A: Use Representative SNP (Max |Global_PC1|)
      
      # 1. Add Global PC1 scores
      data_subset$Global_PC1 <- pca_res_global$x[, 1]
      
      # 2. Select Representative SNP
      global_gene_centroids <- data_subset %>%
          add_count(Cluster_ID, gene, name = "total_snps") %>%
          group_by(Cluster_ID, gene) %>%
          slice_max(order_by = abs(Global_PC1), n = 1, with_ties = FALSE) %>%
          rename(
            mean_aging = b_aging,
            mean_RA = b_RA,
            mean_HZ = b_HZ,
            n_snps = total_snps
          ) %>%
          ungroup()

      if (nrow(global_gene_centroids) > 0) {
           # Prepare vectors
           g_genes_vec <- global_gene_centroids$gene
           g_clus_vec <- global_gene_centroids$Cluster_ID
           g_n_snps_vec <- global_gene_centroids$n_snps
           
           # Pre-allocate result vectors
           g_proj_x <- numeric(nrow(global_gene_centroids))
           g_proj_y <- numeric(nrow(global_gene_centroids))
           g_proj_z <- numeric(nrow(global_gene_centroids))
           g_hover_txts <- character(nrow(global_gene_centroids))
           g_bubble_colors <- character(nrow(global_gene_centroids))
           g_bubble_sizes <- numeric(nrow(global_gene_centroids))
           
           for (k in 1:nrow(global_gene_centroids)) {
                P <- c(global_gene_centroids$mean_aging[k], global_gene_centroids$mean_RA[k], global_gene_centroids$mean_HZ[k])
                P_centered <- P - center_pt_global
                t_proj <- sum(P_centered * direction_vec_global)
                P_proj <- center_pt_global + t_proj * direction_vec_global
                
                g_proj_x[k] <- P_proj[1]
                g_proj_y[k] <- P_proj[2]
                g_proj_z[k] <- P_proj[3]
                
                # Bubble Size
                g_bubble_sizes[k] <- 8 + 3 * log(g_n_snps_vec[k])
                
                # Bubble Color (Gene Color)
                g_name <- g_genes_vec[k]
                g_bubble_colors[k] <- if (g_name %in% names(gene_colors)) gene_colors[g_name] else "grey"
                
                g_hover_txts[k] <- paste0(
                    "<b>Gene: ", g_name, "</b><br>",
                    "Cluster: ", g_clus_vec[k], "<br>",
                    "SNPs: ", g_n_snps_vec[k], "<br>",
                    "Projected Pos: (", round(P_proj[1], 3), ", ", round(P_proj[2], 3), ", ", round(P_proj[3], 3), ")"
                )
           }
           
           fig <- fig %>% add_trace(
                x = g_proj_x, y = g_proj_y, z = g_proj_z,
                type = "scatter3d",
                mode = "markers",
                marker = list(
                    size = g_bubble_sizes, 
                    color = g_bubble_colors, 
                    opacity = 0.9, 
                    line = list(color = 'black', width = 1)
                ),
                text = g_hover_txts,
                hoverinfo = "text",
                name = "Global PCA Bubbles",
                meta = list(type = "pca_gene_bubble", cluster = "Overview"),
                showlegend = FALSE,
                visible = FALSE
           )
           
           trace_registry[[length(trace_registry) + 1]] <- list(
                type = "pca_gene_bubble", cluster = "Overview"
           )
      }
    }
    
    # 5. Controls
    # Only Signal Set Dropdown remains here. Other controls moved to JS.
    
    # Updated: Overview excludes PCA lines (controlled by JS button)
    vis_overview <- sapply(trace_registry, function(x) x$type == "static" || x$type == "overview")
    
    signal_set_buttons <- list()
    
    signal_set_buttons[[1]] <- list(
      label = "Overview (All Signals)", method = "update",
      args = list(list(), list(title = "3D Scatter Plot: Overview", "annotations[0].text" = "<b>Visible SNPs:</b> Calculating..."))
    )
    
    for (clus in sorted_clusters) {
      # Updated: Cluster view excludes PCA lines (controlled by JS button)
      vis_clus <- sapply(trace_registry, function(x) { 
          x$type == "static" || 
          (x$type == "detail" && x$cluster == clus)
      })
      
      signal_set_buttons[[length(signal_set_buttons) + 1]] <- list(
          label = paste("Signal Set:", clus), method = "update",
          args = list(list(), list(title = paste("Signal Set:", clus), "annotations[0].text" = "<b>Visible SNPs:</b> Calculating..."))
      )
    }
    
    updatemenus <- list(
      list(type = "dropdown", direction = "down", active = 0, x = 0.0, xanchor = "left", y = 1.15, yanchor = "top", pad = list(r = 10, t = 10), buttons = signal_set_buttons)
    )
    
    fig <- fig %>% layout(
      title = "3D Scatter Plot: Overview",
      scene = list(xaxis = list(title = '', range = xlim), yaxis = list(title = '', range = ylim), zaxis = list(title = '', range = zlim), aspectmode = 'cube'),
      legend = list(yanchor = "top", y = 0.9, xanchor = "left", x = 1.05, bgcolor = 'rgba(255, 255, 255, 0.5)', font = list(size = 12)),
      annotations = list(
          list(x = 0.5, y = 1.0, xref = "paper", yref = "paper", text = "<b>Visible SNPs:</b> Calculating...", showarrow = FALSE, font = list(size = 14))
      ), 
      updatemenus = updatemenus
    )
    
    # JS Code
    js_code <- "
      function(el, x) {
        var plot = document.getElementById(el.id);

        // --- State Management ---
        var state = {
            pcaVisible: false,
            filterMode: 'all', // 'all' or 'coloc'
            colorMode: 'signal', // 'signal' or 'gene'
            searchActive: false,
            pcaMode: false // New: PCA Mode
        };

        // --- UI Container ---
        var uiContainer = document.createElement('div');
        uiContainer.style.cssText = 'position: absolute; top: 10px; left: 350px; z-index: 1001; display: flex; gap: 8px; align-items: center;';
        el.appendChild(uiContainer);

        // --- Helper: Create Button ---
        function createBtn(id, text, title, onClick) {
            var btn = document.createElement('button');
            btn.id = id;
            btn.innerHTML = text;
            btn.title = title;
            btn.style.cssText = 'padding: 6px 12px; background: #6c757d; color: white; border: none; border-radius: 4px; cursor: pointer; font-family: Arial, sans-serif; font-size: 12px; font-weight: bold; box-shadow: 0 2px 5px rgba(0,0,0,0.2); transition: all 0.2s; min-width: 120px;';
            btn.onclick = onClick;
            uiContainer.appendChild(btn);
            return btn;
        }

        // --- Buttons ---
        // 1. View Filter Button
        var btnFilter = createBtn('btn_filter', 'View: Show All', 'Toggle between All Signals and Colocalization Support Only', function() {
            state.filterMode = (state.filterMode === 'all') ? 'coloc' : 'all';
            updateVisualization();
        });

        // 2. Color Button
        var btnColor = createBtn('btn_color', 'Color: Signal Set', 'Toggle coloring by Signal Set or Gene', function() {
            state.colorMode = (state.colorMode === 'signal') ? 'gene' : 'signal';
            updateVisualization();
        });

        // 3. PCA Lines Button
        var btnPCA = createBtn('btn_pca', 'Show PCA Lines', 'Toggle PCA Regression Lines', function() {
            state.pcaVisible = !state.pcaVisible;
            updateVisualization();
        });
        
        // 4. PCA Mode Button
        var btnPcaMode = createBtn('btn_pca_mode', 'PCA Mode: Off', 'Toggle PCA Dimensionality Reduction View', function() {
            state.pcaMode = !state.pcaMode;
            updateVisualization();
        });

        // --- Update UI Styles ---
        function updateUI() {
            // PCA Button
            if (state.pcaVisible) {
                btnPCA.innerHTML = 'Hide PCA Lines';
                btnPCA.style.background = '#28a745'; // Green
            } else {
                btnPCA.innerHTML = 'Show PCA Lines';
                btnPCA.style.background = '#6c757d'; // Grey
            }

            // Filter Button
            if (state.filterMode === 'coloc') {
                btnFilter.innerHTML = 'View: Coloc Only';
                btnFilter.style.background = '#007bff'; // Blue
            } else {
                btnFilter.innerHTML = 'View: Show All';
                btnFilter.style.background = '#6c757d'; // Grey
            }

            // Color Button
            if (state.colorMode === 'gene') {
                btnColor.innerHTML = 'Color: Gene';
                btnColor.style.background = '#17a2b8'; // Cyan
            } else {
                btnColor.innerHTML = 'Color: Signal Set';
                btnColor.style.background = '#6c757d'; // Grey
            }
            
            // PCA Mode Button
            if (state.pcaMode) {
                btnPcaMode.innerHTML = 'PCA Mode: On';
                btnPcaMode.style.background = '#ffc107'; // Yellow/Orange
                btnPcaMode.style.color = 'black';
                
                // Hide 'Show PCA Lines' simple toggle to avoid confusion
                btnPCA.style.display = 'none'; 
            } else {
                btnPcaMode.innerHTML = 'PCA Mode: Off';
                btnPcaMode.style.background = '#6c757d'; // Grey
                btnPcaMode.style.color = 'white';
                
                // Restore 'Show PCA Lines' simple toggle
                btnPCA.style.display = 'inline-block'; 
            }
        }

        // --- Main Visualization Logic ---
        function updateVisualization() {
            updateUI();
            
            // Get Current Context
            var title = '';
            if (plot.layout.title) {
                 title = typeof plot.layout.title === 'string' ? plot.layout.title : plot.layout.title.text;
            }
            
            var activeCluster = '';
            if (title && title.indexOf('Overview') !== -1) {
                activeCluster = 'Overview';
            } else if (title && title.indexOf('Signal Set:') !== -1) {
                activeCluster = title.split('Signal Set:')[1].trim();
            }

            var update = { visible: [], 'marker.color': [] };
            var visIndices = [];
            var colIndices = [];

            // Iterate Traces
            for (var i = 0; i < plot.data.length; i++) {
                var trace = plot.data[i];
                if (!trace.meta) continue; // Skip unknown traces

                // 1. Color Logic (Always Apply)
                if (trace.meta.type === 'detail') {
                    var newColor = (state.colorMode === 'signal') ? trace.meta.cluster_col : trace.meta.gene_col;
                    // Optimization: Check if color changed? (Skipped for simplicity)
                    update['marker.color'].push(newColor);
                    colIndices.push(i);
                }

                // 2. Visibility Logic
                if (state.searchActive) {
                    // Search Mode: Force hide PCA and Gene Bubbles
                    if (trace.meta.type === 'pca_line' || trace.meta.type === 'pca_gene_bubble') {
                         update.visible.push(false);
                         visIndices.push(i);
                    }
                    // Scatter traces visibility handled by Search function
                } else if (state.pcaMode) {
                    // --- PCA Mode Active ---
                    var shouldBeVisible = false;
                    
                    // Show Static (Axes, Planes)
                    if (trace.meta.type === 'static') shouldBeVisible = true;
                    
                    // Show PCA Lines and Gene Bubbles (matching active cluster)
                    if (trace.meta.type === 'pca_line' || trace.meta.type === 'pca_gene_bubble') {
                         if (activeCluster === 'Overview') {
                             // In Overview, ONLY show Global elements (cluster === 'Overview')
                             if (trace.meta.cluster === 'Overview') shouldBeVisible = true;
                         } else {
                             // Specific Cluster: Show elements for that cluster
                             // (And hide Global elements)
                             if (trace.meta.cluster === activeCluster) shouldBeVisible = true;
                         }
                    }
                    
                    // Hide all Scatter Points (Overview/Detail)
                    if (trace.meta.type === 'overview' || trace.meta.type === 'detail') {
                        shouldBeVisible = false;
                    }

                    update.visible.push(shouldBeVisible);
                    visIndices.push(i);
                } else {
                    // --- Normal Mode ---
                    var shouldBeVisible = false;

                    // A. Base Visibility (Cluster Context)
                    if (activeCluster === 'Overview') {
                        if (trace.meta.type === 'overview' || trace.meta.type === 'static') shouldBeVisible = true;
                        if (trace.meta.type === 'pca_line' && trace.meta.cluster === 'Overview' && state.pcaVisible) shouldBeVisible = true;
                    } else if (activeCluster) {
                        // Specific Cluster
                        if (trace.meta.type === 'detail' && trace.meta.cluster === activeCluster) shouldBeVisible = true;
                        if (trace.meta.type === 'static') shouldBeVisible = true;
                        if (trace.meta.type === 'pca_line' && trace.meta.cluster === activeCluster && state.pcaVisible) shouldBeVisible = true;
                    }

                    // B. Filter Logic (Coloc Support)
                    if (shouldBeVisible && state.filterMode === 'coloc') {
                        // Hide non-assoc traces
                        if (trace.meta.subtype === 'non_assoc') shouldBeVisible = false;
                    }
                    
                    // Ensure PCA Gene Bubbles are HIDDEN in Normal Mode
                    if (trace.meta.type === 'pca_gene_bubble') shouldBeVisible = false;

                    // Push Update
                    update.visible.push(shouldBeVisible);
                    visIndices.push(i);
                }
            }

            // Apply Updates
            if (visIndices.length > 0) {
                Plotly.restyle(el.id, {visible: update.visible}, visIndices).then(function() {
                    updateCount();
                });
            } else {
                updateCount();
            }
            if (colIndices.length > 0) Plotly.restyle(el.id, {'marker.color': update['marker.color']}, colIndices);
        }

        // Hook into Plotly events (e.g., Dropdown changes title -> triggers update)
        el.on('plotly_update', function() {
            // Add delay to let Plotly finish its internal updates
            setTimeout(updateVisualization, 200); 
        });

        // --- Search / Filter Code ---
        var searchContainer = document.createElement('div');
        searchContainer.style.cssText = 'position: absolute; bottom: 20px; left: 20px; width: 220px; background: rgba(255, 255, 255, 0.95); padding: 15px; border-radius: 8px; box-shadow: 0 4px 10px rgba(0,0,0,0.15); z-index: 1000; font-family: Arial, sans-serif; border: 1px solid #ddd;';
        searchContainer.innerHTML = `
              <h4 style=\"margin: 0 0 10px 0; font-size: 14px; color: #333; border-bottom: 1px solid #eee; padding-bottom: 5px;\">Search / Filter</h4>
              <div style=\"margin-bottom: 8px;\"><label style=\"display: block; font-size: 11px; margin-bottom: 2px; color: #555;\">Signal Set:</label><input type=\"text\" id=\"search_signal\" placeholder=\"e.g. bmem_Signal_5\" style=\"width: 100%; padding: 5px; box-sizing: border-box; border: 1px solid #ccc; border-radius: 4px; font-size: 12px;\"></div>
              <div style=\"margin-bottom: 8px;\"><label style=\"display: block; font-size: 11px; margin-bottom: 2px; color: #555;\">Gene:</label><input type=\"text\" id=\"search_gene\" placeholder=\"e.g. TNF\" style=\"width: 100%; padding: 5px; box-sizing: border-box; border: 1px solid #ccc; border-radius: 4px; font-size: 12px;\"></div>
              <div style=\"margin-bottom: 12px;\"><label style=\"display: block; font-size: 11px; margin-bottom: 2px; color: #555;\">SNP:</label><input type=\"text\" id=\"search_snp\" placeholder=\"e.g. rs123\" style=\"width: 100%; padding: 5px; box-sizing: border-box; border: 1px solid #ccc; border-radius: 4px; font-size: 12px;\"></div>
              <div style=\"display: flex; gap: 8px;\"><button id=\"btn_search\" style=\"flex: 1; padding: 6px; background: #007bff; color: white; border: none; border-radius: 4px; cursor: pointer; font-size: 12px; transition: background 0.2s;\">Search</button><button id=\"btn_clear\" style=\"flex: 1; padding: 6px; background: #6c757d; color: white; border: none; border-radius: 4px; cursor: pointer; font-size: 12px; transition: background 0.2s;\">Clear</button></div>
              <div id=\"search_status\" style=\"margin-top: 8px; font-size: 11px; color: #666; min-height: 14px; text-align: center;\"></div>`;
        el.appendChild(searchContainer);
        var btnSearch = searchContainer.querySelector('#btn_search');
        var btnClear = searchContainer.querySelector('#btn_clear');
        btnSearch.onmouseover = function() { this.style.background = '#0056b3'; }; btnSearch.onmouseout = function() { this.style.background = '#007bff'; };
        btnClear.onmouseover = function() { this.style.background = '#5a6268'; }; btnClear.onmouseout = function() { this.style.background = '#6c757d'; };
        
        var originalData = []; var isDataBackedUp = false;
        function backupData() {
          if (isDataBackedUp) return;
          for (var i = 0; i < plot.data.length; i++) {
              var trace = plot.data[i];
              // Backup Scatter3d traces (Overview & Detail)
              if (trace.type === 'scatter3d' && trace.mode && trace.mode.indexOf('markers') !== -1) {
                  var len = trace.x ? trace.x.length : 0;
                  var txt = trace.text;
                  if (txt && !Array.isArray(txt)) txt = Array(len).fill(txt);
                  
                  var htxt = trace.hovertext;
                  if (htxt && !Array.isArray(htxt)) htxt = Array(len).fill(htxt);
                  
                  if (!txt && htxt) txt = htxt;
                  if (!htxt && txt) htxt = txt;
                  
                  originalData[i] = { 
                      x: trace.x ? trace.x.slice() : [], 
                      y: trace.y ? trace.y.slice() : [], 
                      z: trace.z ? trace.z.slice() : [], 
                      text: txt ? txt.slice() : [], 
                      hovertext: htxt ? htxt.slice() : [],
                      customdata: trace.customdata ? trace.customdata.slice() : [], 
                      showlegend: (trace.showlegend === undefined ? true : trace.showlegend) 
                  };
              } else { originalData[i] = null; }
          }
          isDataBackedUp = true;
        }

        function performSearch() {
            backupData();
            // Set Search Active State
            state.searchActive = true;
            state.pcaVisible = false; // Force hide PCA
            updateVisualization(); // Update UI and hide PCA

            var signalQ = document.getElementById('search_signal').value.toLowerCase().trim();
            var signalKeywords = signalQ ? signalQ.split(/\\s+/) : [];
            var geneQ = document.getElementById('search_gene').value.toLowerCase().trim();
            var geneKeywords = geneQ ? geneQ.split(/\\s+/) : [];
            var snpQ = document.getElementById('search_snp').value.toLowerCase().trim();
            var snpKeywords = snpQ ? snpQ.split(/\\s+/) : [];
            
            var statusDiv = document.getElementById('search_status');
            if (signalKeywords.length === 0 && geneKeywords.length === 0 && snpKeywords.length === 0) { statusDiv.innerText = 'Please enter a keyword.'; return; }
            statusDiv.innerText = 'Searching...';
            
            var update = { x: [], y: [], z: [], text: [], hovertext: [], customdata: [], showlegend: [], hoverinfo: [] };
            var indices = []; var visibleIndices = []; var hiddenIndices = []; var matchCount = 0;
            
            for (var i = 0; i < plot.data.length; i++) {
                var trace = plot.data[i]; var orig = originalData[i];
                // Only search in Detail traces (Overview traces don't have detailed customdata usually, or we skip them?)
                // Wait, original code logic:
                // isDetail = scatter3d && customdata !== undefined
                // isOverview = scatter3d && customdata === undefined && !lines
                // isStatic = mesh3d or lines
                
                // Let's adapt original logic
                var isDetail = (trace.type === 'scatter3d' && trace.customdata !== undefined);
                var isOverview = (trace.type === 'scatter3d' && trace.customdata === undefined && (!trace.mode || trace.mode.indexOf('lines') === -1));
                var isStatic = (trace.type === 'mesh3d' || (trace.type === 'scatter3d' && trace.mode && (trace.mode.indexOf('lines') !== -1 || trace.mode.indexOf('text') !== -1)));
                var isPcaLine = (trace.meta && trace.meta.type === 'pca_line');

                // Also respect meta if available
                if (trace.meta) {
                    if (trace.meta.type === 'detail') isDetail = true;
                    if (trace.meta.type === 'overview') isOverview = true;
                    if (trace.meta.type === 'static') isStatic = true;
                }

                if (isPcaLine) {
                    hiddenIndices.push(i);
                } else if (isDetail && orig) {
                    var keep = [];
                    // Search Logic...
                    for (var j = 0; j < orig.customdata.length; j++) {
                        var parts = (orig.customdata[j] || '').split('|');
                        if (parts.length < 3) continue; 
                        var s_set = parts[0].toLowerCase(); var g = parts[1].toLowerCase(); var s = parts[2].toLowerCase();
                        var match = true;
                        if (signalKeywords.length > 0) { var sm = false; for(var k=0; k<signalKeywords.length; k++) if(s_set.indexOf(signalKeywords[k]) !== -1) { sm=true; break; } if(!sm) match=false; }
                        if (match && geneKeywords.length > 0) { var gm = false; for(var k=0; k<geneKeywords.length; k++) if(g.indexOf(geneKeywords[k]) !== -1) { gm=true; break; } if(!gm) match=false; }
                        if (match && snpKeywords.length > 0) { var sm = false; for(var k=0; k<snpKeywords.length; k++) if(s.indexOf(snpKeywords[k]) !== -1) { sm=true; break; } if(!sm) match=false; }
                        if (match) keep.push(j);
                    }
                    if (keep.length > 0) {
                        matchCount += keep.length; visibleIndices.push(i); indices.push(i);
                        update.x.push(keep.map(function(idx){ return orig.x[idx]; }));
                        update.y.push(keep.map(function(idx){ return orig.y[idx]; }));
                        update.z.push(keep.map(function(idx){ return orig.z[idx]; }));
                        update.text.push(keep.map(function(idx){ return orig.text[idx]; }));
                        update.hovertext.push(keep.map(function(idx){ return orig.hovertext[idx]; }));
                        update.customdata.push(keep.map(function(idx){ return orig.customdata[idx]; }));
                        update.showlegend.push(true); // Force show legend if match found
                        update.hoverinfo.push('text');
                    } else { hiddenIndices.push(i); }
                } else if (isStatic) { visibleIndices.push(i); } else if (isOverview) { hiddenIndices.push(i); } // Hide overview in search
            }
            if (indices.length > 0) {
                Plotly.restyle(el.id, update, indices);
            }
            if (visibleIndices.length > 0) Plotly.restyle(el.id, {visible: true}, visibleIndices);
            if (hiddenIndices.length > 0) Plotly.restyle(el.id, {visible: false}, hiddenIndices);
            
            // Handle Legend (simplified)
            // ...
            
            statusDiv.innerText = 'Found ' + matchCount + ' points.'; setTimeout(updateCount, 200);
        }

        function clearSearch() {
            document.getElementById('search_signal').value = ''; document.getElementById('search_gene').value = ''; document.getElementById('search_snp').value = '';
            if (!isDataBackedUp) return;
            var statusDiv = document.getElementById('search_status'); statusDiv.innerText = 'Restoring...';
            
            // Reset State
            state.searchActive = false;
            
            // Restore Data
            var update = { x: [], y: [], z: [], text: [], hovertext: [], customdata: [], showlegend: [], hoverinfo: [] };
            var indices = []; 
            
            for (var i = 0; i < plot.data.length; i++) {
                var trace = plot.data[i]; var orig = originalData[i];
                if (orig) {
                    indices.push(i); 
                    update.x.push(orig.x); update.y.push(orig.y); update.z.push(orig.z); 
                    update.text.push(orig.text); update.hovertext.push(orig.hovertext); 
                    update.customdata.push(orig.customdata); update.showlegend.push(orig.showlegend); 
                    update.hoverinfo.push('text');
                }
            }
            Plotly.restyle(el.id, update, indices).then(function() {
                statusDiv.innerText = 'Cleared.'; 
                // Trigger full visualization update (restore visibility based on Cluster/Filter)
                updateVisualization();
                setTimeout(updateCount, 200);
            });
        }
        document.getElementById('btn_search').onclick = performSearch; document.getElementById('btn_clear').onclick = clearSearch;
        function handleEnter(e) { if (e.key === 'Enter') performSearch(); }
        document.getElementById('search_signal').addEventListener('keypress', handleEnter); document.getElementById('search_gene').addEventListener('keypress', handleEnter); document.getElementById('search_snp').addEventListener('keypress', handleEnter);
        
        function updateCount() {
            var total = 0;
            // Check if PCA Mode is active
            if (state.pcaMode) {
                 // In PCA Mode, count visible Genes (bubbles)
                 for (var i = 0; i < plot.data.length; i++) {
                    var trace = plot.data[i]; 
                    // Use trace.visible explicitly (can be true or undefined)
                    var isVisible = trace.visible === true || trace.visible === undefined;
                    // Count points in 'pca_gene_bubble' traces
                    if (trace.meta && trace.meta.type === 'pca_gene_bubble' && isVisible) { 
                        if (trace.x) total += trace.x.length; 
                    }
                 }
                 var label = '<b>Visible Genes:</b> ';
                 if (plot.layout.annotations && plot.layout.annotations.length > 0) {
                    var newText = label + total;
                    // Always update if content changed
                    if (plot.layout.annotations[0].text !== newText) {
                        Plotly.relayout(el.id, {'annotations[0].text': newText});
                    }
                 }
            } else {
                 // Normal Mode: Count SNPs (scatter points)
                 for (var i = 0; i < plot.data.length; i++) {
                    var trace = plot.data[i]; 
                    var isVisible = trace.visible === true || trace.visible === undefined;
                    // Count points in 'detail' or 'overview' traces (exclude static/pca lines/bubbles)
                    if (trace.type === 'scatter3d' && isVisible && trace.mode && trace.mode.indexOf('markers') !== -1) {
                        // Exclude PCA Bubbles in Normal Mode count (just in case)
                        if (trace.meta && trace.meta.type === 'pca_gene_bubble') continue;
                        if (trace.x) total += trace.x.length; 
                    }
                 }
                 var label = '<b>Visible SNPs:</b> ';
                 if (plot.layout.annotations && plot.layout.annotations.length > 0) {
                    var newText = label + total;
                    if (plot.layout.annotations[0].text !== newText) {
                        Plotly.relayout(el.id, {'annotations[0].text': newText});
                    }
                 }
            }
        }
        
        // Initial Update (Delay to ensure Plotly is ready)
        setTimeout(function() {
            updateVisualization();
            // updateCount will be called inside updateVisualization after restyle
        }, 500);
      }
    "
    
    fig <- htmlwidgets::onRender(fig, js_code)
    html_file <- file.path(plots_interactive_dir, "3D_Scatter_Plot_Interactive_Complete.html")
    tryCatch({ saveWidget(fig, html_file, selfcontained = TRUE); cat(sprintf("Interactive plot saved: %s\n", html_file)) }, 
             error = function(e) { message("Pandoc error, trying non-self-contained..."); saveWidget(fig, html_file, selfcontained = FALSE); cat(sprintf("Interactive plot saved (dependencies): %s\n", html_file)) })
    
  }, error = function(e) { message("Error generating interactive plot: ", e$message) })

  # --- Save Data and Readme ---
  write_csv(data_subset, file.path(data_dir, "merged_intersection_3d_data.csv"))
  readme_content <- c(
    "# 3D Scatter Plot Analysis Report / 3D",
    "",
    paste("Version:", version_name),
    paste("Date / :", Sys.Date),
    paste("Time / :", format(Sys.time, "%H:%M:%S")),
    "",
    "## Overview / ",
    "This analysis visualizes the intersection signals (cell-gene-snp) in a 3D space based on their effect sizes (b values) for Aging, RA, and HZ.",
    "Aging、RAHZ(b), 3D(cell-gene-snp). ",
    "",
    "## Input Data / ",
    paste("- Total Signals:", nrow(data_subset)),
    paste("- Unique Cells:", n_distinct(data_subset$cell)),
    paste("- Unique Genes:", n_distinct(data_subset$gene)),
    "",
    "## Statistics / ",
    "Summary of 'b' values:",
    "Aging:", summary(data_subset$b_aging),
    "RA:", summary(data_subset$b_RA),
    "HZ:", summary(data_subset$b_HZ)
  )
  writeLines(as.character(unlist(readme_content)), file.path(readme_dir, "analysis_report.md"))
  
  cat(sprintf("[%s] Visualization for %s completed.\n", Sys.time(), version_name))
}

# ------------------------------------------------------------------------------
# 5. Execution / 
# ------------------------------------------------------------------------------

# 1. No-Bulk Version (Main Analysis)
# Since we already filtered out bulk data during loading, 'all_data' is the target dataset.
cat(sprintf("\nProcessing No-Bulk Version (Main Analysis) with %d rows.\n", nrow(all_data)))

run_visualization(all_data, "Version_2_No_Bulk_Optimized", FALSE)

cat(sprintf("\n[%s] All analyses completed successfully.\n", Sys.time()))
sink()
