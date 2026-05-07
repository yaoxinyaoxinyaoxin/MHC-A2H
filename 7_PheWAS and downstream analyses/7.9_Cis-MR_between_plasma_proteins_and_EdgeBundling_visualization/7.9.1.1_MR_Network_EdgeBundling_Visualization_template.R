#!/usr/bin/env Rscript
# ==============================================================================
# [Script]: 7.9.1.1_MR_Network_EdgeBundling_Visualization_template.R
# [Method]: Network Analysis & Hierarchical Edge Bundling
# [Step]: Cis-MR between plasma proteins (Validation & Visualization)
# 
# [Function]:
# A generalized template to construct and visualize causal networks derived from 
#       cis-MR results. It employs the Louvain algorithm for community detection 
#       generates both static high-resolution hierarchical edge bundling plots 
#       and interactive 3D Plotly-based HTML visualizations.
#
# [Data Availability]:
# Input requires a standard CSV output from TwoSampleMR with FDR correction.
# ==============================================================================

# 1. Environment & Dependencies
# ------------------------------------------------------------------------------
suppressPackageStartupMessages({
  library(optparse)
  library(dplyr)
  library(readr)
  library(data.table)
  library(openxlsx)
  library(igraph)
  library(ggraph)
  library(ggplot2)
  library(ragg)
  library(jsonlite)
  library(RColorBrewer)
})

# 2. Command Line Arguments
# ------------------------------------------------------------------------------
option_list <- list(
  make_option(c("--input_csv"), type="character", default=NULL,
              help="Path to significant cis-MR results CSV (must contain FDR, method, b, se, exposure, outcome)", metavar="character"),
  make_option(c("--out_dir"), type="character", default="./network_visualization",
              help="Output directory path [default: %default]", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$input_csv)) {
  print_help(opt_parser)
  stop("Missing required argument. Please provide --input_csv.")
}

# 3. Initialization & Directories
# ------------------------------------------------------------------------------
out_dir <- opt$out_dir
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- file.path(out_dir, paste0("Visualization_EdgeBundling_Optimized_", timestamp))

readme_dir <- file.path(output_dir, "readme")
logs_dir <- file.path(output_dir, "logs")
data_dir <- file.path(output_dir, "network_data")
plots_dir <- file.path(output_dir, "plots")
interactive_dir <- file.path(plots_dir, "3D_Interactive")
lib_dir <- file.path(interactive_dir, "libs")

invisible(sapply(c(readme_dir, logs_dir, data_dir, plots_dir, interactive_dir, lib_dir), dir.create, recursive = TRUE, showWarnings = FALSE))

# Logging
log_file <- file.path(logs_dir, paste0("run_log_", timestamp, ".txt"))
sink(log_file, append = TRUE, split = TRUE)

cat("====================================================\n")
cat("Starting MR Network Edge Bundling Visualization\n")
cat("Input File:", opt$input_csv, "\n")
cat("Output Dir:", output_dir, "\n")
cat("====================================================\n")

# Session Info Logging for Reproducibility
sink(file.path(logs_dir, "session_info.txt"))
print(sessionInfo())
sink()

# NPG (Nature Publishing Group) Style Colors
color_pos <- "#E64B35FF" # Muted Red
color_neg <- "#4DBBD5FF" # Muted Blue

# 4. Data Processing
# ------------------------------------------------------------------------------
cat("[INFO] Reading input data...\n")
if (!file.exists(opt$input_csv)) stop("[ERROR] Input file not found.")

df_all <- fread(opt$input_csv)
cat("[INFO] Loaded", nrow(df_all), "rows.\n")

# Filter for main MR methods and FDR < 0.05
df_main <- df_all %>%
  filter(method %in% c("Inverse variance weighted", "Wald ratio"), FDR < 0.05)

cat("[INFO] Extracted IVW and Wald ratio results: ", nrow(df_main), " significant causal pairs.\n")

if(nrow(df_main) == 0) {
  stop("[ERROR] No significant results found after filtering.")
}

# 5. Build Network Data Structures & Clustering
# ------------------------------------------------------------------------------
cat("[INFO] Building nodes and edges, applying Louvain clustering...\n")

edges <- df_main %>%
  select(from = exposure, to = outcome, method, b, se, pval, FDR) %>%
  mutate(
    direction = ifelse(b > 0, "Positive", "Negative"),
    weight = abs(b),
    logFDR = -log10(FDR),
    color = ifelse(b > 0, color_pos, color_neg)
  ) %>%
  group_by(from, to) %>%
  arrange(FDR) %>%
  slice(1) %>%
  ungroup()

unique_proteins <- unique(c(edges$from, edges$to))
nodes <- data.frame(id = unique_proteins, label = unique_proteins, stringsAsFactors = FALSE)

g <- graph_from_data_frame(d = edges, vertices = nodes, directed = TRUE)
g_undirected <- as_undirected(g, mode = "collapse", edge.attr.comb = list(weight = "mean", "ignore"))

# Set seed for reproducibility of stochastic community detection
set.seed(123) 
cluster_res <- cluster_louvain(g_undirected)

mod_score <- modularity(cluster_res)
cat(sprintf("[INFO] Louvain Clustering Modularity Score: %.4f\n", mod_score))

V(g)$degree_total <- igraph::degree(g, mode = "all")
V(g)$degree_in <- igraph::degree(g, mode = "in")
V(g)$degree_out <- igraph::degree(g, mode = "out")
V(g)$betweenness <- igraph::betweenness(g, directed = TRUE)
V(g)$Cluster <- paste0("Cluster_", membership(cluster_res)[V(g)$name])

nodes <- nodes %>%
  mutate(
    Cluster = V(g)$Cluster,
    degree_total = V(g)$degree_total,
    degree_in = V(g)$degree_in,
    degree_out = V(g)$degree_out,
    betweenness = V(g)$betweenness,
    title = paste0("<p><b>", label, "</b><br>Cluster: ", Cluster,
                   "<br>Total Degree: ", degree_total, 
                   "<br>Out (Exposure): ", degree_out, "<br>In (Outcome): ", degree_in, "</p>"),
    value = degree_out + 1 
  )

fwrite(edges, file.path(data_dir, "Network_Edges.csv"))
fwrite(nodes, file.path(data_dir, "Network_Nodes.csv"))
write.xlsx(list(Nodes = nodes, Edges = edges), file.path(data_dir, "Network_Data_Combined.xlsx"))

unique_clusters <- sort(unique(nodes$Cluster))
nb_clusters <- length(unique_clusters)

distinct_palette <- c(
  "#FF7F0E", "#2CA02C", "#9467BD", "#8C564B", "#00A087", 
  "#87CEFA", "#BCBD22", "#17BECF", "#FFBB78", "#98DF8A", 
  "#C5B0D5", "#C49C94", "#F7B6D2", "#DBDB8D", "#9EDAE5",
  "#393B79", "#5254A3", "#6B6ECF", "#9C9EDE", "#637939"
)

if (nb_clusters <= length(distinct_palette)) {
  cluster_colors <- distinct_palette[1:nb_clusters]
} else {
  cluster_colors <- colorRampPalette(distinct_palette)(nb_clusters)
}
names(cluster_colors) <- unique_clusters

# 6. Hierarchical Edge Bundling Plot (Static)
# ------------------------------------------------------------------------------
cat("[INFO] Generating static Hierarchical Edge Bundling plots...\n")

gene_cluster_info <- nodes %>%
  select(gene = id, cluster = Cluster, degree = degree_total) %>%
  arrange(cluster)

hierarchy_edges <- rbind(
  data.frame(from = "Origin", to = unique(gene_cluster_info$cluster)),
  data.frame(from = gene_cluster_info$cluster, to = gene_cluster_info$gene)
)

hierarchy_graph <- graph_from_data_frame(hierarchy_edges)
V(hierarchy_graph)$degree <- 0
V(hierarchy_graph)$is_leaf <- FALSE
V(hierarchy_graph)$cluster <- NA

leaf_indices <- match(gene_cluster_info$gene, V(hierarchy_graph)$name)
V(hierarchy_graph)$degree[leaf_indices] <- gene_cluster_info$degree
V(hierarchy_graph)$is_leaf[leaf_indices] <- TRUE
V(hierarchy_graph)$cluster[leaf_indices] <- as.character(gene_cluster_info$cluster)

from_id <- match(edges$from, V(hierarchy_graph)$name)
to_id <- match(edges$to, V(hierarchy_graph)$name)
valid_edges <- !is.na(from_id) & !is.na(to_id)

layout_obj <- create_layout(hierarchy_graph, layout = 'dendrogram', circular = TRUE)

connect_df <- data.frame(
  from = from_id[valid_edges], 
  to = to_id[valid_edges],
  direction = edges$direction[valid_edges],
  edge_weight = edges$weight[valid_edges],
  source_cluster = layout_obj$cluster[from_id[valid_edges]]
)

# Plot 1: Edge Color = Direction
p_bundling <- ggraph(layout_obj) +
  geom_conn_bundle(data = get_con(from = connect_df$from, to = connect_df$to, 
                                  direction = connect_df$direction, edge_weight = connect_df$edge_weight), 
                   aes(color = direction, edge_width = edge_weight), 
                   alpha = 0.65, tension = 0.85) +
  scale_edge_color_manual(values = c("Positive" = color_pos, "Negative" = color_neg)) +
  scale_edge_width_continuous(range = c(0.1, 1.2)) +
  geom_node_point(aes(filter = is_leaf, x = x*1.05, y = y*1.05, color = cluster, size = degree), alpha = 0.65) +
  scale_color_manual(values = cluster_colors) +
  scale_size_continuous(range = c(1, 6)) +
  geom_node_text(aes(x = x*1.1, y = y*1.1, filter = is_leaf, 
                     label = name, angle = -((-node_angle(x, y) + 90) %% 180) + 90, hjust = 'outward'),
                 size = 2.2, fontface = "bold", family = "sans") +
  coord_fixed() +
  theme_void(base_family = "sans") +
  theme(
    legend.position = "right",
    plot.margin = margin(60, 60, 60, 60), 
    text = element_text(family = "sans")
  ) +
  labs(edge_color = "Effect Direction", edge_width = "Effect Size (|b|)", color = "Protein Cluster", size = "Degree")

ggsave(file.path(plots_dir, "1_Hierarchical_EdgeBundling.png"), plot = p_bundling, device = ragg::agg_png, width = 12, height = 12, dpi = 600)
ggsave(file.path(plots_dir, "1_Hierarchical_EdgeBundling.pdf"), plot = p_bundling, width = 12, height = 12)

# Plot 2: Edge Color = Source Cluster
p_bundling_source_color <- ggraph(layout_obj) +
  geom_conn_bundle(data = get_con(from = connect_df$from, to = connect_df$to, 
                                  source_cluster = connect_df$source_cluster, edge_weight = connect_df$edge_weight), 
                   aes(color = source_cluster, edge_width = edge_weight), 
                   alpha = 0.65, tension = 0.85) +
  scale_edge_color_manual(values = cluster_colors) +
  scale_edge_width_continuous(range = c(0.1, 1.2)) +
  geom_node_point(aes(filter = is_leaf, x = x*1.05, y = y*1.05, color = cluster, size = degree), alpha = 0.65) +
  scale_color_manual(values = cluster_colors) +
  scale_size_continuous(range = c(1, 6)) +
  geom_node_text(aes(x = x*1.1, y = y*1.1, filter = is_leaf, 
                     label = name, angle = -((-node_angle(x, y) + 90) %% 180) + 90, hjust = 'outward'),
                 size = 2.2, fontface = "bold", family = "sans") +
  coord_fixed() +
  theme_void(base_family = "sans") +
  theme(
    legend.position = "right",
    plot.margin = margin(60, 60, 60, 60), 
    text = element_text(family = "sans")
  ) +
  labs(edge_color = "Source Protein Cluster", edge_width = "Effect Size (|b|)", color = "Protein Cluster", size = "Degree")

ggsave(file.path(plots_dir, "1_Hierarchical_EdgeBundling_SourceColor.png"), plot = p_bundling_source_color, device = ragg::agg_png, width = 12, height = 12, dpi = 600)
ggsave(file.path(plots_dir, "1_Hierarchical_EdgeBundling_SourceColor.pdf"), plot = p_bundling_source_color, width = 12, height = 12)


# 7. Interactive Visualization (Plotly 3D HTML)
# ------------------------------------------------------------------------------
cat("[INFO] Generating interactive 3D network plot (HTML/JS)...\n")

leaf_coords <- layout_obj %>% 
  filter(is_leaf == TRUE) %>% 
  select(name, x, y)

scale_factor <- 300 
export_nodes <- nodes %>%
  left_join(leaf_coords, by = c("id" = "name")) %>%
  mutate(
    x = x * scale_factor, 
    y = -y * scale_factor, 
    color = cluster_colors[Cluster]
  ) %>%
  select(id, label, Cluster, color, x, y)

export_edges <- edges %>%
  select(source = from, target = to, b, se, FDR, method, direction) %>%
  mutate(
    OR = exp(b),
    LCI = exp(b - 1.96 * se),
    UCI = exp(b + 1.96 * se),
    color = ifelse(direction == "Positive", 
                   paste0("rgba(", col2rgb(color_pos)[1], ",", col2rgb(color_pos)[2], ",", col2rgb(color_pos)[3], ",0.85)"), 
                   paste0("rgba(", col2rgb(color_neg)[1], ",", col2rgb(color_neg)[2], ",", col2rgb(color_neg)[3], ",0.85)")),
    width = abs(b) * 4
  )

json_nodes <- toJSON(export_nodes, auto_unbox = TRUE)
json_edges <- toJSON(export_edges, auto_unbox = TRUE)

plotly_url <- "https://cdn.plot.ly/plotly-2.32.0.min.js"
jquery_url <- "https://code.jquery.com/jquery-3.6.0.min.js"
select2_js_url <- "https://cdn.jsdelivr.net/npm/select2@4.1.0-rc.0/dist/js/select2.min.js"
select2_css_url <- "https://cdn.jsdelivr.net/npm/select2@4.1.0-rc.0/dist/css/select2.min.css"

safe_download <- function(url, dest) {
  if(!file.exists(dest)) {
    tryCatch({
      download.file(url, dest, quiet = TRUE, method = "libcurl")
    }, error = function(e) {
      system(paste("curl -s -L", shQuote(url), "-o", shQuote(dest)))
    })
  }
}
suppressWarnings({
  safe_download(plotly_url, file.path(lib_dir, "plotly.min.js"))
  safe_download(jquery_url, file.path(lib_dir, "jquery.min.js"))
  safe_download(select2_js_url, file.path(lib_dir, "select2.min.js"))
  safe_download(select2_css_url, file.path(lib_dir, "select2.min.css"))
})

html_content <- paste0('
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="UTF-8">
  <title>Interactive cis-MR Causal Network (3D)</title>
  <style>
    body { margin: 0; padding: 0; overflow: hidden; background-color: #ffffff; font-family: Arial, sans-serif; }
    #graph { width: 100vw; height: 100vh; }
    
    /* Left Search Panel */
    #search-container {
      position: absolute; top: 20px; left: 20px;
      background: rgba(255, 255, 255, 0.95);
      padding: 15px 20px; border-radius: 8px;
      box-shadow: 0 4px 10px rgba(0,0,0,0.15); z-index: 10;
      width: 280px;
    }
    #search-header {
      margin: 0 0 10px 0; font-size: 14px; color: #333; border-bottom: 1px solid #eee; padding-bottom: 5px;
      cursor: move; font-weight: bold; user-select: none;
    }
    .form-group { margin-bottom: 10px; }
    .form-group label { display: block; font-size: 11px; margin-bottom: 4px; color: #555; font-weight: bold; }
    .select2-container { font-size: 12px !important; }
    .btn-row { display: flex; gap: 8px; margin-top: 15px; }
    button { flex: 1; padding: 6px; font-size: 12px; border-radius: 4px; cursor: pointer; border: none; font-weight: bold; color: white; background-color: #007bff; }
    button:hover { background-color: #0056b3; }
    #btn_clear { background-color: #6c757d; }
    #btn_clear:hover { background-color: #5a6268; }
    
    /* Right Sidebar Toggle Button */
    #sidebar-toggle {
      position: absolute; top: 20px; right: 20px; z-index: 20;
      padding: 10px 15px; background: #00A087; box-shadow: 0 4px 10px rgba(0,0,0,0.15);
      width: auto; flex: none;
    }
    #sidebar-toggle:hover { background: #008a74; }
    
    /* Right Edge Data Sidebar */
    #edge-sidebar {
      position: absolute; top: 0; right: -400px; width: 380px; height: 100vh;
      background: rgba(255,255,255,0.98); box-shadow: -2px 0 15px rgba(0,0,0,0.15);
      z-index: 25; transition: right 0.3s ease; display: flex; flex-direction: column;
    }
    #edge-sidebar.open { right: 0; }
    #sidebar-header-box { padding: 15px 20px; border-bottom: 1px solid #ddd; background: #f8f9fa; }
    #sidebar-title { font-size: 15px; font-weight: bold; color: #333; margin: 0; display: inline-block; }
    #close-sidebar { float: right; background: none; border: none; color: #888; font-size: 20px; cursor: pointer; padding: 0; line-height: 1; margin-top: -2px; }
    #close-sidebar:hover { color: #333; }
    #edge-search { width: 100%; padding: 8px; box-sizing: border-box; margin-top: 12px; border: 1px solid #ccc; border-radius: 4px; font-size: 12px; }
    #clear-edge-selection { margin-top: 10px; background: #dc3545; width: 100%; flex: none; }
    #clear-edge-selection:hover { background: #c82333; }
    
    /* Edge List */
    #edge-list { flex: 1; overflow-y: auto; padding: 0; margin: 0; }
    .edge-item {
      font-size: 11px; padding: 12px 20px; border-bottom: 1px solid #eee; display: flex; align-items: flex-start; gap: 10px;
    }
    .edge-item:hover { background: #f1f1f1; }
    .edge-cb { margin-top: 2px; cursor: pointer; }
    .edge-info { flex: 1; line-height: 1.5; color: #444; }
    .edge-info b { color: #111; }
  </style>
  <script src="libs/plotly.min.js"></script>
  <link href="libs/select2.min.css" rel="stylesheet" />
  <script src="libs/jquery.min.js"></script>
  <script src="libs/select2.min.js"></script>
</head>
<body>
  <div id="search-container">
    <div id="search-header">Network Filter</div>
    <div class="form-group">
      <label>Cluster:</label>
      <select id="filter_cluster" class="multi-select" multiple="multiple"></select>
    </div>
    <div class="form-group">
      <label>Exposure Protein:</label>
      <select id="filter_exposure" class="multi-select" multiple="multiple"></select>
    </div>
    <div class="form-group">
      <label>Outcome Protein:</label>
      <select id="filter_outcome" class="multi-select" multiple="multiple"></select>
    </div>
    <div class="btn-row">
      <button id="btn_search">Search</button>
      <button id="btn_clear">Clear / Reset</button>
    </div>
  </div>

  <button id="sidebar-toggle">View Edge Data</button>
  <div id="edge-sidebar">
    <div id="sidebar-header-box">
      <h3 id="sidebar-title">Edge Summary Data</h3>
      <button id="close-sidebar">&times;</button>
      <input type="text" id="edge-search" placeholder="Search exposure, outcome, method...">
      <button id="clear-edge-selection">Clear Selection</button>
    </div>
    <div id="edge-list"></div>
  </div>

  <div id="graph"></div>
  <script>
    const nodes = ', json_nodes, ';
    const edges = ', json_edges, ';
    
    let activeNodeId = null;
    let activeClusters = [];
    let activeExposures = [];
    let activeOutcomes = [];
    let activeEdgeIndices = [];

    // Initialize Draggable for Left Panel
    function makeDraggable(elmnt, headerId) {
      var pos1 = 0, pos2 = 0, pos3 = 0, pos4 = 0;
      var header = document.getElementById(headerId);
      if (header) { header.onmousedown = dragMouseDown; } 
      else { elmnt.onmousedown = dragMouseDown; }

      function dragMouseDown(e) {
        if (["INPUT", "SELECT", "BUTTON", "SPAN"].includes(e.target.tagName)) return;
        e = e || window.event; e.preventDefault();
        pos3 = e.clientX; pos4 = e.clientY;
        document.onmouseup = closeDragElement;
        document.onmousemove = elementDrag;
      }
      function elementDrag(e) {
        e = e || window.event; e.preventDefault();
        pos1 = pos3 - e.clientX; pos2 = pos4 - e.clientY;
        pos3 = e.clientX; pos4 = e.clientY;
        elmnt.style.top = (elmnt.offsetTop - pos2) + "px";
        elmnt.style.left = (elmnt.offsetLeft - pos1) + "px";
        elmnt.style.bottom = "auto"; elmnt.style.right = "auto";
      }
      function closeDragElement() {
        document.onmouseup = null; document.onmousemove = null;
      }
    }
    makeDraggable(document.getElementById("search-container"), "search-header");

    function populateEdgeList() {
        const list = document.getElementById("edge-list");
        let html = "";
        edges.forEach((e, idx) => {
            const sourceNode = nodes.find(n => n.id === e.source);
            const targetNode = nodes.find(n => n.id === e.target);
            if(!sourceNode || !targetNode) return;
            const sLabel = sourceNode.label;
            const tLabel = targetNode.label;
            const searchStr = `${sLabel} ${tLabel} ${e.method}`.toLowerCase();
            html += `
            <div class="edge-item" data-search="${searchStr}">
                <input type="checkbox" class="edge-cb" value="${idx}">
                <div class="edge-info">
                    <b>Exposure:</b> ${sLabel}<br>
                    <b>Outcome:</b> ${tLabel}<br>
                    <b>Method:</b> ${e.method} | <b>FDR:</b> ${e.FDR.toExponential(2)}<br>
                    <b>OR (95% CI):</b> ${e.OR.toFixed(2)} (${e.LCI.toFixed(2)} - ${e.UCI.toFixed(2)})
                </div>
            </div>
            `;
        });
        list.innerHTML = html;

        $(".edge-cb").change(function() {
            activeEdgeIndices = [];
            $(".edge-cb:checked").each(function() {
                activeEdgeIndices.push(parseInt($(this).val()));
            });
            renderGraph(true);
        });
    }

    $(document).ready(function() {
        const clusterSelect = $("#filter_cluster");
        const expSelect = $("#filter_exposure");
        const outSelect = $("#filter_outcome");

        const clusters = [...new Set(nodes.map(n => n.Cluster))].sort();
        clusters.forEach(c => clusterSelect.append(new Option(c, c)));

        const exposures = [...new Set(edges.map(e => {
            const n = nodes.find(x => x.id === e.source); return n ? n.label : null;
        }).filter(Boolean))].sort();
        exposures.forEach(exp => expSelect.append(new Option(exp, exp)));

        const outcomes = [...new Set(edges.map(e => {
            const n = nodes.find(x => x.id === e.target); return n ? n.label : null;
        }).filter(Boolean))].sort();
        outcomes.forEach(out => outSelect.append(new Option(out, out)));

        $(".multi-select").select2({ width: "100%", placeholder: "Select...", allowClear: true });

        // Sidebar interactions
        populateEdgeList();
        $("#sidebar-toggle").click(() => $("#edge-sidebar").addClass("open"));
        $("#close-sidebar").click(() => $("#edge-sidebar").removeClass("open"));
        
        $("#edge-search").on("keyup", function() {
            const val = $(this).val().toLowerCase();
            $(".edge-item").each(function() {
                const text = $(this).attr("data-search");
                $(this).toggle(text.indexOf(val) > -1);
            });
        });

        $("#clear-edge-selection").click(() => {
            $(".edge-cb").prop("checked", false);
            activeEdgeIndices = [];
            renderGraph();
        });
    });

    function getHexColor(hex) {
        if (!hex) return "rgba(100,100,100,0.65)";
        hex = hex.replace("#", "");
        if (hex.length === 3) hex = hex.split("").map(c => c + c).join("");
        let r = parseInt(hex.substring(0, 2), 16) || 0;
        let g = parseInt(hex.substring(2, 4), 16) || 0;
        let b = parseInt(hex.substring(4, 6), 16) || 0;
        return `rgba(${r}, ${g}, ${b}, 0.65)`;
    }

    function checkEdgeActive(e, idx) {
        const hasSidebarFilter = activeEdgeIndices.length > 0;
        const hasLeftFilter = activeClusters.length > 0 || activeExposures.length > 0 || activeOutcomes.length > 0;
        const isAnyFilterActive = activeNodeId || hasLeftFilter || hasSidebarFilter;

        if (!isAnyFilterActive) return true;
        
        // Strict sidebar override
        if (hasSidebarFilter) return activeEdgeIndices.includes(idx);
        
        if (activeNodeId) return (e.source === activeNodeId || e.target === activeNodeId);
        
        if (hasLeftFilter) {
            const s = nodes.find(x => x.id === e.source);
            const t = nodes.find(x => x.id === e.target);
            if (activeClusters.includes(s.Cluster)) return true;
            if (activeExposures.includes(s.label)) return true;
            if (activeOutcomes.includes(t.label)) return true;
            return false;
        }
        return true;
    }

    
    function updateEdgeSidebar() {
        const list = $("#edge-list");
        const items = list.children(".edge-item").get();
        
        const hasLeftFilter = activeClusters.length > 0 || activeExposures.length > 0 || activeOutcomes.length > 0;
        const isAnyFilterActive = activeNodeId || hasLeftFilter;

        items.sort((a, b) => {
            const idxA = parseInt($(a).find(".edge-cb").val());
            const idxB = parseInt($(b).find(".edge-cb").val());
            
            if (!isAnyFilterActive) {
                return idxA - idxB; // default order
            }
            
            const tempSidebar = activeEdgeIndices;
            activeEdgeIndices = [];
            const activeA = checkEdgeActive(edges[idxA], idxA);
            const activeB = checkEdgeActive(edges[idxB], idxB);
            activeEdgeIndices = tempSidebar;
            
            if (activeA && !activeB) return -1;
            if (!activeA && activeB) return 1;
            return idxA - idxB;
        });
        
        $.each(items, function(idx, itm) {
            list.append(itm);
            const edgeIdx = parseInt($(itm).find(".edge-cb").val());
            
            const tempSidebar = activeEdgeIndices;
            activeEdgeIndices = [];
            const isActive = !isAnyFilterActive || checkEdgeActive(edges[edgeIdx], edgeIdx);
            activeEdgeIndices = tempSidebar;
            
            $(itm).css("opacity", isActive ? "1" : "0.3");
            $(itm).css("background-color", isActive ? "" : "#f9f9f9");
        });
        
        list.scrollTop(0);
    }

    function renderGraph(fromSidebar = false) {
        const traces = [];
        const clusters = [...new Set(nodes.map(n => n.Cluster))].sort();
        const clusterSettings = {};
        clusters.forEach((c, idx) => {
            const sign = idx % 2 === 0 ? 1 : -1;
            const height = 150 + Math.floor(idx / 2) * 120;
            clusterSettings[c] = { zOffset: sign * height };
        });

        const hasSidebarFilter = activeEdgeIndices.length > 0;
        const hasLeftFilter = activeClusters.length > 0 || activeExposures.length > 0 || activeOutcomes.length > 0;
        const isAnyFilterActive = activeNodeId || hasLeftFilter || hasSidebarFilter;

        clusters.forEach(clusterName => {
            const cx=[], cy=[], cz=[], ctext=[], csize=[], chover=[], ccolor=[];
            let baseColor = "#333333";

            // 1. Process Edges
            const ex=[], ey=[], ez=[];
            const emid_x=[], emid_y=[], emid_z=[], emid_text=[];

            edges.forEach((e, idx) => {
                const source = nodes.find(n => n.id === e.source);
                const target = nodes.find(n => n.id === e.target);
                if (source && source.Cluster === clusterName && target) {
                    if (checkEdgeActive(e, idx)) {
                        const numPoints = 30;
                        const zOffset = clusterSettings[clusterName].zOffset;
                        for(let i=0; i<=numPoints; i++) {
                            let t = i / numPoints;
                            let x = Math.pow(1-t, 2)*source.x + Math.pow(t, 2)*target.x;
                            let y = Math.pow(1-t, 2)*source.y + Math.pow(t, 2)*target.y;
                            let z = 2*(1-t)*t*zOffset;

                            ex.push(x); ey.push(y); ez.push(z);
                            if (i === 8 || i === 15 || i === 22) {
                                emid_x.push(x); emid_y.push(y); emid_z.push(z);
                                emid_text.push(`<b>Exposure:</b> ${source.label}<br><b>Outcome:</b> ${target.label}<br><b>Method:</b> ${e.method}<br><b>OR (95% CI):</b> ${e.OR.toFixed(2)} (${e.LCI.toFixed(2)} - ${e.UCI.toFixed(2)})<br><b>FDR:</b> ${e.FDR.toExponential(3)}`);
                            }
                        }
                        ex.push(null); ey.push(null); ez.push(null);
                    }
                }
            });

            // 2. Process Nodes
            nodes.forEach(n => {
                if (n.Cluster === clusterName) {
                    baseColor = n.color;
                    let isActive = false;
                    
                    if (!isAnyFilterActive) {
                        isActive = true;
                    } else if (hasSidebarFilter || activeNodeId) {
                        isActive = edges.some((e, idx) => checkEdgeActive(e, idx) && (e.source === n.id || e.target === n.id));
                    } else if (hasLeftFilter) {
                        let selfMatch = false;
                        if (activeClusters.includes(n.Cluster)) selfMatch = true;
                        if (activeExposures.includes(n.label)) selfMatch = true;
                        if (activeOutcomes.includes(n.label)) selfMatch = true;
                        
                        if (selfMatch) isActive = true;
                        else isActive = edges.some((e, idx) => checkEdgeActive(e, idx) && (e.source === n.id || e.target === n.id));
                    }

                    if (isActive) {
                        cx.push(n.x); cy.push(n.y); cz.push(0);
                        ctext.push(n.label);
                        csize.push(n.degree * 0.8 + 12); 
                        chover.push(`<b>${n.label}</b><br>Cluster: ${n.Cluster}<br>Degree: ${n.degree}`);
                        ccolor.push(n.color);
                    }
                }
            });

            // 3. Push Edges Trace
            if (ex.length > 0 || !isAnyFilterActive) {
                const lineWidth = isAnyFilterActive ? 5 : 2;
                const edgeColor = getHexColor(baseColor);
                traces.push({
                    type: "scatter3d", mode: "lines",
                    x: ex, y: ey, z: ez,
                    line: {color: edgeColor, width: lineWidth},
                    hoverinfo: "none",
                    showlegend: true,
                    name: clusterName
                });

                traces.push({
                    type: "scatter3d", mode: "markers",
                    x: emid_x, y: emid_y, z: emid_z, text: emid_text,
                    hoverinfo: "text", marker: {size: 45, color: "rgba(0,0,0,0)"},
                    showlegend: false
                });
            }

            // 4. Push Nodes Trace
            if (cx.length > 0) {
                traces.push({
                    type: "scatter3d", mode: "markers+text",
                    x: cx, y: cy, z: cz, text: ctext,
                    textposition: "top center", hovertext: chover, hoverinfo: "text",
                    marker: {
                        size: csize, color: ccolor, 
                        opacity: 0.65, sizemode: "diameter"
                    },
                    textfont: {size: 11, color: "#333333"},
                    showlegend: false,
                    name: clusterName + " Nodes"
                });
            }
        });

        const layout = {
            uirevision: "true",
            scene: {
                xaxis: {showgrid: false, zeroline: false, showticklabels: false, title: ""},
                yaxis: {showgrid: false, zeroline: false, showticklabels: false, title: ""},
                zaxis: {showgrid: false, zeroline: false, showticklabels: false, title: ""},
                aspectmode: "cube",
                camera: { eye: {x: 0, y: -0.8, z: 0} }
            },
            paper_bgcolor: "#ffffff",
            plot_bgcolor: "#ffffff",
            margin: {l: 0, r: 0, b: 0, t: 0},
            showlegend: true,
            legend: {x: 0.02, y: 0.2, font: {size: 12}, bgcolor: "rgba(255,255,255,0.8)", bordercolor: "#ccc", borderwidth: 1},
            hovermode: "closest"
        };

        Plotly.react("graph", traces, layout);
        if (!fromSidebar) {
            updateEdgeSidebar();
        }
    }

    renderGraph();

    $("#btn_search").click(function() {
        activeNodeId = null;
        activeClusters = $("#filter_cluster").val() || [];
        activeExposures = $("#filter_exposure").val() || [];
        activeOutcomes = $("#filter_outcome").val() || [];
        renderGraph();
    });

    $("#btn_clear").click(function() {
        $("#filter_cluster").val(null).trigger("change");
        $("#filter_exposure").val(null).trigger("change");
        $("#filter_outcome").val(null).trigger("change");
        activeNodeId = null;
        activeClusters = []; activeExposures = []; activeOutcomes = [];
        $(".edge-cb").prop("checked", false);
        activeEdgeIndices = [];
        renderGraph();
    });

    let lastClickTime = 0;
    document.getElementById("graph").on("plotly_click", function(data) {
        const currentTime = new Date().getTime();
        if (currentTime - lastClickTime < 400) { 
            if (data.points[0].text && data.points[0].text.indexOf("<br>") === -1) { 
                const clickedLabel = data.points[0].text;
                const foundNode = nodes.find(n => n.label === clickedLabel);
                if (foundNode) {
                    activeNodeId = foundNode.id;
                    $("#btn_clear").click();
                    activeNodeId = foundNode.id;
                    renderGraph();
                }
            }
        }
        lastClickTime = currentTime;
    });
    
    document.getElementById("graph").on("plotly_doubleclick", function() {
        $("#btn_clear").click();
    });

  </script>
</body>
</html>
')

writeLines(html_content, file.path(interactive_dir, "2_Interactive_Network.html"))

# 8. Generate Readme
# ------------------------------------------------------------------------------
cat("[INFO] Generating Readme file...\n")
readme_file <- file.path(readme_dir, paste0("Visualization_EdgeBundling_Optimized_", timestamp, "_Readme.txt"))

readme_content <- paste0(
  "========================================================================\n",
  "Project: cis-MR Causal Network & Hierarchical Edge Bundling Visualization\n",
  "Date: ", Sys.time(), "\n",
  "========================================================================\n\n",
  "[Task Description]\n",
  "Visualization of significant cis-MR results, combining network graphs with \n",
  "Hierarchical Edge Bundling. Louvain algorithm was used for native clustering.\n",
  "Interactive 3D plot is optimized to inherit the edge-bundling circular layout and supports full-screen.\n\n",
  "[Input Data]\n",
  "MR Results: ", opt$input_csv, "\n",
  "Clustering: Built-in Louvain Algorithm\n\n",
  "[Network Characteristics]\n",
  "Total Proteins (Nodes): ", nrow(nodes), "\n",
  "Total Causal Relationships (Edges): ", nrow(edges), "\n",
  "Positive Effects: ", sum(edges$direction == "Positive"), "\n",
  "Negative Effects: ", sum(edges$direction == "Negative"), "\n\n",
  "[Output Structure]\n",
  "1. plots/\n",
  "   - 1_Hierarchical_EdgeBundling.png/.pdf : Static hierarchical edge bundling plot.\n",
  "   - 3D_Interactive/2_Interactive_Network.html : Interactive web-based 3D network plot.\n",
  "2. network_data/\n",
  "   - Network_Nodes.csv : Node properties (In-degree, Out-degree, Cluster, etc.)\n",
  "   - Network_Edges.csv : Edge properties\n",
  "   - Network_Data_Combined.xlsx : Combined Excel table of nodes and edges\n"
)
writeLines(readme_content, readme_file)

cat("[INFO] Script finished successfully. Check output directory for results.\n")
sink()
