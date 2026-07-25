# ==============================================================================
# Fig1.R — Figure 1: Differential Expression and Pathway Analysis
# ==============================================================================
# Author : Ximing Ran
# Date   : 2026-02-22
#
# Description:
#   Generates all panels of Figure 1 for the ALS plasma proteomics study:
#   Fig1A - Volcano plot of differentially expressed proteins (Miami cohort,
#           linear mixed-effects model; FDR < 0.05 threshold)
#   Fig1B - Pathway enrichment bar plots for up- and down-regulated proteins
#           (g:Profiler GO/KEGG/Reactome enrichment)
#   Fig1C - Protein–protein interaction (PPI) network of significant proteins,
#           coloured by network community / module
#   Fig1E - Expression heatmap of significant proteins across subjects,
#           annotated by PPI module membership (requires individual data)
#
# Input data:
#   data/fig_data/summary_data/Fig1/
#     M1_Mix_effect_model_lmer(2025-02-12).csv  — per-protein lmer results
#     (pathway enrichment CSVs for up- and down-regulated sets)
#     attr_df.csv                               — PPI node attributes
#     interactions.csv                          — PPI edge list
#   data/fig_data/individual_data/Fig1/  (optional — Fig1E skipped if absent)
#     Heatmap_PPI_cluster_mask_with_cluster_matrix.csv
#     Heatmap_PPI_cluster_mask_with_cluster_annotation_col.csv
#     Heatmap_PPI_cluster_mask_with_cluster_annotation_row.csv
#
# Output (./results/):
#   Fig_1A_volcano.{pdf,png}
#   Fig_1B_pathway_up.{pdf,png}, Fig_1B_pathway_down.{pdf,png}
#   Fig_1C_PPI_network.{pdf,png}
#   Fig_1E_heatmap.{pdf,png}  (only if individual data exists)
#
# Prerequisite:
#   Run analysis/01-Differentially_Expressed_Proteins/ Rmd notebooks first.
# ==============================================================================




# ----Setup and Configuration----


# Load required libraries
library(tidyverse)
library(arrow)
library(OlinkAnalyze)
library(gridExtra)
library(ggplot2)
library(stringr)
library(dplyr)
library(tidyr)
library(extrafont)
library(ggsignif)
library(jtools)
library(cowplot)
library(ggpubr)
library(lme4)
library(lmerTest)
library(knitr)
library(kableExtra)
library(progress)
library(ggrepel)
library(gprofiler2)
library(ggraph)
library(tidygraph)
library(scales)
library(ggforce)
library(ggnewscale)
library(igraph)
library(pheatmap)
library(RColorBrewer)
library(ggplotify)
library(ggtext)
library(here)
library(svglite)

# Set random seed for reproducibility
set.seed(2025)

# Set default ggplot theme
theme_set(theme_bw() + theme(legend.position = "bottom"))

# Create output directory if it doesn't exist
output_dir <- here::here("figure", "Fig1", "results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# ----Figure 1A: Volcano Plot of Differential Expression ----

# Load Miami cohort differential expression results
miami_result <- read.csv(
  here::here("data","fig_data", "summary_data", "Fig1", 
             "M1_Mix_effect_model_lmer(2025-02-12).csv"),
  header = TRUE, 
  stringsAsFactors = FALSE
)

# Compute p-value cutoff for FDR ≤ 0.05
p_cutoff <- miami_result %>%
  filter(padj <= 0.05) %>%
  summarise(cutoff = max(Pr_t, na.rm = TRUE)) %>%
  pull(cutoff)

# Prepare volcano plot data
volcano_df <- miami_result %>%
  mutate(
    negLogP = -log10(Pr_t),
    significance = diffexpressed
  )

# Count significant proteins
n_up <- sum(volcano_df$significance == "UP")
n_down <- sum(volcano_df$significance == "DOWN")

# Identify top 10 UP and DOWN proteins by p-value
top_up <- volcano_df %>%
  filter(significance == "UP") %>%
  arrange(Pr_t) %>%
  slice_head(n = 10)

top_down <- volcano_df %>%
  filter(significance == "DOWN") %>%
  arrange(Pr_t) %>%
  slice_head(n = 10)

top_10_genes <- bind_rows(top_up, top_down)

# Create volcano plot
p_a <- ggplot(volcano_df, aes(x = Estimate, y = negLogP, 
                              color = significance)) +
  geom_point(alpha = 0.9, size = 1.5) +
  geom_hline(yintercept = -log10(p_cutoff ), linetype = "dashed") +
  # annotate(
  #   "text",
  #   x = min(volcano_df$Estimate, na.rm = TRUE) + 0.2,
  #   y = max(volcano_df$negLogP, na.rm = TRUE),
  #   label = paste0("DOWN: ", n_down),
  #   hjust = 0, vjust = 1.1,
  #   family = "Arial",
  #   size = 6,
  #   color = "#427ebf"
  # ) +
  # annotate(
  #   "text",
  #   x = max(volcano_df$Estimate, na.rm = TRUE) - 0.2,
  #   y = max(volcano_df$negLogP, na.rm = TRUE),
  #   label = paste0("UP: ", n_up),
  #   hjust = 1, vjust = 1.1,
  #   family = "Arial",
  #   size = 6,
  #   color = "#A70C20"
  # ) +
  scale_color_manual(
    values = c(
      UP = "#A70C20",
      DOWN = "#427ebf",
      NO = "grey50"
    )
  ) +
  labs(
    
    # title = "ALS versus Control (Miami cohort)"
    x = "log2FoldChange",
    y = "-log10(p-value)"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    plot.title = element_text(hjust = 0.5),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  ) +
  geom_label_repel(
    data = top_10_genes,
    aes(label = Protein, color = significance),
    size = 4,
    max.overlaps = 50,
    min.segment.length = 0,
    show.legend = FALSE,
    # family = "Arial",
    box.padding = 0.3,
    point.padding = 0.2
  )

# Save Figure 1A
ggsave(
  file.path(output_dir, "Fig_1A_volcano.pdf"),
  plot = p_a,
  device = cairo_pdf,
  width = 6,
  height = 6
)

# Save Figure 1A
ggsave(
  file.path(output_dir, "Fig_1A_volcano.png"),
  plot = p_a,
  width = 6,
  height = 6,
  dpi = 300
)
ggsave(
  file.path(output_dir, "Fig_1A_volcano.svg"),
  plot = p_a,
  device = svglite::svglite,
  width = 6,
  height = 6
)


# ------Figure 1B: Up-regulated Pathway Enrichment -----------------------------

# Load up-regulated enrichment results
result_up <- read.csv(
  here::here("data","fig_data",  "summary_data", "Fig1", 
             "M1_Mix_effect_model_lmer_up_enrichment(2025-02-12).csv"),
  header = TRUE, 
  stringsAsFactors = FALSE
) %>%
  arrange(p_value) %>%
  select(source, term_name, p_value) %>%
  mutate(term = sprintf("%s - %s", source, term_name))

# Select top 20 pathways
top_30_up <- result_up %>% 
  slice_head(n = 20) %>%
  arrange(desc(p_value))

# Keep original order for axis
term_order <- unique(top_30_up$term)

# Create enrichment plot with highlighted muscle-related terms
p_b <- top_30_up %>%
  mutate(
    highlight = ifelse(grepl("muscle", term, ignore.case = TRUE), 
                       "highlight", "normal"),
    term_md = ifelse(highlight == "highlight", 
                     paste0("<b>", term, "</b>"), term),
    term_md = factor(term_md, levels = ifelse(
      grepl("muscle", term_order, ignore.case = TRUE),
      paste0("<b>", term_order, "</b>"),
      term_order
    ))
  ) %>%
  ggplot(aes(term_md, -log2(p_value), fill = highlight)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  scale_fill_manual(values = c("highlight" = "#A70C20", "normal" = "#d39a8d")) +
  ylab("-log2(FDR)") +
  xlab("") +
  ggtitle("Top 20 Upregulated Pathways") +
  coord_flip() +
  theme_classic(base_size = 14) +
  theme(
    text = element_text(
                        size = 14, 
                        # family = "Arial", 
                        colour = "black"),
    axis.title = element_text(color = "black"),
    axis.text.x = element_text(color = "black"),
    axis.text.y = ggtext::element_markdown(color = "black"),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.title = element_text(
                              size = 14, 
                              # family = "Arial",
                              face = "bold", hjust = 0.5)
  ) +
  geom_hline(yintercept = -log2(0.05), linetype = "dashed", color = "black")

# Save Figure 1B
ggsave(
  file.path(output_dir, "Fig_1B_up_enrichment.pdf"),
  plot = p_b,
  device = cairo_pdf,
  width = 8,
  height = 6
)
ggsave(
  file.path(output_dir, "Fig_1B_up_enrichment.png"),
  plot = p_b,
  width = 8,
  height = 6,
  dpi = 300
)
ggsave(
  file.path(output_dir, "Fig_1B_up_enrichment.svg"),
  plot = p_b,
  device = svglite::svglite,
  width = 8,
  height = 6
)

# -----Figure 1C: Down-regulated Pathway Enrichment ----------------------------

# Load down-regulated enrichment results
result_down <- read.csv(
  here::here("data", "fig_data", "summary_data", "Fig1", 
             "M1_Mix_effect_model_lmer_down_enrichment(2025-02-12).csv"),
  header = TRUE, 
  stringsAsFactors = FALSE
) %>%
  select(source, term_name, p_value) %>%
  mutate(term = sprintf("%s - %s", source, term_name)) %>%
  arrange(p_value)

# Select top 20 pathways
top_30_down <- result_down %>%
  slice_head(n = 20) %>%
  arrange(desc(p_value))

# Create enrichment plot
p_c <- top_30_down %>%
  mutate(
    term = factor(term, levels = unique(term)),
    highlight = ifelse(grepl("muscle", term, ignore.case = TRUE), 
                       "highlight", "normal")
  ) %>%
  ggplot(aes(term, -log2(p_value), fill = highlight)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  scale_fill_manual(values = c("highlight" = "#427ebf", "normal" = "#72b2e0")) +
  geom_hline(yintercept = -log2(0.05), linetype = "dashed", color = "black") +
  ylab("-log2(FDR)") +
  xlab("") +
  ggtitle("Top 20 Downregulated Pathways") +
  coord_flip() +
  theme_classic(base_size = 14) +
  theme(
    text = element_text( size = 14, 
                         # family = "Arial",
                         colour = "black"),
    axis.title = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA),
    plot.title = element_text(
                              size = 14, 
                              # family = "Arial",
                              face = "bold", hjust = 0.5)
  )

# Save Figure 1C
ggsave(
  file.path(output_dir, "Fig_1C_down_enrichment.pdf"),
  plot = p_c,
  device = cairo_pdf,
  width = 8,
  height = 6
)

ggsave(
  file.path(output_dir, "Fig_1C_down_enrichment.png"),
  plot = p_c,
  width = 8,
  height = 6,
  dpi = 300
)
ggsave(
  file.path(output_dir, "Fig_1C_down_enrichment.svg"),
  plot = p_c,
  device = svglite::svglite,
  width = 8,
  height = 6
)


# Create B+C stack without labels for separate export
bc_stack_null <- plot_grid(
  p_b, p_c,
  ncol = 1,
  align = "v", 
  axis = "l"
)

ggsave(
  file.path(output_dir, "Fig_1_BC_stack_null.pdf"),
  plot = bc_stack_null,
  device = cairo_pdf,
  width = 8,
  height = 10
)

ggsave(
  file.path(output_dir, "Fig_1_BC_stack_null.png"),
  plot = bc_stack_null,
  width = 8,
  height = 10,
  dpi = 300
)
ggsave(
  file.path(output_dir, "Fig_1_BC_stack_null.svg"),
  plot = bc_stack_null,
  device = svglite::svglite,
  width = 8,
  height = 10
)


# ----Figure 1D: Protein-Protein Interaction (PPI) Network----------------------

# Load PPI network data
attr_df <- read.csv(
  here::here("data","fig_data",  "summary_data", "Fig1", "attr_df.csv"), 
  header = TRUE, 
  stringsAsFactors = FALSE
)

interactions <- read.csv(
  here::here("data","fig_data",  "summary_data", "Fig1", "interactions.csv"), 
  header = TRUE, 
  stringsAsFactors = FALSE
)

# Build igraph from edges and vertex metadata
ppi_graph <- graph_from_data_frame(
  d = interactions,
  vertices = attr_df %>% distinct(STRING_id, .keep_all = TRUE),
  directed = FALSE
)

# Add vertex attributes
V(ppi_graph)$Protein <- attr_df$Protein[match(V(ppi_graph)$name, attr_df$STRING_id)]
V(ppi_graph)$Eff_ALS <- attr_df$Eff_ALS[match(V(ppi_graph)$name, attr_df$STRING_id)]
V(ppi_graph)$padj_ALS <- attr_df$padj_ALS[match(V(ppi_graph)$name, attr_df$STRING_id)]
V(ppi_graph)$size <- attr_df$log10padj[match(V(ppi_graph)$name, attr_df$STRING_id)]

# Remove isolated nodes
deg <- degree(ppi_graph)
ppi_graph <- delete_vertices(ppi_graph, V(ppi_graph)[which(deg == 0)])

# Add clustering information
V(ppi_graph)$cluster <- attr_df$cluster[match(V(ppi_graph)$name, attr_df$STRING_id)]
V(ppi_graph)$cluster_color <- attr_df$cluster_color[match(V(ppi_graph)$name, attr_df$STRING_id)]

# Add edge attributes
E(ppi_graph)$weight <- E(ppi_graph)$combined_score / 1000

cluster_map <- setNames(V(ppi_graph)$cluster, V(ppi_graph)$name)
E(ppi_graph)$between_cluster <- apply(
  as.data.frame(ends(ppi_graph, E(ppi_graph))),
  1,
  function(pair) cluster_map[pair[1]] != cluster_map[pair[2]]
)

# Create tidygraph and compute layout
ppi_tbl <- as_tbl_graph(ppi_graph)
layout <- create_layout(ppi_tbl, layout = "stress")

# Compute Eff_ALS range for color scale
eff_range <- range(V(ppi_graph)$Eff_ALS, na.rm = TRUE)

# Define module mapping with colors and labels
module_mapping <- data.frame(
  cluster = as.character(c(1, 3, 2, 8, 5, 6)),
  module = c("ModuleA", "ModuleA", "ModuleB", "ModuleC", "ModuleD", "ModuleD"),
  module_color = c("#009E73", "#009E73", "#E69F00", "#56B4E9", "#CC79A7", "#CC79A7"),
  module_label = c(
    "Skeletal Muscle", 
    "Skeletal Muscle",
    "TNF-Mediated Signaling",
    "Regeneration (incl. Neurofilament)",
    "ECM & ECM Interaction",
    "ECM & ECM Interaction"
  ),
  stringsAsFactors = FALSE
)

# Extract nodes for module hulls
hull_df <- as.data.frame(layout) %>%
  left_join(module_mapping, by = "cluster") %>%
  filter(cluster %in% c(1, 3, 2, 8, 5, 6))

# Create PPI network plot
p_d <- ggraph(layout) +
  geom_mark_hull(
    data = hull_df,
    aes(x = x, y = y, group = module, fill = module_color, color = module_color),
    expand = unit(2.5, "mm"),
    alpha = 0.15,
    show.legend = FALSE
  ) +
  scale_fill_identity() +
  scale_color_identity() +
  new_scale_fill() +
  geom_edge_link(
    aes(width = weight),
    color = "grey80",
    linetype = "solid",
    alpha = 0.5
  ) +
  geom_node_point(
    aes(fill = Eff_ALS, size = size),
    shape = 21,
    stroke = 0.1
  ) +
  geom_node_text(
    aes(label = Protein),
    repel = TRUE,
    size = 3
  ) +
  scale_fill_gradientn(
    colors = c("blue", "white", "red"),
    values = rescale(c(eff_range[1], 0, eff_range[2])),
    name = "log2FC"
  ) +
  scale_edge_linetype_manual(
    values = c("TRUE" = "dashed", "FALSE" = "solid"),
    guide = "none"
  ) +
  scale_size_continuous(
    range = c(3, 10),
    name = "-log10(padj)"
  ) +
  scale_edge_width(
    range = c(0.3, 2),
    name = "STRING Score"
  ) +
  theme_void() +
  theme(
    legend.position = "right",
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  ) +
  labs(fill = "ALS")

# Save Figure 1D
# ggsave(
#   file.path(output_dir, "Fig_1D_ppi_network.png"), 
#   plot = p_d, 
#   width = 10, 
#   height = 8, 
#   dpi = 300
# )
ggsave(
  file.path(output_dir, "Fig_1D_ppi_network.pdf"),
  plot = p_d,
  device = cairo_pdf,
  width = 10,
  height = 8
)
ggsave(
  file.path(output_dir, "Fig_1D_ppi_network.png"),
  plot = p_d,
  width = 10,
  height = 8
)
ggsave(
  file.path(output_dir, "Fig_1D_ppi_network.svg"),
  plot = p_d,
  device = svglite::svglite,
  width = 10,
  height = 8
)

# -----Figure 1E: Heatmap of Selected Proteins--------------------------------

# Check for individual data files
fig1e_data_path <- here::here("data", "fig_data", "individual_data", "Fig1",
                               "Heatmap_PPI_cluster_mask_with_cluster_matrix.csv")
if (!file.exists(fig1e_data_path)) {
  cat("  WARNING: Individual data not found for Fig1E. Skipping.\n")
} else {

# Load heatmap data
mat_sel <- read.csv(fig1e_data_path, row.names = 1)

annotation_col <- read.csv(
  here::here("data", "fig_data", "individual_data", "Fig1", 
             "Heatmap_PPI_cluster_mask_with_cluster_annotation_col.csv"),
  row.names = 1
)
colnames(annotation_col) <- "Group"

annotation_row <- read.csv(
  here::here("data","fig_data",  "individual_data", "Fig1", 
             "Heatmap_PPI_cluster_mask_with_cluster_annotation_row.csv"),
  row.names = 1
)

# Align rows between matrix and annotation
common_rows <- intersect(rownames(mat_sel), rownames(annotation_row))
mat_sel <- mat_sel[common_rows, , drop = FALSE]
annotation_row <- annotation_row[common_rows, , drop = FALSE]

stopifnot(all(c("cluster", "cluster_color") %in% colnames(annotation_row)))
annotation_row$cluster <- as.integer(annotation_row$cluster)

# Define module mapping
module_map <- c(
  `1` = "Skeletal Muscle",
  `2` = "TNF-Mediated Signaling",
  `3` = "Skeletal Muscle",
  `5` = "ECM & ECM Interaction",
  `6` = "ECM & ECM Interaction",
  `8` = "Regeneration (incl. Neurofilament)"
)

module_color_map <- c(
  "Skeletal Muscle" = "#009E73",
  "TNF-Mediated Signaling" = "#E69F00",
  "Regeneration (incl. Neurofilament)" = "#56B4E9",
  "ECM & ECM Interaction" = "#CC79A7"
)

# Map clusters to modules and filter
annotation_row$module <- ifelse(
  annotation_row$cluster %in% c(1, 2, 3, 5, 6, 8),
  unname(module_map[as.character(annotation_row$cluster)]),
  NA
)

annotation_row <- annotation_row[!is.na(annotation_row$module), , drop = FALSE]
mat_sel <- mat_sel[rownames(annotation_row), , drop = FALSE]
annotation_row$cluster_color <- unname(module_color_map[annotation_row$module])

# Define annotation colors
ann_colors <- list(
  Group = c("CTRL" = "#63a33b", "ALS-Symptomatic" = "#b32812")
)

# Prepare module annotation for plot
ann_row_for_plot <- annotation_row[, "module", drop = FALSE]
colnames(ann_row_for_plot) <- "Module"
module_levels <- c(
  "Skeletal Muscle",
  "TNF-Mediated Signaling",
  "Regeneration (incl. Neurofilament)",
  "ECM & ECM Interaction"
)
ann_row_for_plot$Module <- factor(ann_row_for_plot$Module, levels = module_levels)

mod_color_tab <- data.frame(
  module = module_levels,
  cluster_color = c("#009E73", "#E69F00", "#56B4E9", "#CC79A7")
)
ann_colors$Module <- setNames(mod_color_tab$cluster_color, mod_color_tab$module)

# Column clustering within groups
col_dist_method <- "correlation"
col_hclust_method <- "ward.D2"

scale_cols <- function(m) scale(m)
col_dist <- function(m, method = "euclidean") {
  if (method == "correlation") {
    as.dist(1 - cor(m, use = "pairwise.complete.obs"))
  } else {
    dist(t(m), method = method)
  }
}

if (!is.factor(annotation_col$Group)) {
  annotation_col$Group <- factor(annotation_col$Group, 
                                 levels = unique(annotation_col$Group))
}

groups_in_order <- levels(annotation_col$Group)
col_idx_list <- split(seq_len(ncol(mat_sel)), annotation_col$Group)
new_col_order <- integer(0)

for (grp in groups_in_order) {
  idx <- col_idx_list[[grp]]
  if (length(idx) == 0) next
  
  submat <- mat_sel[, idx, drop = FALSE]
  submat_scaled <- suppressWarnings(scale_cols(submat))
  
  if (ncol(submat_scaled) > 2) {
    d <- col_dist(submat_scaled, method = col_dist_method)
    hc <- hclust(d, method = col_hclust_method)
    ord <- idx[hc$order]
  } else if (ncol(submat_scaled) == 2) {
    ord <- idx[order(colMeans(submat_scaled, na.rm = TRUE), decreasing = TRUE)]
  } else {
    ord <- idx
  }
  
  new_col_order <- c(new_col_order, ord)
}

mat_ord <- mat_sel[, new_col_order, drop = FALSE]
annotation_col <- annotation_col[new_col_order, , drop = FALSE]

# Row clustering within modules
row_dist_method <- "correlation"
row_hclust_method <- "ward.D2"

scale_rows <- function(m) t(scale(t(m)))
row_dist <- function(m, method = "euclidean") {
  if (method == "correlation") {
    as.dist(1 - cor(t(m), use = "pairwise.complete.obs"))
  } else {
    dist(m, method = method)
  }
}

ann_row_for_plot <- ann_row_for_plot[rownames(mat_ord), , drop = FALSE]
modules_in_order <- levels(ann_row_for_plot$Module)
row_idx_list <- split(seq_len(nrow(mat_ord)), ann_row_for_plot$Module)
new_row_order <- integer(0)

for (mod in modules_in_order) {
  idx <- row_idx_list[[mod]]
  if (length(idx) == 0) next
  
  submat <- mat_ord[idx, , drop = FALSE]
  submat_scaled <- suppressWarnings(scale_rows(submat))
  
  if (nrow(submat_scaled) > 2) {
    d <- row_dist(submat_scaled, method = row_dist_method)
    hc <- hclust(d, method = row_hclust_method)
    ord <- idx[hc$order]
  } else if (nrow(submat_scaled) == 2) {
    ord <- idx[order(rowMeans(submat_scaled, na.rm = TRUE), decreasing = TRUE)]
  } else {
    ord <- idx
  }
  
  new_row_order <- c(new_row_order, ord)
}

mat_ord <- mat_ord[new_row_order, , drop = FALSE]
ann_row_for_plot <- ann_row_for_plot[new_row_order, , drop = FALSE]

# Calculate gaps
sizes_row <- table(ann_row_for_plot$Module)[module_levels]
gaps_row <- cumsum(head(as.integer(sizes_row), -1))

sizes_col <- table(annotation_col$Group)[unique(annotation_col$Group)]
gaps_col <- cumsum(head(as.integer(sizes_col), -1))

# Heatmap color palette
cell_pal <- rev(colorRampPalette(brewer.pal(11, "RdBu"))(256))

# Draw heatmap
ph <- pheatmap(
  mat_ord,
  color = cell_pal,
  scale = "row",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = FALSE,
  annotation_row = ann_row_for_plot,
  annotation_colors = ann_colors,
  gaps_row = gaps_row,
  gaps_col = gaps_col,
  fontsize_row = 10,
  annotation_names_row = FALSE,
  main = ""
)

# Convert to ggplot
gp <- ggplotify::as.ggplot(ph$gtable)
p_e_2 <- gp +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

# Save Figure 1E
ggsave(
  file.path(output_dir, "Fig_1E_heatmap.png"),
  plot = p_e_2,
  width = 12,
  height = 15,
  dpi = 500
)

# save a pdf version
ggsave(
  file.path(output_dir, "Fig_1E_heatmap.pdf"),
  plot = p_e_2,
  device = cairo_pdf,
  width = 12,
  height = 15
)
ggsave(
  file.path(output_dir, "Fig_1E_heatmap.svg"),
  plot = p_e_2,
  device = svglite::svglite,
  width = 12,
  height = 15
)
#
# # ----Combine All Panels into Final Figure 1----------------------------
# 
# # Create B+C stack with labels
# bc_stack <- plot_grid(
#   p_b_2, p_c,
#   ncol = 1,
#   align = "v", 
#   axis = "l",
#   labels = c("B", "C"),
#   label_size = 20, 
#   label_fontface = "bold",
#   hjust = 0, 
#   vjust = 1, 
#   label_x = 0.01, 
#   label_y = 0.99
# )
# 
# # Create A panel with spacing
# p_a_labeled <- plot_grid(
#   NULL, p_a, NULL,
#   ncol = 3, 
#   rel_widths = c(0.15, 0.7, 0.15),
#   hjust = 0, 
#   vjust = 1, 
#   label_x = 0.01, 
#   label_y = 0.99
# )
# 
# # Left column: A + (B,C) stack
# left_col <- plot_grid(
#   p_a_labeled, bc_stack,
#   labels = "A",
#   label_size = 20, 
#   label_fontface = "bold",
#   ncol = 1,
#   rel_heights = c(1, 2)
# )
# 
# # Right column: D + E
# right_top <- plot_grid(p_d, NULL, ncol = 2, rel_widths = c(0.8, 0.2))
# right_bottom <- plot_grid(NULL, p_e_2, nrow = 2, rel_heights = c(0.02, 1.6))
# 
# right_col <- plot_grid(
#   right_top, right_bottom,
#   ncol = 1,
#   rel_heights = c(1, 2),
#   labels = c("D", "E"),
#   label_size = 20, 
#   label_fontface = "bold",
#   hjust = 0, 
#   vjust = 1, 
#   label_x = 0.01, 
#   label_y = 0.99
# )
# 
# # Combine all panels
# combo <- plot_grid(
#   left_col, right_col,
#   ncol = 2,
#   rel_widths = c(1, 1),
#   align = "h"
# )
# 
# # Save combined figure
# ggsave(
#   file.path(output_dir, "Fig_1_combine.png"), 
#   plot = combo, 
#   width = 20, 
#   height = 18, 
#   dpi = 300
# )


} # end Fig1E individual data guard

# =====END OF SCRIPT=========================================================