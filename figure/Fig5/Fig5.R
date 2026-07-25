# ==============================================================================
# Fig5.R — Figure 5: UK Biobank External Validation
# ==============================================================================
# Author : Ximing Ran
# Date   : 2026-02-24
#
# Description:
#   Generates all panels of Figure 5 for the ALS plasma proteomics study,
#   reproducing the UK Biobank (UKB) replication results:
#   Fig5A - Scatter plot correlating per-protein effect sizes (logFC) between
#           the Miami discovery cohort and the UKB replication cohort
#           (95 overlapping proteins); Panel 15 proteins highlighted
#   Fig5B - Temporal heatmap of UKB protein log2FC across 1-year bins
#           relative to ALS diagnosis, grouped by PPI module
#   Fig5C - GAM-smoothed trajectory plots for pre-ALS converters in UKB
#           (up to 2 years before diagnosis)
#   Fig5D - Predicted vs observed ALS onset timing in UKB using the
#           Miami-trained multi-protein prediction model
#
# Input data:
#   data/fig_data/summary_data/Fig5/
#     M1_Mix_effect_model_lmer(2025-02-12).csv           — Miami lmer logFC
#     UKB_limma_als_vs_control_2025-11-17.csv            — UKB limma logFC
#     heatmap_wide_UKB_Convert_—_YrSinceDi.csv           — UKB 1-yr binned
#                                                          mean log2FC
#     (UKB PPI module annotation CSV)
#   data/fig_data/individual_data/Fig5/  (optional — Fig5C/D skipped if absent)
#     confirm_df_Convert.csv              — per-subject UKB trajectory data
#     grid_overall_Convert.csv            — GAM predictions for UKB
#     result_individual_predictions_df.csv — UKB predicted vs observed values
#
# Output (./results/):
#   Fig_5A_Replication_vs_discovery_correlation.{pdf,png}
#   Fig_5B_UKB_heatmap_1yr.{pdf,png}
#   Fig_5C_UKB_pre_ALS_2yr.{pdf,png}     (only if individual data exists)
#   Fig_5D_UKB_DA_model.{pdf,png}         (only if individual data exists)
#
# Prerequisite:
#   Run analysis/06-UKB_analyisis/ Jupyter notebooks first.
# ==============================================================================

# ----------------------#
# 0. Setup & Libraries  #
# ----------------------#

# Core tidyverse and plotting
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(patchwork)

# Paths
library(here)

# Heatmap + helpers
library(pheatmap)
library(RColorBrewer)
library(ggplotify)
library(svglite)

set.seed(2025)

# Global theme for ggplot
theme_set(theme_bw() + theme(legend.position = "bottom"))


############################################################
# Fig5A – UKB vs Miami: Effect-size correlation (95 proteins)
############################################################

# Panel 15 proteins
Panel_15 <- c(
  "NEFL", "DUSP29", "CALCA", "NEB", "NGRN",
  "NOS1", "DTNB", "MEGF10", "APOA4", "ART3",
  "MYL3", "EPHA1", "EDA2R", "TTN", "ACTN2"
)

# 1) Read discovery (Miami) and replication (UKB) results
miami_result <- read.csv(
  here::here("data",  "fig_data", "summary_data", "Fig5", "M1_Mix_effect_model_lmer(2025-02-12).csv"),
  header = TRUE,
  stringsAsFactors = FALSE,
  row.names = 2
)

ukb_result <- read.csv(
  here::here("data", "fig_data",  "summary_data", "Fig5", "UKB_limma_als_vs_control_2025-11-17.csv"),
  header = TRUE,
  stringsAsFactors = FALSE,
  row.names = 1
) %>%
  rownames_to_column("protein") %>% 
  # rename(logFC = Estimate) %>%
  select(protein, logFC)

# 2) Restrict to the 95 DE proteins from discovery and join to UKB
result_combine <- miami_result %>%
  filter(padj < 0.05) %>%
  rownames_to_column("protein") %>%
  select(protein, Estimate) %>%
  left_join(ukb_result, by = "protein") %>%
  filter(!is.na(logFC))   # 95 proteins

# 3) Define direction consistency (Up/Down in both vs others)
result_combine <- result_combine %>%
  mutate(
    color_group = case_when(
      Estimate > 0 & logFC > 0 ~ "Up in both",
      Estimate < 0 & logFC < 0 ~ "Down in both",
      TRUE                     ~ "Opposite/Neutral"
    )
  )

# Proteins for labels (Panel 15)
top_10_genes <- result_combine %>%
  filter(protein %in% Panel_15)

# 4) Correlation statistics
df_plot <- result_combine %>%
  select(protein, Estimate, logFC) %>%
  filter(is.finite(Estimate), is.finite(logFC))

N    <- nrow(df_plot)
ct   <- cor.test(df_plot$Estimate, df_plot$logFC)
corr <- unname(ct$estimate)
pval <- ct$p.value

# 5) Scatter plot: Discovery vs Replication log2FC
p_Fig5A <- ggplot(result_combine,
                  aes(x = Estimate, y = logFC, color = color_group)) +
  geom_smooth(method = "lm", color = "black", se = TRUE) +
  geom_point(alpha = 0.8, size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "black") +
  geom_label_repel(
    data = top_10_genes,
    aes(label = protein),
    size = 4,
    max.overlaps = 50,
    min.segment.length = 0,
    show.legend = FALSE,
    box.padding = 0.3,
    point.padding = 0.2
  ) +
  scale_color_manual(
    values = c(
      "Up in both"        = "#A70C20",
      "Down in both"      = "#427ebf",
      "Opposite/Neutral"  = "grey50"
    )
  ) +
  labs(
    title = "Discovery vs Replication Cohort",
    x     = "log2FC (Discovery)",
    y     = "log2FC (Replication)",
    color = "Direction"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position  = "none",
    plot.title       = element_text(hjust = 0.5, face = "bold"),
    plot.background  = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  ) +
  annotate(
    "text",
    x = -1.1,   # adjust if needed based on data range
    y =  1.2,
    hjust = 0, vjust = 1,
    label = paste0(
      "N = ", N, "\n",
      "Cor = ", round(corr, 2), "\n",
      # "p = ", formatC(pval, format = "e", digits = 2)
      "p < 0.001"
    ),
    size      = 5,
    color     = "black",
    lineheight = 1.1
  )

p_Fig5A

# Save Fig5A
fig5_out <- here::here("figure", "Fig5", "results")
dir.create(fig5_out, recursive = TRUE, showWarnings = FALSE)

ggsave(file.path(fig5_out, "Fig_5A_Replication_vs_discovery_correlation-v2.pdf"),
       plot  = p_Fig5A,
       device = cairo_pdf,
       width  = 6,
       height = 6)
ggsave(file.path(fig5_out, "Fig_5A_Replication_vs_discovery_correlation-v2.png"),
       plot  = p_Fig5A,
       width  = 6,
       height = 6,
       dpi   = 300)
ggsave(file.path(fig5_out, "Fig_5A_Replication_vs_discovery_correlation-v2.svg"),
       plot  = p_Fig5A,
       device = svglite::svglite,
       width  = 6,
       height = 6)

############################################################
# Fig5B – UKB Heatmap of significant proteins over time
############################################################

# 1) Read heatmap data + annotation
df_heat <- read.csv(
  here::here("data", "fig_data","summary_data", "Fig5", "heatmap_wide_UKB_Convert_—_YrSinceDi.csv")
)

# clean column names: X.-13.5 -> -15.5
# colnames(df_heat) <- gsub("^X\\.", "-", colnames(df_heat))
# colnames(df_heat)[2:(ncol(df_heat)-6)] <- c(paste0("-Y", seq(14, 1)), paste0("Y", seq(1, 18))) # Raw
colnames(df_heat)[2:(ncol(df_heat)-6)] <- c(paste0("-Y", seq(12, 1)), paste0("Y", seq(1, 20))) # Raw



annotation_row <- read.csv(
  here::here("data", "fig_data","summary_data", "Fig5",
             "Heatmap_PPI_cluster_mask_with_cluster_annotation_row.csv"),
  row.names = 1
) %>%
  select(cluster) %>%
  rownames_to_column("protein")

# add cluster info to main df
df_heat <- df_heat %>%
  left_join(annotation_row, by = "protein") %>%
  filter(Adj_P_value < 0.05)

# 2) Map cluster -> module names
module_map <- c(
  `1` = "Skeletal Muscle",
  `2` = "TNF-Mediated Signaling",
  `3` = "Skeletal Muscle",
  `5` = "ECM & ECM Interaction",
  `6` = "ECM & ECM Interaction",
  `8` = "Regeneration (incl. Neurofilament)"
)

module_levels <- c(
  "Skeletal Muscle",
  "TNF-Mediated Signaling",
  "Regeneration (incl. Neurofilament)",
  "ECM & ECM Interaction",
  "Others"
)

# Colors per module (for row annotation)
module_color_map <- c(
  "Skeletal Muscle"                        = "#009E73",
  "TNF-Mediated Signaling"                 = "#E69F00",
  "Regeneration (incl. Neurofilament)"      = "#56B4E9",
  "ECM & ECM Interaction"                  = "#CC79A7",
  "Others"                                         = "#D3D3D3"
)

# assign Module based on cluster
df_heat <- df_heat %>%
  mutate(
    Module = ifelse(
      cluster %in% c(1, 2, 3, 5, 6, 8),
      unname(module_map[as.character(cluster)]),
      "Others"
    ),
    Module = factor(Module, levels = module_levels, ordered = TRUE)
  )

# 3) Identify time-bin columns (numeric column names like "-15.5")
time_bin_cols <- colnames(df_heat)[grepl("^-?Y\\d+$", colnames(df_heat))]
# time_bin_cols <- time_bin_cols[order(as.numeric(time_bin_cols))]

# 4) Build base matrix for display and imputation
mat_display <- df_heat %>%
  select(all_of(time_bin_cols)) %>%
  as.matrix()

rownames(mat_display) <- df_heat$protein
mode(mat_display) <- "numeric"

# row-wise mean imputation for NA (for ordering)
impute_row_mean <- function(x) {
  if (all(is.na(x))) return(rep(0, length(x)))
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  x
}
mat_imputed <- t(apply(mat_display, 1, impute_row_mean))

# 5) Row ordering: by Module, then first non-NA position, then row mean
first_non_na_pos <- apply(mat_display, 1, function(r) {
  idx <- which(!is.na(r))
  if (length(idx) == 0) Inf else idx[1]
})
row_mean_imputed <- rowMeans(mat_imputed, na.rm = TRUE)

df_ordered <- df_heat %>%
  mutate(
    first_non_na_pos = first_non_na_pos[match(protein, rownames(mat_display))],
    row_mean_imputed = row_mean_imputed[match(protein, rownames(mat_display))]
  ) %>%
  arrange(Module, first_non_na_pos, desc(row_mean_imputed))

# 6) Reorder matrix accordingly
mat_display_ord <- mat_display[df_ordered$protein, , drop = FALSE]

# 7) Row annotations and row gaps (between modules)
annotation_row_heat <- data.frame(Module = df_ordered$Module)
rownames(annotation_row_heat) <- df_ordered$protein

mod_rle <- rle(as.character(df_ordered$Module))
gap_indices <- cumsum(mod_rle$lengths)
gap_indices <- gap_indices[gap_indices < nrow(df_ordered)]

annotation_colors_UKB <- list(Module = module_color_map)

# 8) Color scale
max_abs <- max(abs(mat_imputed), na.rm = TRUE)
breaks  <- seq(-max_abs, max_abs, length.out = 101)
colors  <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

# numbers to display inside cells (for non-NA entries)
number_matrix <- mat_display_ord
number_matrix[!is.na(number_matrix)] <- round(
  number_matrix[!is.na(number_matrix)], 2
)
number_matrix[is.na(number_matrix)] <- ""

# # 9) Column labels: -12..-1, 1..2 (years from phenoconversion proxy)
# colnames(mat_display_ord) <- c(
#   paste0("-Y", seq(12, 1)),
#   paste0("Y", seq(1, 2))
# )

# 10) Heatmap
ph_Fig5B <- pheatmap::pheatmap(
  mat_display_ord,
  color             = colors,
  breaks            = breaks,
  cluster_rows      = FALSE,
  cluster_cols      = FALSE,
  show_rownames     = TRUE,
  show_colnames     = TRUE,
  annotation_row    = annotation_row_heat,
  annotation_colors = annotation_colors_UKB,
  gaps_row          = gap_indices,
  gaps_col          = c(12,14),
  display_numbers   = number_matrix,
  fontsize_number   = 12,
  fontsize_row      = 12,
  annotation_names_row = FALSE,
  fontsize_col      = 12,
  fontsize          = 12,
  angle_col         = 0,
  cellwidth         = 35,
  border_color      = "black",
  number_color      = "black",
  na_col            = "white"
)

p_Fig5B <- ggplotify::as.ggplot(ph_Fig5B$gtable) +
  theme(
    plot.background  = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )



# Save Fig5B
ggsave(file.path(fig5_out, "Fig_5B_UKB_heatmap_1yr.pdf"),
       plot  = p_Fig5B,
       device = cairo_pdf,
       width = 15,
       height = 3)

ggsave(file.path(fig5_out, "Fig_5B_UKB_heatmap_1yr.png"),
       plot  = p_Fig5B,
       width = 15,
       height = 4,
       dpi   = 300)
ggsave(file.path(fig5_out, "Fig_5B_UKB_heatmap_1yr.svg"),
       plot  = p_Fig5B,
       device = svglite::svglite,
       width = 15,
       height = 3)


############################################################
# Fig5C – UKB pre-ALS longitudinal trajectories (2 yr shift)
############################################################

fig5c_path <- here::here("data", "fig_data", "individual_data", "Fig5",
                          "confirm_df_Convert.csv")
if (!file.exists(fig5c_path)) {
  cat("  WARNING: Individual data not found for Fig5C. Skipping.\n")
} else {

# 1) Load per-individual delta data + GAM grid
confirm_df <- read.csv(fig5c_path)

grid_overall <- read.csv(
  here::here("data", "fig_data",  "individual_data", "Fig5",
             "grid_overall_Convert.csv")
)

# Shift time so that 0 = phenoconversion proxy, negative = pre-ALS
confirm_df$YrSinceOs   <- confirm_df$YrSinceDi + 2
grid_overall$YrSinceOs <- grid_overall$YrSinceDi + 2

# Colors for directions
fill_map <- c(increase = "#D73027", decrease = "#4575B4")

# Proteins to loop over
# protein_list <- unique(confirm_df$protein)

protein_list <- c("CA3", "EDA2R","NEFL")
# Container for per-protein plots
plot_list <- list()

for (prot in protein_list) {
  # y-axis limits
  y_min_val <- -2.0
  y_max_val <-  3.0
  
  # extract per-protein data
  confirm_df_prot <- confirm_df %>% filter(protein == !!prot)
  grid_overall_prot <- grid_overall %>% filter(protein == !!prot)
  
  # points with significant direction at each grid time
  sig_points_df <- grid_overall_prot %>%
    filter(!is.na(dir), YrSinceOs <= 0)
  
  # clip ribbons to plot limits
  grid_overall_prot <- grid_overall_prot %>%
    mutate(
      lower = pmax(lower, y_min_val),
      upper = pmin(upper, y_max_val)
    )
  
  # positions for textual labels
  y_label_val_ind_increase <- -0.4
  y_label_val_ind_decrease <-  0.7
  
  # identify contiguous runs where sig != 0
  r     <- rle(grid_overall_prot$sig)
  ends  <- cumsum(r$lengths)
  starts <- c(1, head(ends, -1) + 1)
  nonzero_runs_idx <- which(r$values != 0)
  
  if (length(nonzero_runs_idx) > 0) {
    run_starts_idx <- starts[nonzero_runs_idx]
    run_ends_idx   <- ends[nonzero_runs_idx]
    run_vals       <- r$values[nonzero_runs_idx]
    
    run_df <- tibble(
      start_idx = run_starts_idx,
      end_idx   = run_ends_idx,
      start_year = grid_overall_prot$YrSinceOs[run_starts_idx],
      end_year   = grid_overall_prot$YrSinceOs[run_ends_idx],
      run_sign   = run_vals
    ) %>%
      mutate(
        dir         = ifelse(run_sign > 0, "increase", "decrease"),
        color_hex   = ifelse(run_sign > 0, "#D73027", "#4575B4"),
        label_start = sprintf("%.1f", start_year),
        label_end   = sprintf("%.1f", end_year),
        y_label_val_ind = ifelse(dir == "increase",
                                 y_label_val_ind_increase,
                                 y_label_val_ind_decrease)
      ) %>%
      # if end year very close to 0, drop the label
      mutate(
        label_end = ifelse(abs(end_year) < 0.1, "", label_end)
      ) %>%
      # do not display segments that start after 0 (we only want pre-ALS)
      filter(start_year <= 0) %>%
      # remove labels for positive times
      mutate(
        label_end = ifelse(end_year >= 0.2, "", label_end)
      )
  } else {
    run_df <- tibble(
      start_idx = integer(), end_idx = integer(),
      start_year = numeric(), end_year = numeric(),
      run_sign = integer(), dir = character(),
      label_start = character(), label_end = character(),
      color_hex = character(), y_label_val_ind = numeric()
    )
  }
  
  # main plot per protein
  p_prot <- ggplot() +
    # Reference line at y=0
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      color = "grey40"
    ) +
    
    # raw individual deltas
    geom_point(
      data = confirm_df_prot,
      aes(x = YrSinceOs, y = delta),
      alpha = 0.6, size = 1, color = "grey10"
    ) +
    # GAM ribbon + fit
    geom_ribbon(
      data = grid_overall_prot,
      aes(x = YrSinceOs, ymin = lower, ymax = upper),
      fill = "grey30", alpha = 0.2
    ) +
    geom_line(
      data = grid_overall_prot,
      aes(x = YrSinceOs, y = fit),
      color = "black", size = 1.6
    ) +
    # vertical line at phenoconversion
    geom_vline(xintercept = 0, color = "black") +
    # Horizontal line at y=0
    # significant time points
    geom_point(
      data = sig_points_df,
      aes(x = YrSinceOs, y = 0, color = dir),
      size = 1.2, alpha = 0.9, show.legend = FALSE
    ) +
    # text annotations for contiguous significant segments
    # geom_text(
    #   data = run_df,
    #   aes(x = start_year, y = y_label_val_ind, label = label_start),
    #   color = run_df$color_hex,
    #   size  = 3.6,
    #   fontface = "bold",
    #   vjust = 0
    # ) +
    # geom_text(
    #   data = run_df,
  #   aes(x = end_year, y = y_label_val_ind, label = label_end),
  #   color = run_df$color_hex,
  #   size  = 3.6,
  #   fontface = "bold",
  #   vjust = 0
  # ) +
  # 
  geom_label(
    data = run_df %>% filter(dir == "increase"),
    aes(x = start_year,
        y = y_label_val_ind - 0.2, label = label_start),
    fill = "white",
    color = (run_df %>% filter(dir == "increase"))$color_hex,
    size = 5,
    fontface = "bold",
    vjust = 0,
    label.size = 0.5
  ) +

    geom_label(
      data = run_df %>% filter(label_end != "", dir == "increase"),
      aes(x = end_year + 1,
          y = y_label_val_ind - 0.2, label = label_end),
      fill = "white",
      color = (run_df %>% filter(label_end != "", dir == "increase"))$color_hex,
      size = 5,
      fontface = "bold",
      vjust = 0,
      label.size = 0.5
    ) +
    
    
    scale_fill_manual(values = fill_map, name = "Direction") +
    scale_color_manual(values = fill_map, guide = FALSE) +
    labs(
      title = prot,
      x = "Years to Phenoconversion Proxy",
      y = "log2FC"
    ) +
    ylim(y_min_val, y_max_val) +
    scale_x_continuous(breaks = seq(-12, 2, by = 2), limits = c(-12, 2))+
    theme_bw() +
    theme(
      plot.title  = element_text(size = 13, face = "bold", hjust = 0.5),
      legend.position = "none",
      legend.title    = element_text(size = 10),
      legend.text     = element_text(size = 9)
    )
  
  plot_list[[prot]] <- p_prot
}

p_prot
# plot_list <- plot_list[protein_order]


# 4) Selected proteins (e.g., NEFL & EDA2R)
# plot_list_select <- plot_list[c("NEFL", "EDA2R")]
# select_plot <- wrap_plots(plot_list, ncol = 3)

library(cowplot)

spacer <- ggplot() + theme_void()

p_c_combined <- plot_grid(
  plot_list[[1]],
  spacer,
  plot_list[[2]],
  spacer,
  plot_list[[3]],
  ncol = 5,
  rel_widths = c(1, 0.05, 1, 0.05, 1),
  align = "h"
)

ggsave(file.path(fig5_out, "Fig_5C_UKB_pre_ALS_2yr.pdf"),
       plot  = p_c_combined,
       device = cairo_pdf,
       width = 12,
       height = 4)
ggsave(file.path(fig5_out, "Fig_5C_UKB_pre_ALS_2yr.png"),
       plot  = p_c_combined,
       width = 12,
       height = 4,
       dpi   = 300)
ggsave(file.path(fig5_out, "Fig_5C_UKB_pre_ALS_2yr.svg"),
       plot  = p_c_combined,
       device = svglite::svglite,
       width = 12,
       height = 4)
############################################################
# Fig5C v2 – With marginal histogram plots (ggExtra)
############################################################
library(ggExtra)

plot_list_ggextra <- list()
for (prot in protein_list) {
  y_min_val <- -2.0
  y_max_val <-  3.0

  confirm_df_prot <- confirm_df %>% filter(protein == !!prot)
  grid_overall_prot <- grid_overall %>%
    filter(protein == !!prot) %>%
    mutate(YrSinceOs = YrSinceDi + 2,
           lower = pmax(lower, y_min_val),
           upper = pmin(upper, y_max_val))

  # Significant direction points (pre-phenoconversion)
  sig_points_df <- grid_overall_prot %>%
    filter(!is.na(dir), YrSinceOs <= 0)

  # Contiguous significant runs for label annotations
  y_label_val_ind_increase <- -0.4
  y_label_val_ind_decrease <-  0.7

  r     <- rle(grid_overall_prot$sig)
  ends  <- cumsum(r$lengths)
  starts <- c(1, head(ends, -1) + 1)
  nonzero_runs_idx <- which(r$values != 0)

  if (length(nonzero_runs_idx) > 0) {
    run_starts_idx <- starts[nonzero_runs_idx]
    run_ends_idx   <- ends[nonzero_runs_idx]
    run_vals       <- r$values[nonzero_runs_idx]

    run_df <- tibble(
      start_idx  = run_starts_idx,
      end_idx    = run_ends_idx,
      start_year = grid_overall_prot$YrSinceOs[run_starts_idx],
      end_year   = grid_overall_prot$YrSinceOs[run_ends_idx],
      run_sign   = run_vals
    ) %>%
      mutate(
        dir         = ifelse(run_sign > 0, "increase", "decrease"),
        color_hex   = ifelse(run_sign > 0, "#D73027", "#4575B4"),
        label_start = sprintf("%.1f", start_year),
        label_end   = sprintf("%.1f", end_year),
        y_label_val_ind = ifelse(dir == "increase",
                                 y_label_val_ind_increase,
                                 y_label_val_ind_decrease)
      ) %>%
      mutate(label_end = ifelse(abs(end_year) < 0.1, "", label_end)) %>%
      filter(start_year <= 0) %>%
      mutate(label_end = ifelse(end_year >= 0.2, "", label_end))
  } else {
    run_df <- tibble(
      start_idx = integer(), end_idx = integer(),
      start_year = numeric(), end_year = numeric(),
      run_sign = integer(), dir = character(),
      label_start = character(), label_end = character(),
      color_hex = character(), y_label_val_ind = numeric()
    )
  }

  fill_map <- c(increase = "#D73027", decrease = "#4575B4")

  # ggMarginal needs aes(x, y) in the base ggplot call,
  # so we rebuild the scatter with all annotation layers
  p_base <- ggplot(confirm_df_prot, aes(x = YrSinceOs, y = delta)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_point(alpha = 0.6, size = 1, color = "grey10") +
    geom_ribbon(
      data = grid_overall_prot,
      aes(x = YrSinceOs, ymin = lower, ymax = upper),
      fill = "grey30", alpha = 0.2, inherit.aes = FALSE
    ) +
    geom_line(
      data = grid_overall_prot,
      aes(x = YrSinceOs, y = fit),
      color = "black", size = 1.6, inherit.aes = FALSE
    ) +
    geom_vline(xintercept = 0, color = "black") +
    # significant time points at y=0
    geom_point(
      data = sig_points_df,
      aes(x = YrSinceOs, y = 0, color = dir),
      size = 1.2, alpha = 0.9, show.legend = FALSE,
      inherit.aes = FALSE
    ) +
    scale_color_manual(values = fill_map, guide = "none") +
    scale_fill_manual(values = fill_map, name = "Direction")

  # Add significant interval labels (increase)
  run_df_inc <- run_df %>% filter(dir == "increase")
  if (nrow(run_df_inc) > 0) {
    p_base <- p_base +
      geom_label(
        data = run_df_inc,
        aes(x = start_year,
            y = y_label_val_ind - 0.2, label = label_start),
        fill = "white",
        color = run_df_inc$color_hex,
        size = 5, fontface = "bold", vjust = 0,
        label.size = 0.5, inherit.aes = FALSE
      )
    run_df_inc_end <- run_df_inc %>% filter(label_end != "")
    if (nrow(run_df_inc_end) > 0) {
      p_base <- p_base +
        geom_label(
          data = run_df_inc_end,
          aes(x = end_year + 1,
              y = y_label_val_ind - 0.2, label = label_end),
          fill = "white",
          color = run_df_inc_end$color_hex,
          size = 5, fontface = "bold", vjust = 0,
          label.size = 0.5, inherit.aes = FALSE
        )
    }
  }

  p_base <- p_base +
    labs(title = prot, x = "Years to Phenoconversion Proxy", y = "log2FC") +
    ylim(y_min_val, y_max_val) +
    scale_x_continuous(breaks = seq(-12, 2, by = 2), limits = c(-12, 2)) +
    theme_bw() +
    theme(
      plot.title = element_text(size = 13, face = "bold", hjust = 0.5),
      legend.position = "none"
    )

  p_marginal <- ggMarginal(
    p_base,
    margins = "y",
    type    = "histogram",
    fill    = "grey60",
    alpha   = 0.5,
    size    = 5
  )

  plot_list_ggextra[[prot]] <- p_marginal
}

p_c_combined_v2_ggextra <- plot_grid(
  plot_list_ggextra[[1]],
  plot_list_ggextra[[2]],
  plot_list_ggextra[[3]],
  ncol = 3
)

ggsave(file.path(fig5_out, "Fig_5C_UKB_pre_ALS_2yr_with_density.pdf"),
       plot  = p_c_combined_v2_ggextra,
       device = cairo_pdf,
       width = 15,
       height = 4)
ggsave(file.path(fig5_out, "Fig_5C_UKB_pre_ALS_2yr_with_density.png"),
       plot  = p_c_combined_v2_ggextra,
       width = 15,
       height = 4,
       dpi   = 300)
ggsave(file.path(fig5_out, "Fig_5C_UKB_pre_ALS_2yr_with_density.svg"),
       plot  = p_c_combined_v2_ggextra,
       device = svglite::svglite,
       width = 15,
       height = 4)

} # end Fig5C individual data guard


############################################################
# Fig5D – UKB phenoconversion prediction (tobit model)
############################################################

# Helper function: predicted vs observed plot with ±1-year band
predict_vs_observed_plot_tobit <- function(
    y_true,
    y_predict,
    vjust   = 0.75,
    title   = NULL,
    MAE,
    RMSE,
    cor_y
) {
  df <- data.frame(
    y_true    = y_true,
    y_predict = y_predict
  )
  
  # Category: earlier vs later onset than predicted
  df$category <- ifelse(
    df$y_predict < df$y_true,
    "Onset Eariler than Predicted",
    "Onset Later than Predicted"
  )
  
  color_map <- c(
    "Onset Later than Predicted"  = "blue",
    "Onset Eariler than Predicted" = "red"
  )
  
  df$residuals <- df$y_true - df$y_predict
  
  # Proportion of |residual| < 1 (not printed, but can be used if needed)
  prop_1 <- sum(abs(df$residuals) < 1, na.rm = TRUE) / nrow(df)
  
  # Common limits for x and y axes
  range_limits <- c(-12, 0)
  
  # Polygon for ±1-year band around y = x
  band_x <- seq(range_limits[1], range_limits[2], length.out = 100)
  band_df <- data.frame(
    x = c(band_x, rev(band_x)),
    y = c(band_x - 1, rev(band_x + 1))
  )
  
  # Clip band to plotting range
  band_df$y[band_df$y < range_limits[1]] <- range_limits[1]
  band_df$y[band_df$y > range_limits[2]] <- range_limits[2]
  
  p <- ggplot(df, aes(x = y_predict, y = y_true, color = category)) +
    geom_polygon(
      data  = band_df,
      aes(x = x, y = y),
      fill  = "grey80",
      alpha = 0.5,
      inherit.aes = FALSE
    ) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(values = color_map) +
    geom_line(
      data  = data.frame(x = band_x, y = band_x),
      aes(x = x, y = y),
      color = "black",
      linetype = "dashed",
      size   = 1,
      inherit.aes = FALSE
    ) +
    labs(
      title = title,
      x     = "Predicted Years to Phenoconversion",
      y     = "Observed Years to Phenoconversion",
      color = "Prediction Type"
    ) +
    theme_classic() +
    guides(color = guide_legend(nrow = 2)) +
    theme(
      plot.title   = element_text(hjust = 0.5, size = 20),
      axis.text.x  = element_text(size = 18),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.text.y  = element_text(size = 18),
      legend.text  = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.position = "none",
      plot.background  = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA)
    ) +
    
    scale_y_continuous(breaks = seq(-15, 0, by = 5), limits = c(-12, 0))+
    scale_x_continuous(breaks = seq(-15, 0, by = 5), limits = c(-12, 0))+
    
    # xlim(range_limits) +
    # ylim(range_limits) +
    annotate(
      "text",
      x = -12,
      y = 0,
      hjust = 0,
      vjust = vjust,
      label =sprintf(
        "MAE = %.3f\nRMSE = %.3f\nCor = %.3f", 
        MAE, RMSE, cor_y, prop_1 * 100
      ),
      size   = 5,
      color  = "black"
    )
  
  p
}

# 1) Read prediction output
fig5d_path <- here::here("data", "fig_data", "individual_data", "Fig5",
                          "result_individual_predictions_df.csv")
if (!file.exists(fig5d_path)) {
  cat("  WARNING: Individual data not found for Fig5D. Skipping.\n")
} else {
df_pred <- read.csv(fig5d_path)

# replace the Panel_1 with NEFL
df_pred <- df_pred %>%
  mutate(Panel = ifelse(Panel == "Panel_1", "NEFL", Panel))

# replace the Panel_data_ukb  with Panel_10
df_pred <- df_pred %>%
  mutate(Panel = ifelse(Panel == "Panel_data_ukb", "Panel_10", Panel))

protein_set_list <- c("Panel_15", "Panel_10","NEFL" )

# Display names for titles
display_names <- c(
  "Panel_15" = "15-Protein Panel\n(UKB subset of 19 proteins)",
  "Panel_10" = "10-Protein Panel\n(UKB subset of 16 proteins)",
  "NEFL"     = "NEFL Only"
)

# shift observed Y to "years to phenoconversion"
df_pred$Y <- df_pred$YrSinceDi + 2
# keep only Y <= 0 (pre-phenoconversion)
df_pred <- df_pred %>% filter(Y <= 0)

plots_d_2yr <- list()

for (protein_set in protein_set_list) {
  
  df_cor <- df_pred %>%
    filter( Panel == protein_set ) %>%
    mutate(y_pred = Y_pred,
           y=Y,
           y_pred  = pmin(y_pred, 0),
           resid   = y - y_pred
    ) %>%
    filter(!is.na(y) & !is.na(y_pred))
  
  MAE   <- mean(abs(df_cor$resid), na.rm = TRUE)
  RMSE  <- sqrt(mean(df_cor$resid^2, na.rm = TRUE))
  cor_y <- cor(df_cor$y, df_cor$y_pred, use = "complete.obs")
  
  p <- predict_vs_observed_plot_tobit(
    y_true  = df_cor$y,
    y_predict = df_cor$y_pred,
    title   = display_names[[protein_set]],
    MAE     = round(MAE, 3),
    RMSE    = round(RMSE, 3),
    cor_y   = round(cor_y, 3)
  ) +
    coord_fixed()
  
  plots_d_2yr[[protein_set]] <- p
}

p_Fig5D <- wrap_plots(plots_d_2yr, ncol = 3)


# Save Fig5D
ggsave(file.path(fig5_out, "Fig_5D_UKB_DA_model.pdf"),
       plot  = p_Fig5D,
       device = cairo_pdf,
       width = 20,
       height = 6.5)
ggsave(file.path(fig5_out, "Fig_5D_UKB_DA_model.png"),
       plot  = p_Fig5D,
       width = 20,
       height = 6.5,
       dpi   = 300)
ggsave(file.path(fig5_out, "Fig_5D_UKB_DA_model.svg"),
       plot  = p_Fig5D,
       device = svglite::svglite,
       width = 20,
       height = 6.5)
} # end Fig5D individual data guard
