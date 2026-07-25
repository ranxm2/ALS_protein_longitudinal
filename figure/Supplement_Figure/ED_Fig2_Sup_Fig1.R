# ==============================================================================
# ED_Fig2_Sup_Fig1.R — Extended Data Figure 2 and Supplementary Figure 1:
#                      Protein Trajectory Atlases and Temporal Heatmaps
# ==============================================================================
# Author : Ximing Ran
# Date   : 2026-02-22
#
# Description:
#   Generates three supplementary figures covering GAM-fitted protein
#   trajectory panels and temporal expression heatmaps:
#
#   ED_Fig2    - Trajectory atlas for 92 proteins with pre-onset trajectory
#               (all time points: converters + affected carriers combined)
#   Sup_Fig1A - Temporal heatmap of mean log2FC across 1-year bins for the
#               73 conversion-period proteins, grouped by PPI module
#   Sup_Fig1B - Trajectory atlas for 73 proteins with pre-onset trajectory
#               (converters only, all time points)
#   Sup_Fig1C - Temporal heatmap for the 52 pre-conversion proteins
#   Sup_Fig1D - Trajectory atlas for 52 proteins with pre-onset trajectory
#               (converters, pre-onset only)
#
#   Each trajectory panel shows individual subject lines (grey), GAM-smoothed
#   mean curve (black), 95% CI ribbon, and coloured rug marks indicating
#   time windows of statistical significance.
#
# Input data:
#   data/fig_data/individual_data/Supplement_Figure/  (optional — trajectory
#     atlases ED_Fig2/Sup_Fig1B/1D skipped if absent)
#     miami_protein_trajectory_gam_delta_confirm_df_sorted_all_proteins_all_time_137.csv
#     miami_protein_trajectory_gam_delta_confirm_df_sorted_all_proteins_convert_137.csv
#     miami_protein_trajectory_gam_delta_confirm_df_sorted_all_proteins_pre_137.csv
#                                        — per-subject delta protein values
#     predictions_Convert_Affected_AllTime.csv
#     predictions_Convert_Only_AllTime.csv
#     predictions_Convert_Only_PreOnset.csv
#                                        — GAM fitted values and CIs
#   data/fig_data/summary_data/Supplement_Figure/
#     combined_summary_all_subsets.csv   — FDR-corrected significance per
#                                          protein and time subset
#     heat_1yr_wide_convert.csv          — 1-yr binned mean log2FC (converters)
#     heat_1yr_wide_pre_convert.csv      — 1-yr binned mean log2FC (pre-onset)
#   data/fig_data/summary_data/Fig2/
#     Heatmap_PPI_cluster_mask_with_cluster_annotation_row.csv
#                                        — PPI module labels for heatmap rows
#
# Output (./results/):
#   ED_Fig2_Miami_ALS_all_92.{pdf,png}
#   Sup_Fig1A_Temporal_Heatmap_Convert.{pdf,png}
#   Sup_Fig1B_Miami_ALS_convert_all_73.{pdf,png}
#   Sup_Fig1C_Temporal_Heatmap_pre_Convert.{pdf,png}
#   Sup_Fig1D_Miami_ALS_convert_pre_52.{pdf,png}
#
# Prerequisite:
#   Run analysis/02-Protein_Trajectories/ Rmd notebooks first.
# ==============================================================================

# Setup ------------------------------------------------------------------------
library(tidyverse)
library(here)
library(patchwork)
library(cowplot)
library(svglite)

output_dir <- here::here("figure", "Supplement_Figure", "results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Set seed for reproducibility
set.seed(2025)

# Define color scheme for increase/decrease
FILL_MAP <- c(increase = "#D73027", decrease = "#4575B4")

# Define plot limits
Y_MIN <- -2.0
Y_MAX <- 3.0

# ==============================================================================
# HELPER FUNCTIONS
# ==============================================================================

#' Load protein trajectory data
#'
#' @param time_period Character: "all_time", "convert", or "pre"
#' @return List containing confirm_df and grid_overall dataframes
load_trajectory_data <- function(time_period) {
  file_suffix <- switch(time_period,
                        "all_time" = "all_time_137",
                        "convert" = "convert_137",
                        "pre" = "pre_137")
  confirm_path <- here("data","fig_data", "individual_data", "Supplement_Figure",
                       paste0("miami_protein_trajectory_gam_delta_confirm_df_sorted_all_proteins_", file_suffix, ".csv"))
  grid_suffix <- switch(time_period,
                        "all_time" = "Convert_Affected_AllTime",
                        "convert" = "Convert_Only_AllTime",
                        "pre" = "Convert_Only_PreOnset")
  grid_path <- here("data", "fig_data","individual_data", "Supplement_Figure",
                    paste0("predictions_", grid_suffix, ".csv"))
  if (!file.exists(confirm_path) || !file.exists(grid_path)) {
    cat(sprintf("  WARNING: Individual data not found for %s. Skipping.\n", time_period))
    return(NULL)
  }
  confirm_df <- read.csv(confirm_path)
  grid_overall <- read.csv(grid_path)
  list(confirm_df = confirm_df, grid_overall = grid_overall)
}

#' Prepare significant region data for labeling
#'
#' @param grid_overall_prot Filtered grid data for single protein
#' @return Tibble with significant region information
prepare_sig_regions <- function(grid_overall_prot) {
  # Identify contiguous significant regions using run-length encoding
  r <- rle(grid_overall_prot$sig)
  ends <- cumsum(r$lengths)
  starts <- c(1, head(ends, -1) + 1)
  nonzero_runs_idx <- which(r$values != 0)
  
  if (length(nonzero_runs_idx) == 0) {
    return(tibble(start_idx = integer(), end_idx = integer(),
                  start_year = numeric(), end_year = numeric(),
                  run_sign = integer(), dir = character(),
                  label_start = character(), label_end = character(),
                  color_hex = character(), y_label = numeric(),
                  y_label_val_ind = numeric(),
                  y_label_end = numeric(), show_label = logical()))
  }
  
  run_starts_idx <- starts[nonzero_runs_idx]
  run_ends_idx <- ends[nonzero_runs_idx]
  run_vals <- r$values[nonzero_runs_idx]
  
  run_df <- tibble(
    y_label = 2.9,
    start_idx = run_starts_idx,
    end_idx = run_ends_idx,
    start_year = grid_overall_prot$YrSinceOs[run_starts_idx],
    end_year = grid_overall_prot$YrSinceOs[run_ends_idx],
    run_sign = run_vals
  ) %>%
    mutate(
      dir = ifelse(run_sign > 0, "increase", "decrease"),
      label_start = sprintf("%.1f", start_year),
      label_end = sprintf("%.1f", end_year),
      color_hex = ifelse(run_sign > 0, "#D73027", "#4575B4"),
      y_label_val_ind = -1.8
    ) %>%
    # Remove labels for endpoints near zero or after conversion
    mutate(
      label_end = ifelse(abs(end_year) < 0.1 | end_year > 0, "", label_end)
    ) %>%
    filter(start_year <= 0)
  
  # Adjust intervals shorter than 3 years
  run_df <- run_df %>%
    mutate(
      end_year = ifelse((end_year - start_year) < 3 & end_year <= 0,
                        start_year + 2, end_year),
      label_end = ifelse((end_year - start_year) < 3 & end_year <= 0,
                         sprintf("%.1f", end_year), label_end)
    ) %>%
    arrange(start_year) %>%
    mutate(
      start_year = ifelse(
        lag((end_year - start_year) < 3 &
              (lead(start_year, default = first(start_year)) - end_year) <= 2,
            default = FALSE),
        lag(start_year) + 5,
        start_year
      )
    ) %>%
    # Shift end label further from y=0 when close to start label
    mutate(
      y_label_end = ifelse(
        abs(end_year - start_year) < 4 & label_end != "",
        ifelse(dir == "increase", y_label_val_ind - 0.5, y_label_val_ind + 0.5),
        y_label_val_ind
      )
    )
  # Only label the interval closest to onset (largest start_year, i.e. closest to 0)
  run_df$show_label <- FALSE
  if (nrow(run_df) > 0) {
    closest_idx <- which.max(run_df$start_year)
    run_df$show_label[closest_idx] <- TRUE
  }
  return(run_df)
}

#' Create trajectory plot for a single protein
#'
#' @param protein Character: protein name
#' @param confirm_df Individual trajectory data
#' @param grid_overall Overall fitted trajectory data
#' @return ggplot object
create_protein_plot <- function(protein, confirm_df, grid_overall) {
  # Filter data for this protein
  confirm_df_prot <- confirm_df %>% filter(protein == !!protein)
  grid_overall_prot <- grid_overall %>% filter(protein == !!protein)
  confirm_df_sorted <- confirm_df_prot %>% arrange(UIDx, YrSinceOs)
  
  # Get significant points (pre-conversion only)
  sig_points_df <- grid_overall_prot %>%
    filter(!is.na(dir), YrSinceOs <= 0)
  
  # Clip confidence intervals to plot limits
  grid_overall_prot <- grid_overall_prot %>%
    mutate(
      lower = ifelse(lower < Y_MIN, Y_MIN, lower),
      upper = ifelse(upper > Y_MAX, Y_MAX, upper)
    )
  
  # Prepare significant region labels
  run_df <- prepare_sig_regions(grid_overall_prot)
  
  # Create plot
  p <- ggplot() +
    # Individual trajectories
    geom_line(data = confirm_df_sorted,
              aes(x = YrSinceOs, y = delta, group = UIDx),
              color = "grey10", alpha = 0.4, size = 0.3) +
    geom_point(data = confirm_df_prot,
               aes(x = YrSinceOs, y = delta),
               alpha = 0.4, size = 1, color = "grey10") +
    
    # Overall fitted trajectory with confidence interval
    geom_ribbon(data = grid_overall_prot,
                aes(x = YrSinceOs, ymin = lower, ymax = upper),
                fill = "grey70", alpha = 0.3) +
    geom_line(data = grid_overall_prot,
              aes(x = YrSinceOs, y = fit),
              color = "black", size = 1.6) +
    
    # Reference lines
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_vline(xintercept = 0, color = "black") +
    
    # Significant points
    geom_point(data = sig_points_df,
               aes(x = YrSinceOs, y = 0, color = dir),
               size = 1.2, alpha = 0.9, show.legend = FALSE) +
    
    # Time period label — only start of interval closest to onset
    geom_label(data = run_df %>% filter(show_label),
               aes(x = start_year, y = y_label_val_ind - 0.2, label = label_start),
               fill = "white", color = (run_df %>% filter(show_label))$color_hex,
               size = 5, fontface = "bold", vjust = 0, label.size = 0.5) +
    
    # Scales and labels
    scale_fill_manual(values = FILL_MAP, name = "Direction") +
    scale_color_manual(values = FILL_MAP, guide = FALSE) +
    labs(title = protein,
         x = "Years to Phenoconversion",
         y = "log2FC") +
    ylim(Y_MIN, Y_MAX) +
    
    # Theme
    theme_bw() +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
      axis.title.x = element_text(size = 18),
      axis.title.y = element_text(size = 18),
      axis.text.x = element_text(size = 16),
      axis.text.y = element_text(size = 16),
      legend.position = "none",
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA)
    )
  
  return(p)
}

#' Generate plots for all proteins
#'
#' @param protein_list Character vector of protein names
#' @param confirm_df Individual trajectory data
#' @param grid_overall Overall fitted trajectory data
#' @return Named list of ggplot objects
generate_all_plots <- function(protein_list, confirm_df, grid_overall) {
  plot_list <- list()
  
  for (protein in protein_list) {
    plot_list[[protein]] <- create_protein_plot(protein, confirm_df, grid_overall)
  }
  
  return(plot_list)
}

# ==============================================================================
# ANALYSIS 1: ALL TIME POINTS (92 proteins with pre-onset trajectory)
# ==============================================================================
message("Generating plots for all time points...")

# Load data
data_all <- load_trajectory_data("all_time")

if (is.null(data_all)) {
  message("  Skipping ED_Fig2 (individual data not available)")
} else {
  # Get list of significant proteins with pre-onset trajectory change
  df_all <- read.csv(here("data", "fig_data","summary_data", "Supplement_Figure",
                          "combined_summary_all_subsets.csv")) %>%
    filter(Convert.Affected_AllTime_fdr_YrSinceOs < 0.05,
           !is.na(Convert.Affected_AllTime_nearest_interval_start) &
             Convert.Affected_AllTime_nearest_interval_start < 0)
  
  protein_list_all <- sort(df_all$protein)
  message(sprintf("  %d proteins pass FDR < 0.05 & nearest_interval_start < 0", length(protein_list_all)))
  
  # Generate plots
  plots_all <- generate_all_plots(protein_list_all,
                                  data_all$confirm_df,
                                  data_all$grid_overall)
  
  # Save combined plot
  final_plot_all <- plot_grid(plotlist = plots_all, ncol = 10, align = "v")
  
  ggsave(file.path(output_dir, "ED_Fig2_Miami_ALS_all_92.pdf"),
         plot = final_plot_all,
         device = cairo_pdf,
         width = 40,
         height = 40)
  
  ggsave(file.path(output_dir, "ED_Fig2_Miami_ALS_all_92.png"),
         plot = final_plot_all,
         width = 40,
         height = 40,
         dpi = 300)
  ggsave(file.path(output_dir, "ED_Fig2_Miami_ALS_all_92.svg"),
         plot = final_plot_all,
         device = svglite::svglite,
         width = 40,
         height = 40)

  message("All time plots complete!")
} # end ED_Fig2 guard

# ==============================================================================
# ANALYSIS 2: CONVERSION PERIOD (73 proteins with pre-onset trajectory)
# ==============================================================================
message("Generating plots for conversion period...")

# Load data
data_convert <- load_trajectory_data("convert")

if (is.null(data_convert)) {
  message("  Skipping Sup_Fig1B (individual data not available)")
} else {
  protein_list_convert <- unique(data_convert$confirm_df$protein)
  
  # Generate plots
  plots_convert <- generate_all_plots(protein_list_convert,
                                      data_convert$confirm_df,
                                      data_convert$grid_overall)
  
  protein_order_73 <- read.csv(here("data", "fig_data","summary_data", "Supplement_Figure",
                                    "combined_summary_all_subsets.csv")) %>%
    filter(Convert_AllTime_fdr_YrSinceOs < 0.05,
           !is.na(Convert_AllTime_nearest_interval_start) &
             Convert_AllTime_nearest_interval_start < 0) %>%
    pull(protein) %>% sort()
  message(sprintf("  %d proteins pass FDR < 0.05 & nearest_interval_start < 0", length(protein_order_73)))
  
  plots_convert_sorted <- plots_convert[protein_order_73]
  
  # Save combined plot
  final_plot_convert <- plot_grid(plotlist = plots_convert_sorted,
                                  ncol = 10, align = "v")
  
  ggsave(file.path(output_dir, "Sup_Fig1B_Miami_ALS_convert_all_73.pdf"),
         plot = final_plot_convert,
         device = cairo_pdf,
         width = 40,
         height = 32)
  
  ggsave(file.path(output_dir, "Sup_Fig1B_Miami_ALS_convert_all_73.png"),
         plot = final_plot_convert,
         width = 40,
         height = 32, dpi = 300)
  ggsave(file.path(output_dir, "Sup_Fig1B_Miami_ALS_convert_all_73.svg"),
         plot = final_plot_convert,
         device = svglite::svglite,
         width = 40,
         height = 32)

  message("Conversion period plots complete!")
} # end Sup_Fig1B guard

# ==============================================================================
# ANALYSIS 3: PRE-CONVERSION PERIOD (proteins with pre-onset trajectory)
# ==============================================================================
message("Generating plots for pre-conversion period...")

# Load data
data_pre <- load_trajectory_data("pre")

if (is.null(data_pre)) {
  message("  Skipping Sup_Fig1D (individual data not available)")
} else {
  protein_list_pre <- unique(data_pre$confirm_df$protein)
  
  # Generate plots
  plots_pre <- generate_all_plots(protein_list_pre,
                                  data_pre$confirm_df,
                                  data_pre$grid_overall)
  
  protein_order_pre <- read.csv(here("data", "fig_data","summary_data", "Supplement_Figure",
                                     "combined_summary_all_subsets.csv")) %>%
    filter(Convert_PreOnset_fdr_YrSinceOs < 0.05,
           !is.na(Convert_PreOnset_nearest_interval_start) &
             Convert_PreOnset_nearest_interval_start < 0) %>%
    pull(protein) %>% sort()
  message(sprintf("  %d proteins pass FDR < 0.05 & nearest_interval_start < 0", length(protein_order_pre)))
  
  # Sort plots
  plots_pre_sorted <- plots_pre[protein_order_pre]
  
  # Save combined plot
  final_plot_pre <- plot_grid(plotlist = plots_pre_sorted,
                              ncol = 10, align = "v")
  
  ggsave(file.path(output_dir, sprintf("Sup_Fig1D_Miami_ALS_convert_pre_%d.pdf", length(protein_order_pre))),
         plot = final_plot_pre,
         device = cairo_pdf,
         width = 40,
         height = 24)
  
  ggsave(file.path(output_dir, sprintf("Sup_Fig1D_Miami_ALS_convert_pre_%d.png", length(protein_order_pre))),
         plot = final_plot_pre,
         width = 40,
         height = 24,
         dpi = 300)
  ggsave(file.path(output_dir, sprintf("Sup_Fig1D_Miami_ALS_convert_pre_%d.svg", length(protein_order_pre))),
         plot = final_plot_pre,
         device = svglite::svglite,
         width = 40,
         height = 24)

  message("Pre-conversion period plots complete!")
} # end Sup_Fig1D guard

# ==============================================================================
# Supplement 3: Temporal Protein Expression Analysis in ALS
# ==============================================================================
# Author: Ximing Ran
# Date: 2026-02-22
# Description:
#   Figure 3A - Temporal heatmap of all convert
#   Figure 3B - Temporal heatmap of pre-convert
# ==============================================================================

# ------------------------------------------------------------------------------
# Setup and Configuration
# ------------------------------------------------------------------------------
# Load required libraries
library(tidyverse)
library(here)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(ggplotify)
library(cowplot)
library(patchwork)

# Set random seed for reproducibility (if applicable)
set.seed(2025)

library(gtable)
library(grid)

# ------------------------------------------------------------------------------
# Helper: stack pheatmap's annotation (Module) legend ABOVE the color legend
# ------------------------------------------------------------------------------
# pheatmap places the continuous color scale and the annotation legend in two
# SEPARATE columns (side by side), each anchored to the top of its cell. To get
# the Module legend stacked ABOVE the color scale, both grobs are merged into a
# single column via a nested 3-row gtable (Module on top, gap, color scale).
#
# Tunables (all in bigpts, independent of the saved figure size):
#   module_legend_height — vertical space reserved for the Module legend on top.
#                          Increase if you add modules so the color scale does
#                          not overlap the last category; decrease to tighten.
#   gap                  — spacing between the Module legend and the color scale.
#   color_legend_height  — height of the continuous color bar (pheatmap default
#                          caps it at 150 bigpts).
stack_legends_module_on_top <- function(ph,
                                        module_legend_height = 150,
                                        gap = 18,
                                        color_legend_height = 150) {
  gt <- ph$gtable
  li <- which(gt$layout$name == "legend")            # continuous color scale
  ai <- which(gt$layout$name == "annotation_legend") # Module legend
  if (length(li) == 0 || length(ai) == 0) return(gt)
  
  leg <- gt$grobs[[li]]
  alg <- gt$grobs[[ai]]
  
  # reuse the column/row span of the two original legend cells
  l_left  <- min(gt$layout$l[c(li, ai)]); r_right <- max(gt$layout$r[c(li, ai)])
  t_top   <- min(gt$layout$t[c(li, ai)]); b_bot   <- max(gt$layout$b[c(li, ai)])
  
  # 3-row stack: Module on top, gap, color scale below
  # (each grob is top-anchored within its own row)
  stack <- gtable(widths  = unit(1, "npc"),
                  heights = unit.c(unit(module_legend_height, "bigpts"),
                                   unit(gap, "bigpts"),
                                   unit(color_legend_height, "bigpts")))
  stack <- gtable_add_grob(stack, alg, t = 1, l = 1, clip = "off", name = "annotation_legend")
  stack <- gtable_add_grob(stack, leg, t = 3, l = 1, clip = "off", name = "legend")
  
  # remove the two originals (drop higher index first) and insert the combined stack
  for (d in sort(c(li, ai), decreasing = TRUE)) {
    gt$grobs[[d]] <- NULL
    gt$layout     <- gt$layout[-d, ]
  }
  gtable_add_grob(gt, stack, t = t_top, l = l_left, b = b_bot, r = r_right,
                  clip = "off", name = "legend_stack")
}

# ==============================================================================
# FIGURE 3A: TEMPORAL HEATMAP (CONVERTERS)
# ==============================================================================

# ------------------------------------------------------------------------------
# Load and Prepare Data for Figure 3A
# ------------------------------------------------------------------------------
# Load temporal expression data
df_3a <- read.csv(
  here::here("data", "fig_data",  "summary_data", "Supplement_Figure", "heat_1yr_wide_convert.csv")
)
cat(sprintf("  Loaded data: %d proteins x %d columns\n", nrow(df_3a), ncol(df_3a)))

# ------------------------------------------------------------------------------
# Fix Column Names for Time Bins
# ------------------------------------------------------------------------------
# Identify time bin columns (encoded by make.names, e.g. X.15.5, X0.5)
# Meta columns are non-numeric; time bin columns encode numeric bin center
cols <- colnames(df_3a)
meta_cols <- c("protein", "start_year", "dir", "P_value", "Adj_P_value",
               "nearest_interval_start",
               "start_idx", "end_idx", "run_sign", "color_hex", "cluster")
time_bin_idx <- which(!cols %in% meta_cols)
time_bin_idx <- time_bin_idx[time_bin_idx > 1 | cols[1] != "protein"]

# Decode make.names encoding back to signed numeric strings
for (i in time_bin_idx) {
  s <- gsub("^X", "", cols[i])
  if (grepl("^\\.", s)) s <- sub("^\\.", "-", s)
  cols[i] <- s
}
colnames(df_3a) <- cols

# Order time bin columns chronologically
time_bin_cols <- cols[time_bin_idx]
time_bin_numeric <- as.numeric(time_bin_cols)
time_bin_cols <- time_bin_cols[order(time_bin_numeric)]

# ------------------------------------------------------------------------------
# Load Module Annotations
# ------------------------------------------------------------------------------
# Load cluster/module annotations for proteins
annotation_row <- read.csv(
  here::here("data", "fig_data",  "summary_data", "Fig2",
             "Heatmap_PPI_cluster_mask_with_cluster_annotation_row.csv"),
  row.names = 1
) %>%
  select(cluster) %>%
  rownames_to_column("protein")

# Join annotations to main data
df_3a <- left_join(df_3a, annotation_row, by = "protein")

# ------------------------------------------------------------------------------
# Define Module Mappings and Colors
# ------------------------------------------------------------------------------
# Map cluster IDs to module names
module_map <- c(
  `1` = "Skeletal Muscle",
  `2` = "TNF-Mediated Signaling",
  `3` = "Skeletal Muscle",
  `5` = "ECM & ECM Interaction",
  `6` = "ECM & ECM Interaction",
  `8` = "Regeneration (incl. Neurofilament)"
)

# Define module colors
module_color_map <- c(
  "Skeletal Muscle" = "#009E73",
  "TNF-Mediated Signaling" = "#E69F00",
  "Regeneration (incl. Neurofilament)" = "#56B4E9",
  "ECM & ECM Interaction" = "#CC79A7",
  "Others" = "#d3d3d3"
)

# Define module order for plotting
module_levels <- c(
  "Skeletal Muscle",
  "TNF-Mediated Signaling",
  "Regeneration (incl. Neurofilament)",
  "ECM & ECM Interaction",
  "Others"
)

# ------------------------------------------------------------------------------
# Filter and Assign Modules
# ------------------------------------------------------------------------------
# Filter to significant proteins with pre-onset trajectory change
# (adjusted p-value < 0.05 AND nearest significant interval starts before phenoconversion)
df_3a <- df_3a %>%
  filter(Adj_P_value < 0.05,
         !is.na(nearest_interval_start) & nearest_interval_start < 0)
cat(sprintf("  Filtered to %d proteins (adj. p < 0.05 & nearest_interval_start < 0)\n", nrow(df_3a)))

# Assign modules based on cluster membership
df_3a$Module <- ifelse(
  df_3a$cluster %in% c(1, 2, 3, 5, 6, 8),
  unname(module_map[as.character(df_3a$cluster)]),
  "Others"
)

# Convert Module to ordered factor
df_3a <- df_3a %>%
  mutate(Module = factor(Module, levels = module_levels, ordered = TRUE))

# ------------------------------------------------------------------------------
# Build Expression Matrix
# ------------------------------------------------------------------------------
# Create numeric matrix of expression values
mat_display <- df_3a %>%
  select(all_of(time_bin_cols)) %>%
  as.matrix()
rownames(mat_display) <- df_3a$protein
mode(mat_display) <- "numeric"

# ------------------------------------------------------------------------------
# Impute Missing Values
# ------------------------------------------------------------------------------
# Function to impute missing values with row mean
impute_row_mean <- function(x) {
  if (all(is.na(x))) return(rep(0, length(x)))
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
}

# Apply imputation row-wise
mat_imputed <- t(apply(mat_display, 1, impute_row_mean))

# ------------------------------------------------------------------------------
# Calculate Row Ordering Metrics
# ------------------------------------------------------------------------------
# Find position of first non-NA value for each protein
first_non_na_pos <- apply(mat_display, 1, function(r) {
  idx <- which(!is.na(r))
  if (length(idx) == 0) Inf else idx[1]
})

# Calculate mean expression for each protein (using imputed values)
row_mean_imputed <- rowMeans(mat_imputed, na.rm = TRUE)

# ------------------------------------------------------------------------------
# Order Data by Module and Expression Pattern
# ------------------------------------------------------------------------------
# Order proteins by:
# 1. Module (factor order)
# 2. First appearance time
# 3. Mean expression (descending)
df_ordered <- df_3a %>%
  mutate(
    first_non_na_pos = first_non_na_pos[match(protein, rownames(mat_display))],
    row_mean_imputed = row_mean_imputed[match(protein, rownames(mat_display))]
  ) %>%
  arrange(Module, first_non_na_pos, desc(row_mean_imputed))

# Reorder display matrix to match protein order
mat_display_ord <- mat_display[df_ordered$protein, , drop = FALSE]

# ------------------------------------------------------------------------------
# Prepare Row Annotations
# ------------------------------------------------------------------------------
# Create row annotation data frame
annotation_row_plot <- data.frame(Module = df_ordered$Module)
rownames(annotation_row_plot) <- df_ordered$protein

# Calculate gap positions between modules
mod_rle <- rle(as.character(df_ordered$Module))
gap_indices <- cumsum(mod_rle$lengths)
gap_indices <- gap_indices[gap_indices < nrow(df_ordered)]

# ------------------------------------------------------------------------------
# Prepare Annotation Colors
# ------------------------------------------------------------------------------
# Create color mapping for modules
module_to_color <- module_color_map
annotation_colors <- list(Module = module_to_color)

# ------------------------------------------------------------------------------
# Prepare Heatmap Color Scale
# ------------------------------------------------------------------------------
# Calculate symmetric color scale based on max absolute value
max_abs <- max(abs(mat_imputed[1:20])+4, na.rm = TRUE)
breaks <- seq(-max_abs, max_abs, length.out = 101)
colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

# ------------------------------------------------------------------------------
# Prepare Display Numbers
# ------------------------------------------------------------------------------
# Create matrix of numbers to display in cells
number_matrix <- mat_display_ord
number_matrix[!is.na(number_matrix)] <-
  sprintf("%.2f", number_matrix[!is.na(number_matrix)])
number_matrix[is.na(number_matrix)] <- ""

# ------------------------------------------------------------------------------
# Format Column Names
# ------------------------------------------------------------------------------
# Create readable time bin labels from numeric column names
bin_vals <- as.numeric(colnames(mat_display_ord))
colnames(mat_display_ord) <- ifelse(
  bin_vals < 0,
  paste0("-Y", abs(bin_vals - 0.5)),
  paste0("Y", bin_vals + 0.5)
)

# ------------------------------------------------------------------------------
# Mark Key Proteins with Asterisks
# ------------------------------------------------------------------------------
# # List of key proteins to highlight (used in Figure 2B)
# key_proteins <- c("CA3", "EDA2R", "NEFL", "CALCA")
#
# # Add ** suffix to key protein names in all relevant objects
# for (obj_name in c("mat_display_ord", "number_matrix", "annotation_row_plot")) {
#   obj <- get(obj_name)
#   rnames <- rownames(obj)
#   mask <- rnames %in% key_proteins
#   rownames(obj)[mask] <- paste0(rnames[mask], "**")
#   assign(obj_name, obj)
# }

# ------------------------------------------------------------------------------
# Create Heatmap (Figure 3A)
# ------------------------------------------------------------------------------
ph_2a <- pheatmap(
  mat_display_ord,
  color = colors,
  breaks = breaks,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_row = annotation_row_plot,
  annotation_colors = annotation_colors,
  gaps_row = gap_indices,
  gaps_col = sum(bin_vals < 0),  # Separate pre- and post-diagnosis periods
  display_numbers = number_matrix,
  fontsize_number = 12,
  fontsize_row = 12,
  annotation_names_row = FALSE,
  fontsize_col = 12,
  fontsize = 12,
  angle_col = 0,
  cellwidth = 35,
  border_color = "black",
  number_color = "black",
  na_col = "white"
)

# Stack the Module legend ABOVE the color legend, then convert to ggplot
gt_3a <- stack_legends_module_on_top(ph_2a)
p_3a <- ggplotify::as.ggplot(gt_3a) +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

# Save Figure 3A in multiple formats
ggsave(
  file.path(output_dir, "Sup_Fig1A_Temporal_Heatmap_Convert.pdf"),
  plot = p_3a,
  device = cairo_pdf,
  width = 24,
  height = 15,
  units = "in"
)

ggsave(
  file.path(output_dir, "Sup_Fig1A_Temporal_Heatmap_Convert.png"),
  plot = p_3a,
  width = 24,
  height = 15,
  units = "in",
  dpi = 300
)
ggsave(
  file.path(output_dir, "Sup_Fig1A_Temporal_Heatmap_Convert.svg"),
  plot = p_3a,
  device = svglite::svglite,
  width = 24,
  height = 15,
  units = "in"
)

# ==============================================================================
# FIGURE 3C: TEMPORAL HEATMAP (PRE-CONVERTERS)
# ==============================================================================

# ------------------------------------------------------------------------------
# Load and Prepare Data for Figure 3C
# ------------------------------------------------------------------------------
# Load temporal expression data
df_3b <- read.csv(
  here::here("data", "fig_data",  "summary_data", "Supplement_Figure", "heat_1yr_wide_pre_convert.csv")
)
cat(sprintf("  Loaded data: %d proteins x %d columns\n", nrow(df_3b), ncol(df_3b)))

# ------------------------------------------------------------------------------
# Fix Column Names for Time Bins
# ------------------------------------------------------------------------------
# Identify time bin columns (same approach as convert section)
cols <- colnames(df_3b)
meta_cols <- c("protein", "start_year", "dir", "P_value", "Adj_P_value",
               "nearest_interval_start",
               "start_idx", "end_idx", "run_sign", "color_hex", "cluster")
time_bin_idx <- which(!cols %in% meta_cols)
time_bin_idx <- time_bin_idx[time_bin_idx > 1 | cols[1] != "protein"]

# Decode make.names encoding back to signed numeric strings
for (i in time_bin_idx) {
  s <- gsub("^X", "", cols[i])
  if (grepl("^\\.", s)) s <- sub("^\\.", "-", s)
  cols[i] <- s
}
colnames(df_3b) <- cols

# Order time bin columns chronologically
time_bin_cols <- cols[time_bin_idx]
time_bin_numeric <- as.numeric(time_bin_cols)
time_bin_cols <- time_bin_cols[order(time_bin_numeric)]

# ------------------------------------------------------------------------------
# Load Module Annotations
# ------------------------------------------------------------------------------
# Load cluster/module annotations for proteins
annotation_row <- read.csv(
  here::here("data", "fig_data",  "summary_data", "Fig2",
             "Heatmap_PPI_cluster_mask_with_cluster_annotation_row.csv"),
  row.names = 1
) %>%
  select(cluster) %>%
  rownames_to_column("protein")

# Join annotations to main data
df_3b <- left_join(df_3b, annotation_row, by = "protein")

# ------------------------------------------------------------------------------
# Define Module Mappings and Colors
# ------------------------------------------------------------------------------
# Map cluster IDs to module names
module_map <- c(
  `1` = "Skeletal Muscle",
  `2` = "TNF-Mediated Signaling",
  `3` = "Skeletal Muscle",
  `5` = "ECM & ECM Interaction",
  `6` = "ECM & ECM Interaction",
  `8` = "Regeneration (incl. Neurofilament)"
)

# Define module colors
module_color_map <- c(
  "Skeletal Muscle" = "#009E73",
  "TNF-Mediated Signaling" = "#E69F00",
  "Regeneration (incl. Neurofilament)" = "#56B4E9",
  "ECM & ECM Interaction" = "#CC79A7",
  "Others" = "#d3d3d3"
)

# Define module order for plotting
module_levels <- c(
  "Skeletal Muscle",
  "TNF-Mediated Signaling",
  "Regeneration (incl. Neurofilament)",
  "ECM & ECM Interaction",
  "Others"
)

# ------------------------------------------------------------------------------
# Filter and Assign Modules
# ------------------------------------------------------------------------------
# Filter to significant proteins with pre-onset trajectory change
df_3b <- df_3b %>%
  filter(Adj_P_value < 0.05,
         !is.na(nearest_interval_start) & nearest_interval_start < 0)
cat(sprintf("  Filtered to %d proteins (adj. p < 0.05 & nearest_interval_start < 0)\n", nrow(df_3b)))

# Assign modules based on cluster membership
df_3b$Module <- ifelse(
  df_3b$cluster %in% c(1, 2, 3, 5, 6, 8),
  unname(module_map[as.character(df_3b$cluster)]),
  "Others"
)

# Convert Module to ordered factor
df_3b <- df_3b %>%
  mutate(Module = factor(Module, levels = module_levels, ordered = TRUE))

# ------------------------------------------------------------------------------
# Build Expression Matrix
# ------------------------------------------------------------------------------
# Create numeric matrix of expression values
mat_display <- df_3b %>%
  select(all_of(time_bin_cols)) %>%
  as.matrix()
rownames(mat_display) <- df_3b$protein
mode(mat_display) <- "numeric"

# ------------------------------------------------------------------------------
# Impute Missing Values
# ------------------------------------------------------------------------------
# Function to impute missing values with row mean
impute_row_mean <- function(x) {
  if (all(is.na(x))) return(rep(0, length(x)))
  x[is.na(x)] <- mean(x, na.rm = TRUE)
  return(x)
}

# Apply imputation row-wise
mat_imputed <- t(apply(mat_display, 1, impute_row_mean))

# ------------------------------------------------------------------------------
# Calculate Row Ordering Metrics
# ------------------------------------------------------------------------------
# Find position of first non-NA value for each protein
first_non_na_pos <- apply(mat_display, 1, function(r) {
  idx <- which(!is.na(r))
  if (length(idx) == 0) Inf else idx[1]
})

# Calculate mean expression for each protein (using imputed values)
row_mean_imputed <- rowMeans(mat_imputed, na.rm = TRUE)

# ------------------------------------------------------------------------------
# Order Data by Module and Expression Pattern
# ------------------------------------------------------------------------------
# Order proteins by:
# 1. Module (factor order)
# 2. First appearance time
# 3. Mean expression (descending)
df_ordered <- df_3b %>%
  mutate(
    first_non_na_pos = first_non_na_pos[match(protein, rownames(mat_display))],
    row_mean_imputed = row_mean_imputed[match(protein, rownames(mat_display))]
  ) %>%
  arrange(Module, first_non_na_pos, desc(row_mean_imputed))

# Reorder display matrix to match protein order
mat_display_ord <- mat_display[df_ordered$protein, , drop = FALSE]

# ------------------------------------------------------------------------------
# Prepare Row Annotations
# ------------------------------------------------------------------------------
# Create row annotation data frame
annotation_row_plot <- data.frame(Module = df_ordered$Module)
rownames(annotation_row_plot) <- df_ordered$protein

# Calculate gap positions between modules
mod_rle <- rle(as.character(df_ordered$Module))
gap_indices <- cumsum(mod_rle$lengths)
gap_indices <- gap_indices[gap_indices < nrow(df_ordered)]

# ------------------------------------------------------------------------------
# Prepare Annotation Colors
# ------------------------------------------------------------------------------
# Create color mapping for modules
module_to_color <- module_color_map
annotation_colors <- list(Module = module_to_color)

# ------------------------------------------------------------------------------
# Prepare Heatmap Color Scale
# ------------------------------------------------------------------------------
# Calculate symmetric color scale based on max absolute value
max_abs <- max(abs(mat_imputed[1:20])+4, na.rm = TRUE)
breaks <- seq(-max_abs, max_abs, length.out = 101)
colors <- colorRampPalette(rev(brewer.pal(11, "RdBu")))(100)

# ------------------------------------------------------------------------------
# Prepare Display Numbers
# ------------------------------------------------------------------------------
# Create matrix of numbers to display in cells
number_matrix <- mat_display_ord
number_matrix[!is.na(number_matrix)] <-
  sprintf("%.2f", number_matrix[!is.na(number_matrix)])
number_matrix[is.na(number_matrix)] <- ""

# ------------------------------------------------------------------------------
# Format Column Names
# ------------------------------------------------------------------------------
# Create readable time bin labels from numeric column names
bin_vals <- as.numeric(colnames(mat_display_ord))
colnames(mat_display_ord) <- ifelse(
  bin_vals < 0,
  paste0("-Y", abs(bin_vals - 0.5)),
  paste0("Y", bin_vals + 0.5)
)

# ------------------------------------------------------------------------------
# Mark Key Proteins with Asterisks
# ------------------------------------------------------------------------------
# # List of key proteins to highlight (used in Figure 2B)
# key_proteins <- c("CA3", "EDA2R", "NEFL", "CALCA")
#
# # Add ** suffix to key protein names in all relevant objects
# for (obj_name in c("mat_display_ord", "number_matrix", "annotation_row_plot")) {
#   obj <- get(obj_name)
#   rnames <- rownames(obj)
#   mask <- rnames %in% key_proteins
#   rownames(obj)[mask] <- paste0(rnames[mask], "**")
#   assign(obj_name, obj)
# }

# ------------------------------------------------------------------------------
# Create Heatmap (Figure 3C)
# ------------------------------------------------------------------------------
ph_2a <- pheatmap(
  mat_display_ord,
  color = colors,
  breaks = breaks,
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_row = annotation_row_plot,
  annotation_colors = annotation_colors,
  gaps_row = gap_indices,
  gaps_col = sum(bin_vals < 0),  # Separate pre- and post-diagnosis periods
  display_numbers = number_matrix,
  fontsize_number = 12,
  fontsize_row = 12,
  annotation_names_row = FALSE,
  fontsize_col = 12,
  fontsize = 12,
  angle_col = 0,
  cellwidth = 35,
  border_color = "black",
  number_color = "black",
  na_col = "white"
)

# Stack the Module legend ABOVE the color legend, then convert to ggplot
gt_3b <- stack_legends_module_on_top(ph_2a)
p_3b <- ggplotify::as.ggplot(gt_3b) +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

# Save Figure 3C in multiple formats
ggsave(
  file.path(output_dir, "Sup_Fig1C_Temporal_Heatmap_pre_Convert.pdf"),
  plot = p_3b,
  device = cairo_pdf,
  width = 24,
  height = 12,
  units = "in"
)

ggsave(
  file.path(output_dir, "Sup_Fig1C_Temporal_Heatmap_pre_Convert.png"),
  plot = p_3b,
  width = 24,
  height = 12,
  units = "in",
  dpi = 300
)
ggsave(
  file.path(output_dir, "Sup_Fig1C_Temporal_Heatmap_pre_Convert.svg"),
  plot = p_3b,
  device = svglite::svglite,
  width = 24,
  height = 12,
  units = "in"
)

# ==============================================================================
# END OF SCRIPT
# ==============================================================================