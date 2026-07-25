# ==============================================================================
# Fig2_genaral.R — Figure 2: Temporal Protein Expression in ALS
# ==============================================================================
# Author : Ximing Ran
# Date   : 2026-02-22
#
# Description:
#   Generates all panels of Figure 2 for the ALS plasma proteomics study:
#   Fig2A - Temporal heatmap of mean protein log2FC across 1-year time bins
#           relative to phenoconversion, with proteins grouped by PPI module
#   Fig2B - Individual spaghetti + GAM-smoothed trajectory plots for key
#           proteins, showing pre-conversion dynamics for converters vs.
#           affected carriers
#
# Input data:
#   data/fig_data/summary_data/Fig2/
#     heat_1yr_wide.csv                                  — 1-yr binned mean
#                                                          log2FC per protein
#     Heatmap_PPI_cluster_mask_with_cluster_annotation_row.csv
#                                                        — PPI module labels
#     predictions_Convert_Affected_AllTime.csv           — GAM fitted values
#   data/fig_data/individual_data/Fig2/  (optional — Fig2B skipped if absent)
#     delta_matrix_wide.csv               — per-subject delta protein values
#     df_visit_info(N=516).csv            — visit metadata
#
# Output (./results/):
#   Fig_2A_Temporal_Heatmap.{pdf,png}
#   Fig_2B_Key_Protein_Trajectories.{pdf,png}  (only if individual data exists)
#
# Prerequisite:
#   Run analysis/02-Protein_Trajectories/ Rmd notebooks first.
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
library(svglite)

# Create output directory if it doesn't exist
output_dir <- here::here("figure", "Fig2", "results")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Set random seed for reproducibility (if applicable)
set.seed(2025)

cat("\n================================================================================\n")
cat("Starting Figure 2 Analysis\n")
cat("================================================================================\n\n")

# ==============================================================================
# FIGURE 2A: TEMPORAL HEATMAP
# ==============================================================================

cat("Generating Figure 2A: Temporal Heatmap...\n")

# ------------------------------------------------------------------------------
# Load and Prepare Data for Figure 2A
# ------------------------------------------------------------------------------

# Load temporal expression data
df_2a <- read.csv(
  here::here("data", "fig_data",  "summary_data", "Fig2", "heat_1yr_wide.csv")
)


cat(sprintf("  Loaded data: %d proteins x %d columns\n", nrow(df_2a), ncol(df_2a)))

# ------------------------------------------------------------------------------
# Fix Column Names for Time Bins
# ------------------------------------------------------------------------------

# Identify time bin columns (encoded by make.names, e.g. X.15.5, X0.5)
# Meta columns are non-numeric; time bin columns encode numeric bin centers
cols <- colnames(df_2a)
meta_cols <- c("protein", "start_year", "dir", "P_value", "Adj_P_value",
               "nearest_interval_start",
               "start_idx", "end_idx", "run_sign", "color_hex", "cluster")
time_bin_idx <- which(!cols %in% meta_cols)
# Exclude the row-name column if present
time_bin_idx <- time_bin_idx[time_bin_idx > 1 | cols[1] != "protein"]

# Decode make.names encoding back to signed numeric strings
for (i in time_bin_idx) {
  s <- gsub("^X", "", cols[i])
  if (grepl("^\\.", s)) s <- sub("^\\.", "-", s)
  cols[i] <- s
}
colnames(df_2a) <- cols

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
df_2a <- left_join(df_2a, annotation_row, by = "protein")

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
df_2a <- df_2a %>%
  filter(Adj_P_value < 0.05,
         !is.na(nearest_interval_start) & nearest_interval_start < 0)

cat(sprintf("  Filtered to %d proteins (adj. p < 0.05 & nearest_interval_start < 0)\n", nrow(df_2a)))

# Assign modules based on cluster membership
df_2a$Module <- ifelse(
  df_2a$cluster %in% c(1, 2, 3, 5, 6, 8),
  unname(module_map[as.character(df_2a$cluster)]),
  "Others"
)

# Convert Module to ordered factor
df_2a <- df_2a %>%
  mutate(Module = factor(Module, levels = module_levels, ordered = TRUE))

# ------------------------------------------------------------------------------
# Build Expression Matrix
# ------------------------------------------------------------------------------

# Create numeric matrix of expression values
mat_display <- df_2a %>% 
  select(all_of(time_bin_cols)) %>% 
  as.matrix()

rownames(mat_display) <- df_2a$protein
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
df_ordered <- df_2a %>%
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
# Create Heatmap (Figure 2A)
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

# Convert pheatmap to ggplot object
p_2a <- ggplotify::as.ggplot(ph_2a$gtable) +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

# Save Figure 2A in multiple formats
ggsave(
  file.path(output_dir, "Fig_2A_Temporal_Heatmap.pdf"),
  plot = p_2a,
  device = cairo_pdf,
  width = 24,
  height = 20,
  units = "in"
)
ggsave(
  file.path(output_dir, "Fig_2A_Temporal_Heatmap.png"),
  plot = p_2a,
  width = 24,
  height = 20,
  units = "in",
  dpi = 300
)
ggsave(
  file.path(output_dir, "Fig_2A_Temporal_Heatmap.svg"),
  plot = p_2a,
  device = svglite::svglite,
  width = 24,
  height = 20,
  units = "in"
)


cat("  Figure 2A saved successfully!\n")
cat(sprintf("    Total proteins displayed: %d\n", nrow(df_ordered)))
cat("    Proteins per module:\n")
module_table <- table(df_ordered$Module)
for (i in seq_along(module_table)) {
  cat(sprintf("      %s: %d\n", names(module_table)[i], module_table[i]))
}

# ==============================================================================
# FIGURE 2B: INDIVIDUAL PROTEIN TRAJECTORIES (UPDATED)
# ==============================================================================

cat("\nGenerating Figure 2B: Individual Protein Trajectories...\n")

# ------------------------------------------------------------------------------
# Check for Data Availability
# ------------------------------------------------------------------------------

# Define file paths
confirm_df_path <- here::here("data", "fig_data", "individual_data", "Fig2",
                              "delta_matrix_wide.csv")
visit_info_path <- here::here("data", "fig_data", "individual_data", "Fig2",
                              "df_visit_info(N=516).csv")
grid_overall_path <- here::here("data", "fig_data", "summary_data", "Fig2",
                                "predictions_Convert_Affected_AllTime.csv")

# Check if individual data files exist
if (!file.exists(confirm_df_path) || !file.exists(grid_overall_path) || !file.exists(visit_info_path)) {
  cat("  WARNING: Individual trajectory data files not found!\n")
  cat("  Expected files:\n")
  cat(sprintf("    - %s\n", confirm_df_path))
  cat(sprintf("    - %s\n", grid_overall_path))
  cat(sprintf("    - %s\n", visit_info_path))
  cat("  Skipping Figure 2B generation.\n\n")
  
} else {
  
  # ----------------------------------------------------------------------------
  # Load Data for Figure 2B
  # ----------------------------------------------------------------------------
  
  cat("  Individual data files found. Loading...\n")
  
  # Load individual protein trajectory data (WIDE format: proteins as columns)
  confirm_df_wide <- read.csv(confirm_df_path)
  
  # Load visit info with subject metadata and time variables
  visit_info <- read.csv(visit_info_path) %>% filter(Group %in% c("Convert", "Affected"))
  
  # Load GAM smoothed trajectory data
  grid_overall <- read.csv(grid_overall_path)
  
  # ----------------------------------------------------------------------------
  # Reshape confirm_df from Wide to Long Format
  # ----------------------------------------------------------------------------
  
  cat("  Reshaping wide-format data to long format...\n")
  
  # Identify protein columns (all columns except SampleID)
  protein_cols <- setdiff(colnames(confirm_df_wide), "SampleID")
  
  # Pivot from wide to long: one row per SampleID x protein
  confirm_df_long <- confirm_df_wide %>%
    tidyr::pivot_longer(
      cols = all_of(protein_cols),
      names_to = "protein",
      values_to = "delta"
    )
  
  # Join with visit_info to get UIDx, YrSinceOs, and other metadata
  # Match on SampleID
  confirm_df <- confirm_df_long %>%
    inner_join(
      visit_info %>% select(SampleID, UIDx, YrSinceOs, CollNum, Group, GenoGroup, Sex),
      by = "SampleID"
    )
  
  cat(sprintf("  Loaded trajectory data for %d proteins\n",
              length(unique(confirm_df$protein))))
  cat(sprintf("  Total observations after join: %d\n", nrow(confirm_df)))
  
  # ----------------------------------------------------------------------------
  # Define Plotting Parameters
  # ----------------------------------------------------------------------------
  
  # Color mapping for direction of change
  fill_map <- c(increase = "#D73027", decrease = "#4575B4")
  
  # Y-axis limits for all plots
  y_min_val <- -2.0
  y_max_val <- 3.0
  
  # Y-coordinates for significance markers and labels
  y_top_val <- 2.6      # For marker row
  y_label_val <- 2.9    # For text labels
  
  # ----------------------------------------------------------------------------
  # Generate Trajectory Plots for All Proteins
  # ----------------------------------------------------------------------------
  
  # Get list of all proteins
  protein_list <- c("CA3", "EDA2R", "NEFL", "CALCA")
  
  # Initialize empty list to store plots
  plot_list <- list()
  
  cat(sprintf("  Generating trajectory plots for %d proteins...\n", length(protein_list)))
  
  # Loop through each protein
  for (prot in protein_list) {
    
    # Filter data for current protein
    confirm_df_prot <- confirm_df %>% filter(protein == prot)
    grid_overall_prot <- grid_overall %>% filter(protein == prot)
    
    # Skip if no GAM data for this protein
    if (nrow(grid_overall_prot) == 0) {
      cat(sprintf("    Skipping %s: no GAM predictions found\n", prot))
      next
    }
    
    # Sort by subject ID and time
    confirm_df_sorted <- confirm_df_prot %>% arrange(UIDx, YrSinceOs)
    
    # Filter to significant points (pre-diagnosis only)
    sig_points_df <- grid_overall_prot %>%
      filter(!is.na(dir), YrSinceOs <= 0)
    
    # Adjust confidence intervals to stay within plot limits
    grid_overall_prot <- grid_overall_prot %>%
      mutate(
        lower = ifelse(lower < y_min_val, y_min_val, lower),
        upper = ifelse(upper > y_max_val, y_max_val, upper)
      )
    
    # ------------------------------------------------------------------------
    # Identify Contiguous Significant Regions
    # ------------------------------------------------------------------------
    
    r <- rle(grid_overall_prot$sig)
    ends <- cumsum(r$lengths)
    starts <- c(1, head(ends, -1) + 1)
    nonzero_runs_idx <- which(r$values != 0)
    
    if (length(nonzero_runs_idx) > 0) {
      run_starts_idx <- starts[nonzero_runs_idx]
      run_ends_idx <- ends[nonzero_runs_idx]
      run_vals <- r$values[nonzero_runs_idx]
      
      run_df <- tibble(
        y_label = y_label_val,
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
          y_label_val_ind = ifelse(dir == "increase", -0.4, 0.2)
        ) %>%
        # Remove end label if very close to diagnosis
        mutate(label_end = ifelse(abs(end_year) < 0.1, "", label_end)) %>%
        # Remove end label for post-diagnosis periods
        mutate(label_end = ifelse(end_year > 0, "", label_end)) %>%
        # Filter to pre-diagnosis periods only
        filter(start_year <= 0)
      
      # Adjust spacing for short intervals
      run_df <- run_df %>%
        mutate(
          end_year = ifelse((end_year - start_year) < 3 & end_year <= 0,
                            start_year + 2, end_year),
          label_end = ifelse((end_year - start_year) < 3 & end_year <= 0,
                             sprintf("%.1f", end_year), label_end)
        ) %>%
        arrange(start_year) %>%
        # Adjust start year if intervals are too close
        mutate(
          start_year = ifelse(
            lag((end_year - start_year) < 3 &
                  (lead(start_year, default = first(start_year)) - end_year) <= 2,
                default = FALSE),
            lag(start_year) + 5,
            start_year
          )
        )
      
    } else {
      # No significant regions found
      run_df <- tibble(
        start_idx = integer(),
        end_idx = integer(),
        start_year = numeric(),
        end_year = numeric(),
        run_sign = integer(),
        dir = character(),
        label_start = character(),
        label_end = character(),
        color_hex = character(),
        y_label = numeric(),
        y_label_val_ind = numeric()
      )
    }
    
    # ------------------------------------------------------------------------
    # Create Individual Protein Plot
    # ------------------------------------------------------------------------
    
    p <- ggplot() +
      # Individual subject trajectories (grey lines)
      geom_line(
        data = confirm_df_sorted,
        aes(x = YrSinceOs, y = delta, group = UIDx),
        color = "grey10",
        alpha = 0.8,
        size = 0.3
      ) +
      # Individual data points
      geom_point(
        data = confirm_df_prot,
        aes(x = YrSinceOs, y = delta),
        alpha = 0.90,
        size = 1,
        color = "grey10"
      ) +
      # GAM confidence interval
      geom_ribbon(
        data = grid_overall_prot,
        aes(x = YrSinceOs, ymin = lower, ymax = upper),
        fill = "grey30",
        alpha = 0.2
      ) +
      # GAM fitted line
      geom_line(
        data = grid_overall_prot,
        aes(x = YrSinceOs, y = fit),
        color = "black",
        size = 1.6
      ) +
      # Reference line at y=0
      geom_hline(
        yintercept = 0,
        linetype = "dashed",
        color = "grey40"
      ) +
      # Vertical line at diagnosis (x=0)
      geom_vline(
        xintercept = 0,
        color = "black"
      ) +
      # Significance markers
      geom_point(
        data = sig_points_df,
        aes(x = YrSinceOs, y = 0, color = dir),
        size = 1.2,
        alpha = 0.9,
        show.legend = FALSE
      ) +
      geom_label(
        data = run_df,
        aes(x = start_year, y = y_label_val_ind-0.2, label = label_start),
        fill = "white",
        color = run_df$color_hex,
        size = 5,
        fontface = "bold",
        vjust = 0,
        label.size = 0.5
      ) +
      # End year labels
      geom_text(
        data = run_df,
        aes(x = end_year, y = y_label_val_ind, label = label_end),
        color = run_df$color_hex,
        size = 5,
        fontface = "bold",
        vjust = 0
      ) +
      # Color scales
      scale_fill_manual(values = fill_map, name = "Direction") +
      scale_color_manual(values = fill_map, guide = FALSE) +
      # Labels and limits
      labs(
        title = prot,
        x = "Years to Phenoconversion",
        y = "log2FC"
      ) +
      ylim(y_min_val, y_max_val) +
      # Theme
      theme_bw() +
      theme(
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.title.x = element_text(size = 18),
        axis.title.y = element_text(size = 18),
        axis.text.x = element_text(size = 16),
        axis.text.y = element_text(size = 16),
        legend.position = "none",
        legend.title = element_text(size = 18),
        legend.text = element_text(size = 18),
        plot.background = element_rect(fill = "transparent", color = NA),
        panel.background = element_rect(fill = "transparent", color = NA)
      )
    
    # Store plot in list
    plot_list[[prot]] <- p
  }
  
  cat("  All trajectory plots generated!\n")
  
  # ----------------------------------------------------------------------------
  # Create Combined Figure 2B (Key Proteins)
  # ----------------------------------------------------------------------------
  
  # Select key proteins for main figure
  plot_list_select <- c(
    plot_list[["CA3"]],
    plot_list[["EDA2R"]],
    plot_list[["NEFL"]],
    plot_list[["CALCA"]]
  )
  
  # Remove any NULLs (proteins not found)
  plot_list_select <- plot_list_select[!sapply(plot_list_select, is.null)]
  
  if (length(plot_list_select) == 0) {
    cat("  WARNING: None of the key_proteins were found in the data.\n")
    cat("  Available proteins (first 20): ",
        paste(head(protein_list, 20), collapse = ", "), "\n")
  } else {
    # Combine plots vertically
    # p_2b_combined <- cowplot::plot_grid(
    #   plotlist = plot_list_select,
    #   ncol = 1,
    #   align = "v"
    # )
    
    
    library(cowplot)
    
    spacer <- ggplot() + theme_void()
    
    p_2b_combined <- plot_grid(
      plot_list_select[[1]],
      spacer,
      plot_list_select[[2]],
      spacer,
      plot_list_select[[3]],
      spacer,
      plot_list_select[[4]],
      ncol = 1,
      rel_heights = c(1, 0.06, 1, 0.06, 1, 0.06, 1),
      align = "v"
    )
    
    
    # Save combined Figure 2B
    ggsave(
      filename = file.path(output_dir, "Fig_2B_Key_Protein_Trajectories.pdf"),
      plot = p_2b_combined,
      device = cairo_pdf,
      width = 5,
      height = 20
    )
    ggsave(
      filename = file.path(output_dir, "Fig_2B_Key_Protein_Trajectories.png"),
      plot = p_2b_combined,
      width = 5,
      height = 20,
      dpi = 300
    )
    ggsave(
      filename = file.path(output_dir, "Fig_2B_Key_Protein_Trajectories.svg"),
      plot = p_2b_combined,
      device = svglite::svglite,
      width = 5,
      height = 20
    )

    cat("  Figure 2B (combined) saved successfully!\n")
  }

  # Save individual NEFL trajectory plot
  if (!is.null(plot_list[["NEFL"]])) {
    p_NEFL <- plot_list[["NEFL"]]
    ggsave(
      filename = file.path(output_dir, "Fig_2B_NEFL_Trajectory.pdf"),
      plot = p_NEFL,
      device = cairo_pdf,
      width = 5,
      height = 5
    )
    ggsave(
      filename = file.path(output_dir, "Fig_2B_NEFL_Trajectory.png"),
      plot = p_NEFL,
      width = 5,
      height = 5,
      dpi = 300
    )
    ggsave(
      filename = file.path(output_dir, "Fig_2B_NEFL_Trajectory.svg"),
      plot = p_NEFL,
      device = svglite::svglite,
      width = 5,
      height = 5
    )
  }
}

# ==============================================================================
# END OF SCRIPT
# ==============================================================================