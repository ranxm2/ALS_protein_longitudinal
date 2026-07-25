# ==============================================================================
# Fig4.R — Figure 4: Phenoconversion Prediction Models
# ==============================================================================
# Author : Ximing Ran
# Date   : 2025-12-01
#
# Description:
#   Generates all panels of Figure 4 for the ALS plasma proteomics study:
#   Fig4A - Heatmap of single-protein AUC across all ML models and time
#           horizons (T = 0.5–10 years before phenoconversion), annotated
#           by PPI module membership
#   Fig4B - ROC curves for logistic regression models at each time horizon,
#           with bootstrapped 95% CI ribbons and AUC annotations
#   Fig4C - UpSet plot showing overlap of proteins selected across time
#           points by logistic regression (77-protein set)
#   Fig4D - Predicted vs observed time-to-phenoconversion scatter plots for
#           all timing-estimation models (Tobit / LME / survival)
#
# Input data:
#   data/fig_data/summary_data/Fig4/
#     LR_single_protein_bootstrap_AUC.csv                — per-protein AUC
#       (from analysis/04-.../2b.Single_protein_AUC_bootstrap.Rmd)
#     Heatmap_PPI_cluster_mask_with_cluster_annotation_row.csv
#                                                        — PPI module labels
#     LR_77_proteins.csv                                 — logistic regression
#                                                          protein-level results
#     ROC_curve_LR.csv                                   — ROC curve points
#       (from analysis/04-.../4.Panel_Evaluation.Rmd)
#     ROC_CI_LR.csv                                      — bootstrapped CI
#       (from analysis/04-.../4.Panel_Evaluation.Rmd)
#     ROC_summary_LR.csv                                 — per-model AUC summary
#       (from analysis/04-.../4.Panel_Evaluation.Rmd)
#   data/fig_data/individual_data/Fig4/  (optional — Fig4D-G skipped if absent)
#     Predicted_Values_All_Models.csv                    — individual predicted
#       (from analysis/05-.../1.Phenoconversion_timing_estimation.Rmd)
#
# Output (./results/):
#   Fig_4A_heatmap_single_protein_LR.{pdf,png}
#   Fig_4B_ROC_LR.{pdf,png}
#   Fig_4C_upset_optimal_model_plot.{pdf,png}
#   Fig_4D-G_<panel>.{pdf,png}           (only if individual data exists)
#
# Prerequisite:
#   Run analysis/04-Phenoconversion_event_prediction/ and
#   analysis/05-Phenoconversion_timing_estimation/ Rmd notebooks first.
# ==============================================================================


# Set seed for reproducibility
set.seed(2025)

# Load required packages -------------------------------------------------------
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
library(purrr)
library(UpSetR)
library(patchwork)
library(ggplotify)
library(svglite)


# Create output directories ----------------------------------------------------
results_dir <- here::here("figure", "Fig4", "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Define file paths ------------------------------------------------------------
data_dir <-  here::here("data","fig_data",  "summary_data", "Fig4")

# Prefer new filenames; fall back to legacy date-suffixed versions
roc_curve_path   <- file.path(data_dir, "ROC_curve_LR.csv")
roc_ci_path      <- file.path(data_dir, "ROC_CI_LR.csv")
roc_summary_path <- file.path(data_dir, "ROC_summary_LR.csv")

if (!file.exists(roc_curve_path))   roc_curve_path   <- file.path(data_dir, "ROC_curve_LR_2025-10-20.csv")
if (!file.exists(roc_ci_path))      roc_ci_path      <- file.path(data_dir, "ROC_CI_LR_2025-10-20.csv")
if (!file.exists(roc_summary_path)) roc_summary_path <- file.path(data_dir, "ROC_summary_LR_2025-10-20.csv")

paths <- list(
  roc_curve = roc_curve_path,
  roc_ci = roc_ci_path,
  roc_summary = roc_summary_path,
  ml_results = file.path(data_dir, "LR_single_protein_bootstrap_AUC.csv"),
  annotation = file.path(data_dir, "Heatmap_PPI_cluster_mask_with_cluster_annotation_row.csv"),
  proteins = file.path(data_dir, "LR_77_proteins.csv")
)





# Figure 4A: HEATMAP ===========================================================

cat("Generating Figure 4B: Heatmap...\n")


# Load data --------------------------------------------------------------------
result_df <- read.csv(paths$ml_results)
annotation_row <- read.csv(paths$annotation, row.names = 1)
protein_all <- read.csv(paths$proteins)

# Filter for proteins of interest
result_df <- result_df %>%
  filter(Protein %in% protein_all$Protein)

# Define color palette ---------------------------------------------------------
cell_pal <- colorRampPalette(brewer.pal(11, "Blues"))(256)

# Build heat matrix ------------------------------------------------------------
result_df_lr <- result_df %>%
  filter(Model == "LR")

heat_mat_raw <- result_df_lr %>%
  select(Protein, T, Average_AUC) %>%
  mutate(T = as.numeric(T)) %>%
  arrange(T) %>%
  mutate(T = as.character(T)) %>%
  pivot_wider(names_from = T, values_from = Average_AUC) %>%
  column_to_rownames("Protein") %>%
  as.matrix()

# Keep T columns in numeric order
heat_mat_raw <- heat_mat_raw[, order(as.numeric(colnames(heat_mat_raw))), drop = FALSE]

# Module annotation ------------------------------------------------------------
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

# Assign modules
annotation_row$module <- ifelse(
  annotation_row$cluster %in% c(1, 2, 3, 5, 6, 8),
  unname(module_map[as.character(annotation_row$cluster)]),
  "Others"
)

# Align rows across matrix and annotation
common_prot <- intersect(rownames(heat_mat_raw), rownames(annotation_row))
heat_mat_raw <- heat_mat_raw[common_prot, , drop = FALSE]
ann_row <- annotation_row[common_prot, , drop = FALSE]

# Define module levels
module_levels <- c(
  "Skeletal Muscle",
  "TNF-Mediated Signaling",
  "Regeneration (incl. Neurofilament)",
  "ECM & ECM Interaction",
  "Others"
)

# Create annotation data frame
ann_row_for_plot <- data.frame(
  Module = factor(ann_row$module, levels = module_levels),
  row.names = rownames(ann_row)
)

# Define module colors (Okabe-Ito with C/D swapped)
module_color_map <- c(
  "Skeletal Muscle"                    = "#009E73",
  "TNF-Mediated Signaling"             = "#E69F00",
  "Regeneration (incl. Neurofilament)" = "#56B4E9",
  "ECM & ECM Interaction"              = "#CC79A7",
  "Others"                                     = "#bdbdbd"
)
ann_colors <- list(Module = module_color_map[module_levels])

# Helper functions for clustering ----------------------------------------------
scale_rows <- function(m) {
  # Ensure matrix is numeric
  if (!is.numeric(m)) {
    stop("Matrix must be numeric for scaling")
  }
  t(scale(t(m)))
}

row_dist <- function(m, method = "euclidean") {
  if (method == "correlation") {
    as.dist(1 - cor(t(m), use = "pairwise.complete.obs"))
  } else {
    dist(m, method = method)
  }
}

# Cluster rows within each module ----------------------------------------------
row_dist_method <- "correlation"
row_hclust_method <- "ward.D2"

modules_in_order <- levels(ann_row_for_plot$Module)
row_idx_list <- split(seq_len(nrow(heat_mat_raw)), ann_row_for_plot$Module)
new_row_order <- integer(0)

# Ensure heat_mat_raw is numeric
heat_mat_raw <- apply(heat_mat_raw, 2, as.numeric)
rownames(heat_mat_raw) <- common_prot

for (mod in modules_in_order) {
  idx <- row_idx_list[[mod]]
  if (length(idx) == 0) next
  
  # Extract submatrix and ensure it's numeric
  submat <- heat_mat_raw[idx, , drop = FALSE]
  
  # Check if submat is numeric
  if (!is.numeric(submat)) {
    cat("Warning: Non-numeric data in module", mod, "- skipping clustering\n")
    new_row_order <- c(new_row_order, idx)
    next
  }
  
  # Scale the submatrix
  submat_scaled <- suppressWarnings(scale_rows(submat))
  
  # Check for valid scaled matrix
  if (all(is.na(submat_scaled)) || !is.numeric(submat_scaled)) {
    cat("Warning: Scaling failed for module", mod, "- using original order\n")
    new_row_order <- c(new_row_order, idx)
    next
  }
  
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

# Apply within-module order and gaps
mat_ord <- heat_mat_raw[new_row_order, , drop = FALSE]
ann_row_for_plot <- ann_row_for_plot[new_row_order, , drop = FALSE]

sizes <- table(ann_row_for_plot$Module)[modules_in_order]
gaps_row <- cumsum(head(as.integer(sizes), -1))

# Define global color scale ----------------------------------------------------
auc_minmax <- c(0.50, 1.00)
breaks <- seq(auc_minmax[1], 0.95, length.out = length(cell_pal) + 1)

# Create display numbers matrix
numbers_mat <- ifelse(is.na(mat_ord), "", sprintf("%.2f", mat_ord))

# Draw heatmap -----------------------------------------------------------------
ph <- pheatmap(
  mat_ord,
  color = cell_pal,
  breaks = breaks,
  scale = "none",
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  show_rownames = TRUE,
  show_colnames = TRUE,
  annotation_row = ann_row_for_plot,
  annotation_colors = ann_colors,
  gaps_row = gaps_row,
  display_numbers = numbers_mat,
  number_color = "black",
  fontsize_number = 8.5,
  border_color = NA,
  angle_col = 0,
  annotation_names_row = FALSE,
  cellwidth = 30,
  legend = TRUE
)

# Convert to ggplot and save ---------------------------------------------------
p_b <- ggplotify::as.ggplot(ph$gtable) +
  theme(
    plot.background = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

ggsave(
  file.path(results_dir, "Fig_4A_heatmap_single_protein_LR.pdf"),
  p_b,
  device = cairo_pdf,
  width = 10, height = 10, units = "in"
)

ggsave(
  file.path(results_dir, "Fig_4A_heatmap_single_protein_LR.png"),
  p_b,
  width = 10, height = 10, units = "in", dpi = 300
)
ggsave(
  file.path(results_dir, "Fig_4A_heatmap_single_protein_LR.svg"),
  p_b,
  device = svglite::svglite,
  width = 10, height = 10, units = "in"
)

cat("Heatmap saved successfully!\n\n")


# Figure 4B: ROC PLOTS =========================================================



# Load data --------------------------------------------------------------------
roc_df <- read.csv(paths$roc_curve) %>%
  filter(Protein_Set %in% c("Panel_1", "Panel_15", "Panel_19")) %>%
  mutate(Protein_Set = ifelse(Protein_Set == "Panel_1", "NEFL", Protein_Set))

roc_ci_df <- read.csv(paths$roc_ci) %>%
  filter(Protein_Set == "Panel_1") %>%
  mutate(Protein_Set = "NEFL")

auc_df <- read.csv(paths$roc_summary) %>%
  filter(Protein_Set %in% c("Panel_1", "Panel_15", "Panel_19")) %>%
  mutate(Protein_Set = ifelse(Protein_Set == "Panel_1", "NEFL", Protein_Set)) %>%
  mutate(Protein_Set = factor(Protein_Set, levels = c("Panel_19", "Panel_15", "NEFL")))

# Helper function: Create ROC plot ---------------------------------------------
create_roc_plot <- function(t_val, roc_data, auc_data, ci_data) {
  # Filter data for this time point
  df_t <- filter(roc_data, T == t_val)
  auc_t <- filter(auc_data, T == t_val)
  ribbon_df <- filter(ci_data, Protein_Set == "NEFL", T == t_val)
  
  # Create labels with AUC values and CI for NEFL
  auc_t <- auc_t %>%
    mutate(label = ifelse(
      Protein_Set == "NEFL",
      sprintf("NEFL only (AUC=%.3f)", CV_AUC),
      sprintf("%s (AUC=%.3f)", Protein_Set, CV_AUC)
    )) %>%
    mutate(Protein_Set = factor(
      Protein_Set,
      levels = rev(c("Panel_19", "Panel_15", "NEFL"))
    ))
  
  # Create named vector for legend
  leg_labs <- auc_t$label
  names(leg_labs) <- auc_t$Protein_Set
  
  # Create plot
  ggplot() +
    geom_line(
      data = df_t,
      aes(x = fpr, y = sensitivity, color = Protein_Set),
      linewidth = 1
    ) +
    geom_abline(
      intercept = 0, slope = 1,
      linetype = "dashed", color = "gray50", linewidth = 0.5
    ) +
    scale_x_continuous(
      labels = percent_format(accuracy = 1),
      breaks = seq(0, 1, 0.2),
      limits = c(0, 1),
      expand = c(0, 0)
    ) +
    scale_y_continuous(
      labels = percent_format(accuracy = 1),
      breaks = seq(0, 1, 0.2),
      limits = c(0, 1),
      expand = c(0, 0)
    ) +
    scale_color_discrete(
      labels = leg_labs,
      breaks = rev(names(leg_labs))
    ) +
    labs(
      title = paste0("T = ", t_val),
      x = "False Positive Rate (1-Specificity)",
      y = "True Positive Rate (Sensitivity)",
      color = NULL
    ) +
    theme_bw() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20, face = "bold"),
      legend.position = c(0.6, 0.15),
      legend.background = element_rect(fill = "white", color = "grey80"),
      legend.text = element_text(size = 12),
      panel.grid.minor = element_blank(),
      axis.text = element_text(size = 14),
      axis.title = element_text(size = 14)
    ) +
    coord_equal() +
    guides(color = guide_legend(ncol = 1))
}

# Generate plots for all time points -------------------------------------------
t_vals <- sort(unique(roc_df$T))

plots <- lapply(t_vals, function(t_val) {
  create_roc_plot(t_val, roc_df, auc_df, roc_ci_df)
})

# Combine plots ----------------------------------------------------------------
combined_plot <- plot_grid(
  plotlist = plots,
  nrow = 1,
  label_size = 20,
  label_fontface = "bold",
  label_x = -0.02,
  label_y = 1
)

# Save plots -------------------------------------------------------------------
ggsave(
  file.path(results_dir, "Fig_4B_ROC_LR.pdf"),
  combined_plot,
  device = cairo_pdf,
  width = 25, height = 5, units = "in"
)

ggsave(
  file.path(results_dir, "Fig_4B_ROC_LR.png"),
  combined_plot,
  width = 25, height = 5, units = "in", dpi = 300
)
ggsave(
  file.path(results_dir, "Fig_4B_ROC_LR.svg"),
  combined_plot,
  device = svglite::svglite,
  width = 25, height = 5, units = "in"
)

cat("ROC plots (all panels) saved successfully!\n")



# Figure 4C: UPSET PLOT ========================================================

cat("Generating Figure 4C: UpSet plot...\n")

# Load data --------------------------------------------------------------------
result_lr <- read.csv(paths$proteins) %>%
  select(Protein, T_0.5, T_1, T_2, T_3, T_5)

# Define time point columns
t_cols <- c("T_0.5", "T_1", "T_2", "T_3", "T_5")

# Convert to logical -----------------------------------------------------------
upset_df_log <- result_lr %>%
  mutate(across(all_of(t_cols), as.logical))

# Count combinations -----------------------------------------------------------
counts_combinations <- upset_df_log %>%
  mutate(
    combination = pmap_chr(
      list(T_0.5, T_1, T_2, T_3, T_5),
      function(lgl1, lgl2, lgl3, lgl4, lgl5) {
        c("T_0.5", "T_1", "T_2", "T_3", "T_5")[c(lgl1, lgl2, lgl3, lgl4, lgl5)] %>%
          paste(collapse = ", ")
      }
    )
  ) %>%
  count(T_0.5, T_1, T_2, T_3, T_5, combination, name = "n")

# Calculate number of sets and assign colors
counts_combinations <- counts_combinations %>%
  mutate(
    n_sets = str_count(combination, "T_"),
    fill_color = case_when(
      combination == "T_1, T_2, T_3" ~ "#226ca7",
      combination == "T_1, T_2, T_3, T_5" ~ "#226ca7",
      combination == "T_0.5, T_1, T_2, T_3, T_5" ~ "#226ca7",
      TRUE ~ "black"
    )
  ) %>%
  arrange(n_sets, desc(n)) %>%
  mutate(combination = fct_inorder(combination))

# Count occurrences by time point ----------------------------------------------
counts_of_subjects <- upset_df_log %>%
  pivot_longer(
    cols = -Protein,
    names_to = "subject",
    values_to = "below"
  ) %>%
  summarize(
    counts = sum(below),
    .by = subject
  )

# Create bar chart for combinations --------------------------------------------
bar_chart <- counts_combinations %>%
  ggplot(aes(x = combination, y = n)) +
  geom_col(aes(fill = fill_color), width = 0.6) +
  scale_fill_identity() +
  geom_text(aes(label = n, color = fill_color), vjust = -0.5, size = 4) +
  scale_color_identity() +
  scale_y_continuous(
    limits = c(0, max(counts_combinations$n) * 1.15),
    expand = c(0, 0)
  ) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),  # ADD THIS
    # 
    # panel.grid.major.x = element_blank(),
    # panel.grid.minor.y = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.y = element_text(size = 14, color = "black")
  ) +
  labs(x = element_blank(), y = "Intersection Size")

# Prepare points data for connection plot --------------------------------------
points_data <- counts_combinations %>%
  mutate(
    subjects = map(
      combination,
      ~ str_split_1(as.character(.x), ",\\s*")
    )
  ) %>%
  unnest(subjects) %>%
  mutate(
    subjects = factor(
      subjects,
      levels = c("T_5", "T_3", "T_2", "T_1", "T_0.5")
    ),
    combination = factor(
      combination,
      levels = levels(counts_combinations$combination)
    )
  ) %>%
  filter(!is.na(subjects))

# Create point chart showing intersections -------------------------------------
point_chart <- points_data %>%
  ggplot(aes(x = combination, y = subjects)) +
  geom_line(aes(group = combination, color = fill_color), linewidth = 2) +
  geom_point(aes(color = fill_color), size = 5) +
  scale_color_identity() +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),  # ADD THIS
    
    # panel.grid.major.x = element_blank(),
    # panel.grid.minor.y = element_blank(),
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.y = element_text(hjust = 0.5, size = 14, color = "black"),
    axis.title.y = element_text(size = 14, color = "black")
  ) +
  labs(y = "Time Points") +
  scale_x_discrete(drop = FALSE)

# Create horizontal bar chart for subject counts -------------------------------
subject_bars <- counts_of_subjects %>%
  mutate(subject = factor(
    subject,
    levels = c("T_5", "T_3", "T_2", "T_1", "T_0.5")
  )) %>%
  ggplot(aes(x = -counts, y = subject)) +
  geom_col(width = 0.6, fill = "black") +
  geom_text(aes(label = counts, x = -counts - 4), hjust = -0.2, size = 4.5, color = "black") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),  # ADD THIS
    axis.text.y = element_text(hjust = 0.5, size = 14, color = "black"),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 14, color = "black")
  ) +
  scale_x_continuous(
    labels = abs,
    limits = c(-max(counts_of_subjects$counts) * 1.15, 0),
    expand = c(0, 0)
  ) +
  labs(x = element_blank(), y = element_blank())

# Combine plots ----------------------------------------------------------------
layout <- "##AAAA
BBCCCC"

p_c <- bar_chart +
  subject_bars +
  point_chart +
  plot_layout(
    design = layout,
    heights = c(0.75, 0.25)
  ) +
  theme(
    # plot.background  = element_rect(fill = "transparent", color = NA),
    # panel.background = element_rect(fill = "transparent", color = NA)
  )

# Save UpSet plot --------------------------------------------------------------
ggsave(
  file.path(results_dir, "Fig_4C_upset_optimal_model_plot.pdf"),
  plot = p_c,
  device = cairo_pdf,
  width = 12,
  height = 6,
  bg = "transparent"
)

ggsave(
  file.path(results_dir, "Fig_4C_upset_optimal_model_plot.png"),
  plot = p_c,
  width = 12,
  height = 6,
  dpi = 300,
  bg = "transparent"
)
ggsave(
  file.path(results_dir, "Fig_4C_upset_optimal_model_plot.svg"),
  plot = p_c,
  device = svglite::svglite,
  width = 12,
  height = 6,
  bg = "transparent"
)



cat("UpSet plot saved successfully!\n\n")


# SETUP ========================================================================

# Set seed for reproducibility
set.seed(2025)

# Load required packages -------------------------------------------------------
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
library(tibble)
library(patchwork)

# Set global theme -------------------------------------------------------------
theme_set(theme_bw() + theme(legend.position = "bottom"))

# Create output directories ----------------------------------------------------
results_dir <- here::here("figure", "Fig4", "results")
dir.create(results_dir, recursive = TRUE, showWarnings = FALSE)

# Define file paths ------------------------------------------------------------
data_file <- here::here(
  "data", "fig_data", "individual_data", "Fig4", 
  "Predicted_Values_All_Models.csv"
)


# HELPER FUNCTION: PREDICTION VS OBSERVED PLOT =================================

predict_vs_observed_plot_tobit <- function(y_true = NULL, 
                                           y_predict = NULL,
                                           vjust = 0.75, 
                                           title = NULL,
                                           MAE = NULL,
                                           RMSE = NULL, 
                                           cor_y = NULL) {
  
  # Create data frame ----------------------------------------------------------
  df <- data.frame(y_true = y_true, y_predict = y_predict)
  
  # Categorize points based on prediction accuracy
  df$category <- ifelse(
    df$y_predict < df$y_true, 
    "Onset Earlier than Predicted", 
    "Onset Later than Predicted"
  )
  
  # Define color mapping
  color_map <- c(
    "Onset Later than Predicted" = "blue", 
    "Onset Earlier than Predicted" = "red"
  )
  
  # Calculate residuals
  df$residuals <- df$y_true - df$y_predict
  
  # Calculate proportion of predictions within ±1 year
  prop_1 <- sum(abs(df$residuals) < 1, na.rm = TRUE) / nrow(df)
  
  # Determine range for axes (same for both x and y) --------------------------
  range_limits <- range(c(y_true, y_predict), na.rm = TRUE)
  
  # Generate coordinates for ±1 year deviation band ---------------------------
  band_x <- seq(range_limits[1], range_limits[2], length.out = 100)
  
  band_df <- data.frame(
    x = c(band_x, rev(band_x)),
    y = c(band_x - 1, rev(band_x + 1))  # Lower (x-1) and upper (x+1) boundaries
  )
  
  # Constrain band to axis limits
  band_df$y[band_df$y < range_limits[1]] <- range_limits[1]
  band_df$y[band_df$y > range_limits[2]] <- range_limits[2]
  
  # Create plot ----------------------------------------------------------------
  p <- ggplot(df, aes(x = y_predict, y = y_true, color = category)) +
    # Add ±1 year deviation band
    geom_polygon(
      data = band_df, 
      aes(x = x, y = y), 
      fill = "grey80", 
      alpha = 0.5, 
      inherit.aes = FALSE
    ) +
    # Add data points
    geom_point(alpha = 0.7, size = 2) +
    # Color scale
    scale_color_manual(values = color_map) +
    # Add y = x reference line
    geom_line(
      data = data.frame(x = band_x, y = band_x), 
      aes(x = x, y = y),
      color = "black", 
      linetype = "dashed", 
      size = 1, 
      inherit.aes = FALSE
    ) +
    # Labels
    labs(
      title = title,
      x = "Predicted Years to Phenoconversion",
      y = "Observed Years to Phenoconversion",
      color = "Prediction Type"
    ) +
    # Set axis limits (equal for both axes)
    xlim(range_limits) +
    ylim(range_limits) +
    coord_fixed() +
    # Theme
    theme_classic() +
    guides(color = guide_legend(nrow = 2)) +
    theme(
      plot.title = element_text(hjust = 0.5, size = 20),
      axis.text.x = element_text(size = 14),
      axis.title.x = element_text(size = 14),
      axis.title.y = element_text(size = 14),
      axis.text.y = element_text(size = 14),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.position = "none",
      plot.background = element_rect(fill = "transparent", color = NA),
      panel.background = element_rect(fill = "transparent", color = NA),
      legend.background = element_rect(fill = "transparent", color = NA)
    ) +
    # Add performance metrics annotation
    annotate(
      "text", 
      x = range_limits[1], 
      y = range_limits[2] - 1, 
      hjust = 0, 
      vjust = vjust,
      label =sprintf(
        "MAE = %.3f\nRMSE = %.3f\nCor = %.3f", 
        MAE, RMSE, cor_y, prop_1 * 100
      ),
      size = 5,
      color = "black"
    )
  
  return(p)
}

# Load prediction data ---------------------------------------------------------
if (!file.exists(data_file)) {
  cat("  WARNING: Individual data not found for Fig4D-G. Skipping.\n")
} else {
df_pred <- read.csv(data_file)

# Define protein sets to analyze
protein_set_list <- c("Panel_19", "Panel_15", "NEFL","Panel_data","Panel_data_ukb")

# Display names for titles
display_names <- c(
  "Panel_19"       = "19-Protein Panel\n(with curation)",
  "Panel_15"       = "15-Protein Panel\n(UKB subset of 19 proteins)",
  "NEFL"           = "NEFL Only\n",
  "Panel_data"     = "16-Protein Panel\n(purely data-driven)",
  "Panel_data_ukb" = "10-Protein Panel\n(UKB subset of 16 proteins)"
)

# Generate plots for each protein set ------------------------------------------
plot_list <- list()

for (protein_set in protein_set_list) {
  
  cat("  Processing:", protein_set, "\n")
  
  # Prepare data
  df_cor <- df_pred %>%
    filter(Protein == protein_set) %>%
    select(y_true, y_predict) %>%
    rename(
      y = y_true,
      y_pred = y_predict
    ) %>%
    mutate(
      y_pred = pmin(y_pred, 0),  # Apply Tobit model censoring at 0
      residual = y - y_pred
    ) %>%
    # Remove rows with missing values
    filter(!is.na(y) & !is.na(y_pred))
  
  # Calculate performance metrics
  MAE <- mean(abs(df_cor$residual), na.rm = TRUE)
  RMSE <- sqrt(mean((df_cor$residual)^2, na.rm = TRUE))
  cor_y <- cor(df_cor$y, df_cor$y_pred, use = "complete.obs")
  
  # Create plot
  p <- predict_vs_observed_plot_tobit(
    y_true = df_cor$y,
    y_predict = df_cor$y_pred,
    title = display_names[[protein_set]],
    MAE = round(MAE, 3),
    RMSE = round(RMSE, 3),
    cor_y = round(cor_y, 3)
  )
  
  # Store plot in list
  plot_list[[protein_set]] <- p
  
  ggsave(
    filename = file.path(results_dir, paste0("Fig_4D-G_", protein_set, ".pdf")),
    plot = p,
    device = cairo_pdf,
    width = 6,
    height = 6.5
  )
  
  ggsave(
    filename = file.path(results_dir, paste0("Fig_4D-G_", protein_set, ".png")),
    plot = p,
    width = 6,
    height = 6.5,
    dpi = 300
  )
  ggsave(
    filename = file.path(results_dir, paste0("Fig_4D-G_", protein_set, ".svg")),
    plot = p,
    device = svglite::svglite,
    width = 6,
    height = 6.5
  )
}

cat("All plots generated successfully!\n\n")


# COMBINE AND SAVE PLOTS =======================================================

cat("Saving combined plot...\n")

# Combine D-G plots in specified order: Panel_19, Panel_15, Panel_data, Panel_data_ukb, NEFL
panel_order <- c("Panel_19", "Panel_15", "Panel_data", "Panel_data_ukb", "NEFL")
combined_DG <- plot_grid(
  plotlist = plot_list[panel_order],
  ncol = 5,
  nrow = 1
)

ggsave(
  filename = file.path(results_dir, "Fig_4D-G_combined.pdf"),
  plot = combined_DG,
  device = cairo_pdf,
  width = 25, height = 5, units = "in"
)

ggsave(
  filename = file.path(results_dir, "Fig_4D-G_combined.png"),
  plot = combined_DG,
  width = 25, height = 5, units = "in", dpi = 300
)
ggsave(
  filename = file.path(results_dir, "Fig_4D-G_combined.svg"),
  plot = combined_DG,
  device = svglite::svglite,
  width = 25, height = 5, units = "in"
)

cat("Combined D-G plot saved successfully!\n\n")


# SAVE INDIVIDUAL PANEL_19 PLOT ================================================
# 
# cat("Saving individual Panel_19 plot...\n")
# 
# # Extract Panel_19 plot
# p1 <- plot_list[["Panel_19"]]
# 
# # Save Panel_19 plot separately
# ggsave(
#   filename = file.path(results_dir, "Fig_4E_Panel_19_all.svg"),
#   plot = p1,
#   width = 6,
#   height = 6.5
# )
# 
# ggsave(
#   filename = file.path(results_dir, "Fig_4E_Panel_19_all.pdf"),
#   plot = p1,
#   width = 6,
#   height = 6.5
# )
# ggsave(
#   filename = file.path(results_dir, "Fig_4E_Panel_19_all.png"),
#   plot = p1,
#   width = 6,
#   height = 6.5,
#   dpi = 300
# )
# 
#
#
# cat("Panel_19 plot saved successfully!\n\n")
} # end Fig4D-G individual data guard

