# ==============================================================================
# Sup_Fig2.R — Supplementary Figure 2: Machine Learning AUC Heatmap
# ==============================================================================
# Author : Ximing Ran
# Date   : 2026-02-22
#
# Description:
#   Generates Supplementary Figure 2 — a heatmap of the best cross-validated
#   AUC achieved by each machine learning model at each prediction time
#   horizon (T = 0.5, 1, 2, 3, 5 years before phenoconversion).
#   Models compared include Logistic Regression, Random Forest, and others.
#   Rows are sorted by descending mean AUC; cell values are annotated.
#
# Input data:
#   data/fig_data/summary_data/Supplement_Figure/
#     ML_result.csv   — columns: Model, T (time horizon), AUC
#                       one row per (model, T, feature set) combination
#
# Output (./results/):
#   Sup_Fig2_AUC.{pdf,png}
#
# Prerequisite:
#   Run analysis/04-Phenoconversion_event_prediction/ Rmd notebooks and/or
#   Python ML scripts first.
# ==============================================================================

library(tidyverse)
library(pheatmap)
library(RColorBrewer)
library(ggplotify)
library(here)
library(svglite)

# ------------------------------------------------------------------------------#
# 1. Load data                                                                  #
# ------------------------------------------------------------------------------#
result_df_ML <- read.csv(
  here("data", "fig_data", "summary_data",
       "Supplement_Figure", "ML_result.csv")
)

# ------------------------------------------------------------------------------#
# 2. Best AUC per (Model, T)                                                    #
# ------------------------------------------------------------------------------#
best_auc <- result_df_ML %>%
  mutate(T = as.numeric(T)) %>%        # ensure T is numeric
  group_by(Model, T) %>%
  summarise(BestAUC = max(AUC, na.rm = TRUE), .groups = "drop")

# ------------------------------------------------------------------------------#
# 3. Wide matrix: rows = Model, cols = T                                        #
# ------------------------------------------------------------------------------#
heat_mat <- best_auc %>%
  mutate(T = as.character(T)) %>%
  pivot_wider(names_from = T, values_from = BestAUC) %>%
  arrange(Model) %>%
  column_to_rownames("Model") %>%
  as.matrix()

# Order columns as T = 0.5, 1, 2, 3, 5 (if present)
heat_mat <- heat_mat[, order(as.numeric(colnames(heat_mat))), drop = FALSE]

# Drop all-NA rows if any
all_na_rows <- apply(heat_mat, 1, function(r) all(is.na(r)))
if (any(all_na_rows)) {
  heat_mat <- heat_mat[!all_na_rows, , drop = FALSE]
}

# Sort rows by average AUC (descending; tie-break by max AUC)
row_order <- order(
  rowMeans(heat_mat, na.rm = TRUE),
  apply(heat_mat, 1, max, na.rm = TRUE),
  decreasing = TRUE
)
heat_mat <- heat_mat[row_order, , drop = FALSE]

# ------------------------------------------------------------------------------#
# 4. Pretty model labels                                                        #
# ------------------------------------------------------------------------------#
rename_map <- c(
  "LR" = "Logistic Regression",
  "RF" = "Random Forest"
)

rn_old <- rownames(heat_mat)
rn_new <- ifelse(rn_old %in% names(rename_map), rename_map[rn_old], rn_old)
rownames(heat_mat) <- rn_new

# ------------------------------------------------------------------------------#
# 5. Color scale and numbers on tiles                                           #
# ------------------------------------------------------------------------------#
# Blues palette: light = low AUC, medium = high AUC (lighter for readability)
cell_pal <- colorRampPalette(brewer.pal(9, "Blues")[1:5])(256)

# Use the observed min and max AUC so the full palette is used
auc_minmax <- range(heat_mat, na.rm = TRUE)
breaks <- seq(auc_minmax[1], auc_minmax[2],
              length.out = length(cell_pal) + 1)

numbers_mat <- ifelse(is.na(heat_mat), "", sprintf("%.3f", heat_mat))

# ------------------------------------------------------------------------------#
# 6. Draw heatmap with pheatmap                                                 #
# ------------------------------------------------------------------------------#
ph <- pheatmap(
  heat_mat,
  color           = cell_pal,
  breaks          = breaks,
  scale           = "none",
  cluster_rows    = FALSE,
  cluster_cols    = FALSE,
  show_rownames   = TRUE,
  show_colnames   = TRUE,
  display_numbers = numbers_mat,
  number_color    = "black",
  border_color    = NA,
  fontsize        = 15,
  fontsize_row    = 15,
  angle_col       = 0
)

# Convert to ggplot object (for combining with other panels if needed)
p <- ggplotify::as.ggplot(ph$gtable) +
  theme(
    plot.background  = element_rect(fill = "transparent", color = NA),
    panel.background = element_rect(fill = "transparent", color = NA)
  )

# ------------------------------------------------------------------------------#
# 7. Save figure                                                                #
# ------------------------------------------------------------------------------#
dir.create("./results", recursive = TRUE, showWarnings = FALSE)
ggsave(
  filename = "./results/Sup_Fig2_AUC.pdf",
  plot     = p,
  width    = 8,
  height   = 6,
  units    = "in"
)

ggsave(
  filename = "./results/Sup_Fig2_AUC.png",
  plot     = p,
  width    = 8,
  height   = 6,
  units    = "in",
  dpi      = 300
)
ggsave(
  filename = "./results/Sup_Fig2_AUC.svg",
  plot     = p,
  device   = svglite::svglite,
  width    = 8,
  height   = 6,
  units    = "in"
)
