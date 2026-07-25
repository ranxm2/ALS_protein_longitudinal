# ==============================================================================
# Fig3.R — Figure 3: Cox Survival Analysis and Kaplan–Meier Curves
# ==============================================================================
# Author : Ximing Ran
# Date   : 2026-02-22
#
# Description:
#   Generates all panels of Figure 3 for the ALS plasma proteomics study:
#   Fig3A - Volcano plot of Cox proportional-hazards results; x-axis = hazard
#           ratio, y-axis = -log10(p-value); proteins significant at FDR < 0.05
#           coloured red (increased risk) or blue (decreased risk)
#   Fig3B - Kaplan–Meier curve stratified by baseline NEFL level (median
#           split: Low NEFL vs High NEFL)
#   Fig3C - Kaplan–Meier curve stratified by LASSO Cox risk score (median
#           split: Low Risk vs High Risk)
#
# Input data:
#   data/fig_data/summary_data/Fig3/
#     fig3a_volcano.csv                  — per-protein Cox HR, p-value, FDR
#   data/fig_data/individual_data/Fig3/  (optional — Fig3B/C skipped if absent)
#     fig4b_km_nefl/raw_data.csv         — individual time_yr, event, nefl_group
#     fig4C_km_lasso/raw_data.csv        — individual time_yr, event, risk_group
#
# Output (./results/):
#   Fig_3A_volcano.{pdf,png}
#   Fig_3B_km_nefl.{pdf,png}             (only if individual data exists)
#   Fig_3C_km_lasso.{pdf,png}            (only if individual data exists)
#   Fig_3C_km_lasso_square.pdf           (only if individual data exists)
#
# Prerequisite:
#   Run analysis/03-Time_to_Event_Analysis/1.Cox_model_baseline.Rmd first.
# ==============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)
library(ggrepel)
library(scales)
library(patchwork)
library(survival)
library(survminer)
library(ggh4x)
library(svglite)

# Output directory for saved plots

out_dir  <- here::here("figure", "Fig3", "results")
data_dir <- here::here("data" ,"fig_data",  "summary_data", "Fig3")
data_ind_dir    <- here::here("data" ,"fig_data",  "individual_data", "Fig3")



dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

theme_set(theme_classic(base_size = 13, base_family = "Arial"))


#!/usr/bin/env Rscript

# Standalone volcano plot extracted from an older Fig3.R version.

library(ggplot2)
library(dplyr)
library(ggrepel)
library(here)

out_dir <- here::here("figure", "Fig3", "results")
data_dir <- here::here("data", "fig_data", "summary_data", "Fig3")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

theme_set(theme_classic(base_size = 13, base_family = "Arial"))

volcano_df <- read.csv(
  file.path(data_dir, "fig3a_volcano.csv"),
  stringsAsFactors = FALSE
)

p_threshold <- unique(volcano_df$p_threshold)
y_max_plot <- unique(volcano_df$y_max_plot)
top_10_genes <- volcano_df %>% filter(is_label)

p_volcano <- ggplot(volcano_df, aes(x = exp(coef), y = negLogP, color = significance)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_hline(
    yintercept = -log10(p_threshold),
    linetype = "dashed",
    color = "gray40",
    linewidth = 0.8
  ) +
  scale_color_manual(
    values = c(UP = "#A70C20", DOWN = "#427ebf", NO = "grey60"),
    labels = c(
      UP = "Increased Risk (FDR < 0.05)",
      DOWN = "Decreased Risk (FDR < 0.05)",
      NO = "Not Significant"
    )
  ) +
  scale_y_continuous(limits = c(0, y_max_plot), expand = c(0, 0)) +
  labs(
    x = "Hazard Ratio",
    y = "-log10(p-value)",
    color = "Association"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "none",
    legend.text = element_text(size = 11),
    plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
    plot.subtitle = element_text(hjust = 0.5, size = 11, color = "gray40")
  ) +
  geom_label_repel(
    data = top_10_genes,
    aes(label = Protein, color = significance),
    size = 3.5,
    max.overlaps = 50,
    min.segment.length = 0,
    show.legend = FALSE,
    box.padding = 0.5,
    point.padding = 0.3,
    ylim = c(NA, y_max_plot)
  )

ggsave(file.path(out_dir, "Fig_3A_volcano.pdf"), p_volcano, device = cairo_pdf, width = 9, height = 6)
ggsave(file.path(out_dir, "Fig_3A_volcano.png"), p_volcano, width = 9, height = 6, dpi = 300)
ggsave(file.path(out_dir, "Fig_3A_volcano.svg"), p_volcano, device = svglite::svglite, width = 9, height = 6)

cat("Volcano plot saved to:", out_dir, "\n")















#----Figure B - Kaplan-Meier by NEFL group  (ggsurvplot from raw data)----
km_nefl_path <- file.path(data_ind_dir, "fig4b_km_nefl", "raw_data.csv")
if (!file.exists(km_nefl_path)) {
  cat("  WARNING: Individual data not found for Fig3B. Skipping.\n")
} else {
km_nefl_raw <- read.csv(km_nefl_path,
                        stringsAsFactors = FALSE) %>%
  mutate(
    nefl_group = factor(nefl_group, levels = c("Low NEFL", "High NEFL"),
                        labels = c("Low", "High"))
  )

km_nefl <- survfit(Surv(time_yr, event) ~ nefl_group, data = km_nefl_raw)

# Version 1: Low/High labels (for main figure)
p_km_nefl <- ggsurvplot(
  km_nefl,
  data              = km_nefl_raw,
  pval              = FALSE,
  conf.int          = TRUE,
  risk.table        = TRUE,
  risk.table.height = 0.25,
  risk.table.y.text.col = TRUE,
  palette           = c("#427ebf", "#A70C20"),
  legend.labs       = c("Low", "High"),
  legend.title      = "",
  xlab              = "Time (years)",
  ylab              = "",
  title             = "",
  xlim              = c(0, 15),
  break.time.by     = 5,
  surv.median.line  = "hv",
  legend            = "none",
  ggtheme           = theme_classic(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.title.y = element_text(margin = margin(r = -2)))
)

png(file.path(out_dir, "Fig_3B_km_nefl.png"),
    width = 9, height = 6, units = "in", res = 300)
print(p_km_nefl)
dev.off()

cairo_pdf(file.path(out_dir, "Fig_3B_km_nefl.pdf"),
    width = 9, height = 6)
print(p_km_nefl, newpage = FALSE)
dev.off()

svglite::svglite(file.path(out_dir, "Fig_3B_km_nefl.svg"),
    width = 9, height = 6)
print(p_km_nefl, newpage = FALSE)
dev.off()

# Version 2: Full legend labels (for legend extraction)
p_km_nefl_legend <- ggsurvplot(
  km_nefl,
  data              = km_nefl_raw,
  pval              = "Log-rank\np = 0.0002",
  pval.coord        = c(0.2, 0.08),
  conf.int          = TRUE,
  risk.table        = FALSE,
  palette           = c("#427ebf", "#A70C20"),
  legend.labs       = c("NEFL only, lower risk score", "NEFL only, higher risk score"),
  legend.title      = "",
  xlab              = "Time (years)",
  ylab              = "",
  title             = "",
  xlim              = c(0, 15),
  break.time.by     = 5,
  surv.median.line  = "hv",
  legend            = c(0.85, 0.95),
  ggtheme           = theme_classic(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.title.y = element_text(margin = margin(r = -2)))
)

png(file.path(out_dir, "Fig_3B_km_nefl_legend.png"),
    width = 9, height = 6, units = "in", res = 300)
print(p_km_nefl_legend)
dev.off()

cairo_pdf(file.path(out_dir, "Fig_3B_km_nefl_legend.pdf"),
    width = 9, height = 6)
print(p_km_nefl_legend, newpage = FALSE)
dev.off()

svglite::svglite(file.path(out_dir, "Fig_3B_km_nefl_legend.svg"),
    width = 9, height = 6)
print(p_km_nefl_legend, newpage = FALSE)
dev.off()

} # end Fig3B individual data guard

#----Figure C - Kaplan-Meier by LASSO risk group  (ggsurvplot from raw data)----
km_lasso_path <- file.path(data_ind_dir, "fig4C_km_lasso", "raw_data.csv")
if (!file.exists(km_lasso_path)) {
  cat("  WARNING: Individual data not found for Fig3C. Skipping.\n")
} else {
km_lasso_raw <- read.csv(km_lasso_path,
                         stringsAsFactors = FALSE) %>%
  mutate(
    risk_group = factor(risk_group, levels = c("Low Risk", "High Risk"),
                        labels = c("Low", "High"))
  )

km_lasso <- survfit(Surv(time_yr, event) ~ risk_group, data = km_lasso_raw)

# Version 1: Low/High labels (for main figure)
p_km_lasso <- ggsurvplot(
  km_lasso,
  data              = km_lasso_raw,
  pval              = FALSE,
  conf.int          = TRUE,
  risk.table        = TRUE,
  risk.table.height = 0.25,
  risk.table.y.text.col = TRUE,
  palette           = c("#427ebf", "#A70C20"),
  legend.labs       = c("Low", "High"),
  legend.title      = "",
  xlab              = "Time (years)",
  ylab              = "",
  title             = "",
  xlim              = c(0, 15),
  break.time.by     = 5,
  surv.median.line  = "hv",
  legend            = "none",
  ggtheme           = theme_classic(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.title.y = element_text(margin = margin(r = -2)))
)

# Version 2: Full legend labels (for legend extraction)
p_km_lasso_legend <- ggsurvplot(
  km_lasso,
  data              = km_lasso_raw,
  pval              = TRUE,
  pval.coord        = c(0.2, 0.08),
  pval.method       = TRUE,
  pval.method.coord = c(0.5, 0.22),
  conf.int          = TRUE,
  risk.table        = FALSE,
  palette           = c("#427ebf", "#A70C20"),
  legend.labs       = c("LASSO (7 proteins), lower risk score", "LASSO (7 proteins), higher risk score"),
  legend.title      = "",
  xlab              = "Time (years)",
  ylab              = "",
  title             = "",
  xlim              = c(0, 15),
  break.time.by     = 5,
  surv.median.line  = "hv",
  legend            = c(0.75, 0.95),
  ggtheme           = theme_classic(base_size = 13) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"),
          axis.title.y = element_text(margin = margin(r = -2)),
          legend.text = element_text(size = 11))
)


pv <- survminer::surv_pvalue(km_lasso, data = km_lasso_raw)$pval
pval_txt <- paste0("Log-rank p = ", formatC(pv, format = "e", digits = 2))

png(file.path(out_dir, "Fig_3C_km_lasso_legend.png"),
    width = 9, height = 6, units = "in", res = 300)
print(p_km_lasso_legend)
dev.off()

cairo_pdf(file.path(out_dir, "Fig_3C_km_lasso_legend.pdf"),
    width = 9, height = 6)
print(p_km_lasso_legend, newpage = FALSE)
dev.off()

svglite::svglite(file.path(out_dir, "Fig_3C_km_lasso_legend.svg"),
    width = 9, height = 6)
print(p_km_lasso_legend, newpage = FALSE)
dev.off()

cairo_pdf(file.path(out_dir, "Fig_3C_km_lasso.pdf"),
    width = 9, height = 6)
print(p_km_lasso, newpage = FALSE)
dev.off()

svglite::svglite(file.path(out_dir, "Fig_3C_km_lasso.svg"),
    width = 9, height = 6)
print(p_km_lasso, newpage = FALSE)
dev.off()


cairo_pdf(file.path(out_dir, "Fig_3C_km_lasso_square.pdf"),
    width = 9, height = 9)
print(p_km_lasso, newpage = FALSE)
dev.off()

svglite::svglite(file.path(out_dir, "Fig_3C_km_lasso_square.svg"),
    width = 9, height = 9)
print(p_km_lasso, newpage = FALSE)
dev.off()

png(file.path(out_dir, "Fig_3C_km_lasso.png"),
    width = 9, height = 5, units = "in", res = 300)
print(p_km_lasso)
dev.off()


} # end Fig3C individual data guard

cat("\n============================================================\n")
cat("All figures saved to:", out_dir, "\n")
cat("============================================================\n")
