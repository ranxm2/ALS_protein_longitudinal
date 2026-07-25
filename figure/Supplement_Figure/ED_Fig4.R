# ==============================================================================
# ED_Fig4.R — Extended Data Figure 4: Kaplan–Meier Curves for Expanded
#             Protein Panels (Panel 15 and Panel 19)
# ==============================================================================
# Author : Ximing Ran
# Date   : 2026-03-23
#
# Description:
#   Generates Extended Data Figure 4 — Kaplan–Meier conversion-free survival
#   curves for two expanded protein panels evaluated in the Cox model panel
#   comparison analysis:
#   ED_Fig4B - KM curve for Panel 15 (15-protein panel); subjects split at
#              the median linear predictor into Low Risk vs High Risk groups
#   ED_Fig4A - KM curve for Panel 19 (19-protein panel); same median split
#   Both panels include log-rank p-value, 95% CI bands, and at-risk tables.
#
# Input data:
#   data/fig_data/individual_data/Fig3/km_panels/  (optional — each panel
#     skipped if its CSV is absent)
#     Panel_15__15P_/raw_data.csv   — individual time_yr, event, risk_group
#     Panel_19__19P_/raw_data.csv   — individual time_yr, event, risk_group
#
# Output (./results/):
#   ED_Fig4B_km_panel15.{pdf,png}  (only if individual data exists)
#   ED_Fig4A_km_panel19.{pdf,png}  (only if individual data exists)
#
# Prerequisite:
#   Run analysis/03-Time_to_Event_Analysis/3.Cox_model_panel_comparison.Rmd
#   first to generate the panel KM CSV files.
# ==============================================================================

library(survival)
library(survminer)
library(dplyr)
library(here)
library(svglite)

out_dir      <- here::here("figure", "Supplement_Figure", "results")
km_panel_dir <- here("data", "fig_data", "individual_data", "Fig3", "km_panels")

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

palette_km <- c("#427ebf", "#A70C20")

# ------------------------------------------------------------------------------
# Helper: read panel CSV, fit KM, plot with ggsurvplot
# ------------------------------------------------------------------------------
# Helper: load data and fit KM
load_km <- function(panel_folder) {
  csv_path <- file.path(km_panel_dir, panel_folder, "raw_data.csv")
  if (!file.exists(csv_path)) {
    cat("  WARNING: Individual data not found:", csv_path, "\n")
    return(NULL)
  }
  d <- read.csv(csv_path, stringsAsFactors = FALSE) %>%
    mutate(risk_group = factor(risk_group,
      levels = c("Low Risk", "High Risk"),
      labels = c("Low", "High")))
  list(data = d, fit = survfit(Surv(time_yr, event) ~ risk_group, data = d))
}

# Simple version: Low/High, no legend
plot_km_simple <- function(km_obj) {
  ggsurvplot(
    km_obj$fit, data = km_obj$data,
    pval = FALSE,
    # pval = fmt_pval(km_obj),
    conf.int = TRUE, risk.table = TRUE,
    risk.table.height = 0.25,
    risk.table.y.text.col = TRUE,
    tables.theme = theme_classic(base_size = 13) +
      theme(plot.title = element_text(hjust = 0.5, face = "bold")),
    palette = palette_km,
    legend.labs = c("Low", "High"),
    legend.title = "", legend = "none",
    xlab = "Time (years)", ylab = "", title = "",
    xlim = c(0, 15), break.time.by = 5,
    surv.median.line = "hv",
    ggtheme = theme_classic(base_size = 13) +
      theme(axis.title.y = element_text(margin = margin(r = -2)))
  )
}

# Helper: formatted exact log-rank p-value
fmt_pval <- function(km_obj) {
  pv <- survminer::surv_pvalue(km_obj$fit, data = km_obj$data)$pval
  if (pv < 1e-4) {
    paste0("Log-rank p = ", formatC(pv, format = "e", digits = 2))
  } else {
    paste0("Log-rank p = ", formatC(pv, format = "f", digits = 3))
  }
}
# Full legend version: for legend extraction
plot_km_legend <- function(km_obj, legend_prefix) {
  ggsurvplot(
    km_obj$fit, data = km_obj$data,
    pval = FALSE,
    conf.int = TRUE, risk.table = FALSE,
    palette = palette_km,
    legend.labs = c(paste0(legend_prefix, ", lower risk score"),
                    paste0(legend_prefix, ", higher risk score")),
    legend.title = "", legend = c(0.75, 0.95),
    xlab = "Time (years)", ylab = "", title = "",
    xlim = c(0, 15), break.time.by = 5,
    surv.median.line = "hv",
    ggtheme = theme_classic(base_size = 13) +
      theme(axis.title.y = element_text(margin = margin(r = -2)),
            legend.text = element_text(size = 11))
  )
}

# ------------------------------------------------------------------------------
# Panel_15 (15P)
# ------------------------------------------------------------------------------
km_15 <- load_km("Panel_15__15P_")

p_15 <- plot_km_simple(km_15)
png(file.path(out_dir, "ED_Fig4B_km_panel15.png"),
    width = 7, height = 6, units = "in", res = 300)
print(p_15); dev.off()
cairo_pdf(file.path(out_dir, "ED_Fig4B_km_panel15.pdf"),
    width = 7, height = 6)
print(p_15, newpage = FALSE); dev.off()

svglite::svglite(file.path(out_dir, "ED_Fig4B_km_panel15.svg"),
    width = 7, height = 6)
print(p_15, newpage = FALSE); dev.off()

p_15_leg <- plot_km_legend(km_15, "15-protein panel")
png(file.path(out_dir, "ED_Fig4B_km_panel15_legend.png"),
    width = 7, height = 6, units = "in", res = 300)
print(p_15_leg); dev.off()
cairo_pdf(file.path(out_dir, "ED_Fig4B_km_panel15_legend.pdf"),
    width = 7, height = 6)
print(p_15_leg, newpage = FALSE); dev.off()

svglite::svglite(file.path(out_dir, "ED_Fig4B_km_panel15_legend.svg"),
    width = 7, height = 6)
print(p_15_leg, newpage = FALSE); dev.off()

# ------------------------------------------------------------------------------
# Panel_19 (19P)
# ------------------------------------------------------------------------------
km_19 <- load_km("Panel_19__19P_")

p_19 <- plot_km_simple(km_19)
png(file.path(out_dir, "ED_Fig4A_km_panel19.png"),
    width = 7, height = 6, units = "in", res = 300)
print(p_19); dev.off()
cairo_pdf(file.path(out_dir, "ED_Fig4A_km_panel19.pdf"),
    width = 7, height = 6)
print(p_19, newpage = FALSE); dev.off()

svglite::svglite(file.path(out_dir, "ED_Fig4A_km_panel19.svg"),
    width = 7, height = 6)
print(p_19, newpage = FALSE); dev.off()

p_19_leg <- plot_km_legend(km_19, "19-protein panel")
png(file.path(out_dir, "ED_Fig4A_km_panel19_legend.png"),
    width = 7, height = 6, units = "in", res = 300)
print(p_19_leg); dev.off()
cairo_pdf(file.path(out_dir, "ED_Fig4A_km_panel19_legend.pdf"),
    width = 7, height = 6)
print(p_19_leg, newpage = FALSE); dev.off()

svglite::svglite(file.path(out_dir, "ED_Fig4A_km_panel19_legend.svg"),
    width = 7, height = 6)
print(p_19_leg, newpage = FALSE); dev.off()

cat("\n============================================================\n")
cat("All figures saved to:", out_dir, "\n")
cat("============================================================\n")


pv <- survminer::surv_pvalue(km_lasso, data = km_lasso_raw)$pval
pval_txt <- paste0("Log-rank p = ", formatC(pv, format = "e", digits = 2))

