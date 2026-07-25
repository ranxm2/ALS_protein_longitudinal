# =============================================================================
# Script 5: Phenoconversion Timing Estimation
# Truncated regression (Tobit) to predict years-to-onset from delta values
#
# Requires: data/intermediate/ (run 01-1 and 02-1 first)
# Inputs:   data/intermediate/delta_matrix_wide.csv
#           data/intermediate/Mix_effect_model_lmer.csv
#           data/df_visit_info.csv
# Outputs:  results/05_phenoconversion_timing/
# =============================================================================

# Set working directory to analysis/discovery/
# In RStudio: setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# From command line: Rscript should be run from this directory
if (requireNamespace("rstudioapi", quietly=TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# ── Libraries ─────────────────────────────────────────────────────────────────
library(tidyverse)
library(AER)       # tobit()
library(patchwork)
library(ggplot2)
library(cowplot)

set.seed(2024)

# ── Output directories ────────────────────────────────────────────────────────
out_dir <- "results/05_phenoconversion_timing"
fig_dir <- file.path(out_dir, "figures")
rds_dir <- file.path(out_dir, "models")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(rds_dir, recursive = TRUE, showWarnings = FALSE)

# ── Decode visit info ─────────────────────────────────────────────────────────
decode_visit_info <- function(vi) {
  vi %>%
    mutate(Sex = ifelse(Male == 1, "Male", "Female")) %>%
    filter(Group != 1)   # exclude Pre-symptomatic (Group == 1)
}

# ── Diagnostic plot helpers ───────────────────────────────────────────────────
Bland_Altman_plot <- function(y_true, y_predict, subtitle = NULL) {
  mean_v <- (y_true + y_predict) / 2
  diff_v <- y_true - y_predict
  md     <- mean(diff_v, na.rm = TRUE)
  sd_d   <- sd(diff_v,   na.rm = TRUE)
  ggplot(data.frame(mean_v, diff_v), aes(x = mean_v, y = diff_v)) +
    geom_point(alpha = 0.7) +
    geom_hline(yintercept = md,              color = "red",  linetype = "dashed") +
    geom_hline(yintercept = md + 1.96 * sd_d, color = "blue", linetype = "dashed") +
    geom_hline(yintercept = md - 1.96 * sd_d, color = "blue", linetype = "dashed") +
    labs(title = "Bland-Altman Plot", subtitle = subtitle,
         x = "Mean of Observed and Predicted", y = "Difference") +
    theme_classic(base_size = 11)
}

residual_hist <- function(y_true, y_predict, subtitle = NULL) {
  ggplot(data.frame(res = y_true - y_predict), aes(x = res)) +
    geom_histogram(binwidth = 0.5, fill = "lightblue", color = "black") +
    labs(title = "Residual Histogram", subtitle = subtitle,
         x = "Residuals (Observed - Predicted)", y = "Frequency") +
    theme_classic(base_size = 11)
}

predict_vs_observed_plot_tobit <- function(y_true, y_predict, subtitle = NULL,
                                            MAE = NA, RMSE = NA,
                                            cor_y = NA, R_squared = NA) {
  df <- data.frame(y_true, y_predict)
  df$category <- ifelse(y_predict < y_true,
                        "Onset Earlier than Predicted",
                        "Onset Later than Predicted")
  rl  <- range(c(y_true, y_predict), na.rm = TRUE)
  bx  <- seq(rl[1], rl[2], length.out = 100)
  bdf <- data.frame(x = c(bx, rev(bx)), y = c(bx - 1, rev(bx + 1)))
  prop1 <- mean(abs(y_true - y_predict) < 1, na.rm = TRUE)

  ggplot(df, aes(x = y_predict, y = y_true, color = category)) +
    geom_polygon(data = bdf, aes(x = x, y = y),
                 fill = "grey80", alpha = 0.5, inherit.aes = FALSE) +
    geom_point(alpha = 0.7, size = 2) +
    scale_color_manual(values = c("Onset Later than Predicted"    = "blue",
                                  "Onset Earlier than Predicted" = "red")) +
    geom_line(data = data.frame(x = bx, y = bx),
              aes(x = x, y = y), color = "black",
              linetype = "dashed", inherit.aes = FALSE) +
    labs(title = "Predicted vs. Observed", subtitle = subtitle,
         x = "Predicted Values", y = "Observed Values") +
    theme_classic(base_size = 11) +
    theme(legend.position = "top") +
    annotate("text", x = rl[1], y = rl[2], hjust = 0, vjust = 1,
             label = sprintf("MAE=%.3f\nRMSE=%.3f\n|Res|<1=%.1f%%\nCor=%.3f\nR²=%.3f",
                             MAE, RMSE, prop1 * 100, cor_y, R_squared),
             size = 4, fontface = "bold")
}

predict_vs_observed_boxplot <- function(y_true, y_predict, subtitle = NULL) {
  df <- data.frame(y_true, y_predict, bin = round(y_predict))
  df$category <- ifelse(y_predict < y_true,
                        "Onset Earlier than Predicted",
                        "Onset Later than Predicted")
  ggplot(df, aes(x = factor(bin), y = y_true)) +
    geom_boxplot(fill = "grey", alpha = 0.6, outlier.shape = NA) +
    geom_jitter(aes(color = category), width = 0.2, alpha = 0.7, size = 1.5) +
    scale_color_manual(values = c("Onset Later than Predicted"    = "blue",
                                  "Onset Earlier than Predicted" = "red")) +
    labs(title = "Predicted vs. Observed (Binned)", subtitle = subtitle,
         x = "Predicted (Rounded)", y = "Observed") +
    theme_classic(base_size = 11) + theme(legend.position = "top")
}

residual_vs_predict <- function(y_true, y_predict, subtitle = NULL) {
  res  <- y_true - y_predict
  df   <- data.frame(y_predict, res,
                     cat = ifelse(res > 0, "Onset Earlier than Predicted",
                                  "Onset Later than Predicted"))
  ir   <- seq(floor(min(y_predict)), ceiling(max(y_predict)), by = 1)
  ggplot(df, aes(x = y_predict, y = res, color = cat)) +
    geom_vline(xintercept = ir, color = "grey70", linetype = "dotted") +
    geom_point(alpha = 0.7, size = 2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c("Onset Later than Predicted"    = "blue",
                                  "Onset Earlier than Predicted" = "red")) +
    labs(title = "Residuals vs. Predicted", subtitle = subtitle,
         x = "Predicted", y = "Residuals") +
    theme_classic(base_size = 11) + theme(legend.position = "top")
}

residual_vs_predict_bin <- function(y_true, y_predict, subtitle = NULL) {
  res <- y_true - y_predict
  df  <- data.frame(y_predict, res, bin = round(y_predict),
                    cat = ifelse(res > 0, "Onset Earlier than Predicted",
                                  "Onset Later than Predicted"))
  ggplot(df, aes(x = factor(bin), y = res)) +
    geom_boxplot(fill = "gray80", alpha = 0.6, outlier.shape = NA) +
    geom_jitter(aes(color = cat), width = 0.2, alpha = 0.7, size = 1.5) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_color_manual(values = c("Onset Later than Predicted"    = "blue",
                                  "Onset Earlier than Predicted" = "red")) +
    labs(title = "Residuals vs. Predicted (Binned)", subtitle = subtitle,
         x = "Predicted (Rounded)", y = "Residuals") +
    theme_classic(base_size = 11) + theme(legend.position = "top")
}

# ── Fit Tobit and collect metrics ────────────────────────────────────────────
fit_tobit_and_eval <- function(df_input, formula, model_name) {
  rg <- tryCatch(
    tobit(formula, left = -Inf, right = 0, data = df_input),
    error = function(e) { message("Tobit failed: ", e$message); NULL }
  )
  if (is.null(rg)) return(NULL)

  y_true    <- df_input$Y
  y_predict <- pmin(predict(rg, newdata = df_input, type = "response"), 0)
  residuals <- y_true - y_predict
  y_mean    <- mean(y_true, na.rm = TRUE)

  MAE       <- mean(abs(residuals),   na.rm = TRUE)
  RMSE      <- sqrt(mean(residuals^2, na.rm = TRUE))
  AIC_val   <- AIC(rg)
  cor_y     <- cor(y_true, y_predict, use = "complete.obs")
  R_squared <- 1 - sum(residuals^2, na.rm = TRUE) /
                   sum((y_true - y_mean)^2, na.rm = TRUE)

  list(model = rg, y_true = y_true, y_predict = y_predict,
       name = model_name, MAE = MAE, RMSE = RMSE, AIC = AIC_val,
       cor = cor_y, R2 = R_squared)
}

# =============================================================================
# 1. Load Data
# =============================================================================
cat("Loading data...\n")

protein_matrix <- read.csv(
  "data/intermediate/delta_matrix_wide.csv"
)

results_protein <- read.csv(
  "data/intermediate/Mix_effect_model_lmer.csv"
)
protein_list_sig <- results_protein %>%
  filter(significant == "Significant") %>% pull(Protein)

# Restrict to proteins in matrix
protein_list_sig <- intersect(protein_list_sig, colnames(protein_matrix))
protein_matrix   <- protein_matrix %>% select(SampleID, all_of(protein_list_sig))

visit_info <- read.csv("data/df_visit_info.csv") %>%
  decode_visit_info() %>%
  mutate(
    Group = factor(Group, levels = c(0, 2, 3)),
    Geno  = factor(Geno,  levels = c(1, 2, 3, 4, 9)),
    Male  = factor(Male,  levels = c(0, 1)),
    Sex   = factor(Sex,   levels = c("Female", "Male"))
  )

# =============================================================================
# 2. Define Protein Panels
# =============================================================================
Panel_1  <- c("NEFL")
Panel_2  <- c("NEFL", "DUSP29")
Panel_3  <- c("NEFL", "DUSP29", "CALCA")
Panel_6  <- c("NEFL", "DUSP29", "CALCA", "NEB", "NOS1", "EDA2R")
Panel_15 <- c("NEFL", "DUSP29", "CALCA", "NEB", "NGRN",
              "NOS1", "DTNB", "MEGF10", "APOA4", "ART3",
              "MYL3", "EPHA1", "EDA2R", "TTN", "ACTN2")
Panel_19 <- c("NEFL", "DUSP29", "CALCA", "NEB", "NGRN",
              "MYL11", "NOS1", "MYH1", "TNNC1", "DTNB",
              "MEGF10", "SYNM", "APOA4", "ART3", "MYL3",
              "EPHA1", "EDA2R", "TTN", "ACTN2")
Panel_data <- c("NEFL", "CALCA", "DUSP29", "NGRN", "CD38",
                "MYH1", "SPINK7", "MSANTD3", "MOCS3", "PPP6R2",
                "CRTAC1", "EPHA1", "FABP2", "NOS1", "SAV1", "MYL1")
Panel_7_cox <- c("NEFL", "DUSP29", "CALCB", "NDE1", "NFAT5",
                 "CLEC4G", "PPFIBP1")

Protein_list_set <- list(
  Panel_1 = Panel_1, Panel_2 = Panel_2, Panel_3 = Panel_3,
  Panel_6 = Panel_6, Panel_15 = Panel_15, Panel_19 = Panel_19,
  Panel_data = Panel_data, Panel_7_cox = Panel_7_cox
)

# Keep only proteins present in delta matrix
Protein_list_set <- lapply(Protein_list_set, function(pl) {
  intersect(pl, colnames(protein_matrix))
})

full_protein_list <- unique(unlist(Protein_list_set))
cat(sprintf("Full protein universe: %d proteins\n", length(full_protein_list)))

# =============================================================================
# 3. Data Preparation
# =============================================================================
visit_info_conv <- visit_info %>% filter(Group == 2)
visit_info_conv_filter <- visit_info_conv %>% filter(TimetoEvent < 0)

# Initial visit protein values per individual
initial_visit_info_conv <- visit_info_conv_filter %>%
  group_by(EID) %>%
  filter(TimetoEvent == min(TimetoEvent)) %>%
  ungroup() %>%
  select(SampleID, EID)

df_init_conv <- protein_matrix %>%
  filter(SampleID %in% initial_visit_info_conv$SampleID) %>%
  left_join(initial_visit_info_conv, by = "SampleID") %>%
  rename_with(~ paste0(.x, "_Init"), !c(SampleID, EID)) %>%
  select(-SampleID) %>%
  relocate(EID)

# Main input dataset (Convert, pre-onset visits)
df_Int_input <- protein_matrix %>%
  filter(SampleID %in% visit_info_conv_filter$SampleID) %>%
  select(SampleID, all_of(full_protein_list)) %>%
  left_join(visit_info_conv_filter, by = "SampleID") %>%
  left_join(df_init_conv, by = "EID") %>%
  mutate(
    Y                    = TimetoEvent,
    X_Male               = ifelse(Male == 1, 1, 0),
    X_Age                = Age,
    X_UID                = as.factor(EID),
    X_Geno_SOD1_A4V      = ifelse(Geno == 1, 1, 0),
    X_Geno_SOD1_non_A4V  = ifelse(Geno == 2, 1, 0),
    X_Geno_C9orf72       = ifelse(Geno == 3, 1, 0)
  )

covariate_list <- c("X_Male", "X_Geno_SOD1_A4V",
                    "X_Geno_SOD1_non_A4V", "X_Geno_C9orf72")

cat(sprintf("Convert pre-onset visits: %d\n", nrow(df_Int_input)))

# =============================================================================
# 4. Single-Protein Models
# =============================================================================
cat("\n--- Single-protein Tobit models ---\n")

result_single <- data.frame(
  Protein     = full_protein_list,
  N_protein   = 1L,
  RMSE        = NA_real_,
  MAE         = NA_real_,
  AIC         = NA_real_,
  Correlation = NA_real_,
  R2_Simple   = NA_real_
)

pdf(file.path(fig_dir, "single_protein_diagnostics.pdf"), width = 14, height = 20)
for (i in seq_along(full_protein_list)) {
  protein       <- full_protein_list[i]
  protein_init  <- paste0(protein, "_Init")
  cat(sprintf("  [%d/%d] %s\n", i, length(full_protein_list), protein))

  # Check required columns exist
  if (!protein %in% colnames(df_Int_input) ||
      !protein_init %in% colnames(df_Int_input)) {
    cat(sprintf("    Skipping (missing column)\n")); next
  }

  formula_str <- paste("Y ~", protein, "+", protein_init, "+",
                       paste(covariate_list, collapse = " + "))
  res <- fit_tobit_and_eval(df_Int_input, as.formula(formula_str),
                             sprintf("Single: %s", protein))
  if (is.null(res)) next

  result_single[i, "RMSE"]        <- res$RMSE
  result_single[i, "MAE"]         <- res$MAE
  result_single[i, "AIC"]         <- res$AIC
  result_single[i, "Correlation"] <- res$cor
  result_single[i, "R2_Simple"]   <- res$R2

  saveRDS(res$model, file.path(rds_dir, sprintf("tobit_single_%s.rds", protein)))

  # Diagnostic plots (6-panel)
  p1 <- Bland_Altman_plot(res$y_true, res$y_predict, protein)
  p2 <- residual_hist(res$y_true, res$y_predict, protein)
  p3 <- predict_vs_observed_plot_tobit(
    res$y_true, res$y_predict, protein,
    MAE = res$MAE, RMSE = res$RMSE, cor_y = res$cor, R_squared = res$R2)
  p4 <- predict_vs_observed_boxplot(res$y_true, res$y_predict, protein)
  p5 <- residual_vs_predict(res$y_true, res$y_predict, protein)
  p6 <- residual_vs_predict_bin(res$y_true, res$y_predict, protein)
  print(((p1 | p2) / (p3 | p4) / (p5 | p6)) + plot_layout(heights = c(1, 1, 1)))
}
dev.off()

result_single <- result_single %>% arrange(MAE)
write.csv(result_single,
          file.path(out_dir, "single_protein_results.csv"),
          row.names = FALSE)
cat("Saved single-protein results.\n")

# =============================================================================
# 5. Multi-Protein Panel Models
# =============================================================================
cat("\n--- Multi-protein Tobit models ---\n")

result_multi <- data.frame(
  Panel       = names(Protein_list_set),
  N_protein   = sapply(Protein_list_set, length),
  RMSE        = NA_real_,
  MAE         = NA_real_,
  AIC         = NA_real_,
  Correlation = NA_real_,
  R2_Simple   = NA_real_
)

pdf(file.path(fig_dir, "multi_protein_diagnostics.pdf"), width = 14, height = 20)
for (i in seq_along(Protein_list_set)) {
  panel_name   <- names(Protein_list_set)[i]
  panel_prots  <- Protein_list_set[[i]]
  panel_inits  <- paste0(panel_prots, "_Init")

  # Filter to available proteins
  avail     <- panel_prots[panel_prots %in% colnames(df_Int_input)]
  avail_ini <- panel_inits[panel_inits %in% colnames(df_Int_input)]
  if (length(avail) == 0) { cat(sprintf("  Skipping %s\n", panel_name)); next }

  cat(sprintf("  [%d/%d] %s (%d proteins)\n", i, length(Protein_list_set),
              panel_name, length(avail)))

  formula_str <- paste(
    "Y ~",
    paste(avail, collapse = " + "), "+",
    if (length(avail_ini) > 0) paste(avail_ini, collapse = " + ") else "1",
    "+", paste(covariate_list, collapse = " + ")
  )

  res <- fit_tobit_and_eval(df_Int_input, as.formula(formula_str), panel_name)
  if (is.null(res)) next

  result_multi[i, "RMSE"]        <- res$RMSE
  result_multi[i, "MAE"]         <- res$MAE
  result_multi[i, "AIC"]         <- res$AIC
  result_multi[i, "Correlation"] <- res$cor
  result_multi[i, "R2_Simple"]   <- res$R2

  saveRDS(res$model, file.path(rds_dir, sprintf("tobit_multi_%s.rds", panel_name)))

  # Diagnostic plots
  p1 <- Bland_Altman_plot(res$y_true, res$y_predict, panel_name)
  p2 <- residual_hist(res$y_true, res$y_predict, panel_name)
  p3 <- predict_vs_observed_plot_tobit(
    res$y_true, res$y_predict, panel_name,
    MAE = res$MAE, RMSE = res$RMSE, cor_y = res$cor, R_squared = res$R2)
  p4 <- predict_vs_observed_boxplot(res$y_true, res$y_predict, panel_name)
  p5 <- residual_vs_predict(res$y_true, res$y_predict, panel_name)
  p6 <- residual_vs_predict_bin(res$y_true, res$y_predict, panel_name)
  print(((p1 | p2) / (p3 | p4) / (p5 | p6)) + plot_layout(heights = c(1, 1, 1)))
}
dev.off()

result_multi <- result_multi %>% arrange(MAE)
write.csv(result_multi,
          file.path(out_dir, "multi_protein_panel_results.csv"),
          row.names = FALSE)

# =============================================================================
# 6. Comparison Summary
# =============================================================================
cat("\n--- Performance comparison ---\n")
comparison <- bind_rows(
  result_single %>% mutate(Type = "Single") %>% rename(Panel = Protein),
  result_multi   %>% mutate(Type = "Multi")
) %>% arrange(MAE)

write.csv(comparison,
          file.path(out_dir, "model_comparison_summary.csv"),
          row.names = FALSE)

p_compare <- comparison %>%
  filter(!is.na(MAE)) %>%
  mutate(Panel = factor(Panel, levels = rev(Panel))) %>%
  ggplot(aes(x = Panel, y = MAE, fill = Type)) +
  geom_bar(stat = "identity", color = "black") +
  geom_text(aes(label = round(MAE, 3)), hjust = -0.1, size = 3.5) +
  coord_flip() +
  scale_fill_manual(values = c(Single = "#bfe5a8", Multi = "#74c476")) +
  labs(title = "Model Comparison: MAE (lower = better)",
       x = "Model", y = "Mean Absolute Error") +
  theme_classic() +
  scale_y_continuous(expand = expansion(mult = c(0, 0.15)))

ggsave(file.path(fig_dir, "model_comparison_MAE.pdf"),
       p_compare, width = 8, height = max(6, nrow(comparison) * 0.3))

cat("\nTop 5 models by MAE:\n")
print(head(comparison %>% select(Panel, Type, MAE, RMSE, Correlation), 5))

cat("\n=== Script 5 complete ===\n")
cat("Outputs in:", out_dir, "\n")
