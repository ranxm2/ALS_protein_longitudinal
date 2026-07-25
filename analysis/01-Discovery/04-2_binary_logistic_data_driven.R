# =============================================================================
# Script 4-2: Binary Logistic Regression + Greedy Search (data-driven)
# Selects proteins maximizing MEAN CV AUROC across all T values simultaneously
#
# Requires: data/intermediate/ (run 01-1 first)
# Inputs:   data/intermediate/final_data_matrix.rds
#           data/intermediate/Mix_effect_model_lmer.csv
#           data/df_visit_info.csv
# Outputs:  results/04-2_binary_logistic_data_driven/
# =============================================================================

# Set working directory to analysis/discovery/
# In RStudio: setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# From command line: Rscript should be run from this directory
if (requireNamespace("rstudioapi", quietly=TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# ── Libraries ─────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(tidyverse)
  library(pROC)
  library(ggplot2)
})

set.seed(2025)

# ── Settings ──────────────────────────────────────────────────────────────────
T_list         <- c(0.5, 1.0, 2.0, 3.0, 5.0)
k              <- 5            # CV folds
greedy_max_steps <- Inf        # Set finite to stop early
min_delta_auc    <- 0          # Minimum AUC improvement to accept protein

# ── Output directories ────────────────────────────────────────────────────────
out_dir <- "results/04-2_binary_logistic_data_driven"
fig_dir <- file.path(out_dir, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ── Decode visit info ─────────────────────────────────────────────────────────
decode_visit_info <- function(vi) {
  vi %>%
    mutate(Sex = ifelse(Male == 1, "Male", "Female")) %>%
    filter(Group != 1)   # exclude Pre-symptomatic (Group == 1)
}

# =============================================================================
# 1. Load Data
# =============================================================================
cat("Loading data...\n")

protein_matrix <- readRDS(
  "data/intermediate/final_data_matrix.rds"
)

results_protein <- read.csv(
  "data/intermediate/Mix_effect_model_lmer.csv"
)
protein_list <- results_protein %>%
  filter(significant == "Significant") %>% pull(Protein)

# Restrict to proteins present in matrix
candidate_proteins <- intersect(protein_list, colnames(protein_matrix))
protein_matrix     <- protein_matrix %>% select(SampleID, all_of(candidate_proteins))

visit_info <- read.csv("data/df_visit_info.csv") %>%
  decode_visit_info() %>%
  mutate(
    Group = factor(Group, levels = c(0, 2, 3)),
    Geno  = factor(Geno,  levels = c(1, 2, 3, 4, 9))
  )

cat(sprintf("Candidate proteins: %d\n", length(candidate_proteins)))
if (length(candidate_proteins) < 2)
  stop("Too few candidate proteins (need >= 2).")

# =============================================================================
# 2. Prepare Binary Datasets (one per T)
# =============================================================================
covariate_list <- c("X_Male", "X_Age", "X_Geno_SOD1_A4V",
                    "X_Geno_SOD1_non_A4V", "X_Geno_C9orf72")

prepare_binary_data <- function(visit_info, protein_matrix, T, k = 5, seed = 2025) {
  cases      <- visit_info %>%
    filter(Group == 2, TimetoEvent >= -T, TimetoEvent < 0) %>%
    mutate(CTRL_Group = "Convert")
  ctrls      <- visit_info %>%
    filter(Group == 0) %>%
    mutate(CTRL_Group = "Control")
  ctrls_conv <- visit_info %>%
    filter(Group == 2, TimetoEvent < -T) %>%
    mutate(CTRL_Group = "Control")

  combined <- bind_rows(cases, ctrls, ctrls_conv)
  set.seed(seed)
  fold_df  <- tibble(EID  = unique(combined$EID),
                     Fold = sample(rep(seq_len(k),
                                       length.out = length(unique(combined$EID)))))

  combined %>%
    left_join(fold_df, by = "EID") %>%
    mutate(
      CASE                = ifelse(CTRL_Group == "Control", 0L, 1L),
      X_Male              = as.integer(Male == 1),
      X_Age               = Age,
      X_UID               = as.factor(EID),
      X_Geno_SOD1_A4V     = as.integer(Geno == 1),
      X_Geno_SOD1_non_A4V = as.integer(Geno == 2),
      X_Geno_C9orf72      = as.integer(Geno == 3)
    ) %>%
    left_join(protein_matrix, by = "SampleID") %>%
    drop_na()
}

cat("Preparing binary datasets for all T values...\n")
cleaned_data_list <- lapply(T_list, function(T) {
  cd <- prepare_binary_data(visit_info, protein_matrix, T, k)
  cat(sprintf("  T=%.1f: N=%d, cases=%d\n", T, nrow(cd), sum(cd$CASE)))
  cd
})
names(cleaned_data_list) <- as.character(T_list)

# =============================================================================
# 3. AUROC Helpers
# =============================================================================
cv_auc_oneT <- function(cleaned_data, protein_set, k = 5) {
  x_cols  <- c(protein_set, covariate_list)
  x_mat   <- as.matrix(cleaned_data[, x_cols, drop = FALSE])
  X       <- cbind(`(Intercept)` = 1, x_mat)
  y       <- cleaned_data$CASE
  folds   <- cleaned_data$Fold

  all_preds <- numeric(0); all_true <- numeric(0)
  for (i in seq_len(k)) {
    ti <- which(folds == i); tri <- which(folds != i)
    fit  <- glm.fit(x = X[tri, , drop=FALSE], y = y[tri], family = binomial())
    beta <- fit$coefficients; beta[is.na(beta)] <- 0
    all_preds <- c(all_preds, plogis(drop(X[ti, , drop=FALSE] %*% beta)))
    all_true  <- c(all_true, y[ti])
  }
  suppressMessages(as.numeric(auc(roc(response   = all_true,
                                      predictor  = all_preds))))
}

cv_auc_avgT <- function(protein_set) {
  aucs <- vapply(as.character(T_list), function(t) {
    cv_auc_oneT(cleaned_data_list[[t]], protein_set, k)
  }, numeric(1))
  mean(aucs, na.rm = TRUE)
}

# =============================================================================
# 4. Step 1: Single-Protein Screening (mean AUC over all T)
# =============================================================================
cat("\n--- Step 1: Single-protein screening (mean AUROC over all T) ---\n")

single_auc_df <- tibble(Protein = candidate_proteins, MeanAUC = NA_real_)

for (idx in seq_along(candidate_proteins)) {
  p <- candidate_proteins[idx]
  single_auc_df$MeanAUC[idx] <- cv_auc_avgT(p)
  cat(sprintf("[%d/%d] %s: MeanAUC=%.4f\n",
              idx, length(candidate_proteins), p,
              single_auc_df$MeanAUC[idx]))
}

single_auc_df <- single_auc_df %>% arrange(desc(MeanAUC))
write.csv(single_auc_df,
          file.path(out_dir, "01-single_protein_mean_auroc.csv"),
          row.names = FALSE)

p_single <- single_auc_df %>% slice_head(n = 30) %>%
  ggplot(aes(x = reorder(Protein, MeanAUC), y = MeanAUC)) +
  geom_bar(stat = "identity", fill = "steelblue", color = "black") +
  geom_text(aes(label = sprintf("%.3f", MeanAUC)), hjust = -0.1, size = 3.5) +
  labs(title = "Top 30 Proteins: Mean AUROC Across All T",
       x = "Protein", y = "Mean CV AUROC") +
  coord_flip() + theme_classic() + ylim(0, 1.05)
ggsave(file.path(fig_dir, "single_protein_mean_auroc.pdf"),
       p_single, width = 6, height = 8)

# =============================================================================
# 5. Step 2: Greedy Forward Selection (maximize mean AUROC)
# =============================================================================
cat("\n--- Step 2: Greedy forward selection (mean AUROC objective) ---\n")

selected_proteins <- single_auc_df$Protein[1]
best_mean_auc     <- single_auc_df$MeanAUC[1]
greedy_log        <- tibble(Step = 1L, Protein = selected_proteins,
                            MeanAUC = best_mean_auc)
remaining_proteins <- setdiff(candidate_proteins, selected_proteins)
step <- 1L

while (length(remaining_proteins) > 0 && step <= greedy_max_steps) {
  step <- step + 1L
  cat(sprintf("  Step %d: searching %d candidates...\n",
              step, length(remaining_proteins)))

  candidate_aucs <- setNames(numeric(length(remaining_proteins)),
                              remaining_proteins)
  for (prot in remaining_proteins)
    candidate_aucs[prot] <- cv_auc_avgT(c(selected_proteins, prot))

  best_prot     <- names(which.max(candidate_aucs))
  best_new_auc  <- candidate_aucs[best_prot]

  if (best_new_auc - best_mean_auc < min_delta_auc) {
    cat(sprintf("  No improvement (delta=%.4f < %.4f). Stopping.\n",
                best_new_auc - best_mean_auc, min_delta_auc))
    break
  }

  selected_proteins  <- c(selected_proteins, best_prot)
  remaining_proteins <- setdiff(remaining_proteins, best_prot)
  best_mean_auc      <- best_new_auc

  greedy_log <- bind_rows(greedy_log,
                           tibble(Step    = step,
                                  Protein = best_prot,
                                  MeanAUC = best_new_auc))
  cat(sprintf("    -> Added: %s (MeanAUC=%.4f)\n", best_prot, best_new_auc))

  write.csv(greedy_log,
            file.path(out_dir, "02-greedy_selection_log.csv"),
            row.names = FALSE)

  rm(candidate_aucs); gc(verbose = FALSE)
}

write.csv(greedy_log,
          file.path(out_dir, "02-greedy_selection_log.csv"),
          row.names = FALSE)

# Greedy curve plot
max_step <- which.max(greedy_log$MeanAUC)
p_greedy <- ggplot(greedy_log, aes(x = Step, y = MeanAUC)) +
  geom_point() + geom_line(linetype = "dotted") +
  geom_point(data = greedy_log[max_step, ], color = "red", size = 3) +
  geom_text(data = greedy_log[max_step, ],
            aes(label = sprintf("%.3f\n%s", MeanAUC, Protein)),
            vjust = -0.5, size = 3.5) +
  labs(title = "Greedy Selection: Mean AUROC Across All T",
       x = "Number of Proteins Added", y = "Mean CV AUROC") +
  theme_classic()
ggsave(file.path(fig_dir, "greedy_mean_auroc_curve.pdf"),
       p_greedy, width = 12, height = 4)

# =============================================================================
# 6. Evaluate Final Optimal Panel Per T
# =============================================================================
cat("\n--- Evaluating optimal panel per T ---\n")

optimal_panel <- greedy_log$Protein[1:max_step]
cat(sprintf("Optimal panel: %d proteins\n  %s\n",
            length(optimal_panel), paste(optimal_panel, collapse = ", ")))

per_T_auc <- tibble(T = T_list, AUC = NA_real_)
for (tt in seq_along(T_list)) {
  T_val <- T_list[tt]
  per_T_auc$AUC[tt] <- cv_auc_oneT(cleaned_data_list[[as.character(T_val)]],
                                    optimal_panel, k)
  cat(sprintf("  T=%.1f: AUC=%.4f\n", T_val, per_T_auc$AUC[tt]))
}

write.csv(per_T_auc,
          file.path(out_dir, "03-optimal_panel_per_T_auroc.csv"),
          row.names = FALSE)

p_perT <- ggplot(per_T_auc, aes(x = factor(T), y = AUC)) +
  geom_bar(stat = "identity", fill = "#2ca25f", color = "black", width = 0.6) +
  geom_text(aes(label = sprintf("%.3f", AUC)), vjust = -0.3, size = 4) +
  labs(title = sprintf("Optimal Panel (%d proteins) AUROC per T",
                       length(optimal_panel)),
       x = "Time Horizon T (years)", y = "CV AUROC") +
  ylim(0, 1.05) + theme_classic()
ggsave(file.path(fig_dir, "optimal_panel_per_T.pdf"), p_perT, width = 6, height = 5)

# Summary report
writeLines(
  c(
    sprintf("Optimal protein panel (%d proteins):", length(optimal_panel)),
    paste("  ", optimal_panel, collapse = "\n"),
    "",
    "Per-T AUROC:",
    apply(per_T_auc, 1, function(r)
      sprintf("  T=%.1f yr: AUC=%.4f", as.numeric(r[1]), as.numeric(r[2])))
  ),
  file.path(out_dir, "04-summary_report.txt")
)

cat("\n=== Script 4-2 complete ===\n")
cat("Outputs in:", out_dir, "\n")
