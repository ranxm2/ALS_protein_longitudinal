# =============================================================================
# Script 4-1: Binary Logistic Regression + Greedy Search (per-T optimization)
# For each time horizon T, greedily select proteins maximizing CV AUROC at T
#
# Requires: data/intermediate/ (run 01-1 first)
# Inputs:   data/intermediate/final_data_matrix.rds
#           results/01_differential_expression/Mix_effect_model_lmer.csv
#           data/df_visit_info.csv
# Outputs:
#   results/04-1_binary_logistic_local/  (summary, greedy results)
#   data/intermediate/                   (individual-level binary data)
# =============================================================================

if (requireNamespace("rstudioapi", quietly=TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# ── Libraries ─────────────────────────────────────────────────────────────────
library(tidyverse)
library(progress)

set.seed(2025)

# ── Output directories ────────────────────────────────────────────────────────
out_dir  <- "results/04-1_binary_logistic_local"
int_dir  <- "data/intermediate"
fold_dir <- "data/Binary"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(int_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Load Data
# =============================================================================
cat("Loading data...\n")

protein_matrix <- readRDS("data/intermediate/final_data_matrix.rds")

results_protein <- read.csv(
  "results/01_differential_expression/Mix_effect_model_lmer.csv"
)
protein_list <- results_protein %>%
  filter(significant == "Significant") %>% pull(Protein)

protein_list <- intersect(protein_list, colnames(protein_matrix))
protein_matrix <- protein_matrix %>% select(SampleID, all_of(protein_list))

visit_info <- read.csv("data/df_visit_info.csv") %>%
  filter(Group != 1) %>%
  mutate(
    Group = factor(Group, levels = c(0, 2, 3)),
    Geno  = factor(Geno,  levels = c(1, 2, 3, 4, 9))
  )

cat(sprintf("Significant proteins in matrix: %d\n", length(protein_list)))

# =============================================================================
# 2. Binary Data Construction
# =============================================================================
prepare_binary_data <- function(visit_info, protein_matrix, T_val, k = 5,
                                 fold_dir = "data/Binary") {
  cases <- visit_info %>%
    filter(Group == 2, TimetoEvent >= -T_val, TimetoEvent < 0) %>%
    mutate(CTRL_Group = "Convert")

  ctrls <- visit_info %>%
    filter(Group == 0) %>%
    mutate(CTRL_Group = "Control")

  ctrls_conv <- visit_info %>%
    filter(Group == 2, TimetoEvent < -T_val) %>%
    mutate(CTRL_Group = "Control")

  combined <- bind_rows(cases, ctrls, ctrls_conv)

  fold_file <- file.path(fold_dir,
                          sprintf("5_Fold_mannual_first_onset_%.1f.csv", T_val))
  fold_df <- read.csv(fold_file) %>% rename(EID = Row.Labels)
  cat(sprintf("  Loaded fold partition: %s (%d individuals)\n",
              basename(fold_file), nrow(fold_df)))

  combined <- combined %>% left_join(fold_df, by = "EID")

  unassigned_eids <- unique(combined$EID[is.na(combined$Fold)])
  if (length(unassigned_eids) > 0) {
    set.seed(2025)
    extra_folds <- tibble(
      EID  = unassigned_eids,
      Fold = sample(rep(seq_len(k), length.out = length(unassigned_eids)))
    )
    combined <- combined %>%
      mutate(Fold = ifelse(is.na(Fold),
                           extra_folds$Fold[match(EID, extra_folds$EID)],
                           Fold))
  }

  combined <- combined %>%
    mutate(
      CASE                = as.integer(CTRL_Group != "Control"),
      X_Male              = as.integer(Male == 1),
      X_Age               = Age,
      X_Geno_SOD1_A4V     = as.integer(Geno == 1),
      X_Geno_SOD1_non_A4V = as.integer(Geno == 2),
      X_Geno_C9orf72      = as.integer(Geno == 3)
    )

  cleaned <- combined %>%
    left_join(protein_matrix, by = "SampleID") %>%
    drop_na(any_of(c("CASE", "Fold", "X_Male", "X_Age",
                     setdiff(colnames(protein_matrix), "SampleID"))))

  cat(sprintf("  T=%.1f | N=%d | cases=%d | controls=%d | folds=%d\n",
              T_val, nrow(cleaned), sum(cleaned$CASE),
              sum(cleaned$CASE == 0), k))

  cleaned
}

# =============================================================================
# 3. Fast AUC & CV Utilities
# =============================================================================
fast_auc <- function(y, p) {
  r <- rank(p, ties.method = "average")
  n1 <- sum(y == 1L); n0 <- length(y) - n1
  if (n1 == 0L || n0 == 0L) return(NA_real_)
  (sum(r[y == 1L]) - n1 * (n1 + 1) / 2) / (n1 * n0)
}

cv_auc_with_design <- function(X, y, folds) {
  # folds is a list of list(train, test) — precomputed
  pred <- numeric(length(y))
  for (f in folds) {
    fit  <- glm.fit(x = X[f$train, , drop = FALSE],
                    y = y[f$train], family = binomial())
    beta <- fit$coefficients
    beta[is.na(beta)] <- 0
    pred[f$test] <- plogis(drop(X[f$test, , drop = FALSE] %*% beta))
  }
  fast_auc(y, pred)
}

# =============================================================================
# 4. Main Loop: Greedy Selection Per T
# =============================================================================
T_list         <- c(0.5, 1.0, 2.0, 3.0, 5.0)
covariate_list <- c("X_Male", "X_Age", "X_Geno_SOD1_A4V",
                    "X_Geno_SOD1_non_A4V", "X_Geno_C9orf72")
k <- 5

summary_df <- expand.grid(T = T_list, AUC = NA, N_protein = NA,
                           Protein_list = NA, stringsAsFactors = FALSE)

for (T_val in T_list) {
  cat(sprintf("\n=== T = %.1f ===\n", T_val))

  cleaned <- prepare_binary_data(visit_info, protein_matrix, T_val, k, fold_dir)

  # Precompute fold train/test indices once
  fold_vec <- cleaned$Fold
  n <- nrow(cleaned)
  all_idx <- seq_len(n)
  folds <- lapply(seq_len(k), function(i) {
    test <- which(fold_vec == i)
    list(train = setdiff(all_idx, test), test = test)
  })

  write.csv(cleaned,
            file.path(int_dir, sprintf("binary_data_T_%.1f.csv", T_val)),
            row.names = FALSE)

  # Build full numeric matrix once
  x_full <- as.matrix(cleaned[, c(protein_list, covariate_list), drop = FALSE])
  y      <- as.integer(cleaned$CASE)

  # Precompute covariate-only design (intercept + covariates)
  cov_design <- cbind(1, x_full[, covariate_list, drop = FALSE])
  colnames(cov_design)[1] <- "(Intercept)"

  # ── Step 1: Single-protein screening ──────────────────────────────────────
  cat("  Evaluating single proteins...\n")
  single_aucs <- numeric(length(protein_list))

  pb <- progress_bar$new(total = length(protein_list), clear = FALSE,
                          format = "  [:bar] :percent")
  for (i in seq_along(protein_list)) {
    X <- cbind(cov_design, x_full[, protein_list[i]])
    single_aucs[i] <- cv_auc_with_design(X, y, folds)
    pb$tick()
  }
  names(single_aucs) <- protein_list

  # Sort by AUC descending
  ord <- order(single_aucs, decreasing = TRUE)
  protein_auc <- tibble(Protein = protein_list[ord], AUC = single_aucs[ord])

  write.csv(protein_auc,
            file.path(out_dir, sprintf("Single_protein_T_%.1f.csv", T_val)),
            row.names = FALSE)

  # ── Step 2: Greedy forward selection ──────────────────────────────────────
  cat("  Greedy forward selection...\n")
  n_prot <- length(protein_list)
  greedy_proteins <- character(n_prot)
  greedy_aucs     <- numeric(n_prot)

  greedy_proteins[1] <- protein_auc$Protein[1]
  greedy_aucs[1]     <- protein_auc$AUC[1]
  remaining          <- setdiff(protein_list, greedy_proteins[1])

  # Base design: intercept + covariates + selected proteins so far
  base_design <- cbind(cov_design, x_full[, greedy_proteins[1], drop = FALSE])

  pb2 <- progress_bar$new(total = n_prot, clear = FALSE,
                           format = "  Greedy step :current/:total [:bar]")
  pb2$tick()

  for (step in 2:n_prot) {
    pb2$tick()
    if (length(remaining) == 0) break

    best_auc  <- -Inf
    best_prot <- remaining[1]

    for (prot in remaining) {
      X   <- cbind(base_design, x_full[, prot])
      auc <- cv_auc_with_design(X, y, folds)
      if (!is.na(auc) && auc > best_auc) {
        best_auc  <- auc
        best_prot <- prot
      }
    }

    greedy_proteins[step] <- best_prot
    greedy_aucs[step]     <- best_auc
    remaining <- setdiff(remaining, best_prot)

    # Extend base design with the newly selected protein
    base_design <- cbind(base_design, x_full[, best_prot, drop = FALSE])
  }

  greedy_results <- tibble(Protein = greedy_proteins, AUC = greedy_aucs)

  write.csv(greedy_results,
            file.path(out_dir, sprintf("Greedy_protein_T_%.1f.csv", T_val)),
            row.names = FALSE)

  best_n     <- which.max(greedy_results$AUC)
  best_panel <- greedy_results$Protein[1:best_n]
  summary_df[summary_df$T == T_val, "AUC"]          <- greedy_results$AUC[best_n]
  summary_df[summary_df$T == T_val, "N_protein"]    <- best_n
  summary_df[summary_df$T == T_val, "Protein_list"] <- paste(best_panel, collapse = ", ")

  cat(sprintf("  Best: %d proteins, AUROC=%.4f\n", best_n, greedy_results$AUC[best_n]))
}

write.csv(summary_df, file.path(out_dir, "summary.csv"), row.names = FALSE)
print(summary_df %>% select(T, AUC, N_protein))

cat("\n=== Script 4-1 complete ===\n")
cat("Outputs in:", out_dir, "\n")
