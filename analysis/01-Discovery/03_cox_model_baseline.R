# =============================================================================
# Script 3: Time-to-Event Analysis — Cox Proportional Hazards Model
# Univariate Cox + LASSO Cox for ALS onset prediction in people at risk
#
# Requires: data/intermediate/ (run 01-1 and 02-1 first)
# Inputs:   data/intermediate/delta_matrix_wide.csv
#           results/01_differential_expression/Mix_effect_model_lmer.csv
#           data/df_visit_info.csv
# Outputs:
#   results/03_cox_model_baseline/  (summary: Cox results, LASSO results, figures)
#   data/intermediate/              (individual-level survival data)
# =============================================================================

# Set working directory to analysis/discovery/
# In RStudio: setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# From command line: Rscript should be run from this directory
if (requireNamespace("rstudioapi", quietly=TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# ── Libraries ─────────────────────────────────────────────────────────────────
library(tidyverse)
library(survival)
library(glmnet)
library(ggplot2)
library(ggrepel)

set.seed(2025)

# ── Output directories ────────────────────────────────────────────────────────
out_dir  <- "results/03_cox_model_baseline"
fig_dir  <- file.path(out_dir, "figures")
int_dir  <- "data/intermediate"
dir.create(out_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir,  recursive = TRUE, showWarnings = FALSE)
dir.create(int_dir,  recursive = TRUE, showWarnings = FALSE)

# ── Significance stars ────────────────────────────────────────────────────────
map_significance_stars <- function(p) {
  ifelse(p < 0.001, "***",
  ifelse(p < 0.01,  "**",
  ifelse(p < 0.05,  "*",
  ifelse(p < 0.1,   ".", ""))))
}

# ── Decode visit info ─────────────────────────────────────────────────────────
decode_visit_info <- function(vi) {
  vi %>%
    mutate(Sex = ifelse(Male == 1, "Male", "Female")) %>%
    filter(Group %in% c(1, 2))   # People at risk (onset and not onset)
}

# =============================================================================
# 1. Load Data
# =============================================================================
cat("Loading data...\n")

protein.mat.miami <- read.csv(
  "data/intermediate/delta_matrix_wide.csv"
)

protein.de.miami <- read.csv(
  "results/01_differential_expression/Mix_effect_model_lmer.csv"
)

protein.sign.miami <- protein.de.miami %>%
  arrange(padj) %>%
  filter(significant == "Significant") %>%
  pull(Protein)

cat(sprintf("Significant proteins: %d\n", length(protein.sign.miami)))

# Load visit info
visit.info.miami <- read.csv("data/df_visit_info.csv") %>%
  decode_visit_info() %>%
  mutate(
    Group = factor(Group, levels = c(1, 2)),
    Geno  = factor(Geno,  levels = c(1, 2, 3, 4, 9)),
    Male  = factor(Male,  levels = c(0, 1)),
    Sex   = factor(Sex,   levels = c("Female", "Male"))
  )

# =============================================================================
# 2. Prepare Survival Data
# =============================================================================
cat("Preparing survival data...\n")

# People at risk with onset (Group 2): event=1, time from baseline to onset
vis_onset <- visit.info.miami %>%
  filter(Group == 2) %>%
  arrange(EID, VisitIndex) %>%
  group_by(EID) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    event   = 1L,
    time_yr = -as.numeric(TimetoEvent)
  )

# People at risk without onset (Group 1): event=0, time from baseline to censoring
vis_no_onset <- visit.info.miami %>%
  filter(Group == 1) %>%
  arrange(EID, VisitIndex) %>%
  group_by(EID) %>%
  slice(1) %>%
  ungroup() %>%
  mutate(
    event   = 0L,
    time_yr = -as.numeric(TimetoCen)
  )

vis_sub <- bind_rows(vis_onset, vis_no_onset) %>%
  mutate(Geno = factor(Geno, levels = c(4, 1, 2, 3))) %>%
  filter(is.finite(time_yr), !is.na(time_yr), time_yr > 0)

prot_df <- protein.mat.miami %>% as.data.frame()
dat0    <- vis_sub %>% left_join(prot_df, by = "SampleID")

cat(sprintf("Survival data: N=%d, events=%d\n",
            nrow(dat0), sum(dat0$event)))

covar <- c("Age", "Sex", "Geno")

# =============================================================================
# 3. Univariate Cox Regression
# =============================================================================
cat("\nFitting univariate Cox models...\n")

results_list <- vector("list", length(protein.sign.miami))

for (i in seq_along(protein.sign.miami)) {
  protein_name <- protein.sign.miami[i]

  tryCatch({
    dat_p <- dat0 %>%
      select(EID, SampleID, time_yr, event,
             all_of(protein_name), all_of(covar)) %>%
      drop_na()

    if (nrow(dat_p) < 10 || sum(dat_p$event) < 5) {
      results_list[[i]] <- data.frame(
        Protein = protein_name, N = nrow(dat_p), N_events = sum(dat_p$event),
        coef = NA, exp_coef = NA, se_coef = NA, z = NA, p_value = NA,
        Status = "Insufficient data"
      )
    } else {
      formula_str <- paste0("Surv(time_yr, event) ~ ", protein_name,
                            " + ", paste(covar, collapse = " + "))
      fit <- coxph(as.formula(formula_str), data = dat_p)
      s   <- summary(fit)$coefficients
      pc  <- s[1, , drop = FALSE]

      results_list[[i]] <- data.frame(
        Protein  = protein_name,
        N        = nrow(dat_p),
        N_events = sum(dat_p$event),
        coef     = pc[1, "coef"],
        exp_coef = pc[1, "exp(coef)"],
        se_coef  = pc[1, "se(coef)"],
        z        = pc[1, "z"],
        p_value  = pc[1, "Pr(>|z|)"],
        Status   = "Success"
      )
    }
  }, error = function(e) {
    results_list[[i]] <<- data.frame(
      Protein = protein_name, N = NA, N_events = NA,
      coef = NA, exp_coef = NA, se_coef = NA, z = NA, p_value = NA,
      Status = paste("Error:", e$message)
    )
  })
}

# NEFL score: linear predictor from univariate Cox model
nefl_fit <- coxph(Surv(time_yr, event) ~ NEFL + Age + Sex + Geno, data = dat0)
dat0$NEFL_score <- predict(nefl_fit, newdata = dat0, type = "lp")

cox_results <- bind_rows(results_list) %>%
  mutate(fdr = ifelse(!is.na(p_value), p.adjust(p_value, method = "fdr"), NA)) %>%
  arrange(p_value)

write.csv(cox_results,
          file.path(out_dir, "cox_results_univariate.csv"),
          row.names = FALSE)
cat(sprintf("Saved univariate Cox results: %d proteins\n", nrow(cox_results)))

# =============================================================================
# 4. Volcano Plot
# =============================================================================
volcano_df <- cox_results %>%
  filter(!is.na(p_value), !is.na(coef)) %>%
  mutate(
    negLogP      = -log10(p_value),
    significance = case_when(
      fdr < 0.05 & coef > 0 ~ "UP",
      fdr < 0.05 & coef < 0 ~ "DOWN",
      TRUE                   ~ "NO"
    )
  )

bh_threshold_val <- {
  m  <- length(volcano_df$p_value)
  ps <- sort(volcano_df$p_value)
  ki <- which(ps <= seq_len(m) / m * 0.05)
  if (length(ki) == 0) 0.05 / m else ps[max(ki)]
}

n_up   <- sum(volcano_df$significance == "UP")
n_down <- sum(volcano_df$significance == "DOWN")

top_up   <- volcano_df %>% filter(significance == "UP")   %>% arrange(p_value) %>% slice_head(n = 15)
top_down <- volcano_df %>% filter(significance == "DOWN")  %>% arrange(p_value) %>% slice_head(n = 10)
top_genes <- bind_rows(top_up, top_down)

y_max_data <- max(volcano_df$negLogP, na.rm = TRUE)

p_volcano <- ggplot(volcano_df, aes(x = coef, y = negLogP, color = significance)) +
  geom_point(alpha = 0.6, size = 2) +
  geom_hline(yintercept = -log10(bh_threshold_val),
             linetype = "dashed", color = "gray40") +
  annotate("text", x = min(volcano_df$coef, na.rm=TRUE) + 0.1,
           y = y_max_data * 1.1, label = paste0("Protective: ", n_down),
           hjust = 0, color = "#427ebf", size = 5, fontface = "bold") +
  annotate("text", x = max(volcano_df$coef, na.rm=TRUE) - 0.1,
           y = y_max_data * 1.1, label = paste0("Risk-Associated: ", n_up),
           hjust = 1, color = "#A70C20", size = 5, fontface = "bold") +
  scale_color_manual(values = c(UP = "#A70C20", DOWN = "#427ebf", NO = "grey60")) +
  scale_y_continuous(limits = c(0, y_max_data * 1.3), expand = c(0, 0)) +
  labs(x = "Cox Coefficient (log Hazard Ratio)",
       y = expression(-log[10]("P-value")),
       title = "Proteomic Biomarkers of ALS Conversion Risk") +
  theme_classic(base_size = 14) +
  geom_label_repel(data = top_genes, aes(label = Protein, color = significance),
                   size = 3.5, max.overlaps = 50, show.legend = FALSE,
                   fontface = "bold")

ggsave(file.path(fig_dir, "volcano_cox.pdf"), p_volcano, width = 10, height = 7)
cat("Saved Cox volcano plot.\n")

# =============================================================================
# 5. LASSO Cox Regression
# =============================================================================
cat("\nFitting LASSO Cox model...\n")

dat_lasso <- dat0 %>%
  select(EID, SampleID, time_yr, event,
         all_of(protein.sign.miami), all_of(covar)) %>%
  drop_na()

cat(sprintf("LASSO data: N=%d, events=%d, proteins=%d\n",
            nrow(dat_lasso), sum(dat_lasso$event),
            length(protein.sign.miami)))

# Save individual-level LASSO input data
write.csv(dat_lasso,
          file.path(int_dir, "cox_lasso_data.csv"),
          row.names = FALSE)

# Build feature matrices
X_proteins   <- as.matrix(dat_lasso[, protein.sign.miami])
X_sex        <- as.numeric(dat_lasso$Sex == "Male")
X_geno       <- model.matrix(~ Geno - 1, data = dat_lasso)
colnames(X_geno) <- paste0("Geno_", gsub("^Geno", "", colnames(X_geno)))
X_covariates <- cbind(Age = dat_lasso$Age, Sex_Male = X_sex, X_geno)

# Remove reference level (Geno == 4, Other Genotype) if present
ref_col <- grep("Geno_4$", colnames(X_covariates))
if (length(ref_col) > 0) X_covariates <- X_covariates[, -ref_col, drop = FALSE]

X_full          <- cbind(X_proteins, X_covariates)
penalty_factors <- c(rep(1, ncol(X_proteins)), rep(0, ncol(X_covariates)))
y_surv          <- Surv(dat_lasso$time_yr, dat_lasso$event)

set.seed(2025)
cat("Running LOOCV LASSO (may take several minutes)...\n")
cv_fit <- cv.glmnet(
  x              = X_full,
  y              = y_surv,
  family         = "cox",
  penalty.factor = penalty_factors,
  nfolds         = nrow(X_full),  # LOOCV
  type.measure   = "deviance",
  alpha          = 1
)

lambda_min <- cv_fit$lambda.min
lambda_1se <- cv_fit$lambda.1se

cat(sprintf("Lambda.min = %.4e | Lambda.1se = %.4e\n", lambda_min, lambda_1se))

# Extract selected proteins at lambda.min
fit_min   <- glmnet(X_full, y_surv, family = "cox",
                    penalty.factor = penalty_factors,
                    lambda = lambda_min, alpha = 1)
coef_min  <- coef(fit_min, s = lambda_min)
nz_min    <- which(coef_min != 0)
sel_names <- rownames(coef_min)[nz_min]
sel_coefs <- as.vector(coef_min[nz_min])

# Filter to proteins only
prot_mask <- sel_names %in% protein.sign.miami
lasso_min_df <- data.frame(
  Protein   = sel_names[prot_mask],
  LASSO_Coef = sel_coefs[prot_mask],
  LASSO_HR   = exp(sel_coefs[prot_mask])
) %>% arrange(desc(abs(LASSO_Coef)))

write.csv(lasso_min_df,
          file.path(out_dir, "lasso_selected_proteins_min.csv"),
          row.names = FALSE)
cat(sprintf("LASSO selected %d proteins at lambda.min\n", nrow(lasso_min_df)))

# Extract at lambda.1se
fit_1se   <- glmnet(X_full, y_surv, family = "cox",
                    penalty.factor = penalty_factors,
                    lambda = lambda_1se, alpha = 1)
coef_1se  <- coef(fit_1se, s = lambda_1se)
nz_1se    <- which(coef_1se != 0)
sel_1se   <- rownames(coef_1se)[nz_1se]
sel_c_1se <- as.vector(coef_1se[nz_1se])
prot_1se  <- sel_1se %in% protein.sign.miami

lasso_1se_df <- data.frame(
  Protein    = sel_1se[prot_1se],
  LASSO_Coef = sel_c_1se[prot_1se],
  LASSO_HR   = exp(sel_c_1se[prot_1se])
) %>% arrange(desc(abs(LASSO_Coef)))

write.csv(lasso_1se_df,
          file.path(out_dir, "lasso_selected_proteins_1se.csv"),
          row.names = FALSE)
cat(sprintf("LASSO selected %d proteins at lambda.1se\n", nrow(lasso_1se_df)))

# Save CV curve data
cv_curve_df <- data.frame(
  log_lambda = log(cv_fit$lambda),
  lambda     = cv_fit$lambda,
  cvm        = cv_fit$cvm,
  cvsd       = cv_fit$cvsd,
  nzero      = cv_fit$nzero
)
write.csv(cv_curve_df,
          file.path(out_dir, "lasso_cv_curve.csv"),
          row.names = FALSE)

# CV curve plot
pdf(file.path(fig_dir, "lasso_cv_curve.pdf"), width = 10, height = 5)
par(mfrow = c(1, 2))
plot(cv_fit, main = "LASSO Cox: Cross-Validation Curve")
abline(v = log(lambda_min), col = "red",  lty = 2, lwd = 2)
abline(v = log(lambda_1se), col = "blue", lty = 2, lwd = 2)
legend("topleft", legend = c("Lambda.min", "Lambda.1se"),
       col = c("red", "blue"), lty = 2, lwd = 2, bty = "n")
plot(cv_fit$glmnet.fit, xvar = "lambda", main = "LASSO Coefficient Paths")
abline(v = log(lambda_min), col = "red",  lty = 2, lwd = 2)
abline(v = log(lambda_1se), col = "blue", lty = 2, lwd = 2)
par(mfrow = c(1, 1))
dev.off()
cat("Saved LASSO CV curve.\n")

# =============================================================================
# 6. Panel Cox Models and Scores
# =============================================================================
# LASSO score = linear predictor from LASSO Cox at lambda.min
dat0$LASSO_score <- as.numeric(predict(fit_min, newx = X_full, s = lambda_min, type = "link"))

# Panel definitions
Panel_15 <- c("NEFL", "DUSP29", "CALCA", "NEB", "NGRN",
              "NOS1", "DTNB", "MEGF10", "APOA4", "ART3",
              "MYL3", "EPHA1", "EDA2R", "TTN", "ACTN2")

Panel_19 <- c("NEFL", "CALCA", "DUSP29", "EDA2R", "ACTN2",
              "TTN", "DTNB", "NOS1", "NGRN", "EPHA1",
              "NEB", "MYL3", "MEGF10", "APOA4", "ART3",
              "MYL11", "MYH1", "TNNC1", "SYNM")

# Panel-15 Cox model
panel15_formula <- as.formula(paste0("Surv(time_yr, event) ~ ",
                                      paste(Panel_15, collapse = " + "),
                                      " + ", paste(covar, collapse = " + ")))
panel15_fit <- coxph(panel15_formula, data = dat0)
dat0$Panel_15_score <- predict(panel15_fit, newdata = dat0, type = "lp")
cat(sprintf("Panel-15 Cox model fitted (%d proteins)\n", length(Panel_15)))

# Panel-19 Cox model
panel19_formula <- as.formula(paste0("Surv(time_yr, event) ~ ",
                                      paste(Panel_19, collapse = " + "),
                                      " + ", paste(covar, collapse = " + ")))
panel19_fit <- coxph(panel19_formula, data = dat0)
dat0$Panel_19_score <- predict(panel19_fit, newdata = dat0, type = "lp")
cat(sprintf("Panel-19 Cox model fitted (%d proteins)\n", length(Panel_19)))

# Save individual-level survival data with all scores
write.csv(dat0,
          file.path(int_dir, "cox_survival_data.csv"),
          row.names = FALSE)

# Also update cox_lasso_data with scores
dat_lasso$NEFL_score     <- dat0$NEFL_score[match(dat_lasso$SampleID, dat0$SampleID)]
dat_lasso$LASSO_score    <- dat0$LASSO_score[match(dat_lasso$SampleID, dat0$SampleID)]
dat_lasso$Panel_15_score <- dat0$Panel_15_score[match(dat_lasso$SampleID, dat0$SampleID)]
dat_lasso$Panel_19_score <- dat0$Panel_19_score[match(dat_lasso$SampleID, dat0$SampleID)]
write.csv(dat_lasso,
          file.path(int_dir, "cox_lasso_data.csv"),
          row.names = FALSE)

cat(sprintf("Saved survival data with 4 scores: N=%d\n", nrow(dat0)))

# =============================================================================
# 7. Kaplan-Meier Curves (median-split by each score)
# =============================================================================
library(survminer)

score_info <- list(
  list(name = "NEFL",      col = "NEFL_score"),
  list(name = "LASSO",     col = "LASSO_score"),
  list(name = "Panel-15",  col = "Panel_15_score"),
  list(name = "Panel-19",  col = "Panel_19_score")
)

for (si in score_info) {
  dat0$risk_group <- ifelse(dat0[[si$col]] >= median(dat0[[si$col]], na.rm = TRUE),
                            "High", "Low")
  dat0$risk_group <- factor(dat0$risk_group, levels = c("Low", "High"))

  km_fit <- survfit(Surv(time_yr, event) ~ risk_group, data = dat0)
  logrank <- survdiff(Surv(time_yr, event) ~ risk_group, data = dat0)
  pval_exact <- 1 - pchisq(logrank$chisq, df = 1)

  p <- ggsurvplot(
    km_fit, data = dat0,
    pval = FALSE, conf.int = TRUE,
    risk.table = TRUE, risk.table.col = "strata",
    palette = c("#427ebf", "#A70C20"),
    title = paste0("Kaplan-Meier: ", si$name, " Score (Median Split)"),
    subtitle = sprintf("Log-rank p = %.4e", pval_exact),
    xlab = "Time (years)", ylab = "Survival Probability",
    legend.labs = c("Low Risk", "High Risk"),
    ggtheme = theme_bw()
  )

  pdf(file.path(out_dir, paste0("KM_", gsub("-", "_", si$name), ".pdf")),
      width = 8, height = 7, onefile = FALSE)
  print(p)
  dev.off()
  cat(sprintf("Saved KM curve: %s\n", si$name))
}

cat("\n=== Script 3 complete ===\n")
cat("Outputs in:", out_dir, "\n")
