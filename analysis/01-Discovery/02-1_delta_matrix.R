# =============================================================================
# Script 2-1: Delta Protein Matrix Calculation
# Fit GAMM on Controls, predict baselines, compute delta = observed - expected
#
# Requires: data/intermediate/ (run 01-1_differential_expression.R first)
# Inputs:   data/intermediate/final_data_matrix.rds
#           results/01_differential_expression/Mix_effect_model_lmer.csv
#           data/df_visit_info.csv
# Outputs (all in data/intermediate/):
#   delta_matrix_long.csv
#   delta_matrix_wide.csv
# =============================================================================

# Set working directory to analysis/discovery/
# In RStudio: setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# From command line: Rscript should be run from this directory
if (requireNamespace("rstudioapi", quietly=TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# ── Libraries ─────────────────────────────────────────────────────────────────
library(tidyverse)
library(mgcv)
library(progress)
set.seed(2025)

# ── Output directories ────────────────────────────────────────────────────────
int_dir <- "data/intermediate"
dir.create(int_dir, recursive = TRUE, showWarnings = FALSE)

# ── Group / genotype decode helper ───────────────────────────────────────────
decode_visit_info <- function(vi) {
  vi %>%
    mutate(Sex = ifelse(Male == 1, "Male", "Female"))
}

# =============================================================================
# 1. Load Data
# =============================================================================
cat("Loading data...\n")

protein.mat.miami <- readRDS(
  "data/intermediate/final_data_matrix.rds"
)

protein.de.miami <- read.csv(
  "results/01_differential_expression/Mix_effect_model_lmer.csv"
)

visit.info.miami <- read.csv("data/df_visit_info.csv") %>%
  decode_visit_info() %>%
  mutate(
    Group = factor(Group, levels = c(0, 1, 2, 3)),
    Geno  = factor(Geno,  levels = c(1, 2, 3, 4, 9)),
    Male  = factor(Male,  levels = c(0, 1)),
    Sex   = factor(Sex,   levels = c("Female", "Male"))
  )

protein.all.miami  <- colnames(protein.mat.miami)[-1]
protein.sign.miami <- protein.de.miami %>%
  filter(significant == "Significant") %>% pull(Protein)

cat(sprintf("Total proteins: %d | Significant: %d\n",
            length(protein.all.miami), length(protein.sign.miami)))

# =============================================================================
# 2. Define Analysis Groups
# =============================================================================
visit.info.miami.ctrl <- visit.info.miami %>% filter(Group == 0)

cat(sprintf("Control samples (for GAMM fit): %d | Total samples (for delta): %d\n",
            nrow(visit.info.miami.ctrl), nrow(visit.info.miami)))

# =============================================================================
# 3. Main Loop: Delta Calculation
# =============================================================================
protein_list <- protein.sign.miami

delta_df_long <- data.frame(
  SampleID     = character(),
  Delta_Value  = numeric(),
  Protein_Name = character(),
  EID          = character(),
  Age          = numeric(),
  TimetoEvent  = numeric(),
  Sex          = character(),
  Geno         = integer(),
  stringsAsFactors = FALSE
)

pb <- progress_bar$new(
  format = "  Processing [:bar] :percent | :current/:total | ETA: :eta",
  total  = length(protein_list), clear = FALSE
)

for (protein in protein_list) {
  pb$tick()

  # Control data for fitting baseline
  ctrl_df <- protein.mat.miami %>%
    inner_join(visit.info.miami.ctrl, by = "SampleID") %>%
    select(SampleID, EID, Age, Sex, Geno,
           protein_val = !!sym(protein)) %>%
    filter(!is.na(protein_val)) %>%
    mutate(EID  = as.factor(EID),
           Age  = as.numeric(Age),
           Sex  = as.factor(Sex),
           Geno = as.factor(Geno))

  # Confirmed data for delta prediction
  confirm_df <- protein.mat.miami %>%
    inner_join(visit.info.miami, by = "SampleID") %>%
    select(SampleID, EID, Age, TimetoEvent, Sex, Geno,
           protein_val = !!sym(protein)) %>%
    filter(!is.na(protein_val)) %>%
    mutate(EID         = as.factor(EID),
           Age         = as.numeric(Age),
           TimetoEvent = as.numeric(TimetoEvent),
           Sex         = factor(Sex,  levels = levels(ctrl_df$Sex)),
           Geno        = factor(Geno, levels = levels(ctrl_df$Geno)))

  if (nrow(ctrl_df) < 20 || nrow(confirm_df) < 5) {
    warning(paste("Skipping", protein, "- insufficient data"))
    next
  }

  # Fit GAMM on controls
  mod_ctrl <- tryCatch(
    mgcv::gamm(
      formula = protein_val ~ Sex + s(Age, k = 10),
      random  = list(EID = ~1),
      data    = ctrl_df,
      method  = "REML"
    ),
    error = function(e) NULL
  )

  if (is.null(mod_ctrl)) next

  gam_ctrl <- mod_ctrl$gam

  # Predict expected values and compute delta
  confirm_df$pred_expected <- predict(gam_ctrl, newdata = confirm_df,
                                      type = "response")
  confirm_df$delta <- confirm_df$protein_val - confirm_df$pred_expected

  delta_df_long <- rbind(delta_df_long, data.frame(
    SampleID     = confirm_df$SampleID,
    Delta_Value  = confirm_df$delta,
    Protein_Name = protein,
    EID          = as.character(confirm_df$EID),
    Age          = confirm_df$Age,
    TimetoEvent  = confirm_df$TimetoEvent,
    Sex          = as.character(confirm_df$Sex),
    Geno         = as.integer(as.character(confirm_df$Geno)),
    stringsAsFactors = FALSE
  ))
}

# =============================================================================
# 4. Save Delta Matrix (all samples)
# =============================================================================

cat(sprintf("\nDelta matrix: %d observations, %d proteins, %d subjects\n",
            nrow(delta_df_long),
            length(unique(delta_df_long$Protein_Name)),
            length(unique(delta_df_long$EID))))

write.csv(delta_df_long,
          file.path(int_dir, "delta_matrix_long.csv"), row.names = FALSE)

delta_df_wide <- delta_df_long %>%
  select(SampleID, Protein_Name, Delta_Value) %>%
  pivot_wider(names_from = Protein_Name, values_from = Delta_Value)

write.csv(delta_df_wide,
          file.path(int_dir, "delta_matrix_wide.csv"), row.names = FALSE)

cat("\n=== Script 2-1 complete ===\n")
cat("Outputs in:", int_dir, "\n")
