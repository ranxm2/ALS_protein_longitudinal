# =============================================================================
# Script 2-2: Protein Trajectory Analysis Using Delta Matrix
# Fit GAMM: delta ~ s(TimetoEvent) + (1|individual)
#
# Requires: data/intermediate/ (run 01-1 and 02-1 first)
# Inputs:   data/intermediate/delta_matrix_wide.csv
#           results/01_differential_expression/Mix_effect_model_lmer.csv
#           data/df_visit_info.csv
# Outputs:
#   results/02-Protein_trajectories/  (summary: predictions, summary CSVs, PDFs)
#   data/intermediate/                (individual-level data)
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

set.seed(42)

# ── Output directories ────────────────────────────────────────────────────────
out_dir <- "results/02-Protein_trajectories"
int_dir <- "data/intermediate"
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(int_dir, recursive = TRUE, showWarnings = FALSE)

# ── Decode visit info ─────────────────────────────────────────────────────────
decode_visit_info <- function(vi) {
  vi %>%
    mutate(Sex = ifelse(Male == 1, "Male", "Female")) %>%
    filter(Group %in% c(2, 3))   # only Convert and Affected
}

# =============================================================================
# 1. Load Data
# =============================================================================
cat("Loading data...\n")

delta_matrix_wide <- read.csv(
  "data/intermediate/delta_matrix_wide.csv"
)

visit_info <- read.csv("data/df_visit_info.csv") %>%
  decode_visit_info() %>%
  mutate(
    Group = factor(Group, levels = c(2, 3)),
    Geno  = factor(Geno,  levels = c(9, 1, 2, 3, 4)),
    Male  = factor(Male,  levels = c(0, 1)),
    Sex   = factor(Sex,   levels = c("Female", "Male"))
  )

protein.de.miami <- read.csv(
  "results/01_differential_expression/Mix_effect_model_lmer.csv"
)
protein.sign.miami <- protein.de.miami %>%
  filter(significant == "Significant") %>%
  arrange(padj) %>% pull(Protein)

cat(sprintf("Delta matrix: %d rows x %d cols\n",
            nrow(delta_matrix_wide), ncol(delta_matrix_wide)))

# =============================================================================
# 2. Create Analysis Subsets
# =============================================================================
visit_convert_affected <- visit_info %>%
  filter(Group %in% c(2, 3)) %>%
  mutate(X_UID = as.factor(EID))

visit_convert_all <- visit_info %>%
  filter(Group == 2) %>%
  mutate(X_UID = as.factor(EID))

visit_convert_preonset <- visit_info %>%
  filter(Group == 2, TimetoEvent < 0) %>%
  mutate(X_UID = as.factor(EID))

cat(sprintf("Subset 1 (Convert+Affected): %d | Subset 2 (Convert): %d | Subset 3 (Pre-onset): %d\n",
            nrow(visit_convert_affected), nrow(visit_convert_all),
            nrow(visit_convert_preonset)))

# =============================================================================
# 3. Helper Functions
# =============================================================================
get_gam_object <- function(mod) {
  if (is.null(mod))           return(NULL)
  if (inherits(mod, "gam"))   return(mod)
  if (is.list(mod) && "gam" %in% names(mod)) return(mod$gam)
  NULL
}

get_smooth_pvalue <- function(gam_obj) {
  if (is.null(gam_obj)) return(NA)
  tryCatch({
    st <- summary(gam_obj)$s.table
    if (!is.null(st) && "s(TimetoEvent)" %in% rownames(st))
      return(st["s(TimetoEvent)", "p-value"])
    NA
  }, error = function(e) NA)
}

get_significance_start_time <- function(predictions_df) {
  sig_points <- predictions_df %>% filter(sig != 0)
  if (nrow(sig_points) == 0)
    return(list(start_time_increase = NA, start_time_decrease = NA,
                start_time_any = NA, nearest_interval = NA,
                nearest_interval_start = NA))

  r      <- rle(predictions_df$sig)
  ends   <- cumsum(r$lengths)
  starts <- c(1, head(ends, -1) + 1)
  nzi    <- which(r$values != 0)

  intervals <- tibble(
    interval_start = predictions_df$TimetoEvent[starts[nzi]],
    interval_end   = predictions_df$TimetoEvent[ends[nzi]],
    run_sign       = r$values[nzi]
  )

  spanning <- intervals %>%
    filter(interval_start < 0 & interval_end > 0) %>%
    arrange(abs(interval_start))
  pre_only <- intervals %>%
    filter(interval_end <= 0) %>%
    arrange(desc(interval_end))

  nearest <- if (nrow(spanning) > 0) spanning %>% slice(1) else pre_only %>% slice(1)
  nearest_label <- if (nrow(nearest) > 0)
    sprintf("[%.1f, %.1f]", nearest$interval_start, nearest$interval_end)
  else NA_character_
  nearest_start <- if (nrow(nearest) > 0) nearest$interval_start else NA_real_

  inc <- predictions_df %>% filter(sig ==  1)
  dec <- predictions_df %>% filter(sig == -1)
  list(
    start_time_increase    = if (nrow(inc) > 0) min(inc$TimetoEvent) else NA,
    start_time_decrease    = if (nrow(dec) > 0) min(dec$TimetoEvent) else NA,
    start_time_any         = min(sig_points$TimetoEvent),
    nearest_interval       = nearest_label,
    nearest_interval_start = nearest_start
  )
}

# ── Trajectory fit ────────────────────────────────────────────────────────────
fit_delta_trajectory <- function(delta_df, protein_name,
                                  subset_name, min_obs = 10) {
  if (nrow(delta_df) < min_obs) return(NULL)
  df <- delta_df %>% filter(!is.na(delta)) %>%
    mutate(X_UID = as.factor(X_UID))
  if (nrow(df) < min_obs) return(NULL)

  formula_delta <- delta ~ s(TimetoEvent, k = 10)
  mod_delta <- NULL
  tryCatch({
    mod_delta <- mgcv::gamm(formula_delta, random = list(X_UID = ~1),
                             data = df, method = "REML")
  }, error = function(e) message("  Mixed GAM failed: ", e$message))

  if (is.null(mod_delta)) {
    tryCatch({
      mod_delta <- list(gam = mgcv::gam(formula_delta, data = df, method = "REML"))
    }, error = function(e) message("  Fixed GAM also failed"))
  }
  if (is.null(mod_delta)) return(NULL)

  gam_obj <- get_gam_object(mod_delta)
  pval    <- get_smooth_pvalue(gam_obj)

  yrs <- seq(min(df$TimetoEvent, na.rm = TRUE),
             max(df$TimetoEvent, na.rm = TRUE), length.out = 300)
  gp  <- data.frame(TimetoEvent = yrs)
  pr  <- predict(gam_obj, newdata = gp, se.fit = TRUE)
  gp$fit <- as.numeric(pr$fit); gp$se <- as.numeric(pr$se.fit)
  gp$lower <- gp$fit - 1.96 * gp$se
  gp$upper <- gp$fit + 1.96 * gp$se
  gp <- gp %>% mutate(
    sig = case_when(lower > 0 ~ 1L, upper < 0 ~ -1L, TRUE ~ 0L),
    dir = case_when(sig > 0 ~ "increase", sig < 0 ~ "decrease",
                    TRUE ~ NA_character_)
  )

  st <- get_significance_start_time(gp)
  list(protein = protein_name, model = mod_delta, gam = gam_obj,
       data = df, predictions = gp, formula = formula_delta,
       n_obs = nrow(df), n_subjects = length(unique(df$X_UID)),
       pvalue_TimetoEvent = pval,
       start_time_increase = st$start_time_increase,
       start_time_decrease = st$start_time_decrease,
       start_time_any = st$start_time_any,
       nearest_interval = st$nearest_interval,
       nearest_interval_start = st$nearest_interval_start)
}

# ── Batch analysis ────────────────────────────────────────────────────────────
run_batch_analysis <- function(delta_wide, visit_design, subset_name, outdir, datadir) {
  protein_list_sub <- setdiff(colnames(delta_wide), "SampleID")
  cat(sprintf("\n=== Running analysis: %s | Proteins: %d ===\n",
              subset_name, length(protein_list_sub)))

  results_list <- list()
  grid_all <- data.frame(); data_all <- data.frame()

  for (protein in protein_list_sub) {
    delta_df <- delta_wide %>%
      select(SampleID, all_of(protein)) %>%
      rename(delta = !!sym(protein)) %>%
      inner_join(visit_design, by = "SampleID")

    result <- fit_delta_trajectory(delta_df, protein, subset_name)
    if (!is.null(result)) {
      results_list[[protein]] <- result
      grid_all <- bind_rows(grid_all,
                            result$predictions %>% mutate(protein = protein))
      data_all <- bind_rows(data_all,
                            result$data %>%
                              select(SampleID, X_UID, TimetoEvent, delta) %>%
                              mutate(protein = protein))
    }
  }

  subdir <- file.path(outdir, gsub(" ", "_", subset_name))
  dir.create(subdir, recursive = TRUE, showWarnings = FALSE)

  # Summary-level (trajectories, 95% CrI) → results
  write.csv(grid_all,
            file.path(subdir, paste0("predictions_",
                                     gsub(" ", "_", subset_name), ".csv")),
            row.names = FALSE)

  # Individual-level data → data/intermediate
  write.csv(data_all,
            file.path(datadir, paste0("trajectory_data_",
                                      gsub(" ", "_", subset_name), ".csv")),
            row.names = FALSE)

  summary_df <- tibble(
    protein                = names(results_list),
    n_obs                  = sapply(results_list, `[[`, "n_obs"),
    n_subjects             = sapply(results_list, `[[`, "n_subjects"),
    pvalue_TimetoEvent     = sapply(results_list, `[[`, "pvalue_TimetoEvent"),
    start_time_increase    = sapply(results_list, `[[`, "start_time_increase"),
    start_time_decrease    = sapply(results_list, `[[`, "start_time_decrease"),
    start_time_any         = sapply(results_list, `[[`, "start_time_any"),
    nearest_interval       = sapply(results_list, `[[`, "nearest_interval"),
    nearest_interval_start = sapply(results_list, `[[`, "nearest_interval_start")
  ) %>%
    mutate(fdr_TimetoEvent     = p.adjust(pvalue_TimetoEvent, method = "BH"),
           significant_nominal = pvalue_TimetoEvent < 0.05,
           significant_fdr     = fdr_TimetoEvent    < 0.05) %>%
    arrange(pvalue_TimetoEvent)

  write.csv(summary_df,
            file.path(subdir, paste0("summary_",
                                     gsub(" ", "_", subset_name), ".csv")),
            row.names = FALSE)

  list(results = results_list, predictions = grid_all,
       data = data_all, summary = summary_df, subset_name = subset_name)
}

# =============================================================================
# 4. Run Three Subsets
# =============================================================================
analysis_1 <- run_batch_analysis(delta_matrix_wide, visit_convert_affected,
                                  "Convert_Affected_AllTime", out_dir, int_dir)
analysis_2 <- run_batch_analysis(delta_matrix_wide, visit_convert_all,
                                  "Convert_Only_AllTime", out_dir, int_dir)
analysis_3 <- run_batch_analysis(delta_matrix_wide, visit_convert_preonset,
                                  "Convert_Only_PreOnset", out_dir, int_dir)

# =============================================================================
# 5. Cross-Subset Comparison
# =============================================================================
all_summaries <- bind_rows(
  analysis_1$summary %>% mutate(subset = "Convert+Affected_AllTime"),
  analysis_2$summary %>% mutate(subset = "Convert_AllTime"),
  analysis_3$summary %>% mutate(subset = "Convert_PreOnset")
)

all_summaries <- all_summaries %>%
  mutate(nearest_interval_start = ifelse(
    !is.na(fdr_TimetoEvent) & fdr_TimetoEvent < 0.05,
    nearest_interval_start, NA_real_
  ))

write.csv(all_summaries,
          file.path(out_dir, "combined_summary_all_subsets_long.csv"),
          row.names = FALSE)

summary_wide <- all_summaries %>%
  select(protein, subset, n_obs, pvalue_TimetoEvent, fdr_TimetoEvent,
         start_time_any, nearest_interval_start) %>%
  pivot_wider(names_from  = subset,
              values_from = c(n_obs, pvalue_TimetoEvent, fdr_TimetoEvent,
                              start_time_any, nearest_interval_start),
              names_glue  = "{subset}_{.value}")

write.csv(summary_wide,
          file.path(out_dir, "combined_summary_all_subsets.csv"),
          row.names = FALSE)

# Report significant counts
for (a in list(analysis_1, analysis_2, analysis_3)) {
  cat(sprintf("\n%s: nominal=%d, FDR=%d\n",
              a$subset_name,
              sum(a$summary$significant_nominal, na.rm = TRUE),
              sum(a$summary$significant_fdr,     na.rm = TRUE)))
}

cat("\n=== Script 2-2 complete ===\n")
cat("Outputs in:", out_dir, "\n")
