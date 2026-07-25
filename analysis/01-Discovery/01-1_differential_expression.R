# =============================================================================
# Script 01-1: Differentially Expressed Proteins
# ALS vs. Healthy Controls — Mixed-Effects Models
#
# Inputs:  data/discovery/protein_matrix.csv
#          data/discovery/df_visit_info.csv
# Outputs: results/01_differential_expression/
#   Mix_effect_model_lmer.csv
#   figures/volcano_plot.pdf
# =============================================================================

# Set working directory to the script's folder
# Set working directory to analysis/discovery/
# In RStudio: setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
# From command line: Rscript should be run from this directory
if (requireNamespace("rstudioapi", quietly=TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}
# If not in RStudio, set manually: setwd("/path/to/analysis/discovery")

# ── Libraries ─────────────────────────────────────────────────────────────────
library(here)
library(tidyverse)
library(lme4)
library(lmerTest)
library(progress)
library(ggrepel)

set.seed(2024)

# ── Output directories ────────────────────────────────────────────────────────
out_dir  <- here::here("analysis/01-Discovery/results/01_differential_expression")
fig_dir  <- file.path(out_dir, "figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ── Group / genotype decode helper ───────────────────────────────────────────
decode_visit_info <- function(vi) {
  vi %>%
    mutate(Sex = ifelse(Male == 1, "Male", "Female")) %>%
    filter(OlinkGrp != 1)   # exclude Pre-symptomatic (OlinkGrp == 1)
}

# =============================================================================
# 1. Load Data
# =============================================================================
cat("Loading protein matrix...\n")
protein_matrix <- read.csv(here::here("data/discovery/protein_matrix.csv"), check.names = FALSE)

cat("Loading visit info...\n")
visit_info <- read.csv(here::here("data/discovery/df_visit_info.csv")) %>% decode_visit_info()

protein_list <- colnames(protein_matrix)[-1]

df_joined <- protein_matrix %>% left_join(visit_info, by = "SampleID")

cat(sprintf("Proteins: %d\n", length(protein_list)))
cat(sprintf("Samples in joined data: %d\n", nrow(df_joined)))

# =============================================================================
# 2. Mixed-Effects Model Loop (all proteins)
# =============================================================================
df_model <- df_joined %>% filter(OlinkGrp %in% c(0, 3))

cat("\nFitting mixed-effects models for all proteins...\n")
n_proteins <- length(protein_list)

pb <- progress_bar$new(
  format = "Processing [:bar] :percent (:current/:total) :elapsed",
  total = n_proteins, width = 60
)

results_protein <- data.frame(
  Protein   = protein_list,
  Estimate  = NA_real_,
  Std_Error = NA_real_,
  df        = NA_real_,
  t_value   = NA_real_,
  Pr_t      = NA_real_
)

for (i in seq_len(n_proteins)) {
  protein <- protein_list[i]
  df_p <- df_model %>%
    select(SampleID, NPX = all_of(protein), OlinkGrp, GeneGrp4, Male, CollAge, UIDx) %>%
    mutate(
      X_Genotype = fct_relevel(as.factor(GeneGrp4), "9"),
      X_ALS      = ifelse(OlinkGrp == 3, 1, 0),
      X_Male     = ifelse(Male == 1, 1, 0),
      X_Age      = CollAge,
      X_UID      = as.factor(UIDx)
    )

  tryCatch({
    m <- lmerTest::lmer(
      NPX ~ X_ALS + X_Male + X_Genotype + X_Age + (1 | X_UID),
      data    = df_p,
      control = lme4::lmerControl(optimizer = "Nelder_Mead",
                                  check.conv.singular = "ignore")
    )
    coefs <- summary(m)$coefficients
    results_protein[i, 2:6] <- coefs["X_ALS", 1:5]
  }, error = function(e) {
    message(sprintf("  Failed for %s: %s", protein, e$message))
  })

  pb$tick()
}

results_protein <- results_protein %>%
  mutate(
    padj        = p.adjust(Pr_t, method = "BH"),
    significant = ifelse(padj < 0.05, "Significant", "Not-Significant"),
    diffexpressed = ifelse(Estimate > 0 & padj < 0.05, "UP",
                           ifelse(Estimate < 0 & padj < 0.05, "DOWN", "NO"))
  )

write.csv(results_protein,
          file.path(out_dir, "Mix_effect_model_lmer.csv"),
          row.names = FALSE)
cat(sprintf("\nSaved lmer results: %d proteins, %d significant\n",
            nrow(results_protein), sum(results_protein$significant == "Significant")))

# =============================================================================
# 3. Volcano Plot
# =============================================================================
mix_effect_protein <- results_protein %>%
  filter(significant == "Significant") %>% pull(Protein)

results_protein_plot <- results_protein %>% arrange(desc(padj))
num_up   <- sum(results_protein_plot$diffexpressed == "UP")
num_down <- sum(results_protein_plot$diffexpressed == "DOWN")
max_fc   <- max(abs(results_protein_plot$Estimate), na.rm = TRUE)
max_p    <- max(-log2(results_protein_plot$padj),   na.rm = TRUE)

top_10_up <- results_protein_plot %>%
  filter(Estimate > 0) %>% arrange(padj) %>% slice_head(n = 10) %>%
  mutate(label_color = "red")
top_10_down <- results_protein_plot %>%
  filter(Estimate < 0) %>% arrange(padj) %>% slice_head(n = 10) %>%
  mutate(label_color = "blue")
top_10_genes <- bind_rows(top_10_up, top_10_down)

p_volcano <- ggplot(results_protein_plot,
                    aes(Estimate, -log2(padj), color = diffexpressed)) +
  geom_point(show.legend = FALSE) +
  scale_color_manual(values = c("DOWN" = "blue", "NO" = "grey", "UP" = "red")) +
  scale_x_continuous(limits = c(-max_fc, max_fc), name = "Effect Size (Beta1)") +
  scale_y_continuous(limits = c(0, max_p), name = "-log2(Adjusted P-value)") +
  geom_hline(yintercept = -log2(0.05), linetype = "dashed", color = "black") +
  geom_label_repel(data = top_10_genes,
                   aes(label = Protein, color = diffexpressed),
                   size = 4, max.overlaps = 50, show.legend = FALSE) +
  annotate("text", x =  max_fc * 0.6, y = 0.85 * max_p,
           label = paste("UP:", num_up),   color = "red",  size = 6) +
  annotate("text", x = -max_fc * 0.6, y = 0.85 * max_p,
           label = paste("DOWN:", num_down), color = "blue", size = 6) +
  theme_classic() +
  ggtitle("Differential Protein Expression: ALS vs. Control")

ggsave(file.path(fig_dir, "volcano_plot.pdf"), p_volcano, width = 8, height = 5)
cat("Saved volcano plot.\n")

cat("\n=== Script 01-1 complete ===\n")
cat("Outputs in:", out_dir, "\n")
