# =============================================================================
# Script 01-2: Pathway Enrichment Analysis (GO)
# Up- and down-regulated proteins from mixed-effects model
#
# If enrichment CSVs already exist, skips the API call and plots directly.
#
# Inputs:  data/intermediate/Mix_effect_model_lmer.csv
# Outputs: data/intermediate/
#   Mix_effect_model_lmer_up_enrichment.csv
#   Mix_effect_model_lmer_down_enrichment.csv
#   figures/up_enrichment_plot.pdf
#   figures/down_enrichment_plot.pdf
# =============================================================================

# Set working directory to analysis/discovery/
if (requireNamespace("rstudioapi", quietly=TRUE) && rstudioapi::isAvailable()) {
  setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
}

# ── Libraries ─────────────────────────────────────────────────────────────────
library(tidyverse)

# ── Output directories ────────────────────────────────────────────────────────
int_dir  <- "data/intermediate"
out_dir  <- "results/01_differential_expression"
fig_dir  <- file.path(out_dir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ── File paths ────────────────────────────────────────────────────────────────
input_file       <- file.path(out_dir, "Mix_effect_model_lmer.csv")
up_enrich_file   <- file.path(out_dir, "Mix_effect_model_lmer_up_enrichment.csv")
down_enrich_file <- file.path(out_dir, "Mix_effect_model_lmer_down_enrichment.csv")

# =============================================================================
# 1. Load differential expression results
# =============================================================================
cat("Loading differential expression results from:", input_file, "\n")
results_protein <- read.csv(input_file)

up_proteins   <- results_protein %>% filter(diffexpressed == "UP")   %>% pull(Protein)
down_proteins <- results_protein %>% filter(diffexpressed == "DOWN") %>% pull(Protein)

# =============================================================================
# 2. Up-regulated enrichment
# =============================================================================
if (file.exists(up_enrich_file)) {
  cat("Up-regulated enrichment CSV found, skipping API call.\n")
  result_up <- read.csv(up_enrich_file)
} else {
  cat("Running enrichment analysis for up-regulated proteins...\n")
  library(gprofiler2)
  if (length(up_proteins) > 0) {
    tryCatch({
      up.enrichment <- gost(query = up_proteins,
                            organism = "hsapiens", correction_method = "fdr",
                            evcodes = TRUE)
      result_up <- subset(up.enrichment$result, select = -c(parents))
      write.csv(result_up, up_enrich_file, row.names = FALSE)
      write.csv(result_up, file.path(out_dir, "Mix_effect_model_lmer_up_enrichment.csv"), row.names = FALSE)
    }, error = function(e) {
      warning("Up-regulated enrichment failed (network?): ", conditionMessage(e))
      result_up <- NULL
    })
  }
}

if (exists("result_up") && !is.null(result_up) && nrow(result_up) > 0) {
  write.csv(result_up, file.path(out_dir, "Mix_effect_model_lmer_up_enrichment.csv"), row.names = FALSE)
  top_30_up <- result_up %>%
    select(source, term_name, p_value) %>%
    mutate(term = sprintf("%s - %s", source, term_name)) %>%
    arrange(p_value) %>% slice_head(n = 30) %>% arrange(desc(p_value))

  p_up <- top_30_up %>%
    mutate(term = factor(term, levels = unique(term)),
           highlight = ifelse(grepl("muscle", term, ignore.case = TRUE),
                              "highlight", "normal")) %>%
    ggplot(aes(term, -log2(p_value), fill = highlight)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    scale_fill_manual(values = c("highlight" = "#e31a1c", "normal" = "#fc9272")) +
    ylab("-log2(FDR)") + xlab("") + ggtitle("Top 30 Up-regulated Pathways") +
    coord_flip() + theme_classic(base_size = 12) +
    geom_hline(yintercept = -log2(0.05), linetype = "dashed")
  ggsave(file.path(fig_dir, "up_enrichment_plot.pdf"), p_up, width = 10, height = 7)
  cat("Saved up_enrichment_plot.pdf\n")
}

# =============================================================================
# 3. Down-regulated enrichment
# =============================================================================
if (file.exists(down_enrich_file)) {
  cat("Down-regulated enrichment CSV found, skipping API call.\n")
  result_down <- read.csv(down_enrich_file)
} else {
  cat("Running enrichment analysis for down-regulated proteins...\n")
  if (!requireNamespace("gprofiler2", quietly = TRUE)) library(gprofiler2)
  if (length(down_proteins) > 0) {
    tryCatch({
      down.enrichment <- gost(query = down_proteins,
                              organism = "hsapiens", correction_method = "fdr",
                              evcodes = TRUE)
      result_down <- subset(down.enrichment$result, select = -c(parents))
      write.csv(result_down, down_enrich_file, row.names = FALSE)
      write.csv(result_down, file.path(out_dir, "Mix_effect_model_lmer_down_enrichment.csv"), row.names = FALSE)
    }, error = function(e) {
      warning("Down-regulated enrichment failed (network?): ", conditionMessage(e))
      result_down <- NULL
    })
  }
}

if (exists("result_down") && !is.null(result_down) && nrow(result_down) > 0) {
  write.csv(result_down, file.path(out_dir, "Mix_effect_model_lmer_down_enrichment.csv"), row.names = FALSE)
  top_30_down <- result_down %>%
    select(source, term_name, p_value) %>%
    mutate(term = sprintf("%s - %s", source, term_name)) %>%
    arrange(p_value) %>% slice_head(n = 30) %>% arrange(desc(p_value))

  p_down <- top_30_down %>%
    mutate(term = factor(term, levels = unique(term)),
           highlight = ifelse(grepl("muscle", term, ignore.case = TRUE),
                              "highlight", "normal")) %>%
    ggplot(aes(term, -log2(p_value), fill = highlight)) +
    geom_bar(stat = "identity", show.legend = FALSE) +
    scale_fill_manual(values = c("highlight" = "#2171b5", "normal" = "#9ecae1")) +
    ylab("-log2(FDR)") + xlab("") + ggtitle("Top 30 Down-regulated Pathways") +
    coord_flip() + theme_classic(base_size = 12) +
    geom_hline(yintercept = -log2(0.05), linetype = "dashed")
  ggsave(file.path(fig_dir, "down_enrichment_plot.pdf"), p_down, width = 10, height = 7)
  cat("Saved down_enrichment_plot.pdf\n")
}

cat("\n=== Script 01-2 complete ===\n")
cat("Outputs in:", out_dir, "\n")
