#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

input_dir <- if (length(args) >= 1) args[1] else "."
head_lines <- if (length(args) >= 2) as.integer(args[2]) else 12L
heads_output_name <- if (length(args) >= 3) args[3] else "all_raw_heads_1_12.txt"
summary_output_name <- if (length(args) >= 4) args[4] else "mutation_summary_table.tsv"

if (is.na(head_lines) || head_lines <= 0) {
  stop("head_lines must be a positive integer.")
}

input_dir <- normalizePath(input_dir, winslash = "/", mustWork = TRUE)
raw_files <- sort(list.files(input_dir, pattern = "\\.raw$", full.names = TRUE))

if (length(raw_files) == 0) {
  stop(sprintf("No .raw files found in: %s", input_dir))
}

heads_output <- file.path(input_dir, heads_output_name)
summary_output <- file.path(input_dir, summary_output_name)

extract_gene <- function(file_name) {
  m <- regexec("^UKB_(.+?)_pathogenic_corrected\\.raw$", file_name)
  parts <- regmatches(file_name, m)[[1]]
  if (length(parts) >= 2) parts[2] else tools::file_path_sans_ext(file_name)
}

# 1) Write combined heads
con_out <- file(heads_output, open = "wt", encoding = "UTF-8")
on.exit(close(con_out), add = TRUE)

for (raw_file in raw_files) {
  name <- basename(raw_file)
  writeLines(sprintf("@%s (1-%d)", name, head_lines), con_out)
  head_txt <- readLines(raw_file, n = head_lines, warn = FALSE, encoding = "UTF-8")
  if (length(head_txt) > 0) writeLines(head_txt, con_out)
  writeLines("", con_out)
}

# 2) Build mutation summary table
summary_list <- vector("list", length(raw_files))
row_count <- 0L

for (i in seq_along(raw_files)) {
  raw_file <- raw_files[i]
  gene <- extract_gene(basename(raw_file))

  # Read header line only
  header <- readLines(raw_file, n = 1, warn = FALSE, encoding = "UTF-8")
  if (length(header) == 0) {
    summary_list[[i]] <- data.frame(
      gene = character(0),
      mutation = character(0),
      carried_count = integer(0),
      stringsAsFactors = FALSE
    )
    next
  }

  cols <- strsplit(header, "\t", fixed = TRUE)[[1]]
  if (length(cols) <= 6) {
    summary_list[[i]] <- data.frame(
      gene = character(0),
      mutation = character(0),
      carried_count = integer(0),
      stringsAsFactors = FALSE
    )
    next
  }

  mutation_cols <- cols[7:length(cols)]

  df <- utils::read.table(
    raw_file,
    sep = "\t",
    header = TRUE,
    quote = "",
    comment.char = "",
    fill = TRUE,
    check.names = FALSE,
    stringsAsFactors = FALSE
  )

  # genotype starts from 7th column in PLINK .raw
  geno <- df[, 7:ncol(df), drop = FALSE]

  carried_counts <- vapply(
    geno,
    FUN = function(x) {
      suppressWarnings(x_num <- as.numeric(x))
      sum(!is.na(x_num) & x_num > 0)
    },
    FUN.VALUE = integer(1)
  )

  this_df <- data.frame(
    gene = rep(gene, length(mutation_cols)),
    mutation = mutation_cols,
    carried_count = as.integer(carried_counts),
    stringsAsFactors = FALSE
  )

  summary_list[[i]] <- this_df
  row_count <- row_count + nrow(this_df)
}

summary_df <- do.call(rbind, summary_list)
utils::write.table(
  summary_df,
  file = summary_output,
  sep = "\t",
  row.names = FALSE,
  quote = FALSE
)

cat(sprintf("Processed raw files: %d\n", length(raw_files)))
cat(sprintf("Heads output: %s\n", heads_output))
cat(sprintf("Summary output: %s\n", summary_output))
cat(sprintf("Summary rows: %d\n", row_count))
