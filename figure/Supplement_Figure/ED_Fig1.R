# ==============================================================================
# ED_Fig1.R — Extended Data Figure 1: Patient Visit Timeline
# ==============================================================================
# Author : Ximing Ran
# Date   : 2026-02-22
#
# Description:
#   Generates Extended Data Figure 1 — spaghetti line plots of individual
#   patient visit schedules across the four clinical groups in the Miami
#   familial ALS cohort (Control, Pre-symptomatic, Convert, Affected).
#   Each horizontal line represents one subject; points mark individual
#   sample collection visits. Subjects are ordered by follow-up duration
#   within each group.
#
# Input data:
#   data/fig_data/individual_data/Supplement_Figure/
#     df_visit_info(N=516).csv   — one row per visit with columns:
#                                  UIDx, Group, FUdur_yr, YrSinceOs, etc.
#
# Output (./results/):
#   ED_Fig1A_All_group_visit_followup_time.{pdf,png}
#   (skipped entirely if individual data is absent)
#
# Prerequisite:
#   Run analysis/01-Differentially_Expressed_Proteins/ Rmd notebooks first
#   (visit info is exported during data preparation).
# ==============================================================================

# Load Required Libraries ======================================================
library(tidyverse)
library(cowplot)
library(dplyr)
library(here)
library(svglite)

# Define Functions =============================================================

#' Create Line Plot for Patient Visits
#'
#' @param data Data frame containing visit information
#' @param patient_label Character string for plot title (optional)
#' @param group_label Character string for grouping variable (e.g., "GenoGroup", "Sex")
#' @param x_val Character string specifying the x-axis variable name
#' @param x_label Character string for x-axis label
#' @param last_visit Logical or character; if TRUE, marks last visit; 
#'                   "FIRST_AFTER_ONSET" for first visit after onset
#' @param first_visit Logical; if TRUE, marks first visit with blue triangle
#' @param Onset Logical; if TRUE, orders patients by onset time
#' @param last_year Numeric; marks visits within last N years
#' @param after_onset Logical or "FIRST"; marks visits after onset
#' @param censor Logical; if TRUE, uses censoring time reference
#'
#' @return A ggplot object
#'
#' @examples
#' visit_line_plot(data = df_visit_info, x_val = "FUdur_yr", 
#'                 x_label = "Follow-up Duration (years)")
visit_line_plot <- function(data, 
                            patient_label = NULL,
                            group_label = NULL,
                            x_val = NULL,
                            x_label = NULL,
                            last_visit = NULL,
                            first_visit = NULL,
                            Onset = NULL,
                            last_year = NULL,
                            after_onset = NULL,
                            censor = NULL) {
  
  # Set default group label
  if (is.null(group_label)) {
    group_label <- "Group"
  }
  
  # Data Processing ============================================================
  
  # Determine ordering based on onset/censor settings
  if (!is.null(Onset) || !is.null(censor)) {
    # Order by minimum X value (earliest time point)
    order_df <- data %>%
      mutate(X_value = .data[[x_val]]) %>%
      select(UIDx, X_value, all_of(group_label)) %>%
      group_by(UIDx, .data[[group_label]]) %>%
      summarise(order_value = min(X_value, na.rm = TRUE), .groups = "drop") %>%
      arrange(.data[[group_label]], order_value) %>%
      mutate(order_index = row_number()) %>%
      select(-all_of(group_label))
  } else {
    # Order by maximum X value (latest time point)
    order_df <- data %>%
      mutate(X_value = .data[[x_val]]) %>%
      select(UIDx, X_value, all_of(group_label)) %>%
      group_by(UIDx, .data[[group_label]]) %>%
      summarise(order_value = max(X_value, na.rm = TRUE), .groups = "drop") %>%
      arrange(.data[[group_label]], desc(order_value)) %>%
      mutate(order_index = row_number()) %>%
      select(-all_of(group_label))
  }
  
  # Prepare main data frame with ordering
  df <- data %>%
    mutate(X_value = .data[[x_val]]) %>%
    select(UIDx, CollNum, X_value, all_of(group_label)) %>%
    left_join(order_df, by = "UIDx") %>%
    arrange(order_index) %>%
    mutate(UIDx = factor(UIDx, levels = rev(unique(.$UIDx))))
  
  # Define color mapping based on group label
  if (group_label == "GenoGroup") {
    colors_map <- c(
      "SOD1 A4V" = "red",
      "SOD1 nonA4V" = "blue",
      "C9orf72" = "orange",
      "Other Genotype" = "green",
      "Unknown" = "darkgreen"
    )
  } else if (group_label == "Sex") {
    colors_map <- c(
      "Male" = "blue",
      "Female" = "red"
    )
  } else {
    colors_map <- "black"
  }
  
  # Handle censoring time adjustments
  if (!is.null(censor)) {
    data_last_visit <- data %>%
      group_by(UIDx) %>%
      filter(CollNum == max(CollNum)) %>%
      ungroup()
    x_label <- "Follow-up Duration (years)"
    data$X_value <- data$YrSinceCen
    data_last_visit$X_value <- data_last_visit$YrSinceCen
  } else {
    data_last_visit <- NULL
    x_label <- ifelse(
      group_label %in% c("Control", "Pre-symptomatic"),
      "Follow-up Duration (years)",
      "Follow-up Duration (years)"
    )
  }
  
  # Base Plot Construction =====================================================
  
  p <- ggplot(df, aes(x = X_value, y = UIDx, group = UIDx, 
                      color = .data[[group_label]]))
  
  # Add vertical reference line (except for age plots)
  if (x_val != "CollAge") {
    p <- p + geom_vline(xintercept = 0, linetype = "dashed", color = "grey")
  }
  
  # Add main plot elements
  p <- p +
    geom_line(alpha = 1, size = 1.2) +
    geom_point(alpha = 1, size = 1.5) +
    scale_color_manual(values = colors_map) +
    labs(
      title = patient_label,
      x = x_label,
      y = "Participant",
      color = group_label
    ) +
    theme_classic(base_size = 14) +
    theme(
      strip.text = element_text(size = 14, face = "bold"),
      axis.text.y = element_blank(),
      plot.title = element_text(hjust = 0.5, size = 14),
      panel.grid.minor = element_blank(),
      legend.position = if (!is.null(group_label) && group_label == "Group") {
        "none"
      } else {
        "right"
      }
    )
  
  # Add Special Markers ========================================================
  
  # Mark last visits
  if (!is.null(last_visit)) {
    if (isTRUE(last_visit) && isTRUE(Onset)) {
      # Last visit before onset for Convert group
      data_last_visit <- data %>%
        filter(YrSinceOs < 0, Group == "Convert") %>%
        group_by(UIDx) %>%
        filter(YrSinceOs == max(YrSinceOs, na.rm = TRUE)) %>%
        ungroup() %>%
        mutate(X_value = .data[[x_val]]) %>%
        select(UIDx, X_value)
    } else if (last_visit == "FIRST_AFTER_ONSET") {
      # First visit after onset
      data_last_visit <- data %>%
        filter(YrSinceOs >= 0) %>%
        group_by(UIDx) %>%
        filter(CollNum == min(CollNum)) %>%
        ungroup() %>%
        mutate(X_value = .data[[x_val]]) %>%
        select(UIDx, X_value)
    } else {
      # Overall last visit
      data_last_visit <- data %>%
        group_by(UIDx) %>%
        filter(CollNum == max(CollNum)) %>%
        ungroup()
    }
    
    p <- p + geom_point(
      data = data_last_visit,
      aes(x = X_value, y = UIDx),
      shape = 2, color = "red", size = 3
    )
  }
  
  # Mark first visits
  if (!is.null(first_visit)) {
    data_first_visit <- data %>%
      group_by(UIDx) %>%
      filter(CollNum == min(CollNum)) %>%
      ungroup() %>%
      mutate(X_value = .data[[x_val]]) %>%
      select(UIDx, X_value)
    
    p <- p + geom_point(
      data = data_first_visit,
      aes(x = X_value, y = UIDx),
      shape = 2, color = "blue", size = 3
    )
  }
  
  # Add censoring indicators
  if (!is.null(censor) && x_val == "YrSinceCen") {
    data_last_visit <- data %>%
      group_by(UIDx) %>%
      filter(CollNum == max(CollNum)) %>%
      ungroup()
    
    p <- p +
      geom_segment(
        data = data_last_visit,
        aes(x = 0, y = UIDx, xend = YrSinceCen, yend = UIDx),
        linetype = "dashed", color = "red", size = 0.5, alpha = 0.5
      ) +
      geom_point(
        data = data_last_visit,
        aes(x = 0, y = UIDx),
        shape = 21, fill = "white", color = "red", size = 2
      )
    
    if (!is.null(last_year)) {
      p <- p + geom_vline(
        xintercept = -last_year,
        linetype = "dashed",
        color = "grey"
      )
    }
  }
  
  # Mark visits in last N years
  if (!is.null(last_year)) {
    data_last_yrs <- data %>%
      filter(YrSinceCen > -last_year) %>%
      ungroup()
    
    p <- p + geom_point(
      data = data_last_yrs,
      aes(x = X_value, y = UIDx),
      shape = 4, color = "black", size = 3
    )
  }
  
  # Mark visits after onset
  if (!is.null(after_onset)) {
    if (isTRUE(after_onset)) {
      data_after_onset <- data %>%
        filter(YrSinceOs > 0) %>%
        mutate(X_value = .data[[x_val]]) %>%
        ungroup()
    } else if (after_onset == "FIRST") {
      data_after_onset <- data %>%
        filter(YrSinceOs >= 0) %>%
        group_by(UIDx) %>%
        filter(row_number() > 1) %>%
        mutate(X_value = .data[[x_val]]) %>%
        ungroup()
    }
    
    p <- p + geom_point(
      data = data_after_onset,
      aes(x = X_value, y = UIDx),
      shape = 4, color = "black", size = 3
    )
  }
  
  return(p)
}

# Main Analysis ================================================================

# Load Data
sup_fig1_path <- here("data", "fig_data", "individual_data", "Supplement_Figure",
                       "df_visit_info(N=516).csv")
if (!file.exists(sup_fig1_path)) {
  cat("  WARNING: Individual data not found for ED_Fig1. Skipping.\n")
} else {
df_visit_info <- read.csv(sup_fig1_path)

# Create Individual Plots ------------------------------------------------------

# Control group
p1 <- visit_line_plot(
  data = filter(df_visit_info, Group == "Control"),
  x_val = "FUdur_yr",
  x_label = "Follow-up Duration (years)"
) +
  theme(aspect.ratio = 1)

# Pre-symptomatic group (with spacing)
p_temp <- visit_line_plot(
  data = filter(df_visit_info, Group == "Pre-symptomatic"),
  x_val = "FUdur_yr",
  x_label = "Follow-up Duration (years)"
)

# Add empty rows for visual spacing
uids <- levels(p_temp$data$UIDx)
n_empty <- 10
empty_levels <- paste0(" ", seq_len(n_empty))
new_levels <- c(empty_levels, uids)
p_temp$data$UIDx <- factor(p_temp$data$UIDx, levels = new_levels)

p2 <- p_temp +
  scale_y_discrete(drop = FALSE) +
  theme(aspect.ratio = 1)

# Convert group
p3 <- visit_line_plot(
  data = filter(df_visit_info, Group == "Convert"),
  x_val = "FUdur_yr",
  x_label = "Follow-up Duration (years)"
) +
  theme(aspect.ratio = 1)

# Affected group
p4 <- visit_line_plot(
  data = filter(df_visit_info, Group == "Affected"),
  x_val = "FUdur_yr",
  x_label = "Follow-up Duration (years)"
) +
  theme(aspect.ratio = 1)

# Combine and Save Plots -------------------------------------------------------

combination_plot <- plot_grid(p1, p2, p3, p4, nrow = 1)

# Ensure output directory exists
dir.create("results", showWarnings = FALSE, recursive = TRUE)

# Save plot
ggsave(
  filename = file.path("results", 
                       "ED_Fig1A_All_group_visit_followup_time.pdf"),
  plot = combination_plot,
  width = 28,
  height = 7
)
# Save plot
ggsave(
  filename = file.path("results",
                       "ED_Fig1A_All_group_visit_followup_time.png"),
  plot = combination_plot,
  width = 28,
  height = 7,
  dpi = 300
)
ggsave(
  filename = file.path("results",
                       "ED_Fig1A_All_group_visit_followup_time.svg"),
  plot = combination_plot,
  device = svglite::svglite,
  width = 28,
  height = 7
)

} # end ED_Fig1 individual data guard

# ==============================================================================
# End of Script
# ==============================================================================