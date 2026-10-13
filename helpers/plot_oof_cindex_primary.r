# ============================================================================
# Dumbbell Plot for OOF C-Index Comparisons
# ============================================================================
# Creates dumbbell plots comparing SOC vs NMR models across different
# covariate sets (BASE, CLIN, MARKER, etc.)
#
# Usage:
#   # Your data format:
#   df <- data.frame(
#     model = c("base", "base", "clin", "clin", "marker", "marker"),
#     type = c("SOC", "NMR", "SOC", "NMR", "SOC", "NMR"),
#     c_index = c(0.52, 0.67, 0.62, 0.67, 0.65, 0.69),
#     ci_lower = c(0.52, 0.66, 0.61, 0.66, 0.64, 0.68),
#     ci_upper = c(0.53, 0.67, 0.63, 0.68, 0.66, 0.70)
#   )
#   
#   # Create plot
#   p <- plot_cindex_dumbbell(df, title = "KPRB Internal: T2D")
#   print(p)
# ============================================================================

library(ggplot2)
library(dplyr)
library(tidyr)

#' Plot C-index comparison as dumbbell plot
#'
#' @param df Data frame with columns: model, type, c_index, ci_lower, ci_upper
#' @param title Plot title
#' @param subtitle Optional subtitle
#' @param xlim Optional x-axis limits (e.g., c(0.5, 0.85))
#' @param model_labels Optional named vector to relabel models (e.g., c("base" = "BASE", "clin" = "CLIN"))
#' @param type_labels Optional named vector to relabel types (e.g., c("SOC" = "No NMR", "NMR" = "NMR"))
#' @param use_ggsci Use ggsci::scale_color_npg() for colors (requires ggsci package)
#' 
#' @return ggplot object
plot_cindex_dumbbell <- function(
  df,
  title = NULL,
  subtitle = NULL,
  xlim = NULL,
  model_labels = NULL,
  type_labels = NULL,
  use_ggsci = TRUE
) {
  
  # Validate input
  required_cols <- c("model", "type", "c_index", "ci_lower", "ci_upper")
  missing_cols <- setdiff(required_cols, colnames(df))
  
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }
  
  # Apply custom labels if provided
  df_plot <- df
  
  if (!is.null(model_labels)) {
    df_plot <- df_plot %>%
      mutate(model = recode(model, !!!model_labels))
  }
  
  if (!is.null(type_labels)) {
    df_plot <- df_plot %>%
      mutate(type = recode(type, !!!type_labels))
  }
  
  # Convert to factors with proper ordering
  unique_models <- unique(df_plot$model)
    df_plot <- df_plot %>%
    mutate(
        model = factor(model, levels = rev(unique_models)),
        type  = factor(type, levels = unique(type))
    )
  
  # Create wide format for dumbbell segments
  df_wide <- df_plot %>%
    select(model, type, c_index) %>%
    pivot_wider(names_from = type, values_from = c_index)
  
  # Get the two type names dynamically
  type_names <- colnames(df_wide)[colnames(df_wide) != "model"]
  
  if (length(type_names) != 2) {
    stop("Data must have exactly 2 types per model for dumbbell plot")
  }
  
  type1 <- type_names[1]
  type2 <- type_names[2]
  
  # Create plot
  p <- ggplot() +
    # Dumbbell connector lines
    geom_segment(
      data = df_wide,
      aes(y = model, yend = model, 
          x = .data[[type1]], xend = .data[[type2]]),
      color = "grey70",
      linewidth = 0.7,
      linetype = "dotted"
    ) +
    # Confidence intervals (thicker/longer CIs)
    geom_errorbarh(
      data = df_plot,
      aes(
        y = model,
        x = c_index,
        xmin = ci_lower,
        xmax = ci_upper,
        color = I("#000000")
   
      ),
      linewidth = 0.9,
      height = 0.22,
      alpha = 0.9,
      show.legend = FALSE
    ) + 
    # Points (increase dot size)
    geom_point(
      data = df_plot,
      aes(y = model, x = c_index, color = type),
      size = 1.25       # INCREASED: was 1
    ) +
    # Styling
    labs(
      x = "C-index (95% CI)",
      y = NULL,
      color = NULL,
      title = title,
      subtitle = subtitle
    ) +
    theme_bw() +
    theme(
      panel.grid.major.y = element_blank(),
      plot.title = element_text(hjust = 0.5, face = 'bold', size = 12),
      plot.subtitle = element_text(hjust = 0.5, size = 10),
      legend.position = "right"
    )
  
  # Add color scale
  if (use_ggsci && requireNamespace("ggsci", quietly = TRUE)) {
    p <- p + ggsci::scale_color_npg()
  } else {
    # Fallback to default colors
    p <- p + scale_color_manual(
      values = c("#E64B35", "#4DBBD5")  # Red and blue
    )
  }
  
  # Apply x-axis limits if specified
  if (!is.null(xlim)) {
    p <- p + coord_cartesian(xlim = xlim)
  }
  
  return(p)
}

#' Convenience wrapper for OOF results with standard naming
#'
#' Automatically applies BASE/CLIN/MARKER labels
#' 
#' @param df Data frame with columns: model, type, c_index, ci_lower, ci_upper
#' @param title Plot title
#' @param disease Disease name (e.g., "T2D", "CHD", "NAFLD", "CKD")
#' @param xlim Optional x-axis limits
#' @param marker_label Optional y-axis label for the `marker` model row.
#'   Defaults to `"MARKER"`. Use e.g. `"SoC"` or `"PANEL"` to override.
#' @param type_labels Optional named vector to relabel types in the legend
#'   (e.g., c("SOC" = "No NMR", "NMR" = "NMR")). Defaults to NULL (unchanged).
#' 
#' @return ggplot object
plot_oof_cindex <- function(
  df,
  disease = NULL,
  title = NULL,
  subtitle = NULL,
  xlim = NULL,
  marker_label = "MARKER",
  type_labels = NULL
) {
  
  # Generate title if not provided
  if (is.null(title) && !is.null(disease)) {
    title <- paste0("KPRB Internal (OOF): ", disease)
  }
  
  # Standard model labels
  model_labels <- c(
    "base" = "BASE",
    "clin" = "CLIN", 
    "marker" = marker_label
  )
  
  # Filter to only models that exist in the data
  existing_models <- intersect(names(model_labels), unique(df$model))
  model_labels_filtered <- model_labels[existing_models]
  
  # Filter to only types that exist in the data
  type_labels_filtered <- NULL
  if (!is.null(type_labels)) {
    existing_types <- intersect(names(type_labels), unique(df$type))
    type_labels_filtered <- type_labels[existing_types]
  }

  # Call main plotting function
  plot_cindex_dumbbell(
    df = df,
    title = title,
    subtitle = subtitle,
    xlim = xlim,
    model_labels = model_labels_filtered,
    type_labels = type_labels_filtered,
    use_ggsci = TRUE
  )
}

#' Create multi-panel plot for multiple diseases
#'
#' @param results_list Named list of data frames (one per disease)
#' @param xlim Optional x-axis limits (applied to all panels)
#' 
#' @return patchwork object with combined plots
plot_oof_multi_disease <- function(
  results_list,
  xlim = NULL
) {
  
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("Package 'patchwork' is required for multi-panel plots. Install with: install.packages('patchwork')")
  }
  
  # Create individual plots
  plots <- lapply(names(results_list), function(disease) {
    plot_oof_cindex(
      df = results_list[[disease]],
      disease = disease,
      xlim = xlim
    )
  })
  
  # Combine using patchwork
  combined <- patchwork::wrap_plots(plots, ncol = 2)
  
  return(combined)
}

#' Reshape a transport/OOF summary table into dumbbell-plot format
#'
#' Accepts a data frame with a `spec` column (e.g. base, clin, marker,
#' base_nmr, clin_nmr, marker_nmr) and returns columns:
#' model, type, c_index, ci_lower, ci_upper.
#'
#' @param summary_df  summary table with spec + C-index columns
#' @param spec_col    name of the spec column (default "spec")
#' @param endpoint_col optional endpoint column; if present, returns one
#'                     dumbbell-ready data frame per endpoint as a named list
#' @param endpoint_order optional character vector specifying facet order
#'                       (e.g. c("T2D", "CKD", "MASLD", ...)); defaults to the
#'                       order of first appearance in `summary_df`
#'
#' @return data.frame, or named list of data.frames when `endpoint_col` is set
prep_cindex_dumbbell_from_specs <- function(
  summary_df,
  spec_col = "spec",
  endpoint_col = NULL,
  endpoint_order = NULL
) {
  required <- c(spec_col, "c_index", "ci_lower", "ci_upper")
  missing_cols <- setdiff(required, colnames(summary_df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  reshape_one <- function(df) {
    df %>%
      mutate(
        model = sub("_nmr$", "", .data[[spec_col]]),
        type  = ifelse(grepl("_nmr$", .data[[spec_col]]), "NMR", "SOC")
      ) %>%
      select(model, type, c_index, ci_lower, ci_upper)
  }

  if (!is.null(endpoint_col)) {
    if (!endpoint_col %in% colnames(summary_df)) {
      stop("endpoint_col '", endpoint_col, "' not found in summary_df.")
    }

    endpoints <- unique(summary_df[[endpoint_col]])
    if (is.null(endpoint_order)) {
      endpoint_order <- endpoints
    } else {
      missing_eps <- setdiff(endpoint_order, endpoints)
      if (length(missing_eps) > 0) {
        stop("endpoint_order contains labels not in summary_df: ",
             paste(missing_eps, collapse = ", "))
      }
      # keep only endpoints present in the data, in the requested order
      endpoint_order <- intersect(endpoint_order, endpoints)
    }

    out <- split(summary_df, summary_df[[endpoint_col]])
    out <- lapply(out, reshape_one)
    out[endpoint_order]
  } else {
    reshape_one(summary_df)
  }
}

#' Build a named list of dumbbell-ready data frames from transport results
#'
#' @param transport_results named list returned by `transport_all_specs()` /
#'   `run_kprb_transport()` (each element must have a `$summary` slot)
#' @param endpoint_labels optional named vector mapping result names to plot
#'                        labels (e.g. c(t2d = "T2D", ckd = "CKD"))
#'
#' @return named list of dumbbell-ready data frames
prep_transport_results_list <- function(
  transport_results,
  endpoint_labels = NULL
) {
  if (is.null(names(transport_results))) {
    stop("transport_results must be a named list.")
  }

  out <- lapply(transport_results, function(res) {
    if (is.null(res$summary)) stop("Each transport result must contain `$summary`.")
    prep_cindex_dumbbell_from_specs(res$summary)
  })

  if (!is.null(endpoint_labels)) {
    nm <- intersect(names(endpoint_labels), names(out))
    names(out)[match(nm, names(out))] <- endpoint_labels[nm]
  }

  out
}

#' Faceted dumbbell plot for multiple diseases/endpoints
#'
#' Same visual style as `plot_cindex_dumbbell()`, but one panel per disease.
#'
#' @param results_list Named list of dumbbell-ready data frames (model, type,
#'                     c_index, ci_lower, ci_upper)
#' @param xlim Optional x-axis limits (applied to all panels)
#' @param ncol Number of columns in the facet grid
#' @param title Optional overall plot title
#' @param marker_label Optional y-axis label for the `marker` model row.
#'   Defaults to `"MARKER"`. Use e.g. `"SoC"` or `"PANEL"` to override.
#'
#' @return ggplot object
plot_oof_cindex_faceted <- function(
  results_list,
  xlim = c(0.5, 0.85),
  ncol = 3,
  title = NULL,
  marker_label = "MARKER"
) {
  model_labels <- c(
    "base" = "BASE",
    "clin" = "CLIN",
    "marker" = marker_label
  )

  df_plot <- bind_rows(results_list, .id = "disease") %>%
    mutate(
      disease = factor(disease, levels = names(results_list)),
      model = recode(model, !!!model_labels),
      model = factor(model, levels = rev(c("BASE", "CLIN", marker_label))),
      type = factor(type, levels = c("SOC", "NMR")),
      y_base = as.numeric(model)
    )

  df_wide <- df_plot %>%
    select(disease, model, y_base, type, c_index) %>%
    pivot_wider(names_from = type, values_from = c_index)

  p <- ggplot() +
    geom_segment(
      data = df_wide,
      aes(y = y_base, yend = y_base, x = SOC, xend = NMR),
      color = "grey85",
      linewidth = 0.35,
      linetype = "dotted"
    ) +
    geom_errorbar(
      data = df_plot,
      aes(y = y_base, xmin = ci_lower, xmax = ci_upper),
      color = "black",
      orientation = "y",
      width = 0.2,
      linewidth = 1.5,
      show.legend = FALSE
    ) +
    geom_point(
      data = df_plot,
      aes(y = y_base, x = c_index, color = type),
      size = 2
    ) +
    facet_wrap(~ disease, ncol = ncol) +
    scale_y_continuous(
      breaks = seq_along(levels(df_plot$model)),
      labels = levels(df_plot$model),
      limits = c(0.5, 3.5),
      expand = expansion(mult = 0)
    ) +
    scale_x_continuous(
      breaks = seq(xlim[1], xlim[2], by = 0.1)
    ) +
    coord_cartesian(xlim = xlim, clip = "off") +
    scale_color_manual(
      values = c("SOC" = "#E64B35", "NMR" = "#4DBBD5")
    ) +
    labs(
      x = "C-index (95% CI)",
      y = NULL,
      color = NULL,
      title = title
    ) +
    theme_bw(base_size = 10) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 12),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      legend.position = "bottom",
      legend.direction = "horizontal",
      panel.spacing = unit(0.8, "lines"),
      plot.margin = margin(4, 4, 4, 4),
      plot.title = element_text(hjust = 0.5, face = "bold", size = 12)
    )

  p
}

