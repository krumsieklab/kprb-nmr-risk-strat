# ============================================================================
# plot_within_stratum_cindex_dumbbell()
# ============================================================================
# Horizontal dumbbell plot: marker_SOC vs marker_SOC+NMR per stratum, with
# model-specific global C-index reference lines. Panels are stacked vertically
# by strat_name.
#
# Input: long-format output from subgroup_cindex(), e.g. t2d_cindex_df
#   Required columns: strat_name, stratum, model, c_index
#   Optional: ci_lower, ci_upper, n, n_events, stratum
#
# Usage:
#   source(here::here("helpers", "plot_within_stratum_cindex_dumbbell.r"))
#   p <- plot_within_stratum_cindex_dumbbell(
#     cidx       = t2d_cindex_df,
#     soc_model  = "marker_SOC",
#     nmr_model  = "marker_SOC+NMR",
#     soc_global = 0.7406,
#     nmr_global = 0.7915,
#     title      = "Within-stratum discrimination",
#     x_by       = 0.10          # tick every 0.1 (0.60, 0.70, 0.80, ...)
#     # x_breaks = c(0.60, 0.70, 0.80)  # or pass explicit breaks
#     # model_labels = c("SOC", "SOC+NMR")
#     # legend_position = "none"        # hide legend
#   )
#   print(p)
# ============================================================================

if (!exists("%>%", mode = "function")) {
  `%>%` <- dplyr::`%>%`
}

plot_within_stratum_cindex_dumbbell <- function(
  cidx,
  soc_model = "marker_SOC",
  nmr_model = "marker_SOC+NMR",
  soc_global,
  nmr_global,
  title = "Within-stratum discrimination",
  subtitle = NULL,
  xlim = c(0.60, 0.88),
  x_by = 0.10,
  x_breaks = NULL,
  sort_by_soc = TRUE,
  show_ci = TRUE,
  show_global_labels = FALSE,
  col_soc = "#E64B35",
  col_nmr = "#4DBBD5",
  model_labels = c("SOC", "SOC+NMR"),
  legend_position = "bottom",
  point_size = 3.2,
  segment_linewidth = 0.5,
  base_size = 11,
  strip_switch = "y"
) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  if (!requireNamespace("stringr", quietly = TRUE)) stop("Package 'stringr' is required for y-axis label formatting.")

  required_cols <- c("strat_name", "stratum", "model", "c_index")
  missing_cols <- setdiff(required_cols, names(cidx))
  if (length(missing_cols) > 0) {
    stop(sprintf(
      "cidx is missing required columns: %s",
      paste(missing_cols, collapse = ", ")
    ))
  }

  if (!is.numeric(soc_global) || length(soc_global) != 1 || is.na(soc_global)) {
    stop("soc_global must be a single non-missing numeric value.")
  }
  if (!is.numeric(nmr_global) || length(nmr_global) != 1 || is.na(nmr_global)) {
    stop("nmr_global must be a single non-missing numeric value.")
  }
  if (length(model_labels) != 2) {
    stop("model_labels must be length 2, e.g. c('SOC', 'SOC+NMR').")
  }

  has_ci <- all(c("ci_lower", "ci_upper") %in% names(cidx))
  if (show_ci && !has_ci) {
    warning("ci_lower/ci_upper not found in cidx; plotting without error bars.")
    show_ci <- FALSE
  }

  soc_raw <- cidx %>%
    dplyr::filter(as.character(.data$model) == soc_model)

  if (!"stratum" %in% names(soc_raw)) {
    if ("n" %in% names(soc_raw)) {
      soc_raw$stratum <- paste0(soc_raw$stratum, " (n=", soc_raw$n, ")")
    } else {
      soc_raw$stratum <- soc_raw$stratum
    }
  }

  soc_df <- soc_raw %>%
    dplyr::transmute(
      strat_name = .data$strat_name,
      stratum = .data$stratum,
      stratum = .data$stratum,
      c_soc = .data$c_index,
      ci_soc_lower = .data$ci_lower,
      ci_soc_upper = .data$ci_upper
    )

  nmr_df <- cidx %>%
    dplyr::filter(as.character(.data$model) == nmr_model) %>%
    dplyr::transmute(
      strat_name = .data$strat_name,
      stratum = .data$stratum,
      c_nmr = .data$c_index,
      ci_nmr_lower = .data$ci_lower,
      ci_nmr_upper = .data$ci_upper
    )

  if (nrow(soc_df) == 0 || nrow(nmr_df) == 0) {
    stop(sprintf(
      "No rows found for soc_model='%s' and/or nmr_model='%s'.",
      soc_model, nmr_model
    ))
  }

  plot_df <- soc_df %>%
    dplyr::inner_join(nmr_df, by = c("strat_name", "stratum"))

  if (nrow(plot_df) == 0) {
    stop("No matched strata after joining SOC and NMR rows.")
  }

  strat_name_levels <- plot_df %>%
    dplyr::distinct(.data$strat_name) %>%
    dplyr::pull(.data$strat_name) %>%
    as.character()

  plot_df <- plot_df %>% dplyr::group_by(.data$strat_name)

  if (sort_by_soc) {
    plot_df <- plot_df %>%
      dplyr::arrange(dplyr::desc(.data$c_soc), .by_group = TRUE)
  }

  plot_df <- plot_df %>%
    dplyr::mutate(
      stratum_plot = factor(
        .data$stratum,
        levels = unique(.data$stratum)
      )
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      strat_name = factor(as.character(.data$strat_name), levels = strat_name_levels)
    )

  pts <- dplyr::bind_rows(
    plot_df %>%
      dplyr::transmute(
        strat_name = .data$strat_name,
        stratum = .data$stratum,
        stratum_plot = .data$stratum_plot,
        model = soc_model,
        c_index = .data$c_soc,
        ci_lower = .data$ci_soc_lower,
        ci_upper = .data$ci_soc_upper
      ),
    plot_df %>%
      dplyr::transmute(
        strat_name = .data$strat_name,
        stratum = .data$stratum,
        stratum_plot = .data$stratum_plot,
        model = nmr_model,
        c_index = .data$c_nmr,
        ci_lower = .data$ci_nmr_lower,
        ci_upper = .data$ci_nmr_upper
      )
  ) %>%
    dplyr::mutate(model = factor(
      .data$model,
      levels = c(soc_model, nmr_model),
      labels = model_labels
    ))

  p <- ggplot2::ggplot(plot_df, ggplot2::aes(y = .data$stratum_plot)) +
    ggplot2::geom_vline(
      xintercept = soc_global,
      linetype = "dashed",
      colour = col_soc,
      linewidth = 0.5,
      alpha = 0.7
    ) +
    ggplot2::geom_vline(
      xintercept = nmr_global,
      linetype = "dashed",
      colour = col_nmr,
      linewidth = 0.5,
      alpha = 0.7
    ) +
    ggplot2::geom_segment(
      ggplot2::aes(x = .data$c_soc, xend = .data$c_nmr, yend = .data$stratum_plot),
      colour = "grey60",
      linewidth = segment_linewidth,
      linetype = "dotted"
    )

  if (show_ci) {
    p <- p +
      ggplot2::geom_errorbarh(
        data = pts,
        ggplot2::aes(
          xmin = .data$ci_lower,
          xmax = .data$ci_upper,
          colour = .data$model
        ),
        height = 0,
        linewidth = 0.45,
        alpha = 0.45
      )
  }

  p <- p +
    ggplot2::geom_point(
      data = pts,
      ggplot2::aes(x = .data$c_index, colour = .data$model),
      size = point_size
    ) +
    ggplot2::scale_colour_manual(
      values = stats::setNames(c(col_soc, col_nmr), model_labels),
      name = NULL
    )

  if (!is.null(xlim)) {
    p <- p + ggplot2::coord_cartesian(xlim = xlim)
  }

  if (is.null(x_breaks) && !is.null(x_by) && !is.null(xlim)) {
    x_breaks <- seq(xlim[[1]], xlim[[2]], by = x_by)
  }

  if (!is.null(x_breaks)) {
    p <- p + ggplot2::scale_x_continuous(
      breaks = x_breaks,
      labels = function(x) sprintf("%.2f", x)
    )
  }

  if (show_global_labels) {
    p <- p +
      ggplot2::annotate(
        "text",
        x = soc_global,
        y = Inf,
        label = "SoC global",
        colour = col_soc,
        vjust = -0.4,
        size = 3
      ) +
      ggplot2::annotate(
        "text",
        x = nmr_global,
        y = Inf,
        label = "NMR global",
        colour = col_nmr,
        vjust = -0.4,
        size = 3
      )
  }

  p <- p +
    ggplot2::facet_grid(
      rows = ggplot2::vars(strat_name),
      cols = ggplot2::vars(),
      scales = "free_y",
      space = "free_y",
      switch = strip_switch
    ) +
    ggplot2::scale_y_discrete(
      labels = function(x) {
        x %>%
          stringr::str_replace(" \\(n=", "\nn=") %>%
          stringr::str_replace_all("\\(|\\)", "")
      }
    ) +
    ggplot2::labs(
      x = "C-index",
      y = NULL,
      title = title,
      subtitle = subtitle
    ) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      legend.position = legend_position,
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major.y = ggplot2::element_blank(),
      panel.spacing.y = ggplot2::unit(0.35, "lines"),
      strip.text.y = ggplot2::element_text(face = "bold", angle = 0, hjust = 0.5),
      strip.background = ggplot2::element_rect(fill = "gray90"),
      plot.title = ggplot2::element_text(hjust = 0.65),
      plot.subtitle = ggplot2::element_text(hjust = 0.5),
      plot.title.position = "plot"
    )

  p
}
