# ============================================================================
# plot_calibration_adequacy_heatmap()
# ============================================================================
# Calibration-only subgroup figure. Shows ABSOLUTE per-model calibration
# adequacy (NOT a delta) on a continuous colour scale, with SoC and NMR placed
# ADJACENT within each disease so the SoC-vs-NMR comparison reads directly as
# a colour shift. Two value schemes:
#
#   metric = "absolute"  (default, recommended)
#       raw ICI on a continuous scale (legend in percentage points).
#
#   metric = "relative"  (cross-disease comparability, but use with care)
#       ICI / observed incidence. Comparable across diseases with very
#       different base rates, BUT penalises rare-outcome strata: a tiny absolute
#       error becomes a large fraction of a small base rate, so well-calibrated
#       rare strata can read high. Not a published standard (ICI: Austin &
#       Steyerberg 2019) -- interpret with care.
#
# WHY NOT A DELTA: calibration is not an incremental-value question, and a
# delta is blind to level -- a stratum bad under BOTH models has delta ~0 and
# looks benign. Per-model levels keep that visible.
#
# INPUT (one row per disease x stratum x model):
#   required: disease, strat_name, stratum, model
#   value:    ici  OR  calibration_score (read as 1 - ICI)
#   relative: incidence column (pass incidence_col=); if absent, falls back to
#             crude n_events/n with a warning (prefer KM cumulative incidence
#             at t_star).
#
# NB: model strings may be whitespace-padded ("marker_SOC    "); trimmed here.
#
# Usage:
#   source(here::here("helpers", "plot_adjacency_calibration.r"))
#   p <- plot_calibration_adequacy_heatmap(
#     calibration_df = primary_calibration_df,        # must have a 'disease' col
#     diseases       = c("T2D", "NAFLD", "CKD"),
#     disease_labels = c("NAFLD" = "MASLD"),
#     metric         = "absolute"
#   )
#   p
#
# Optional global / overall reference rows (same convention as
# plot_adjacency_cindex / plot_adjacency_calibration_score):
#   bind rows with strat_name = "All", stratum = NA, ici, and model.
#   These appear as a top "All" band.
#
#   p <- plot_calibration_adequacy_heatmap(
#     dplyr::bind_rows(primary_calibration_df, global_calibration_df),
#     metric = "absolute"
#   )
# ============================================================================

if (!exists("%>%", mode = "function")) {
  `%>%` <- dplyr::`%>%`
}

.default_calibration_stratum_order <- function() {
  list(
    Age  = c("Age < 55", "Age > 55"),
    BMI  = c("BMI < 30", "BMI > 30"),
    Sex  = c("Female", "Male"),
    Ethnicity = c("Asian", "Black", "Hispanic", "White"),
    PRS  = c("Low PRS", "High PRS")
  )
}

.calibration_fill_defaults <- function(metric) {
  if (metric == "absolute") {
    list(
      colours = c("#8FCBB0", "#F0C56B", "#D9594C"),
      anchors = c(0, 5, 10)   # percentage points; upper anchor may extend to fill_limits
    )
  } else {
    list(
      colours = c("#3E9B6E", "#B7DBA8", "#F0C56B", "#D9594C"),
      anchors = c(0, 0.10, 0.25, 0.50)
    )
  }
}

.calibration_fill_spec <- function(
  metric,
  fill_values,
  fill_limits = NULL,
  fill_colours = NULL,
  fill_anchors = NULL
) {
  fd <- .calibration_fill_defaults(metric)
  if (is.null(fill_colours)) fill_colours <- fd$colours
  if (is.null(fill_anchors)) fill_anchors <- fd$anchors

  if (is.null(fill_limits)) {
    hi <- max(fill_values, na.rm = TRUE)
    fill_limits <- c(
      0,
      if (metric == "absolute") max(hi * 1.05, fd$anchors[length(fd$anchors)])
      else max(hi * 1.05, tail(fd$anchors, 1))
    )
  }

  anchors <- unique(c(fill_anchors, fill_limits[2]))
  anchors <- anchors[anchors >= fill_limits[1] & anchors <= fill_limits[2]]
  if (length(anchors) < 2) anchors <- fill_limits

  if (length(fill_colours) == 1) {
    values <- c(0, 1)
    fill_colours <- rep(fill_colours, 2)
  } else if (length(anchors) == length(fill_colours)) {
    values <- scales::rescale(anchors, from = fill_limits, to = c(0, 1))
  } else {
    values <- seq(0, 1, length.out = length(fill_colours))
  }

  list(
    limits = fill_limits,
    colours = fill_colours,
    values = values
  )
}

.calibration_prepare_df <- function(
  calibration_df,
  metric,
  diseases = NULL,
  disease_labels = NULL,
  drop_strata = NULL,
  strat_name_levels = names(stratum_order),
  stratum_order = .default_calibration_stratum_order(),
  value_col = NULL,
  incidence_col = NULL,
  global_strat_name = "All",
  global_stratum_label = "All",
  global_facet_label = "",
  model_levels = c("marker_SOC", "marker_SOC+NMR"),
  model_labels = c("marker_SOC" = "SoC", "marker_SOC+NMR" = "NMR")
) {
  required <- c("disease", "strat_name", "stratum", "model")
  missing  <- setdiff(required, names(calibration_df))
  if (length(missing) > 0) {
    stop(sprintf(paste0(
      "calibration_df is missing required columns: %s.\n",
      "  If your per-disease tables are separate, rbind them with a 'disease' column first, e.g.\n",
      "  dplyr::bind_rows(T2D = cal_t2d, NAFLD = cal_nafld, CKD = cal_ckd, .id = 'disease')"),
      paste(missing, collapse = ", ")))
  }

  df <- calibration_df
  df$model <- trimws(as.character(df$model))   # source strings can be padded

  # Remap global reference rows BEFORE as.character(stratum) turns NA into "NA".
  # Same convention as plot_adjacency_cindex(): strat_name == "All", stratum = NA.
  is_global_row <- as.character(df$strat_name) == global_strat_name & is.na(df$stratum)
  df$stratum[is_global_row] <- global_stratum_label
  df$strat_name[is_global_row] <- global_facet_label

  # --- ICI ---
  if (is.null(value_col)) {
    value_col <- if ("ici" %in% names(df)) "ici"
                 else if ("calibration_score" %in% names(df)) "calibration_score"
                 else stop("Need an 'ici' or 'calibration_score' column.")
  }
  if (!(value_col %in% names(df))) {
    stop(sprintf("value_col '%s' not found in calibration_df.", value_col))
  }
  df$.ici <- if (identical(value_col, "calibration_score")) 1 - df[[value_col]] else df[[value_col]]

  # --- filters ---
  if (!is.null(drop_strata))
    df <- df[!(as.character(df$strat_name) %in% drop_strata), , drop = FALSE]
  if (!is.null(diseases))
    df <- df[as.character(df$disease) %in% diseases, , drop = FALSE]
  df <- df[df$model %in% model_levels, , drop = FALSE]
  if (nrow(df) == 0) stop("No calibration rows after filtering.")

  # --- continuous fill value ---
  if (metric == "absolute") {
    df$.value <- df$.ici
    df$.fill  <- df$.ici * 100   # legend / tile labels in percentage points
  } else {
    if (!is.null(incidence_col)) {
      if (!(incidence_col %in% names(df)))
        stop(sprintf("incidence_col '%s' not found in calibration_df.", incidence_col))
      denom <- df[[incidence_col]]
    } else if (all(c("n_events", "n") %in% names(df))) {
      warning("metric='relative' without incidence_col: using crude n_events/n as the denominator. Prefer observed (KM) cumulative incidence at t_star.")
      denom <- df$n_events / df$n
    } else {
      stop("metric='relative' needs an incidence_col, or n_events and n to fall back on.")
    }
    denom[!is.finite(denom) | denom <= 0] <- NA_real_
    df$.value <- df$.ici / denom
    df$.fill  <- df$.value
  }

  # --- optional tile value label ---
  df$.val_label <- if (metric == "absolute") {
    formatC(df$.fill, format = "f", digits = 1)
  } else {
    formatC(df$.fill, format = "f", digits = 2)
  }

  # --- disease labels / ordering ---
  if (!is.null(disease_labels))
    df$disease <- dplyr::recode(as.character(df$disease), !!!disease_labels)
  disease_levels <- if (!is.null(diseases)) {
    if (!is.null(disease_labels)) unname(dplyr::recode(diseases, !!!disease_labels)) else diseases
  } else unique(df$disease)

  # --- factor ordering ---
  stratum_levels <- unlist(
    stratum_order[intersect(strat_name_levels, names(stratum_order))],
    use.names = FALSE)
  if (any(as.character(df$stratum) == global_stratum_label, na.rm = TRUE)) {
    stratum_levels <- c(global_stratum_label, stratum_levels)
  }
  stratum_levels <- c(stratum_levels, setdiff(unique(as.character(df$stratum)), stratum_levels))
  present_strata <- intersect(strat_name_levels, unique(as.character(df$strat_name)))
  if (any(as.character(df$stratum) == global_stratum_label, na.rm = TRUE)) {
    present_strata <- c(global_facet_label, present_strata)
  }

  df %>%
    dplyr::mutate(
      disease    = factor(.data$disease, levels = disease_levels),
      strat_name = factor(as.character(.data$strat_name), levels = present_strata),
      stratum    = factor(as.character(.data$stratum), levels = stratum_levels),
      model      = factor(dplyr::recode(.data$model, !!!model_labels),
                          levels = unname(model_labels[model_levels]))
    ) %>%
    dplyr::filter(
      !is.na(.data$strat_name),
      !is.na(.data$stratum),
      !is.na(.data$disease),
      !is.na(.data$.fill)
    )
}

plot_calibration_adequacy_heatmap <- function(
  calibration_df,
  diseases = c("T2D", "NAFLD", "CKD"),
  disease_labels = c("NAFLD" = "MASLD"),
  metric = c("absolute", "relative"),
  drop_strata = c("PRS 80%/20%"),
  value_col = NULL,
  incidence_col = NULL,
  fill_limits = NULL,
  fill_colours = NULL,
  fill_anchors = NULL,
  model_levels = c("marker_SOC", "marker_SOC+NMR"),
  model_labels = c("marker_SOC" = "SoC", "marker_SOC+NMR" = "NMR"),
  strat_name_levels = names(stratum_order),
  stratum_order = .default_calibration_stratum_order(),
  global_strat_name = "All",
  global_stratum_label = "All",
  global_facet_label = "",
  title = NULL,
  subtitle = NULL,
  panel_tag = NULL,
  show_values = FALSE,
  value_size = 2.5,
  base_size = 11,
  tile_border = "white",
  legend_position = "right",
  return_data = FALSE
) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")
  if (!requireNamespace("scales", quietly = TRUE)) stop("Package 'scales' is required.")
  metric <- match.arg(metric)

  cal <- .calibration_prepare_df(
    calibration_df, metric = metric, diseases = diseases,
    disease_labels = disease_labels, drop_strata = drop_strata,
    strat_name_levels = strat_name_levels, stratum_order = stratum_order,
    value_col = value_col, incidence_col = incidence_col,
    global_strat_name = global_strat_name,
    global_stratum_label = global_stratum_label,
    global_facet_label = global_facet_label,
    model_levels = model_levels, model_labels = model_labels)

  fill_spec <- .calibration_fill_spec(
    metric = metric,
    fill_values = cal$.fill,
    fill_limits = fill_limits,
    fill_colours = fill_colours,
    fill_anchors = fill_anchors
  )

  fill_name <- if (metric == "absolute") "ICI (pp)" else "Relative ICI"
  fill_labels <- if (metric == "absolute") {
    scales::label_number(accuracy = 1, suffix = " pp")
  } else {
    scales::label_number(accuracy = 0.01)
  }

  if (!is.null(panel_tag) && nzchar(panel_tag) && !is.null(title)) {
    title <- paste0(panel_tag, "  ", title)
  }

  p <- ggplot2::ggplot(
    cal, ggplot2::aes(x = .data$model, y = .data$stratum, fill = .data$.fill)
  ) +
    ggplot2::geom_tile(colour = tile_border, linewidth = 0.35) +
    ggplot2::facet_grid(rows = ggplot2::vars(strat_name),
                        cols = ggplot2::vars(disease),
                        scales = "free_y", space = "free_y", switch = "y") +
    ggplot2::scale_fill_gradientn(
      colours = fill_spec$colours,
      values = fill_spec$values,
      limits = fill_spec$limits,
      oob = scales::squish,
      name = fill_name,
      labels = fill_labels,
      na.value = "grey90") +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.subtitle = ggplot2::element_text(colour = "grey35", size = base_size - 2),
      axis.text.x   = ggplot2::element_text(
        face = "bold", size = base_size - 1,
        angle = 45, hjust = 0.5, vjust = 0.5),
      axis.text.y   = ggplot2::element_text(size = base_size - 1),
      strip.text.x  = ggplot2::element_text(
        face = "bold", angle = 0, hjust = 0.5, vjust = 0.3),
      strip.text.y.left = ggplot2::element_text(face = "bold", angle = 0, hjust = 0),
      strip.background  = ggplot2::element_blank(),
      strip.placement   = "outside",
      panel.grid      = ggplot2::element_blank(),
      panel.spacing.y = ggplot2::unit(0.15, "lines"),
      panel.spacing.x = ggplot2::unit(0.6, "lines"),
      legend.position = legend_position,
      legend.key.width  = ggplot2::unit(0.5, "cm"),
      legend.key.height = ggplot2::unit(0.35, "cm"),
      plot.margin = ggplot2::margin(6, 8, 4, 8)
    )

  if (show_values) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = .data$.val_label),
      size = value_size, colour = "grey15")
  }

  if (return_data) attr(p, "data") <- cal
  p
}
