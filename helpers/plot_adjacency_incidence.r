# ============================================================================
# plot_adjacency_incidence()
# ============================================================================
# Enrichment subgroup figure. Shows ABSOLUTE per-model observed incidence
# within a top-risk percentile (NOT a delta), with SoC and NMR placed
# ADJACENT within each disease so the SoC-vs-NMR comparison reads directly as
# a colour shift.
#
# By default the fill is DISCRETIZED into interpretation bands (mirrors
# plot_adjacency_cindex()'s banded C-index scale) rather than a continuous
# gradient, since raw incidence spans a wide range across diseases (roughly
# 0-90% across strata: CKD/T2D top-1% incidence reaches ~80%+, while
# AcuteMI/CHD/COPD/TIA-Stroke top-1% incidence mostly sits under ~20%).
# Discrete bands make within-disease SoC-vs-NMR shifts (and cross-disease
# comparisons) easier to read at a glance than a continuous gradient.
# Pass discrete = FALSE to fall back to the original continuous gradient.
#
# INPUT (one row per disease x stratum x model x top_prop), e.g. stacked output
# from subgroup_top_risk_incidence_calculator() with a disease column:
#   required: disease, strat_name, stratum, model, top_label (or top_prop)
#   value:    obs_incidence  (or pass value_col=)
#
# Filter to one percentile via top_label (default "Top 1%") or top_prop.
#
# Usage:
#   source(here::here("helpers", "plot_adjacency_incidence.r"))
#   primary_top_df <- purrr::map_dfr(disease_ids, function(id) {
#     top_dfs[[id]] %>% dplyr::mutate(disease = dplyr::recode(id, !!!labels))
#   })
#   plot_adjacency_incidence(
#     incidence_df = primary_top_df,
#     top_label      = "Top 1%",
#     diseases       = c("T2D", "NAFLD", "CKD"),
#     disease_labels = c("NAFLD" = "MASLD")
#   )
#
# Optional global / overall reference rows (same pattern as plot_adjacency_cindex):
#   bind rows with strat_name = "All", stratum = NA, obs_incidence, and matching
#   model / top_label columns. These appear as a top "All" band.
#
#   global_incidence_df <- data.frame(
#     disease = rep(c("T2D", "CKD"), each = 2),
#     stratum = NA, strat_name = "All",
#     model = rep(c("marker_SOC", "marker_SOC+NMR"), 2),
#     top_label = "Top 1%",
#     obs_incidence = c(0.45, 0.72, 0.50, 0.80)
#   )
#   plot_adjacency_incidence(
#     dplyr::bind_rows(primary_top_df, global_incidence_df),
#     top_label = "Top 1%"
#   )
#
#   # Customize the bands (breaks must have one more element than labels/colours):
#   plot_adjacency_incidence(
#     incidence_df = primary_top_df,
#     band_breaks  = c(-Inf, 0.10, 0.25, 0.50, Inf),
#     band_labels  = c("< 10%", "10-25%", "25-50%", ">= 50%"),
#     band_colours = c("#EDEAE4", "#B7DBA8", "#8FCBB0", "#3E9B6E")
#   )
#
#   # Original continuous gradient:
#   plot_adjacency_incidence(incidence_df = primary_top_df, discrete = FALSE)
# ============================================================================

if (!exists("%>%", mode = "function")) {
  `%>%` <- dplyr::`%>%`
}

.default_incidence_stratum_order <- function() {
  list(
    Age  = c("Age < 55", "Age > 55"),
    BMI  = c("BMI < 30", "BMI > 30"),
    Sex  = c("Female", "Male"),
    Ethnicity = c("Asian", "Black", "Hispanic", "White"),
    PRS  = c("Low PRS", "High PRS")
  )
}

.incidence_fill_defaults <- function() {
  list(
    limits  = c(0, 1),
    colours = c("#EDEAE4", "#E0DDD6", "#B7DBA8", "#8FCBB0", "#3E9B6E"),
    anchors = c(0, 0.25, 0.50, 0.75, 1.0)
  )
}

.incidence_fill_spec <- function(
  fill_values,
  fill_limits = NULL,
  fill_colours = NULL,
  fill_anchors = NULL
) {
  if (!requireNamespace("scales", quietly = TRUE)) {
    stop("Package 'scales' is required.")
  }

  fd <- .incidence_fill_defaults()
  if (is.null(fill_colours)) fill_colours <- fd$colours
  if (is.null(fill_anchors)) fill_anchors <- fd$anchors
  if (is.null(fill_limits)) fill_limits <- fd$limits

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

# ---- Discrete band defaults -----------------------------------------------
# Chosen from the observed spread of top-1% incidence across the 9 primary
# diseases (min ~0, max ~0.89): low-incidence diseases (AcuteMI, CHD, COPD,
# TIA/Stroke) mostly fall in the first two bands; mid-incidence diseases
# (AFib, NAFLD, HF/CM) mostly land in the middle bands; CKD/T2D reach into
# the top bands. Override via band_breaks / band_labels / band_colours.
.incidence_band_defaults <- function() {
  list(
    breaks = c(-Inf, 0.05, 0.15, 0.30, 0.45, 0.65, 0.80, Inf),
    labels = c(
      "< 5%",
      "5-15%",
      "15-30%",
      "30-45%",
      "45-65%",
      "65-80%",
      ">= 80%"
    ),
    colours = c(
      "#F7F4EE", "#EDEAE4", "#D9E8C9", "#B7DBA8", "#8FCBB0",
      "#3E9B6E", "#145A32"
    )
  )
}

.incidence_band_spec <- function(
  band_breaks = NULL,
  band_labels = NULL,
  band_colours = NULL
) {
  bd <- .incidence_band_defaults()
  if (is.null(band_breaks))  band_breaks  <- bd$breaks
  if (is.null(band_labels))  band_labels  <- bd$labels
  if (is.null(band_colours)) band_colours <- bd$colours

  if (length(band_breaks) != length(band_labels) + 1) {
    stop("band_breaks must have exactly one more element than band_labels (n bins = n breaks - 1).")
  }
  if (length(band_colours) != length(band_labels)) {
    stop("band_colours must have the same length as band_labels, one colour per incidence band.")
  }

  list(
    breaks  = band_breaks,
    labels  = band_labels,
    colours = stats::setNames(band_colours, band_labels)
  )
}

.incidence_prepare_df <- function(
  incidence_df,
  top_label = "Top 1%",
  top_prop = NULL,
  diseases = NULL,
  disease_labels = NULL,
  drop_strata = NULL,
  strat_name_levels = names(stratum_order),
  stratum_order = .default_incidence_stratum_order(),
  value_col = NULL,
  global_strat_name = "All",
  global_stratum_label = "All",
  global_facet_label = "",
  model_levels = c("marker_SOC", "marker_SOC+NMR"),
  model_labels = c("marker_SOC" = "SoC", "marker_SOC+NMR" = "NMR")
) {
  required <- c("disease", "strat_name", "stratum", "model")
  missing  <- setdiff(required, names(incidence_df))
  if (length(missing) > 0) {
    stop(sprintf(paste0(
      "incidence_df is missing required columns: %s.\n",
      "  Stack per-disease tables with a 'disease' column first, e.g.\n",
      "  purrr::map_dfr(disease_ids, ~ top_dfs[[.x]] %>% mutate(disease = ...))"),
      paste(missing, collapse = ", ")))
  }

  df <- incidence_df
  df$model <- trimws(as.character(df$model))

  # Remap global reference rows BEFORE as.character(stratum) turns NA into "NA".
  # Same convention as plot_adjacency_cindex(): strat_name == "All", stratum = NA.
  is_global_row <- as.character(df$strat_name) == global_strat_name & is.na(df$stratum)
  df$stratum[is_global_row] <- global_stratum_label
  df$strat_name[is_global_row] <- global_facet_label

  if (!is.null(top_prop)) {
    if (!("top_prop" %in% names(df))) {
      stop("top_prop filter requested but 'top_prop' column is missing.")
    }
    keep <- abs(df$top_prop - top_prop) < sqrt(.Machine$double.eps)
    keep[is.na(keep)] <- FALSE
    # Global rows with a missing top_prop are treated as an overlay for this plot.
    # Do NOT keep every global percentile — that overplots the All band.
    missing_prop <- is_global_row & is.na(df$top_prop)
    df <- df[keep | missing_prop, , drop = FALSE]
  } else if (!is.null(top_label)) {
    if (!("top_label" %in% names(df))) {
      stop("top_label filter requested but 'top_label' column is missing.")
    }
    keep <- trimws(as.character(df$top_label)) == trimws(top_label)
    keep[is.na(keep)] <- FALSE
    # Keep unlabeled global rows as an overlay; still require a matching
    # top_label when it is present (e.g. Top 1% / 5% / 10% All rows).
    missing_label <- is_global_row & (
      is.na(df$top_label) | !nzchar(trimws(as.character(df$top_label)))
    )
    df <- df[keep | missing_label, , drop = FALSE]
  }

  if (nrow(df) == 0) {
    stop(sprintf(
      "No rows match top_label='%s'%s.",
      top_label,
      if (!is.null(top_prop)) sprintf(" / top_prop=%s", top_prop) else ""
    ))
  }

  # Recompute after filtering (row subset may have changed)
  is_global_row <- as.character(df$strat_name) == global_facet_label &
    as.character(df$stratum) == global_stratum_label

  if (is.null(value_col)) {
    value_col <- if ("obs_incidence" %in% names(df)) {
      "obs_incidence"
    } else {
      stop("Need an 'obs_incidence' column in incidence_df.")
    }
  }
  if (!(value_col %in% names(df))) {
    stop(sprintf("value_col '%s' not found in incidence_df.", value_col))
  }

  df$.value <- df[[value_col]]
  df$.fill  <- df$.value
  df$.val_label <- formatC(100 * df$.fill, format = "f", digits = 1)

  if (!is.null(drop_strata)) {
    df <- df[!(as.character(df$strat_name) %in% drop_strata), , drop = FALSE]
  }
  if (!is.null(diseases)) {
    df <- df[as.character(df$disease) %in% diseases, , drop = FALSE]
  }
  df <- df[df$model %in% model_levels, , drop = FALSE]
  if (nrow(df) == 0) stop("No incidence rows after filtering.")

  if (!is.null(disease_labels)) {
    df$disease <- dplyr::recode(as.character(df$disease), !!!disease_labels)
  }
  disease_levels <- if (!is.null(diseases)) {
    if (!is.null(disease_labels)) {
      unname(dplyr::recode(diseases, !!!disease_labels))
    } else {
      diseases
    }
  } else {
    unique(df$disease)
  }

  stratum_levels <- unlist(
    stratum_order[intersect(strat_name_levels, names(stratum_order))],
    use.names = FALSE
  )
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
      model      = factor(
        dplyr::recode(.data$model, !!!model_labels),
        levels = unname(model_labels[model_levels])
      )
    ) %>%
    dplyr::filter(
      !is.na(.data$strat_name),
      !is.na(.data$stratum),
      !is.na(.data$disease),
      !is.na(.data$.fill)
    )
}

plot_adjacency_incidence <- function(
  incidence_df,
  top_label = "Top 1%",
  top_prop = NULL,
  diseases = c("T2D", "NAFLD", "CKD"),
  disease_labels = c("NAFLD" = "MASLD"),
  drop_strata = c("PRS 80%/20%"),
  value_col = NULL,
  discrete = TRUE,
  band_breaks = NULL,
  band_labels = NULL,
  band_colours = NULL,
  dark_band_threshold = NULL,
  fill_limits = c(0, 1),
  fill_colours = NULL,
  fill_anchors = NULL,
  model_levels = c("marker_SOC", "marker_SOC+NMR"),
  model_labels = c("marker_SOC" = "SoC", "marker_SOC+NMR" = "NMR"),
  strat_name_levels = names(stratum_order),
  stratum_order = .default_incidence_stratum_order(),
  global_strat_name = "All",
  global_stratum_label = "All",
  global_facet_label = "",
  title = NULL,
  subtitle = NULL,
  panel_tag = NULL,
  no_title = TRUE,
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

  plot_label <- if (!is.null(top_label)) {
    top_label
  } else {
    sprintf("Top %.0f%%", 100 * top_prop)
  }

  plot_df <- .incidence_prepare_df(
    incidence_df = incidence_df,
    top_label = top_label,
    top_prop = top_prop,
    diseases = diseases,
    disease_labels = disease_labels,
    drop_strata = drop_strata,
    strat_name_levels = strat_name_levels,
    stratum_order = stratum_order,
    value_col = value_col,
    global_strat_name = global_strat_name,
    global_stratum_label = global_stratum_label,
    global_facet_label = global_facet_label,
    model_levels = model_levels,
    model_labels = model_labels
  )

  if (isTRUE(no_title)) {
    title <- NULL
    subtitle <- NULL
  } else {
    if (is.null(title)) {
      title <- sprintf("%s observed incidence by model", plot_label)
    }
    if (is.null(subtitle)) {
      subtitle <- if (isTRUE(discrete)) {
        sprintf(
          "Absolute %s observed incidence per model, in bands \u00b7 darker is higher",
          plot_label
        )
      } else {
        sprintf(
          "Absolute %s observed incidence per model (not differenced) \u00b7 higher is better",
          plot_label
        )
      }
    }
    if (!is.null(panel_tag) && nzchar(panel_tag)) {
      title <- paste0(panel_tag, "  ", title)
    }
  }

  if (isTRUE(discrete)) {
    band_spec <- .incidence_band_spec(
      band_breaks  = band_breaks,
      band_labels  = band_labels,
      band_colours = band_colours
    )
    plot_df$.fill_plot <- cut(
      plot_df$.fill,
      breaks = band_spec$breaks,
      labels = band_spec$labels,
      include.lowest = TRUE,
      right = FALSE
    )
    if (is.null(dark_band_threshold)) {
      # White labels for the top two green bands (65-80% and >=80% by default).
      finite_breaks <- band_spec$breaks[is.finite(band_spec$breaks)]
      dark_band_threshold <- if (length(finite_breaks) >= 2) {
        utils::tail(finite_breaks, 2)[1]
      } else {
        utils::tail(finite_breaks, 1)
      }
    }
    plot_df$.text_colour <- ifelse(plot_df$.fill >= dark_band_threshold, "white", "grey15")
  } else {
    fill_spec <- .incidence_fill_spec(
      fill_values = plot_df$.fill,
      fill_limits = fill_limits,
      fill_colours = fill_colours,
      fill_anchors = fill_anchors
    )
    plot_df$.fill_plot <- plot_df$.fill
    plot_df$.text_colour <- "grey15"
  }

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$model, y = .data$stratum, fill = .data$.fill_plot)
  ) +
    ggplot2::geom_tile(colour = tile_border, linewidth = 0.35) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(strat_name),
      cols = ggplot2::vars(disease),
      scales = "free_y",
      space = "free_y",
      switch = "y"
    )

  if (isTRUE(discrete)) {
    p <- p + ggplot2::scale_fill_manual(
      values = band_spec$colours,
      drop = FALSE,
      name = "Obs. Incident\n(Top 5%)",
      na.value = "grey90"
    )
  } else {
    p <- p + ggplot2::scale_fill_gradientn(
      colours = fill_spec$colours,
      values = fill_spec$values,
      limits = fill_spec$limits,
      oob = scales::squish,
      name = "Obs. Incident\n(Top 5%)",
      labels = scales::label_percent(accuracy = 1),
      na.value = "grey90"
    )
  }

  p <- p +
    ggplot2::labs(title = title, subtitle = subtitle, x = NULL, y = NULL) +
    ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(
      plot.title    = ggplot2::element_text(face = "bold", size = base_size + 1),
      plot.subtitle = ggplot2::element_text(colour = "grey35", size = base_size - 2),
      axis.text.x   = ggplot2::element_text(
        face = "bold", size = base_size - 1,
        angle = 45, hjust = 0.5, vjust = 0.5
      ),
      axis.text.y   = ggplot2::element_text(size = base_size - 1),
      strip.text.x  = ggplot2::element_text(
        face = "bold", angle = 0, hjust = 0.5, vjust = 0.5
      ),
      strip.text.y.left = ggplot2::element_text(face = "bold", angle = 0, hjust = 0),
      strip.background  = ggplot2::element_blank(),
      strip.placement   = "outside",
      panel.grid        = ggplot2::element_blank(),
      panel.spacing.y   = ggplot2::unit(0.15, "lines"),
      panel.spacing.x   = ggplot2::unit(0.6, "lines"),
      legend.position   = legend_position,
      legend.key.width  = ggplot2::unit(0.5, "cm"),
      legend.key.height = ggplot2::unit(0.35, "cm"),
      plot.margin = ggplot2::margin(6, 8, 4, 8)
    )

  if (show_values) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = .data$.val_label, colour = .data$.text_colour),
      size = value_size
    ) +
      ggplot2::scale_colour_identity(guide = "none")
  }

  if (return_data) attr(p, "data") <- plot_df
  p
}
