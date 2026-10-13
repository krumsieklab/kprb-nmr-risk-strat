# ============================================================================
# plot_adjacency_cindex()
# ============================================================================
# Discrimination-only subgroup figure. Shows ABSOLUTE per-model C-index
# (NOT a delta) in interpretation bands, with SoC and NMR placed ADJACENT
# within each disease so the SoC-vs-NMR comparison reads directly as
# a colour shift. Default fill is sequential: darker = stronger discrimination.
#
# WHY NOT A DELTA: per-model levels keep poor discrimination under BOTH models
# visible; a delta can look benign when both models are weak.
#
# INPUT (one row per disease x stratum x model), e.g. stacked output from
# subgroup_cindex() with a disease column:
#   required: disease, strat_name, stratum, model
#   value:    c_index  (or pass value_col=)
#
# NB: model strings may be whitespace-padded ("marker_SOC    "); trimmed here.
#
# Usage:
#   source(here::here("helpers", "plot_adjacency_cindex.r"))
#   primary_cindex_df <- purrr::map_dfr(disease_ids, function(id) {
#     cindex_dfs[[id]] %>% dplyr::mutate(disease = dplyr::recode(id, !!!labels))
#   })
#   p <- plot_adjacency_cindex(
#     cindex_df      = primary_cindex_df,
#     diseases       = c("T2D", "NAFLD", "CKD"),
#     disease_labels = c("NAFLD" = "MASLD"),
#     no_title       = TRUE  # <---- Add this argument to suppress title space
#   )
#   p
# ============================================================================

if (!exists("%>%", mode = "function")) {
  `%>%` <- dplyr::`%>%`
}

.default_cindex_stratum_order <- function() {
  list(
    Age  = c("Age < 55", "Age > 55"),
    BMI  = c("BMI < 30", "BMI > 30"),
    Sex  = c("Female", "Male"),
    Ethnicity = c("Asian", "Black", "Hispanic", "White"),
    PRS  = c("Low PRS", "High PRS")
  )
}

.cindex_fill_defaults <- function() {
  list(
    breaks = c(-Inf, 0.60, 0.65, 0.70, 0.75, 0.80, Inf),
    labels = c(
      "< 0.60 (negligible)",
      "0.60-0.65 (near-uninformative)",
      "0.65-0.70 (weak)",
      "0.70-0.75 (moderate)",
      "0.75-0.80 (good)",
      ">= 0.80 (strong)"
    ),
    colours = c("#FBE7E7", "#F3F0E4", "#DDECCB", "#A8DDB5", "#5AB4AC", "#01665E")
  )
}

.cindex_fill_spec <- function(
  fill_values,
  fill_limits = NULL,
  fill_colours = NULL,
  fill_anchors = NULL
) {
  fd <- .cindex_fill_defaults()
  if (is.null(fill_colours)) fill_colours <- fd$colours
  if (length(fill_colours) != length(fd$labels)) {
    stop(sprintf(
      "fill_colours must have length %s, one colour per C-index band.",
      length(fd$labels)
    ))
  }

  list(
    breaks = fd$breaks,
    labels = fd$labels,
    colours = stats::setNames(fill_colours, fd$labels)
  )
}

.cindex_prepare_df <- function(
  cindex_df,
  diseases = NULL,
  disease_labels = NULL,
  drop_strata = NULL,
  strat_name_levels = names(stratum_order),
  stratum_order = .default_cindex_stratum_order(),
  value_col = NULL,
  global_strat_name = "All",
  global_stratum_label = "All",
  global_facet_label = "",
  model_levels = c("marker_SOC", "marker_SOC+NMR"),
  model_labels = c("marker_SOC" = "SOC", "marker_SOC+NMR" = "SOC+NMR")
) {
  required <- c("disease", "strat_name", "stratum", "model")
  missing  <- setdiff(required, names(cindex_df))
  if (length(missing) > 0) {
    stop(sprintf(paste0(
      "cindex_df is missing required columns: %s.\n",
      "  Stack per-disease tables with a 'disease' column first, e.g.\n",
      "  purrr::map_dfr(disease_ids, ~ cindex_dfs[[.x]] %>% mutate(disease = ...))"),
      paste(missing, collapse = ", ")))
  }

  df <- cindex_df
  df$model <- trimws(as.character(df$model))

  if (is.null(value_col)) {
    has_c_index <- "c_index" %in% names(df)
    has_cindex <- "cindex" %in% names(df)
    if (has_c_index && has_cindex) {
      df$.value <- dplyr::coalesce(df$c_index, df$cindex)
    } else if (has_c_index) {
      df$.value <- df$c_index
    } else if (has_cindex) {
      df$.value <- df$cindex
    } else {
      stop("Need a 'c_index' (or 'cindex') column in cindex_df.")
    }
  } else {
    if (!(value_col %in% names(df))) {
      stop(sprintf("value_col '%s' not found in cindex_df.", value_col))
    }
    df$.value <- df[[value_col]]
  }

  df$.fill  <- df$.value
  df$.val_label <- formatC(df$.fill, format = "f", digits = 3)

  if (!is.null(drop_strata)) {
    df <- df[!(as.character(df$strat_name) %in% drop_strata), , drop = FALSE]
  }
  if (!is.null(diseases)) {
    df <- df[as.character(df$disease) %in% diseases, , drop = FALSE]
  }
  df <- df[df$model %in% model_levels, , drop = FALSE]
  if (nrow(df) == 0) stop("No C-index rows after filtering.")

  is_global_row <- as.character(df$strat_name) == global_strat_name & is.na(df$stratum)
  df$stratum[is_global_row] <- global_stratum_label
  df$strat_name[is_global_row] <- global_facet_label

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
  if (any(as.character(df$stratum) == global_stratum_label)) {
    stratum_levels <- c(global_stratum_label, stratum_levels)
  }
  stratum_levels <- c(stratum_levels, setdiff(unique(df$stratum), stratum_levels))
  present_strata <- intersect(strat_name_levels, unique(as.character(df$strat_name)))
  if (any(as.character(df$stratum) == global_stratum_label)) {
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

plot_adjacency_cindex <- function(
  cindex_df,
  diseases = c("T2D", "NAFLD", "CKD"),
  disease_labels = c("NAFLD" = "MASLD"),
  drop_strata = c("PRS 80%/20%"),
  value_col = NULL,
  fill_limits = c(0.50, 0.90),
  fill_colours = NULL,
  fill_anchors = NULL,
  model_levels = c("marker_SOC", "marker_SOC+NMR"),
  model_labels = c("marker_SOC" = "SOC", "marker_SOC+NMR" = "NMR"),
  strat_name_levels = names(stratum_order),
  stratum_order = .default_cindex_stratum_order(),
  global_strat_name = "All",
  global_stratum_label = "All",
  global_facet_label = "",
  title = "C-index by model",
  subtitle = NULL,
  panel_tag = "A",
  show_values = TRUE,
  value_size = 2.5,
  base_size = 11,
  tile_border = "white",
  legend_position = "right",
  return_data = FALSE,
  no_title = FALSE   # <-- NEW ARGUMENT. Use TRUE to hide both title and subtitle.
) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Package 'dplyr' is required.")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Package 'ggplot2' is required.")

  plot_df <- .cindex_prepare_df(
    cindex_df = cindex_df,
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

  fill_spec <- .cindex_fill_spec(
    fill_values = plot_df$.fill,
    fill_limits = fill_limits,
    fill_colours = fill_colours,
    fill_anchors = fill_anchors
  )
  plot_df$.fill_band <- cut(
    plot_df$.fill,
    breaks = fill_spec$breaks,
    labels = fill_spec$labels,
    include.lowest = TRUE,
    right = FALSE
  )
  plot_df$.text_colour <- ifelse(plot_df$.fill >= 0.80, "white", "grey15")

  # Make both title and subtitle NULL if no_title = TRUE
  if (no_title) {
    title <- NULL
    subtitle <- NULL
  } else {
    if (is.null(subtitle)) {
      subtitle <- "Absolute C-index per model in interpretation bands \u00b7 darker is better"
    }
    if (!is.null(panel_tag) && nzchar(panel_tag)) {
      title <- paste0(panel_tag, "  ", title)
    }
  }

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = .data$model, y = .data$stratum, fill = .data$.fill_band)
  ) +
    ggplot2::geom_tile(colour = tile_border, linewidth = 0.35) +
    ggplot2::facet_grid(
      rows = ggplot2::vars(strat_name),
      cols = ggplot2::vars(disease),
      scales = "free_y",
      space = "free_y",
      switch = "y"
    ) +
    ggplot2::scale_fill_manual(
      values = fill_spec$colours,
      drop = FALSE,
      name = "C-index",
      na.value = "grey90"
    ) +
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
        face = "bold", angle = 0, hjust = 0.5, vjust = 0.3
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
      size = value_size,
    ) +
      ggplot2::scale_colour_identity(guide = "none")
  }

  if (return_data) attr(p, "data") <- plot_df
  p
}
