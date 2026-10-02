# Observed incidence from precomputed out-of-fold absolute risks.
#
# Percentile curves bin a supplied risk column and estimate Kaplan-Meier
# incidence in each bin. Top-risk incidence ranks that same column and
# estimates incidence in the top predicted-risk groups.

.enrichment_oof <- function(results) {
  if (is.data.frame(results)) {
    results
  } else if (is.list(results) && !is.data.frame(results) &&
             "oof_predictions" %in% names(results)) {
    results$oof_predictions
  } else {
    stop("`results` must be an OOF data.frame or a list containing `$oof_predictions`.")
  }
}

.enrichment_prepare <- function(results, t_star, risk_col, time_col, status_col) {
  oof <- .enrichment_oof(results)
  required_cols <- c(time_col, status_col, risk_col)
  missing_cols <- setdiff(required_cols, names(oof))
  if (length(missing_cols) > 0L) {
    stop("Missing OOF columns: ", paste(missing_cols, collapse = ", "))
  }

  risk <- oof[[risk_col]]
  time <- oof[[time_col]]
  status <- as.numeric(oof[[status_col]])
  if (!is.numeric(risk) || anyNA(risk) || any(!is.finite(risk)) || any(risk < 0 | risk > 1)) {
    stop("`", risk_col, "` must contain finite absolute risks in [0, 1].")
  }
  if (anyNA(time) || any(!is.finite(time)) || any(time < 0)) {
    stop("Follow-up times must be finite, non-negative values.")
  }
  if (anyNA(status) || any(!status %in% c(0, 1))) {
    stop("Event status must be coded 0/1 without missing values.")
  }
  if (t_star > max(time)) {
    stop(sprintf("`t_star` (%.3f) exceeds maximum follow-up (%.3f).", t_star, max(time)))
  }

  oof$.risk <- risk
  oof$.time <- time
  oof$.status <- status
  oof
}

.enrichment_km_incidence <- function(time, status, t_star) {
  fit <- survival::survfit(survival::Surv(time, status) ~ 1)
  at_horizon <- summary(fit, times = t_star, extend = TRUE)
  if (length(at_horizon$surv) == 0L || all(is.na(at_horizon$surv))) {
    return(list(obs_incidence = NA_real_, obs_lower = NA_real_, obs_upper = NA_real_))
  }
  list(
    obs_incidence = 1 - at_horizon$surv[1],
    obs_lower = 1 - at_horizon$upper[1],
    obs_upper = 1 - at_horizon$lower[1]
  )
}

.enrichment_top_label <- function(p) {
  vapply(p, function(x) {
    pct <- x * 100
    if (abs(pct - round(pct)) < 1e-8) {
      paste0("Top ", round(pct), "%")
    } else {
      paste0("Top ", signif(pct, 3), "%")
    }
  }, character(1))
}

compute_percentile_incidence_from_oof_risk <- function(
  results,
  t_star,
  risk_col,
  n_percentiles = 100,
  time_col = "time",
  status_col = "status"
) {
  if (!is.numeric(n_percentiles) || n_percentiles < 2 || n_percentiles != as.integer(n_percentiles)) {
    stop("`n_percentiles` must be an integer of at least 2.")
  }
  oof <- .enrichment_prepare(results, t_star, risk_col, time_col, status_col)
  if (n_percentiles > nrow(oof)) {
    stop("`n_percentiles` cannot exceed the number of OOF participants.")
  }

  oof$.percentile <- dplyr::ntile(oof$.risk, as.integer(n_percentiles))
  rows <- lapply(seq_len(n_percentiles), function(percentile) {
    group <- oof[oof$.percentile == percentile, , drop = FALSE]
    km <- .enrichment_km_incidence(group$.time, group$.status, t_star)
    data.frame(
      risk_percentile = percentile,
      n = nrow(group),
      n_events = sum(group$.status),
      pred_risk_mean = mean(group$.risk),
      obs_incidence = km$obs_incidence,
      obs_lower = km$obs_lower,
      obs_upper = km$obs_upper,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, rows)
  attr(out, "t_star") <- t_star
  attr(out, "risk_col") <- risk_col
  out
}

plot_percentile_incidence <- function(
  percentile_df,
  title = "Observed Incidence by Predicted Risk Percentile",
  t_star = attr(percentile_df, "t_star"),
  y_limits = NULL,
  x_limits = c(1, 100),
  point_size = 1.8,
  palette = c("SOC" = "#E64B35", "SOC+NMR" = "#4DBBD5")
) {
  plot_df <- percentile_df[!is.na(percentile_df$obs_incidence), , drop = FALSE]
  plot_df$model <- factor(plot_df$model, levels = intersect(names(palette), unique(plot_df$model)))
  y_label <- if (is.null(t_star)) {
    "Observed Incidence"
  } else {
    sprintf("Observed %d-yr Incidence", as.integer(round(t_star)))
  }

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = risk_percentile, y = obs_incidence, color = model, shape = model)
  ) +
    ggplot2::geom_point(size = point_size, alpha = 0.8) +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::scale_shape_manual(values = c("SOC" = 16, "SOC+NMR" = 17)) +
    ggplot2::scale_x_continuous(breaks = seq(0, 100, by = 10), limits = x_limits) +
    ggplot2::scale_y_continuous(
      labels = function(x) sprintf("%.0f%%", x * 100),
      limits = y_limits
    ) +
    ggplot2::labs(
      title = title,
      x = "Predicted Risk Percentile",
      y = y_label,
      color = NULL,
      shape = NULL
    ) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom"
    )
}

plot_percentile_incidence_faceted <- function(
  percentile_df,
  y_limits,
  t_star = attr(percentile_df, "t_star"),
  ncol = 3,
  x_limits = c(-5, 105),
  point_size = 1.25,
  alpha_range = c(0.05, 0.80),
  palette = c("SOC" = "#E64B35", "SOC+NMR" = "#4DBBD5")
) {
  y_label <- if (is.null(t_star)) {
    "Observed Incidence"
  } else {
    sprintf("Observed %d-yr Incidence", as.integer(round(t_star)))
  }
  disease_levels <- names(y_limits)
  plot_df <- percentile_df
  plot_df$disease <- factor(plot_df$disease, levels = disease_levels)
  plot_df$model <- factor(plot_df$model, levels = names(palette))

  y_bounds <- do.call(rbind, lapply(disease_levels, function(disease) {
    data.frame(
      disease = disease,
      risk_percentile = x_limits,
      obs_incidence = c(0, y_limits[[disease]]),
      stringsAsFactors = FALSE
    )
  }))
  y_bounds$disease <- factor(y_bounds$disease, levels = disease_levels)

  ggplot2::ggplot(
    plot_df,
    ggplot2::aes(x = risk_percentile, y = obs_incidence, color = model, group = model)
  ) +
    ggplot2::geom_blank(
      data = y_bounds,
      ggplot2::aes(x = risk_percentile, y = obs_incidence),
      inherit.aes = FALSE
    ) +
    ggplot2::geom_point(ggplot2::aes(alpha = risk_percentile), size = point_size) +
    ggplot2::facet_wrap(~disease, ncol = ncol, scales = "free") +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::scale_alpha_continuous(range = alpha_range, guide = "none") +
    ggplot2::scale_x_continuous(breaks = seq(0, 100, by = 20), limits = x_limits) +
    ggplot2::scale_y_continuous(
      labels = function(x) sprintf("%.0f%%", x * 100),
      expand = ggplot2::expansion(mult = c(0, 0.03))
    ) +
    ggplot2::coord_cartesian(xlim = x_limits, clip = "off") +
    ggplot2::labs(
      x = "Predicted Risk Percentile",
      y = y_label,
      color = NULL
    ) +
    ggplot2::theme_bw(base_size = 10) +
    ggplot2::theme(
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = 12),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position = "bottom",
      panel.spacing = grid::unit(0.8, "lines")
    )
}

.enrichment_ipcw_weights <- function(time, status, t_star) {
  censored <- survival::survfit(survival::Surv(time, 1 - status) ~ 1)
  surv_steps <- c(1, censored$surv)
  g_before <- function(x) {
    idx <- findInterval(x, censored$time, left.open = TRUE)
    pmax(surv_steps[idx + 1], .Machine$double.eps)
  }
  event_by_tstar <- time <= t_star & status == 1
  weight <- rep(0, length(time))
  if (any(event_by_tstar)) {
    weight[event_by_tstar] <- 1 / g_before(time[event_by_tstar])
  }
  list(event_by_tstar = event_by_tstar, weight = weight)
}

top_risk_incidence_from_oof_risk <- function(
  disease_models,
  t_star,
  risk_col,
  top_props = c(0.01, 0.05, 0.10),
  reference_model = "SOC",
  time_col = "time",
  status_col = "status",
  verbose = FALSE
) {
  if (is.null(names(disease_models)) || any(names(disease_models) == "")) {
    stop("`disease_models` must be a named list of diseases.")
  }
  if (!reference_model %in% unique(unlist(lapply(disease_models, names)))) {
    stop("`reference_model` was not found among the model names.")
  }

  rows <- list()
  for (disease in names(disease_models)) {
    models <- disease_models[[disease]]
    for (model in names(models)) {
      oof <- .enrichment_prepare(models[[model]], t_star, risk_col, time_col, status_col)
      ipcw <- .enrichment_ipcw_weights(oof$.time, oof$.status, t_star)
      oof$event_by_tstar <- ipcw$event_by_tstar
      oof$ipcw_event_w <- ipcw$weight

      total_events <- sum(oof$event_by_tstar)
      total_ipcw <- sum(oof$ipcw_event_w)
      oof <- oof[order(oof$.risk, decreasing = TRUE), , drop = FALSE]

      for (prop in top_props) {
        n_top <- max(1L, min(nrow(oof), as.integer(ceiling(nrow(oof) * prop))))
        top <- oof[seq_len(n_top), , drop = FALSE]
        km <- .enrichment_km_incidence(top$.time, top$.status, t_star)
        top_events <- sum(top$event_by_tstar)
        top_ipcw <- sum(top$ipcw_event_w)
        rows[[length(rows) + 1L]] <- data.frame(
          disease = disease,
          model = model,
          t_star = t_star,
          top_prop = prop,
          top_label = .enrichment_top_label(prop),
          n_total = nrow(oof),
          n_top = n_top,
          n_events_top_by_tstar = top_events,
          obs_incidence = km$obs_incidence,
          obs_lower = km$obs_lower,
          obs_upper = km$obs_upper,
          event_capture_raw = if (total_events > 0) top_events / total_events else NA_real_,
          event_capture_ipcw = if (total_ipcw > 0) top_ipcw / total_ipcw else NA_real_,
          stringsAsFactors = FALSE
        )
      }
    }
  }

  out <- dplyr::bind_rows(rows)
  reference <- out[out$model == reference_model, c("disease", "top_prop", "obs_incidence", "event_capture_raw", "event_capture_ipcw")]
  names(reference)[3:5] <- c("ref_obs_incidence", "ref_event_capture_raw", "ref_event_capture_ipcw")
  out <- dplyr::left_join(out, reference, by = c("disease", "top_prop"))
  out$delta_obs_incidence_vs_ref <- out$obs_incidence - out$ref_obs_incidence
  out$delta_event_capture_raw_vs_ref <- out$event_capture_raw - out$ref_event_capture_raw
  out$delta_event_capture_ipcw_vs_ref <- out$event_capture_ipcw - out$ref_event_capture_ipcw
  out$ref_obs_incidence <- NULL
  out$ref_event_capture_raw <- NULL
  out$ref_event_capture_ipcw <- NULL

  out$disease <- factor(out$disease, levels = names(disease_models))
  out$model <- factor(out$model, levels = unique(unlist(lapply(disease_models, names))))
  out$top_label <- factor(out$top_label, levels = .enrichment_top_label(top_props))
  attr(out, "t_star") <- t_star
  attr(out, "risk_col") <- risk_col

  if (verbose) {
    print(out[, c("disease", "model", "top_label", "n_top", "obs_incidence", "delta_obs_incidence_vs_ref")], row.names = FALSE)
  }
  out
}
