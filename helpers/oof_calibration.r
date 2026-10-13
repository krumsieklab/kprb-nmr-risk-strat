# Grouped calibration from precomputed out-of-fold absolute risks.
#
# compute_calibration_df_from_oof_risk() bins a supplied risk column and
# compares mean predicted risk with Kaplan-Meier observed risk.
# plot_calibration_compare() overlays those group summaries.

compute_calibration_df_from_oof_risk <- function(
  results,
  t_star,
  risk_col,
  n_groups = 10,
  time_col = "time",
  status_col = "status",
  verbose = TRUE
) {
  is_results_object <- is.list(results) && !is.data.frame(results) &&
    "oof_predictions" %in% names(results)
  oof <- if (is.data.frame(results)) {
    results
  } else if (is_results_object) {
    results$oof_predictions
  } else {
    stop("`results` must be an OOF data.frame or a list containing `$oof_predictions`.")
  }

  required_cols <- c(time_col, status_col, risk_col)
  missing_cols <- setdiff(required_cols, names(oof))
  if (length(missing_cols) > 0L) {
    stop("Missing OOF columns: ", paste(missing_cols, collapse = ", "))
  }
  if (!is.numeric(t_star) || length(t_star) != 1L || is.na(t_star) ||
      !is.finite(t_star) || t_star <= 0) {
    stop("`t_star` must be one finite, positive number.")
  }
  if (!is.numeric(n_groups) || length(n_groups) != 1L ||
      is.na(n_groups) || n_groups < 2 || n_groups != as.integer(n_groups)) {
    stop("`n_groups` must be one integer of at least 2.")
  }
  n_groups <- as.integer(n_groups)
  if (n_groups > nrow(oof)) {
    stop("`n_groups` cannot exceed the number of OOF participants.")
  }

  risk <- oof[[risk_col]]
  if (!is.numeric(risk) || anyNA(risk) || any(!is.finite(risk)) ||
      any(risk < 0 | risk > 1)) {
    stop("`", risk_col, "` must contain finite absolute risks in [0, 1].")
  }

  time <- oof[[time_col]]
  status <- as.numeric(oof[[status_col]])
  if (anyNA(time) || any(!is.finite(time)) || any(time < 0)) {
    stop("Follow-up times must be finite, non-negative values.")
  }
  if (anyNA(status) || any(!status %in% c(0, 1))) {
    stop("Event status must be coded 0/1 without missing values.")
  }
  if (t_star > max(time)) {
    stop(sprintf(
      "`t_star` (%.3f) exceeds maximum OOF follow-up (%.3f).",
      t_star, max(time)
    ))
  }
  followup_90 <- as.numeric(stats::quantile(time, 0.90))
  if (t_star > followup_90) {
    warning(sprintf(
      "`t_star` (%.3f) exceeds the 90th percentile of follow-up (%.3f); some group KM estimates may be unstable.",
      t_star, followup_90
    ))
  }

  working <- oof
  working$.precomputed_oof_risk <- risk
  working$.calibration_group <- dplyr::ntile(working$.precomputed_oof_risk, n_groups)

  km_formula <- stats::as.formula(
    sprintf("survival::Surv(%s, %s) ~ 1", time_col, status_col)
  )

  calibration_groups <- lapply(seq_len(n_groups), function(group) {
    group_data <- working[working$.calibration_group == group, , drop = FALSE]

    predicted_mean <- mean(group_data$.precomputed_oof_risk)
    predicted_lower <- as.numeric(stats::quantile(group_data$.precomputed_oof_risk, 0.10))
    predicted_upper <- as.numeric(stats::quantile(group_data$.precomputed_oof_risk, 0.90))

    km_fit <- survival::survfit(km_formula, data = group_data)
    km_at_horizon <- summary(km_fit, times = t_star, extend = TRUE)

    if (length(km_at_horizon$surv) == 0L || all(is.na(km_at_horizon$surv))) {
      observed_risk <- NA_real_
      observed_lower <- NA_real_
      observed_upper <- NA_real_
    } else {
      observed_risk <- 1 - km_at_horizon$surv[1]
      observed_lower <- 1 - km_at_horizon$upper[1]
      observed_upper <- 1 - km_at_horizon$lower[1]
    }

    data.frame(
      group = group,
      n = nrow(group_data),
      n_events = sum(group_data[[status_col]]),
      pred_mean = predicted_mean,
      pred_lower = predicted_lower,
      pred_upper = predicted_upper,
      obs_risk = observed_risk,
      obs_lower = observed_lower,
      obs_upper = observed_upper,
      stringsAsFactors = FALSE
    )
  })

  calibration <- do.call(rbind, calibration_groups)
  calibration$calib_error <- calibration$obs_risk - calibration$pred_mean
  calibration$abs_calib_error <- abs(calibration$calib_error)

  valid <- !is.na(calibration$abs_calib_error)
  if (!any(valid)) {
    mean_absolute_error <- NA_real_
    grouped_ici <- NA_real_
    calibration_score <- NA_real_
    maximum_absolute_error <- NA_real_
  } else {
    mean_absolute_error <- mean(calibration$abs_calib_error[valid])
    grouped_ici <- stats::weighted.mean(
      calibration$abs_calib_error[valid],
      calibration$n[valid]
    )
    calibration_score <- 1 - grouped_ici
    maximum_absolute_error <- max(calibration$abs_calib_error[valid])
  }

  if (verbose) {
    cat("\nCalibration from precomputed fully OOF absolute risks\n")
    cat(sprintf(
      "Horizon: %.3f | groups: %d | risk column: %s\n",
      t_star, n_groups, risk_col
    ))
    print(
      calibration[, c("group", "n", "n_events", "pred_mean", "obs_risk", "calib_error")],
      row.names = FALSE
    )
    cat(sprintf("Mean absolute calibration error: %.4f\n", mean_absolute_error))
    cat(sprintf("Grouped ICI: %.4f\n", grouped_ici))
    cat(sprintf("Calibration score (1 - grouped ICI): %.4f\n", calibration_score))
    cat(sprintf("Maximum absolute calibration error: %.4f\n\n", maximum_absolute_error))
  }

  attr(calibration, "t_star") <- t_star
  attr(calibration, "n_groups") <- n_groups
  attr(calibration, "risk_col") <- risk_col
  attr(calibration, "n") <- nrow(oof)
  attr(calibration, "n_events") <- sum(status)
  attr(calibration, "mean_abs_calib_error") <- mean_absolute_error
  attr(calibration, "ici") <- grouped_ici
  attr(calibration, "calibration_score") <- calibration_score
  attr(calibration, "max_abs_calib_error") <- maximum_absolute_error

  if (is_results_object && "c_index" %in% names(results)) {
    attr(calibration, "c_index") <- results$c_index
  }

  calibration
}


plot_calibration_compare <- function(
  calib_list,
  title = "Calibration Plot Comparison",
  t_star = NULL,
  axis_limits = NULL,
  show_error_bars = TRUE,
  show_lines = TRUE,
  label_groups = FALSE,
  label_top_n = NULL,
  palette = NULL
) {
  if (!is.list(calib_list) || is.data.frame(calib_list)) {
    stop("`calib_list` must be a named list of calibration data frames.")
  }
  if (is.null(names(calib_list)) || any(names(calib_list) == "")) {
    stop("`calib_list` must be named, e.g. list(SOC = df1, `SOC+NMR` = df2).")
  }

  calib_both <- dplyr::bind_rows(lapply(seq_along(calib_list), function(i) {
    df <- calib_list[[i]]
    df$model <- names(calib_list)[i]
    df
  }))

  required_cols <- c("group", "pred_mean", "obs_risk", "model")
  missing_cols <- setdiff(required_cols, names(calib_both))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (is.null(t_star)) {
    t_star <- attr(calib_list[[1]], "t_star")
  }
  t_label <- if (!is.null(t_star)) sprintf("%d-yr", as.integer(round(t_star))) else "t*"

  plot_df <- calib_both[!is.na(calib_both$obs_risk), , drop = FALSE]
  plot_df$model <- factor(plot_df$model, levels = unique(plot_df$model))
  model_levels <- levels(plot_df$model)
  n_models <- length(model_levels)

  if (is.null(axis_limits)) {
    all_vals <- c(
      plot_df$pred_mean,
      plot_df$obs_risk,
      if (show_error_bars) c(plot_df$obs_lower, plot_df$obs_upper) else NULL
    )
    ax_min <- max(0, floor(min(all_vals, na.rm = TRUE) * 20) / 20 - 0.02)
    ax_max <- min(1, ceiling(max(all_vals, na.rm = TRUE) * 20) / 20 + 0.02)
  } else {
    ax_min <- axis_limits[1]
    ax_max <- axis_limits[2]
  }

  if (is.null(palette)) {
    palette <- scales::hue_pal()(n_models)
    names(palette) <- model_levels
  }
  shape_values <- setNames(c(16, 17, 15, 18, 8, 3, 7, 9)[seq_len(n_models)], model_levels)

  p <- ggplot2::ggplot(
    plot_df,
    ggplot2::aes(
      x = pred_mean,
      y = obs_risk,
      color = model,
      shape = model,
      group = model
    )
  ) +
    ggplot2::geom_abline(
      intercept = 0, slope = 1,
      linetype = "dashed", color = "gray40", linewidth = 0.7
    )

  if (isTRUE(show_lines)) {
    p <- p + ggplot2::geom_line(linewidth = 0.7, alpha = 0.55)
  }
  p <- p + ggplot2::geom_point(size = 3, alpha = 0.95)

  if (show_error_bars && all(c("obs_lower", "obs_upper") %in% names(plot_df))) {
    bar_width <- (ax_max - ax_min) * 0.015
    p <- p + ggplot2::geom_errorbar(
      ggplot2::aes(ymin = obs_lower, ymax = obs_upper),
      width = bar_width, linewidth = 0.56, alpha = 0.55
    )
  }

  if (isTRUE(label_groups)) {
    label_df <- plot_df
    if (!is.null(label_top_n)) {
      max_group <- max(plot_df$group, na.rm = TRUE)
      label_df <- plot_df[plot_df$group >= (max_group - label_top_n + 1), , drop = FALSE]
    }
    p <- p + ggplot2::geom_text(
      data = label_df,
      ggplot2::aes(label = group),
      show.legend = FALSE, size = 3, vjust = -0.7
    )
  }

  p +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::scale_shape_manual(values = shape_values) +
    ggplot2::scale_x_continuous(
      limits = c(ax_min, ax_max),
      labels = function(x) sprintf("%.0f%%", x * 100)
    ) +
    ggplot2::scale_y_continuous(
      limits = c(ax_min, ax_max),
      labels = function(x) sprintf("%.0f%%", x * 100)
    ) +
    ggplot2::coord_fixed(ratio = 1, xlim = c(ax_min, ax_max), ylim = c(ax_min, ax_max)) +
    ggplot2::labs(
      title = title,
      x = sprintf("Predicted %s Risk", t_label),
      y = sprintf("Observed %s Risk (KM)", t_label),
      color = "Model",
      shape = "Model"
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      aspect.ratio = 1,
      legend.position = "bottom"
    )
}
