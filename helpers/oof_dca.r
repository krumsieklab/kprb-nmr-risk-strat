# Decision curves from precomputed fold-specific absolute risks.
# Uses pred_risk_{t}yr. Does not refit a Cox model or convert linear predictors.

.dca_oof <- function(results) {
  if (is.data.frame(results)) return(results)
  if (is.list(results) && "oof_predictions" %in% names(results)) {
    return(results$oof_predictions)
  }
  stop("`results` must be an OOF data.frame or a list containing `$oof_predictions`.")
}

# IPCW net benefit per 1,000 patients across probability thresholds.
compute_dca_df_from_oof_risk <- function(
  results,
  t_star,
  risk_col        = NULL,
  threshold_range = seq(0.001, 0.50, by = 0.001),
  time_col        = "time",
  status_col      = "status",
  verbose         = FALSE
) {
  oof <- .dca_oof(results)
  if (is.null(risk_col)) risk_col <- sprintf("pred_risk_%dyr", as.integer(round(t_star)))

  missing_cols <- setdiff(c(time_col, status_col, risk_col), names(oof))
  if (length(missing_cols) > 0) {
    stop("Missing columns in oof_predictions: ", paste(missing_cols, collapse = ", "))
  }
  if (any(threshold_range <= 0) || any(threshold_range >= 1)) {
    stop("threshold_range values must be strictly between 0 and 1.")
  }

  oof$pred_risk <- as.numeric(oof[[risk_col]])
  if (anyNA(oof$pred_risk) || any(!is.finite(oof$pred_risk)) ||
      any(oof$pred_risk < 0 | oof$pred_risk > 1)) {
    stop("`", risk_col, "` must contain finite absolute risks in [0, 1].")
  }

  max_time   <- max(oof[[time_col]], na.rm = TRUE)
  pct90_time <- as.numeric(stats::quantile(oof[[time_col]], 0.90, na.rm = TRUE))
  if (t_star > max_time) {
    stop(sprintf("t_star (%.2f) exceeds maximum observed follow-up (%.2f).", t_star, max_time))
  }
  if (t_star > pct90_time) {
    warning(sprintf(
      "t_star (%.2f) exceeds the 90th percentile of follow-up (%.2f). IPCW weights may be unstable.",
      t_star, pct90_time
    ))
  }

  time   <- oof[[time_col]]
  status <- oof[[status_col]]
  cens_km <- survival::survfit(survival::Surv(time, 1 - status) ~ 1)
  G_func  <- stats::stepfun(cens_km$time, c(1, cens_km$surv), right = TRUE)
  G_tstar <- G_func(t_star)
  if (G_tstar <= 0) {
    stop(sprintf("Censoring survival G(t_star = %.2f) = 0.", t_star))
  }

  event_by_tstar        <- oof[[time_col]] <= t_star & oof[[status_col]] == 1
  observed_past_tstar   <- oof[[time_col]] > t_star
  censored_before_tstar <- oof[[time_col]] <= t_star & oof[[status_col]] == 0

  oof$Y           <- as.integer(event_by_tstar)
  oof$informative <- !censored_before_tstar
  oof$ipcw_w      <- NA_real_
  if (any(event_by_tstar)) {
    g_ev <- G_func(oof[[time_col]][event_by_tstar])
    g_ev[g_ev <= 0] <- NA
    oof$ipcw_w[event_by_tstar] <- 1 / g_ev
  }
  oof$ipcw_w[observed_past_tstar] <- 1 / G_tstar

  n_total   <- nrow(oof)
  inf_mask  <- oof$informative & !is.na(oof$ipcw_w)
  w_all     <- oof$ipcw_w[inf_mask]
  Y_all     <- oof$Y[inf_mask]
  sum_wY    <- sum(w_all * Y_all, na.rm = TRUE)
  sum_w1mY  <- sum(w_all * (1 - Y_all), na.rm = TRUE)

  nb_treat_all  <- sum_wY / n_total -
    (threshold_range / (1 - threshold_range)) * (sum_w1mY / n_total)
  nb_treat_none <- rep(0, length(threshold_range))

  pos_base <- inf_mask & !is.na(oof$pred_risk)
  nb_model <- vapply(threshold_range, function(pt) {
    test_pos <- pos_base & (oof$pred_risk >= pt)
    if (!any(test_pos)) return(0)
    w_pos  <- oof$ipcw_w[test_pos]
    Y_pos  <- oof$Y[test_pos]
    n_TP_w <- sum(w_pos * Y_pos, na.rm = TRUE)
    n_FP_w <- sum(w_pos * (1 - Y_pos), na.rm = TRUE)
    (n_TP_w / n_total) - (pt / (1 - pt)) * (n_FP_w / n_total)
  }, numeric(1))

  dca_df <- data.frame(
    threshold     = threshold_range,
    nb_model      = nb_model * 1000,
    nb_treat_all  = nb_treat_all * 1000,
    nb_treat_none = nb_treat_none * 1000,
    stringsAsFactors = FALSE
  )
  attr(dca_df, "t_star")       <- t_star
  attr(dca_df, "risk_col")     <- risk_col
  attr(dca_df, "n")            <- n_total
  attr(dca_df, "n_events")     <- sum(oof[[status_col]], na.rm = TRUE)
  attr(dca_df, "overall_risk") <- sum_wY / n_total

  if (isTRUE(verbose)) {
    cat(sprintf(
      "t* = %g | %s | N = %d | IPCW risk = %.3f | max NB = %.2f\n",
      t_star, risk_col, n_total, sum_wY / n_total, max(dca_df$nb_model)
    ))
  }
  dca_df
}

# Net benefit (per 1,000) at one threshold. Uses the nearest computed threshold.
nb_at_threshold <- function(dca_df, threshold) {
  idx <- which.min(abs(dca_df$threshold - threshold))
  dca_df$nb_model[[idx]]
}

# Named list of model curves (SOC, SOC+NMR) plus treat-all / treat-none.
plot_dca_compare <- function(
  dca_list,
  title          = "Decision Curve Analysis",
  x_limits       = NULL,
  show_treat_all = TRUE,
  show_treat_none = TRUE,
  annotate_refs  = TRUE,
  palette        = c("SOC" = "#E64B35", "SOC+NMR" = "#4DBBD5"),
  t_star         = NULL,
  line_size      = 0.8
) {
  if (is.data.frame(dca_list)) {
    dca_list <- list(Model = dca_list)
  }
  if (is.null(names(dca_list)) || any(names(dca_list) == "")) {
    stop("dca_list must be a named list of DCA data.frames.")
  }
  if (is.null(t_star)) t_star <- attr(dca_list[[1]], "t_star")
  t_label <- if (!is.null(t_star)) sprintf("%d-yr", as.integer(round(t_star))) else "t*"

  model_df <- dplyr::bind_rows(lapply(names(dca_list), function(lab) {
    df <- dca_list[[lab]]
    data.frame(threshold = df$threshold, nb = df$nb_model, model = lab, stringsAsFactors = FALSE)
  }))
  model_df$model <- factor(model_df$model, levels = names(dca_list))
  ref_df <- dca_list[[1]]

  if (!is.null(x_limits)) {
    ref_df   <- ref_df[ref_df$threshold >= x_limits[1] & ref_df$threshold <= x_limits[2], , drop = FALSE]
    model_df <- model_df[model_df$threshold >= x_limits[1] & model_df$threshold <= x_limits[2], , drop = FALSE]
  }

  ax_min_x <- min(ref_df$threshold)
  ax_max_x <- max(ref_df$threshold)
  all_nb   <- c(model_df$nb, if (show_treat_all) ref_df$nb_treat_all else NULL)
  nb_upper <- max(all_nb, na.rm = TRUE) * 1.08
  nb_lower <- max(-5, min(all_nb, na.rm = TRUE) - 1)

  p <- ggplot2::ggplot()
  if (show_treat_none) {
    p <- p + ggplot2::geom_hline(
      yintercept = 0, linetype = "dotted", color = "gray50", linewidth = line_size * 0.75
    )
  }
  if (show_treat_all) {
    p <- p + ggplot2::geom_line(
      data = ref_df,
      ggplot2::aes(x = threshold, y = nb_treat_all),
      linetype = "dashed", color = "gray50", linewidth = line_size * 0.75
    )
  }
  p <- p + ggplot2::geom_line(
    data = model_df,
    ggplot2::aes(x = threshold, y = nb, color = model, linetype = model),
    linewidth = line_size
  )

  if (annotate_refs) {
    annot_x <- ax_min_x + (ax_max_x - ax_min_x) * 0.02
    if (show_treat_all) {
      nb_all <- stats::approx(ref_df$threshold, ref_df$nb_treat_all, xout = annot_x, rule = 2)$y
      if (!is.na(nb_all) && nb_all >= nb_lower) {
        p <- p + ggplot2::annotate(
          "text", x = annot_x, y = nb_all, label = "Treat all",
          hjust = 0, vjust = -0.5, color = "gray50", size = 3.2, fontface = "italic"
        )
      }
    }
    if (show_treat_none && nb_lower <= 0) {
      p <- p + ggplot2::annotate(
        "text", x = annot_x, y = 0, label = "Treat none",
        hjust = 0, vjust = -0.5, color = "gray50", size = 3.2, fontface = "italic"
      )
    }
  }

  p +
    ggplot2::scale_color_manual(values = palette) +
    ggplot2::scale_linetype_manual(values = stats::setNames(rep("solid", nlevels(model_df$model)), levels(model_df$model))) +
    ggplot2::scale_x_continuous(
      limits = c(ax_min_x, ax_max_x),
      labels = function(x) sprintf("%.0f%%", x * 100),
      expand = c(0.01, 0.01)
    ) +
    ggplot2::scale_y_continuous(expand = c(0.02, 0.02)) +
    ggplot2::coord_cartesian(ylim = c(nb_lower, nb_upper)) +
    ggplot2::labs(
      title = title,
      x = sprintf("Threshold Probability (%s)", t_label),
      y = "Net Benefit (per 1000 patients)",
      color = "Model",
      linetype = "Model"
    ) +
    ggplot2::theme_bw(base_size = 13) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(face = "bold", hjust = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      legend.position  = "bottom"
    )
}
