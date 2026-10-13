# ============================================================================
# top_risk_incidence_calculator()
# ============================================================================
# Computes observed incidence and event capture in the top predicted-risk groups.
#
# Supports:
#   1) A single disease list:
#        top_risk_incidence_calculator(t2d_list, ...)
#
#   2) A named config list:
#        top_risk_incidence_calculator(risk_percentile_configs, ...)
#
#      where each disease config is:
#        list(input_list = t2d_list, y_limits = c(0, 0.8))
#
# Method:
#   - Convert OOF Cox LPs to absolute risk at t_star.
#   - Rank participants by predicted risk within each model.
#   - Select top 1%, top 5%, top 10%, etc.
#   - Estimate observed incidence using KM within each top-risk group.
#   - Compute deltas versus a reference model.
# ============================================================================

top_risk_incidence_calculator <- function(
  input,
  disease_name          = NULL,
  model_keys            = c("marker_SOC", "marker_SOC+NMR"),
  model_labels          = c("MARKER", "MARKER+NMR"),
  reference_model_key   = model_keys[1],
  top_props             = c(0.01, 0.05, 0.10),
  t_star                = 5,
  recalibrate           = FALSE,
  time_col              = "time",
  status_col            = "status",
  lp_col                = "lp",
  return_wide           = FALSE,
  verbose               = TRUE
) {

  # ---- Utility --------------------------------------------------------------
  null_if_missing <- function(x, default) {
    if (is.null(x)) default else x
  }

    make_top_label <- function(p) {
    vapply(p, function(x) {
        pct <- x * 100

        if (abs(pct - round(pct)) < 1e-8) {
        paste0("Top ", round(pct), "%")
        } else {
        paste0("Top ", signif(pct, 3), "%")
        }
    }, character(1))
    }

  get_oof_df <- function(x, model_key) {
    if (is.data.frame(x)) {
      oof <- x
    } else if (is.list(x) && is.data.frame(x$oof_predictions)) {
      oof <- x$oof_predictions
    } else {
      stop(
        "Element '", model_key, "' must be either a data.frame or a list ",
        "containing $oof_predictions."
      )
    }

    required_cols <- c(time_col, status_col, lp_col)
    missing_cols  <- setdiff(required_cols, names(oof))

    if (length(missing_cols) > 0) {
      stop(
        "Element '", model_key, "' is missing required column(s): ",
        paste(missing_cols, collapse = ", ")
      )
    }

    oof <- oof[stats::complete.cases(oof[, required_cols, drop = FALSE]), , drop = FALSE]

    if (nrow(oof) == 0) {
      stop("Element '", model_key, "' has zero complete rows after filtering.")
    }

    oof
  }

  # ---- Detect input type ----------------------------------------------------
  is_single_disease <- all(model_keys %in% names(input))

  is_config_list <- is.list(input) &&
    !is.null(names(input)) &&
    all(vapply(input, function(x) is.list(x) && !is.null(x$input_list), logical(1)))

  if (is_single_disease) {
    disease_name <- null_if_missing(disease_name, "Disease")

    disease_configs <- list()
    disease_configs[[disease_name]] <- list(input_list = input)

  } else if (is_config_list) {
    disease_configs <- input

  } else {
    stop(
      "input must be either a single disease list containing model_keys, ",
      "or a named config list where each element has $input_list."
    )
  }

  # ---- Validate model/reference setup --------------------------------------
  if (length(model_keys) != length(model_labels)) {
    stop("model_keys and model_labels must have the same length.")
  }

  if (!reference_model_key %in% model_keys) {
    stop("reference_model_key must be one of model_keys.")
  }

  if (any(top_props <= 0) || any(top_props > 1)) {
    stop("top_props must be probabilities in (0, 1].")
  }

  reference_model_label <- model_labels[match(reference_model_key, model_keys)]

  # ---- LP -> absolute risk --------------------------------------------------
  add_predicted_risk <- function(oof, model_key) {

    max_time   <- max(oof[[time_col]], na.rm = TRUE)
    pct90_time <- as.numeric(stats::quantile(oof[[time_col]], 0.90, na.rm = TRUE))

    if (t_star > max_time) {
      stop(sprintf(
        "For model '%s', t_star (%.2f) exceeds maximum follow-up time (%.2f).",
        model_key, t_star, max_time
      ))
    }

    if (t_star > pct90_time) {
      warning(sprintf(
        "For model '%s', t_star (%.2f) exceeds 90th percentile follow-up (%.2f). KM/IPCW estimates may be unstable.",
        model_key, t_star, pct90_time
      ))
    }

    if (isTRUE(recalibrate)) {
      calib_formula <- stats::as.formula(
        sprintf("survival::Surv(%s, %s) ~ %s", time_col, status_col, lp_col)
      )

      calib_fit <- survival::coxph(
        calib_formula,
        data = oof,
        x = FALSE,
        y = FALSE
      )

      calibration_slope <- as.numeric(stats::coef(calib_fit)[1])

    } else {
      calib_formula <- stats::as.formula(
        sprintf("survival::Surv(%s, %s) ~ offset(%s)", time_col, status_col, lp_col)
      )

      calib_fit <- survival::coxph(
        calib_formula,
        data = oof,
        x = FALSE,
        y = FALSE
      )

      calibration_slope <- 1
    }

    bh <- survival::basehaz(calib_fit, centered = FALSE)

    H0_tstar <- stats::approx(
      x      = bh$time,
      y      = bh$hazard,
      xout   = t_star,
      method = "constant",
      rule   = 2,
      f      = 0
    )$y

    oof$lp_eff    <- calibration_slope * oof[[lp_col]]
    oof$pred_surv <- exp(-H0_tstar * exp(oof$lp_eff))
    oof$pred_risk <- 1 - oof$pred_surv

    list(
      oof               = oof,
      calibration_slope = calibration_slope,
      H0_tstar          = H0_tstar,
      baseline_risk     = 1 - exp(-H0_tstar)
    )
  }

  # ---- IPCW event capture support ------------------------------------------
  add_ipcw_event_weights <- function(oof) {

    cens_formula <- stats::as.formula(
      sprintf("survival::Surv(%s, 1 - %s) ~ 1", time_col, status_col)
    )

    cens_km <- survival::survfit(cens_formula, data = oof)

    surv_steps <- c(1, cens_km$surv)
    time_steps <- cens_km$time

    G_before <- function(x) {
      idx <- findInterval(x, time_steps, left.open = TRUE)
      out <- surv_steps[idx + 1]
      pmax(out, .Machine$double.eps)
    }

    event_by_tstar <- oof[[time_col]] <= t_star & oof[[status_col]] == 1

    oof$event_by_tstar <- event_by_tstar
    oof$ipcw_event_w   <- 0

    if (any(event_by_tstar)) {
      oof$ipcw_event_w[event_by_tstar] <- 1 / G_before(oof[[time_col]][event_by_tstar])
    }

    oof
  }

  # ---- Observed incidence in selected group --------------------------------
  km_incidence <- function(df) {
    km_formula <- stats::as.formula(
      sprintf("survival::Surv(%s, %s) ~ 1", time_col, status_col)
    )

    km  <- survival::survfit(km_formula, data = df)
    kms <- summary(km, times = t_star, extend = TRUE)

    if (length(kms$surv) == 0 || all(is.na(kms$surv))) {
      return(list(
        obs_incidence = NA_real_,
        obs_lower     = NA_real_,
        obs_upper     = NA_real_
      ))
    }

    list(
      obs_incidence = 1 - kms$surv[1],
      obs_lower     = 1 - kms$upper[1],
      obs_upper     = 1 - kms$lower[1]
    )
  }

  # ---- Main loop ------------------------------------------------------------
  if (verbose) {
    cat("\n============================================\n")
    cat("Computing Top-Risk Incidence / Event Capture\n")
    cat("============================================\n")
    cat(sprintf("Time horizon (t_star)  : %.2f\n", t_star))
    cat(sprintf("Top groups             : %s\n", paste(make_top_label(top_props), collapse = ", ")))
    cat(sprintf("Models                 : %s\n", paste(model_labels, collapse = " vs ")))
    cat(sprintf("Reference model        : %s\n", reference_model_label))
    cat(sprintf("Recalibration          : %s\n",
                ifelse(recalibrate,
                       "TRUE  (estimate slope + baseline)",
                       "FALSE (baseline only; slope fixed at 1)")))
  }

  all_rows <- list()

  row_idx <- 1

  for (disease in names(disease_configs)) {

    input_list <- disease_configs[[disease]]$input_list

    missing_keys <- setdiff(model_keys, names(input_list))
    if (length(missing_keys) > 0) {
      stop(
        "For disease '", disease, "', missing model_keys: ",
        paste(missing_keys, collapse = ", ")
      )
    }

    if (verbose) {
      cat("\n", disease, "\n", sep = "")
    }

    for (m in seq_along(model_keys)) {

      model_key   <- model_keys[m]
      model_label <- model_labels[m]

      oof_raw  <- get_oof_df(input_list[[model_key]], model_key)
      risk_obj <- add_predicted_risk(oof_raw, model_key)

      oof <- risk_obj$oof
      oof <- add_ipcw_event_weights(oof)

      n_total <- nrow(oof)

      total_events_by_tstar_raw <- sum(oof$event_by_tstar, na.rm = TRUE)
      total_events_by_tstar_ipcw <- sum(oof$ipcw_event_w, na.rm = TRUE)

      # Rank descending by predicted risk.
      oof <- oof[order(oof$pred_risk, decreasing = TRUE), , drop = FALSE]
      oof$risk_rank <- seq_len(nrow(oof))

      for (p in top_props) {

        n_top <- ceiling(n_total * p)
        n_top <- max(1, min(n_top, n_total))

        top_df <- oof[seq_len(n_top), , drop = FALSE]

        km_obj <- km_incidence(top_df)

        top_events_raw <- sum(top_df$event_by_tstar, na.rm = TRUE)
        top_events_ipcw <- sum(top_df$ipcw_event_w, na.rm = TRUE)

        event_capture_raw <- ifelse(
          total_events_by_tstar_raw > 0,
          top_events_raw / total_events_by_tstar_raw,
          NA_real_
        )

        event_capture_ipcw <- ifelse(
          total_events_by_tstar_ipcw > 0,
          top_events_ipcw / total_events_by_tstar_ipcw,
          NA_real_
        )

        all_rows[[row_idx]] <- data.frame(
          disease                    = disease,
          model_key                  = model_key,
          model                      = model_label,
          reference_model_key        = reference_model_key,
          reference_model            = reference_model_label,
          t_star                     = t_star,
          top_prop                   = p,
          top_label                  = make_top_label(p),
          n_total                    = n_total,
          n_top                      = n_top,
          n_events_total_by_tstar    = total_events_by_tstar_raw,
          n_events_top_by_tstar      = top_events_raw,
          event_capture_raw          = event_capture_raw,
          event_capture_ipcw         = event_capture_ipcw,
          pred_risk_mean_top         = mean(top_df$pred_risk, na.rm = TRUE),
          pred_risk_median_top       = stats::median(top_df$pred_risk, na.rm = TRUE),
          pred_risk_min_top          = min(top_df$pred_risk, na.rm = TRUE),
          pred_risk_max_top          = max(top_df$pred_risk, na.rm = TRUE),
          obs_incidence              = km_obj$obs_incidence,
          obs_lower                  = km_obj$obs_lower,
          obs_upper                  = km_obj$obs_upper,
          calibration_slope          = risk_obj$calibration_slope,
          baseline_risk_tstar        = risk_obj$baseline_risk,
          recalibrate                = recalibrate,
          stringsAsFactors           = FALSE
        )

        row_idx <- row_idx + 1
      }

      if (verbose) {
        cat(sprintf(
          "  %-18s done | N = %d | events by t* = %d\n",
          model_label,
          n_total,
          total_events_by_tstar_raw
        ))
      }
    }
  }

  out <- dplyr::bind_rows(all_rows)

  out <- out %>%
    dplyr::mutate(
      disease   = factor(disease, levels = names(disease_configs)),
      model     = factor(model, levels = model_labels),
      top_label = factor(top_label, levels = make_top_label(top_props))
    )

  # ---- Add deltas versus reference model -----------------------------------
  ref_df <- out %>%
    dplyr::filter(model_key == reference_model_key) %>%
    dplyr::select(
      disease,
      top_prop,
      ref_obs_incidence        = obs_incidence,
      ref_event_capture_raw    = event_capture_raw,
      ref_event_capture_ipcw   = event_capture_ipcw,
      ref_pred_risk_mean_top   = pred_risk_mean_top,
      ref_n_events_top_by_tstar = n_events_top_by_tstar
    )

  out <- out %>%
    dplyr::left_join(ref_df, by = c("disease", "top_prop")) %>%
    dplyr::mutate(
      delta_obs_incidence_vs_ref =
        obs_incidence - ref_obs_incidence,

      delta_event_capture_raw_vs_ref =
        event_capture_raw - ref_event_capture_raw,

      delta_event_capture_ipcw_vs_ref =
        event_capture_ipcw - ref_event_capture_ipcw,

      delta_pred_risk_mean_top_vs_ref =
        pred_risk_mean_top - ref_pred_risk_mean_top,

      delta_n_events_top_by_tstar_vs_ref =
        n_events_top_by_tstar - ref_n_events_top_by_tstar
    )

  # ---- Optional wide output -------------------------------------------------
  if (isTRUE(return_wide)) {
    out_wide <- out %>%
      dplyr::select(
        disease,
        model_key,
        model,
        top_label,
        n_top,
        n_events_top_by_tstar,
        event_capture_raw,
        event_capture_ipcw,
        obs_incidence,
        obs_lower,
        obs_upper,
        delta_obs_incidence_vs_ref,
        delta_event_capture_raw_vs_ref,
        delta_event_capture_ipcw_vs_ref,
        delta_n_events_top_by_tstar_vs_ref
      ) %>%
      tidyr::pivot_wider(
        names_from = top_label,
        values_from = c(
          n_top,
          n_events_top_by_tstar,
          event_capture_raw,
          event_capture_ipcw,
          obs_incidence,
          obs_lower,
          obs_upper,
          delta_obs_incidence_vs_ref,
          delta_event_capture_raw_vs_ref,
          delta_event_capture_ipcw_vs_ref,
          delta_n_events_top_by_tstar_vs_ref
        )
      )

    attr(out_wide, "long") <- out

    if (verbose) {
      cat("============================================\n\n")
    }

    return(out_wide)
  }

  if (verbose) {
    cat("============================================\n\n")
  }

  return(out)
}