# ============================================================================
# Calibration summaries from precomputed fully OOF absolute risks
# ============================================================================
#
# Unlike compute_calibration_df() in calibration_oof_recal.r, this helper does
# not fit coxph() and does not estimate a baseline hazard. It only groups a
# supplied absolute-risk column and compares mean predicted risk with
# Kaplan-Meier observed risk at the requested horizon.
#
# The returned columns match plot_calibration_compare(), so the existing
# plotter can be reused without modification.
# ============================================================================

#' Compute grouped calibration from a precomputed OOF risk column.
#'
#' @param results Either a data.frame of OOF predictions or a results list
#'   containing `$oof_predictions`.
#' @param t_star Calibration horizon in the same units as `time_col`.
#' @param risk_col Name of the precomputed absolute-risk column.
#' @param n_groups Number of equal-size predicted-risk groups.
#' @param time_col,status_col Follow-up and event-indicator column names.
#' @param verbose Print group and error summaries.
#'
#' @return A data.frame compatible with `plot_calibration_compare()`.
compute_calibration_df_from_oof_risk <- function(
  results,
  t_star,
  risk_col,
  n_groups = 10,
  time_col = "time",
  status_col = "status",
  verbose = TRUE
) {
  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package `survival` is required for Kaplan-Meier estimates.")
  }
  if (!requireNamespace("dplyr", quietly = TRUE)) {
    stop("Package `dplyr` is required for equal-size risk grouping.")
  }

  is_results_object <- is.list(results) && !is.data.frame(results) &&
    "oof_predictions" %in% names(results)
  oof <- if (is.data.frame(results)) {
    results
  } else if (is_results_object) {
    results$oof_predictions
  } else {
    stop(
      "`results` must be an OOF data.frame or a list containing ",
      "`$oof_predictions`."
    )
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
      t_star,
      max(time)
    ))
  }
  followup_90 <- as.numeric(stats::quantile(time, 0.90))
  if (t_star > followup_90) {
    warning(sprintf(
      paste0(
        "`t_star` (%.3f) exceeds the 90th percentile of follow-up (%.3f); ",
        "some group KM estimates may be unstable."
      ),
      t_star,
      followup_90
    ))
  }

  working <- oof
  working$.precomputed_oof_risk <- risk
  working$.calibration_group <- dplyr::ntile(
    working$.precomputed_oof_risk,
    n_groups
  )

  km_formula <- stats::as.formula(
    sprintf("survival::Surv(%s, %s) ~ 1", time_col, status_col)
  )

  calibration_groups <- lapply(seq_len(n_groups), function(group) {
    group_data <- working[
      working$.calibration_group == group,
      ,
      drop = FALSE
    ]

    predicted_mean <- mean(group_data$.precomputed_oof_risk)
    predicted_lower <- as.numeric(stats::quantile(
      group_data$.precomputed_oof_risk,
      0.10
    ))
    predicted_upper <- as.numeric(stats::quantile(
      group_data$.precomputed_oof_risk,
      0.90
    ))

    km_fit <- survival::survfit(km_formula, data = group_data)
    km_at_horizon <- summary(km_fit, times = t_star, extend = TRUE)

    if (length(km_at_horizon$surv) == 0L ||
        all(is.na(km_at_horizon$surv))) {
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
  calibration$calib_error <-
    calibration$obs_risk - calibration$pred_mean
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
      t_star,
      n_groups,
      risk_col
    ))
    print(
      calibration[
        ,
        c("group", "n", "n_events", "pred_mean", "obs_risk", "calib_error")
      ],
      row.names = FALSE
    )
    cat(sprintf("Mean absolute calibration error: %.4f\n", mean_absolute_error))
    cat(sprintf("Grouped ICI: %.4f\n", grouped_ici))
    cat(sprintf("Calibration score (1 - grouped ICI): %.4f\n", calibration_score))
    cat(sprintf("Maximum absolute calibration error: %.4f\n\n",
                maximum_absolute_error))
  }

  attr(calibration, "t_star") <- t_star
  attr(calibration, "n_groups") <- n_groups
  attr(calibration, "risk_col") <- risk_col
  attr(calibration, "absolute_risk_method") <-
    "precomputed fold-specific training-derived baseline"
  attr(calibration, "n") <- nrow(oof)
  attr(calibration, "n_events") <- sum(status)
  attr(calibration, "mean_abs_calib_error") <- mean_absolute_error
  attr(calibration, "ici") <- grouped_ici
  attr(calibration, "calibration_score") <- calibration_score
  attr(calibration, "max_abs_calib_error") <- maximum_absolute_error
  attr(calibration, "ici_method") <-
    "grouped weighted mean absolute calibration error"
  attr(calibration, "calibration_score_method") <- "1 - grouped ICI"

  if (is_results_object && "c_index" %in% names(results)) {
    attr(calibration, "c_index") <- results$c_index
  }

  calibration
}
