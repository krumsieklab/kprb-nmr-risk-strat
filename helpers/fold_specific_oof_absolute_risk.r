# ============================================================================
# Fold-specific training-derived baseline hazards and fully OOF absolute risks
# ============================================================================
#
# This file intentionally does not replace or modify cv_glmnet_cox_oof().
# It post-processes an existing glmnet OOF result using:
#   * the saved outer-fold assignments;
#   * the saved fold-specific cv.glmnet fits and selected lambdas; and
#   * the original model data plus the exact vars_to_scale vector.
#
# Existing RDS files do not store raw outer-training rows or fold-specific
# scaling parameters. Therefore, an RDS alone is insufficient: `df` and
# `vars_to_scale` must describe the exact data/preprocessing used for the fit.
# Reproduced test LPs are checked against the saved OOF LPs before risks are
# calculated. A mismatch fails rather than silently producing invalid risks.
#
# glmnet 4.1.8 uses the Breslow approximation for Cox ties. The baseline below
# is estimated explicitly from each outer-training fold:
#
#   dH0(t_j) = d_j / sum_{i in R(t_j)} exp(lp_i)
#
# This avoids coxph offset-centering ambiguity and ensures that held-out
# outcomes do not contribute to the baseline used for their predictions.
# ============================================================================

.fold_oof_logspace_add <- function(log_x, log_y) {
  if (is.infinite(log_x) && log_x < 0) return(log_y)
  if (is.infinite(log_y) && log_y < 0) return(log_x)

  max_log <- max(log_x, log_y)
  max_log + log(exp(log_x - max_log) + exp(log_y - max_log))
}

.fold_oof_logsumexp <- function(x) {
  max_x <- max(x)
  max_x + log(sum(exp(x - max_x)))
}

.fold_oof_risk_column <- function(prediction_time) {
  if (prediction_time == round(prediction_time)) {
    label <- format(round(prediction_time), trim = TRUE, scientific = FALSE)
  } else {
    label <- format(
      prediction_time,
      trim = TRUE,
      scientific = FALSE,
      digits = 15
    )
    label <- sub("0+$", "", label)
    label <- sub("\\.$", "", label)
  }
  label <- gsub("\\.", "p", label, fixed = FALSE)
  paste0("pred_risk_", label, "yr")
}

.fold_oof_compute_scale_params <- function(data, vars_to_scale) {
  params <- vector("list", length(vars_to_scale))
  names(params) <- vars_to_scale

  for (var in vars_to_scale) {
    params[[var]] <- list(
      center = mean(data[[var]], na.rm = TRUE),
      scale = stats::sd(data[[var]], na.rm = TRUE)
    )
  }

  params
}

.fold_oof_apply_scaling <- function(data, scale_params) {
  scaled <- data

  for (var in names(scale_params)) {
    center <- scale_params[[var]]$center
    scale <- scale_params[[var]]$scale

    # Mirrors apply_scaling() in cv_glmnet_oof.r
    if (is.finite(scale) && !is.na(scale) && scale > 0) {
      scaled[[var]] <- (scaled[[var]] - center) / scale
    }
  }

  scaled
}

.fold_oof_make_x_matrix <- function(data, predictors) {
  xdf <- data[, predictors, drop = FALSE]

  # Mirrors make_x_matrix() in cv_glmnet_oof.r
  for (name in names(xdf)) {
    if (is.factor(xdf[[name]])) {
      xdf[[name]] <- as.numeric(as.character(xdf[[name]]))
    } else if (is.logical(xdf[[name]])) {
      xdf[[name]] <- as.integer(xdf[[name]])
    }
  }

  x <- data.matrix(xdf)
  storage.mode(x) <- "double"
  colnames(x) <- predictors

  if (anyNA(x)) {
    stop("Predictor matrix contains NA after conversion.")
  }

  x
}

#' Estimate an LP=0 Breslow cumulative baseline hazard.
#'
#' @param time Numeric follow-up times.
#' @param status Event indicators coded 0/1.
#' @param lp Numeric Cox linear predictors from one fitted model.
#'
#' @return A data.frame with one row per event time and columns
#'   `event_time`, `n_events`, `risk_set_size`, `hazard_increment`,
#'   `log_hazard_increment`, `baseline_cumhaz`, and
#'   `log_baseline_cumhaz`.
estimate_breslow_baseline_from_lp <- function(time, status, lp) {
  if (!is.numeric(time) || !is.numeric(lp)) {
    stop("`time` and `lp` must be numeric.")
  }
  if (length(time) != length(status) || length(time) != length(lp)) {
    stop("`time`, `status`, and `lp` must have equal lengths.")
  }
  if (length(time) == 0L) {
    stop("Cannot estimate a baseline hazard from zero observations.")
  }
  if (anyNA(time) || anyNA(status) || anyNA(lp)) {
    stop("`time`, `status`, and `lp` cannot contain missing values.")
  }
  if (any(!is.finite(time)) || any(time < 0)) {
    stop("`time` must contain finite, non-negative values.")
  }
  if (any(!is.finite(lp))) {
    stop("`lp` must contain finite values; LPs are not truncated or winsorized.")
  }

  status_numeric <- as.numeric(status)
  if (any(!status_numeric %in% c(0, 1))) {
    stop("`status` must be coded 0/1.")
  }

  order_desc <- order(time, decreasing = TRUE)
  time_desc <- time[order_desc]
  status_desc <- status_numeric[order_desc]
  lp_desc <- lp[order_desc]

  runs <- rle(time_desc)
  group_end <- cumsum(runs$lengths)
  group_start <- c(1L, head(group_end, -1L) + 1L)

  risk_set_log_weight <- -Inf
  risk_set_size <- 0L
  event_rows <- vector("list", length(runs$values))
  event_count <- 0L

  for (group in seq_along(runs$values)) {
    idx <- group_start[group]:group_end[group]
    group_log_weight <- .fold_oof_logsumexp(lp_desc[idx])
    risk_set_log_weight <- .fold_oof_logspace_add(
      risk_set_log_weight,
      group_log_weight
    )
    risk_set_size <- risk_set_size + length(idx)

    deaths <- sum(status_desc[idx] == 1)
    if (deaths > 0L) {
      event_count <- event_count + 1L
      log_increment <- log(deaths) - risk_set_log_weight

      event_rows[[event_count]] <- data.frame(
        event_time = runs$values[group],
        n_events = deaths,
        risk_set_size = risk_set_size,
        hazard_increment = exp(log_increment),
        log_hazard_increment = log_increment,
        stringsAsFactors = FALSE
      )
    }
  }

  if (event_count == 0L) {
    return(data.frame(
      event_time = numeric(),
      n_events = integer(),
      risk_set_size = integer(),
      hazard_increment = numeric(),
      log_hazard_increment = numeric(),
      baseline_cumhaz = numeric(),
      log_baseline_cumhaz = numeric(),
      stringsAsFactors = FALSE
    ))
  }

  baseline <- do.call(rbind, event_rows[seq_len(event_count)])
  baseline <- baseline[order(baseline$event_time), , drop = FALSE]
  rownames(baseline) <- NULL

  log_cumulative <- numeric(nrow(baseline))
  running_log_cumulative <- -Inf
  for (row in seq_len(nrow(baseline))) {
    running_log_cumulative <- .fold_oof_logspace_add(
      running_log_cumulative,
      baseline$log_hazard_increment[row]
    )
    log_cumulative[row] <- running_log_cumulative
  }

  baseline$baseline_cumhaz <- exp(log_cumulative)
  baseline$log_baseline_cumhaz <- log_cumulative
  baseline
}

#' Evaluate a Breslow baseline curve at fixed prediction times.
#'
#' The cumulative hazard is a right-continuous step function. It is zero before
#' the first event and remains at its last value beyond the final event.
evaluate_breslow_baseline <- function(baseline, prediction_times) {
  if (!is.data.frame(baseline) ||
      !all(c("event_time", "baseline_cumhaz", "log_baseline_cumhaz") %in%
           names(baseline))) {
    stop("`baseline` must be returned by estimate_breslow_baseline_from_lp().")
  }
  if (!is.numeric(prediction_times) || length(prediction_times) == 0L ||
      anyNA(prediction_times) || any(!is.finite(prediction_times)) ||
      any(prediction_times <= 0)) {
    stop("`prediction_times` must contain finite, positive numeric values.")
  }

  prediction_times <- as.numeric(prediction_times)
  index <- findInterval(prediction_times, baseline$event_time)
  cumulative <- numeric(length(prediction_times))
  log_cumulative <- rep(-Inf, length(prediction_times))

  has_event <- index > 0L
  cumulative[has_event] <- baseline$baseline_cumhaz[index[has_event]]
  log_cumulative[has_event] <-
    baseline$log_baseline_cumhaz[index[has_event]]

  data.frame(
    prediction_time = prediction_times,
    baseline_cumhaz = cumulative,
    log_baseline_cumhaz = log_cumulative,
    stringsAsFactors = FALSE
  )
}

#' Convert Cox LPs and an LP=0 baseline cumulative hazard to absolute risks.
#'
#' Uses a stable form of `1 - exp(-H0 * exp(lp))`. LPs are never removed,
#' truncated, or winsorized. Very large cumulative hazards correctly approach
#' a risk of 1.
absolute_risk_from_lp <- function(
  lp,
  baseline_cumhaz,
  log_baseline_cumhaz = NULL
) {
  if (!is.numeric(lp) || anyNA(lp) || any(!is.finite(lp))) {
    stop("`lp` must contain finite numeric values.")
  }
  if (!is.numeric(baseline_cumhaz) || length(baseline_cumhaz) != 1L ||
      is.na(baseline_cumhaz) || baseline_cumhaz < 0) {
    stop("`baseline_cumhaz` must be one non-negative numeric value.")
  }

  if (is.null(log_baseline_cumhaz)) {
    log_baseline_cumhaz <- if (baseline_cumhaz == 0) {
      -Inf
    } else {
      log(baseline_cumhaz)
    }
  }

  if (!is.numeric(log_baseline_cumhaz) ||
      length(log_baseline_cumhaz) != 1L ||
      is.na(log_baseline_cumhaz)) {
    stop("`log_baseline_cumhaz` must be one numeric value.")
  }

  log_subject_cumhaz <- lp + log_baseline_cumhaz
  risk <- numeric(length(lp))

  finite_subject_hazard <- log_subject_cumhaz <= log(.Machine$double.xmax)
  subject_cumhaz <- exp(log_subject_cumhaz[finite_subject_hazard])
  risk[finite_subject_hazard] <- -expm1(-subject_cumhaz)
  risk[!finite_subject_hazard] <- 1

  if (any(!is.finite(risk)) || any(risk < 0 | risk > 1)) {
    stop("Computed absolute risks are not finite values in [0, 1].")
  }

  risk
}

#' Add fully OOF absolute risks using fold-specific training-derived baselines.
#'
#' @param results Existing result from the repository's
#'   `cv_glmnet_cox_oof()` glmnet path.
#' @param df The exact original model dataframe supplied to the CV workflow
#'   (additional columns are allowed).
#' @param vars_to_scale The exact variables scaled in each outer fold. Supply
#'   `character()` when no explicit fold scaling was used.
#' @param prediction_times Positive time horizon(s), e.g. `c(1, 3, 5)`.
#' @param time_col,status_col,studyid_col Column names in `df`.
#' @param fold_col,lp_col Column names in `results$oof_predictions`.
#' @param lp_tolerance Relative tolerance for checking reproduced held-out LPs
#'   against the saved OOF LPs.
#' @param retain_cv_objects Keep the large saved `fold_models`/`full_model`
#'   objects in the returned list. Defaults to `FALSE` because they are no
#'   longer needed after risks are generated.
#' @param verbose Print fold summaries.
#'
#' @return A copy of `results` whose `oof_predictions` contains one
#'   `pred_risk_<time>yr` column per requested horizon. Additional elements are
#'   `fold_baseline_cumhaz`, `prediction_times`, `risk_column_map`, and
#'   `absolute_risk_method`.
compute_fold_specific_oof_absolute_risk <- function(
  results,
  df,
  vars_to_scale,
  prediction_times = c(1, 3, 5),
  time_col = "time",
  status_col = "status",
  studyid_col = "StudyID",
  fold_col = "fold",
  lp_col = "lp",
  lp_tolerance = 1e-8,
  retain_cv_objects = FALSE,
  verbose = TRUE
) {
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("Package `glmnet` is required to predict from saved fold models.")
  }
  if (!is.list(results) ||
      !all(c("oof_predictions", "fold_models", "fold_info", "predictors") %in%
           names(results))) {
    stop(
      "`results` must contain oof_predictions, fold_models, fold_info, ",
      "and predictors."
    )
  }
  if (!is.data.frame(df)) {
    stop("`df` must be a data.frame.")
  }

  model_engine <- results$model_engine
  if (is.null(model_engine)) {
    model_engine <- if (inherits(results$fold_models[[1]], "cv.glmnet")) {
      "glmnet"
    } else {
      NA_character_
    }
  }
  if (!identical(model_engine, "glmnet")) {
    stop(
      "This post-processing helper currently supports the repository's ",
      "glmnet CV path only; its explicit baseline uses glmnet's Breslow ties."
    )
  }

  predictors <- unique(results$predictors)
  vars_to_scale <- unique(vars_to_scale)
  invalid_scale <- setdiff(vars_to_scale, predictors)
  if (length(invalid_scale) > 0L) {
    stop(
      "`vars_to_scale` contains variables not used by the saved model: ",
      paste(invalid_scale, collapse = ", ")
    )
  }

  if (!is.numeric(prediction_times) || length(prediction_times) == 0L ||
      anyNA(prediction_times) || any(!is.finite(prediction_times)) ||
      any(prediction_times <= 0)) {
    stop("`prediction_times` must contain finite, positive numeric values.")
  }
  prediction_times <- sort(unique(as.numeric(prediction_times)))
  risk_columns <- vapply(
    prediction_times,
    .fold_oof_risk_column,
    character(1)
  )
  if (anyDuplicated(risk_columns)) {
    stop("`prediction_times` produce duplicate risk column names.")
  }

  if (!is.numeric(lp_tolerance) || length(lp_tolerance) != 1L ||
      is.na(lp_tolerance) || lp_tolerance < 0) {
    stop("`lp_tolerance` must be one non-negative number.")
  }
  if (!is.logical(retain_cv_objects) || length(retain_cv_objects) != 1L ||
      is.na(retain_cv_objects)) {
    stop("`retain_cv_objects` must be TRUE or FALSE.")
  }

  oof <- results$oof_predictions
  required_oof <- c(studyid_col, fold_col, lp_col, time_col, status_col)
  missing_oof <- setdiff(required_oof, names(oof))
  if (length(missing_oof) > 0L) {
    stop(
      "Missing columns in results$oof_predictions: ",
      paste(missing_oof, collapse = ", ")
    )
  }
  if (anyNA(oof[[fold_col]]) || anyNA(oof[[lp_col]]) ||
      any(!is.finite(oof[[lp_col]]))) {
    stop("Saved OOF folds and LPs must be complete and finite.")
  }

  required_df <- unique(c(studyid_col, time_col, status_col, predictors))
  missing_df <- setdiff(required_df, names(df))
  if (length(missing_df) > 0L) {
    stop("Missing columns in `df`: ", paste(missing_df, collapse = ", "))
  }

  complete <- stats::complete.cases(df[, required_df, drop = FALSE])
  complete_df <- df[complete, required_df, drop = FALSE]

  oof_ids <- as.character(oof[[studyid_col]])
  complete_ids <- as.character(complete_df[[studyid_col]])
  if (anyDuplicated(oof_ids) || anyDuplicated(complete_ids)) {
    stop("Participant identifiers must be unique in both `df` and OOF output.")
  }
  if (length(oof_ids) != length(complete_ids) ||
      !setequal(oof_ids, complete_ids)) {
    stop(
      "Complete-case participant identities in `df` do not exactly match ",
      "the saved OOF output. Supply the exact dataframe used for this model."
    )
  }

  matched <- match(oof_ids, complete_ids)
  model_df <- complete_df[matched, , drop = FALSE]
  rownames(model_df) <- NULL

  same_time <- isTRUE(all.equal(
    as.numeric(model_df[[time_col]]),
    as.numeric(oof[[time_col]]),
    tolerance = 1e-12,
    check.attributes = FALSE
  ))
  same_status <- identical(
    as.numeric(model_df[[status_col]]),
    as.numeric(oof[[status_col]])
  )
  if (!same_time || !same_status) {
    stop(
      "`df` outcomes do not match the saved OOF output. Supply the exact ",
      "model dataframe used by the CV run."
    )
  }

  folds <- sort(unique(oof[[fold_col]]))
  fold_info_rows <- match(folds, results$fold_info$fold)
  if (anyNA(fold_info_rows) || !"lambda_used" %in% names(results$fold_info)) {
    stop("Saved fold_info must contain one `lambda_used` value per fold.")
  }

  for (risk_col in risk_columns) {
    oof[[risk_col]] <- NA_real_
  }

  fold_baselines <- vector("list", length(folds))

  for (fold_index in seq_along(folds)) {
    fold <- folds[fold_index]
    test_index <- which(oof[[fold_col]] == fold)
    train_index <- which(oof[[fold_col]] != fold)

    train_data <- model_df[train_index, , drop = FALSE]
    test_data <- model_df[test_index, , drop = FALSE]

    scale_params <- .fold_oof_compute_scale_params(
      train_data,
      vars_to_scale
    )
    train_scaled <- .fold_oof_apply_scaling(train_data, scale_params)
    test_scaled <- .fold_oof_apply_scaling(test_data, scale_params)
    x_train <- .fold_oof_make_x_matrix(train_scaled, predictors)
    x_test <- .fold_oof_make_x_matrix(test_scaled, predictors)

    fold_name <- paste0("fold_", fold)
    fold_model <- results$fold_models[[fold_name]]
    if (is.null(fold_model) && fold_index <= length(results$fold_models)) {
      fold_model <- results$fold_models[[fold_index]]
    }
    if (!inherits(fold_model, "cv.glmnet")) {
      stop("Saved model for fold ", fold, " is not a cv.glmnet object.")
    }

    info_row <- fold_info_rows[fold_index]
    selected_lambda <- results$fold_info$lambda_used[info_row]
    if (length(selected_lambda) != 1L || !is.finite(selected_lambda)) {
      stop("Fold ", fold, " does not have one finite selected lambda.")
    }

    train_lp <- as.numeric(stats::predict(
      fold_model,
      newx = x_train,
      s = selected_lambda,
      type = "link"
    ))
    reproduced_test_lp <- as.numeric(stats::predict(
      fold_model,
      newx = x_test,
      s = selected_lambda,
      type = "link"
    ))

    saved_test_lp <- oof[[lp_col]][test_index]
    absolute_lp_difference <- abs(reproduced_test_lp - saved_test_lp)
    relative_lp_difference <- absolute_lp_difference /
      pmax(1, abs(saved_test_lp))
    max_relative_lp_difference <- max(relative_lp_difference)

    if (any(relative_lp_difference > lp_tolerance)) {
      stop(
        sprintf(
          paste0(
            "Reproduced held-out LPs do not match saved OOF LPs in fold %s ",
            "(max relative difference %.3g; tolerance %.3g). Check `df` and ",
            "`vars_to_scale`."
          ),
          as.character(fold),
          max_relative_lp_difference,
          lp_tolerance
        )
      )
    }

    baseline_curve <- estimate_breslow_baseline_from_lp(
      time = train_data[[time_col]],
      status = train_data[[status_col]],
      lp = train_lp
    )
    baseline_at_times <- evaluate_breslow_baseline(
      baseline_curve,
      prediction_times
    )

    max_train_followup <- max(train_data[[time_col]])
    beyond_followup <- prediction_times > max_train_followup
    if (any(beyond_followup)) {
      warning(
        sprintf(
          "Fold %s prediction time(s) exceed training follow-up (%.3f): %s",
          as.character(fold),
          max_train_followup,
          paste(prediction_times[beyond_followup], collapse = ", ")
        )
      )
    }

    for (time_index in seq_along(prediction_times)) {
      oof[[risk_columns[time_index]]][test_index] <- absolute_risk_from_lp(
        lp = saved_test_lp,
        baseline_cumhaz =
          baseline_at_times$baseline_cumhaz[time_index],
        log_baseline_cumhaz =
          baseline_at_times$log_baseline_cumhaz[time_index]
      )
    }

    fold_baselines[[fold_index]] <- data.frame(
      fold = fold,
      prediction_time = prediction_times,
      baseline_cumhaz = baseline_at_times$baseline_cumhaz,
      log_baseline_cumhaz = baseline_at_times$log_baseline_cumhaz,
      n_train = length(train_index),
      n_train_events = sum(train_data[[status_col]]),
      n_test = length(test_index),
      selected_lambda = selected_lambda,
      max_train_followup = max_train_followup,
      prediction_time_beyond_train_followup = beyond_followup,
      max_relative_test_lp_difference = max_relative_lp_difference,
      stringsAsFactors = FALSE
    )

    if (verbose) {
      summaries <- paste(
        sprintf(
          "H0(%g)=%.6g",
          prediction_times,
          baseline_at_times$baseline_cumhaz
        ),
        collapse = " | "
      )
      cat(sprintf(
        "Fold %s: n_train=%d, events=%d, lambda=%.6g, %s\n",
        as.character(fold),
        length(train_index),
        sum(train_data[[status_col]]),
        selected_lambda,
        summaries
      ))
    }
  }

  for (risk_col in risk_columns) {
    if (anyNA(oof[[risk_col]]) || any(!is.finite(oof[[risk_col]])) ||
        any(oof[[risk_col]] < 0 | oof[[risk_col]] > 1)) {
      stop("Column `", risk_col, "` must contain finite risks in [0, 1].")
    }
  }

  out <- results
  if (!retain_cv_objects) {
    out$fold_models <- NULL
    out$full_model <- NULL
  }
  out$oof_predictions <- oof
  out$fold_baseline_cumhaz <- do.call(rbind, fold_baselines)
  rownames(out$fold_baseline_cumhaz) <- NULL
  out$prediction_times <- prediction_times
  out$risk_column_map <- data.frame(
    prediction_time = prediction_times,
    risk_column = risk_columns,
    stringsAsFactors = FALSE
  )
  out$absolute_risk_method <-
    "fold-specific training-derived Breslow baseline; fully OOF absolute risk"
  out$vars_to_scale <- vars_to_scale

  out
}
