# ============================================================================
# Final-Model Cox Fit (no cross-validation, no held-out predictions)
#
# Trains a single Cox model on the FULL complete-case dataset, mirroring
# the regularisation/preprocessing pipeline used by `cv_glmnet_oof.r`,
# but without outer folds.  Intended for "ready-to-deploy" models that
# will be applied to an external cohort.
#
# Engines (parameter `model_engine`):
#   - "glmnet" (default) : regularised Cox via glmnet::cv.glmnet, with the
#                          *inner* CV used to select lambda.min.  This is the
#                          same scheme each fold of cv_glmnet_oof used,
#                          just applied to all observations at once.
#   - "coxph"            : standard Cox PH via survival::coxph (no penalty).
#
# Returned object is self-contained: it stores the scaling parameters
# computed on the training set, the predictor list, the lambda actually
# used (glmnet), and the fitted model.  Use `predict_lp_full()` to score
# new data — it reapplies the saved scaling and produces linear predictors
# on the same scale as the training-time LP.
#
# DEPENDENCIES
#   This file reuses small utilities defined in `cv_glmnet_oof.r`
#   (compute_scale_params, apply_scaling, make_x_matrix,
#    build_penalty_factor, safe_lambda_for_prediction, get_cvm_at_lambda,
#    coerce_predictors_numeric, build_coxph_formula).
#   Make sure that file is sourced before this one.
# ============================================================================

library(survival)
library(glmnet)

# Fail fast if the shared helpers are not loaded.
.required_helpers <- c(
  "compute_scale_params", "apply_scaling", "make_x_matrix",
  "build_penalty_factor", "safe_lambda_for_prediction", "get_cvm_at_lambda",
  "coerce_predictors_numeric", "build_coxph_formula"
)
.missing_helpers <- .required_helpers[!vapply(.required_helpers, exists, logical(1))]
if (length(.missing_helpers) > 0) {
  stop(
    "fit_full_cox.r requires helpers from cv_glmnet_oof.r. ",
    "Missing: ", paste(.missing_helpers, collapse = ", "),
    ". Please source('helpers/cv_glmnet_oof.r') first."
  )
}
rm(.required_helpers, .missing_helpers)

# ----------------------------------------------------------------------------
# Main entry point
# ----------------------------------------------------------------------------
#' Fit a final Cox model on the full dataset
#'
#' @param df              data.frame with `studyid_col`, `time_col`,
#'                        `status_col`, and all `predictors`
#' @param predictors      character vector of predictor column names
#' @param time_col        time column (default "time")
#' @param status_col      status column (default "status")
#' @param studyid_col     subject ID column (default "StudyID")
#' @param vars_to_scale   subset of `predictors` to z-score using mean/sd
#'                        computed from the full training set
#' @param inner_nfolds    folds for cv.glmnet's lambda search (glmnet only)
#' @param alpha           glmnet alpha (0 = ridge, default)
#' @param type.measure    cv.glmnet type.measure (default "C")
#' @param seed            RNG seed for reproducibility
#' @param standardize     glmnet standardize flag (forced FALSE if vars_to_scale
#'                        is supplied)
#' @param penalty.factor  named numeric, applied to a subset of `predictors`;
#'                        same convention as `cv_glmnet_cox_oof()`
#' @param nlambda         glmnet nlambda
#' @param verbose         print progress
#' @param model_engine    "glmnet" or "coxph"
#'
#' @return list with the fitted model and everything needed to score a
#'   new cohort.  See `predict_lp_full()`.
fit_full_cox <- function(
  df,
  predictors,
  time_col = "time",
  status_col = "status",
  studyid_col = "StudyID",
  vars_to_scale = NULL,
  inner_nfolds = 10,
  alpha = 0,
  type.measure = "C",
  seed = 123,
  standardize = FALSE,
  penalty.factor = NULL,
  nlambda = 100,
  verbose = TRUE,
  model_engine = "glmnet"
) {
  model_engine <- match.arg(model_engine, choices = c("glmnet", "coxph"))

  predictors <- unique(predictors)

  if (!is.data.frame(df)) {
    stop("df must be a data.frame")
  }

  required_cols <- c(studyid_col, time_col, status_col, predictors)
  missing_cols <- setdiff(required_cols, names(df))
  if (length(missing_cols) > 0) {
    stop("Missing required columns in df: ", paste(missing_cols, collapse = ", "))
  }

  if (!is.null(vars_to_scale)) {
    vars_to_scale <- unique(vars_to_scale)
    invalid_scale <- setdiff(vars_to_scale, predictors)
    if (length(invalid_scale) > 0) {
      stop("vars_to_scale contains variables not in predictors: ",
           paste(invalid_scale, collapse = ", "))
    }
  }

  dat <- df[, required_cols, drop = FALSE]
  dat <- dat[stats::complete.cases(dat), , drop = FALSE]

  n_complete <- nrow(dat)
  n_total <- nrow(df)
  if (n_complete == 0L) {
    stop("No complete cases remain after filtering required columns.")
  }

  # Engine-specific validation -----------------------------------------------
  if (model_engine == "glmnet") {
    if (!is.null(vars_to_scale) && standardize) {
      warning("Both vars_to_scale and standardize=TRUE specified. ",
              "Using vars_to_scale and forcing standardize=FALSE.")
      standardize <- FALSE
    }
    if (!is.null(penalty.factor)) {
      invalid_pf <- setdiff(names(penalty.factor), predictors)
      if (length(invalid_pf) > 0) {
        stop("penalty.factor contains variables not in predictors: ",
             paste(invalid_pf, collapse = ", "))
      }
    }
    if (inner_nfolds < 2L) {
      stop("inner_nfolds must be >= 2")
    }
  }

  if (verbose) {
    cat("\n============================================\n")
    cat("Final Cox Model (no fold splits)\n")
    cat("============================================\n")
    cat("Model engine:", model_engine, "\n")
    cat("Predictors:", length(predictors), "\n")
    if (!is.null(vars_to_scale)) {
      cat("Variables to scale:", length(vars_to_scale), "\n")
      cat("  (", paste(vars_to_scale, collapse = ", "), ")\n", sep = "")
    }
    if (model_engine == "glmnet") {
      cat("Inner folds (lambda tuning):", inner_nfolds, "\n")
      cat("Alpha:", alpha, "\n")
      cat("Type measure:", type.measure, "\n")
    }
    cat("Complete cases:", n_complete, "of", n_total,
        sprintf("(%.1f%%)\n", 100 * n_complete / n_total))
    cat("Seed:", seed, "\n\n")
  }

  # Scaling on the FULL training set (saved for downstream scoring) ---------
  scale_params <- NULL
  if (!is.null(vars_to_scale)) {
    scale_params <- compute_scale_params(dat, vars_to_scale)
    dat <- apply_scaling(dat, scale_params)
  }

  # Engine dispatch ----------------------------------------------------------
  if (model_engine == "glmnet") {
    fit_result <- .fit_full_glmnet(
      dat = dat,
      predictors = predictors,
      time_col = time_col,
      status_col = status_col,
      inner_nfolds = inner_nfolds,
      alpha = alpha,
      type.measure = type.measure,
      seed = seed,
      standardize = standardize,
      penalty.factor = penalty.factor,
      nlambda = nlambda,
      verbose = verbose
    )
  } else {
    fit_result <- .fit_full_coxph(
      dat = dat,
      predictors = predictors,
      time_col = time_col,
      status_col = status_col,
      seed = seed,
      verbose = verbose
    )
  }

  # In-sample (apparent) C-index --------------------------------------------
  surv_obj <- survival::Surv(dat[[time_col]], dat[[status_col]])
  c_result <- survival::concordance(surv_obj ~ I(-fit_result$lp_train))
  c_index_apparent <- c_result$concordance
  se_apparent      <- sqrt(c_result$var)

  if (verbose) {
    cat("\n--- In-sample (apparent) performance ---\n")
    cat("  C-index:", round(c_index_apparent, 4),
        sprintf("[SE %.4f]\n", se_apparent))
    cat("  (NOTE: optimistic; use OOF results from cv_glmnet_cox_oof()\n",
        "         for honest performance estimates.)\n", sep = "")
    cat("============================================\n")
    cat("Final model fit complete!\n")
    cat("============================================\n\n")
  }

  out <- list(
    model               = fit_result$model,
    model_engine        = model_engine,
    predictors          = predictors,
    vars_to_scale       = vars_to_scale,
    scale_params        = scale_params,
    time_col            = time_col,
    status_col          = status_col,
    studyid_col         = studyid_col,
    n                   = n_complete,
    n_events            = sum(dat[[status_col]]),
    seed                = seed,
    c_index_apparent    = c_index_apparent,
    se_apparent         = se_apparent,
    train_lp            = fit_result$lp_train,
    train_studyids      = dat[[studyid_col]]
  )

  if (model_engine == "glmnet") {
    out$alpha          <- alpha
    out$inner_nfolds   <- inner_nfolds
    out$type.measure   <- type.measure
    out$nlambda        <- nlambda
    out$standardize    <- standardize
    out$penalty.factor <- penalty.factor
    out$lambda_min     <- fit_result$lambda_min
    out$lambda_used    <- fit_result$lambda_used
    out$cvm_selected   <- fit_result$cvm_selected
  } else {
    out$alpha          <- NA_real_
    out$inner_nfolds   <- NA_integer_
    out$lambda_min     <- NA_real_
    out$lambda_used    <- NA_real_
  }

  out
}

# ----------------------------------------------------------------------------
# Engine: glmnet
# ----------------------------------------------------------------------------
.fit_full_glmnet <- function(
  dat,
  predictors,
  time_col,
  status_col,
  inner_nfolds,
  alpha,
  type.measure,
  seed,
  standardize,
  penalty.factor,
  nlambda,
  verbose
) {
  x_train <- make_x_matrix(dat, predictors)
  y_train <- survival::Surv(dat[[time_col]], dat[[status_col]])

  pf <- build_penalty_factor(predictors, penalty.factor)

  set.seed(seed)

  glmnet_args <- list(
    x = x_train,
    y = y_train,
    family = "cox",
    type.measure = type.measure,
    alpha = alpha,
    nfolds = inner_nfolds,
    parallel = FALSE,
    standardize = standardize,
    nlambda = nlambda
  )
  if (!is.null(pf)) {
    glmnet_args$penalty.factor <- pf
  }

  model <- do.call(glmnet::cv.glmnet, glmnet_args)

  lambda_used  <- safe_lambda_for_prediction(model, verbose = verbose)
  cvm_selected <- get_cvm_at_lambda(model, lambda_used)

  lp_train <- as.numeric(
    predict(model, newx = x_train, s = lambda_used, type = "link")
  )

  if (verbose) {
    cat(sprintf(
      "Final glmnet fit: lambda.min=%s | lambda.used=%.6f | CV score=%.4f\n",
      ifelse(is.finite(model$lambda.min),
             sprintf("%.6f", model$lambda.min), "Inf"),
      lambda_used,
      cvm_selected
    ))
  }

  list(
    model        = model,
    lp_train     = lp_train,
    lambda_min   = model$lambda.min,
    lambda_used  = lambda_used,
    cvm_selected = cvm_selected
  )
}

# ----------------------------------------------------------------------------
# Engine: coxph
# ----------------------------------------------------------------------------
.fit_full_coxph <- function(
  dat,
  predictors,
  time_col,
  status_col,
  seed,
  verbose
) {
  dat <- coerce_predictors_numeric(dat, predictors)
  fml <- build_coxph_formula(time_col, status_col, predictors)

  set.seed(seed)
  model <- survival::coxph(fml, data = dat)

  lp_train <- as.numeric(predict(model, newdata = dat, type = "lp"))

  list(
    model    = model,
    lp_train = lp_train
  )
}

# ----------------------------------------------------------------------------
# Score a new cohort with a saved fit
# ----------------------------------------------------------------------------
#' Apply a fitted full Cox model to new data
#'
#' Predicts the linear predictor (`lp`) for `newdata`.  By default the function
#' will fail on rows missing any required predictor; pass `drop_incomplete = TRUE`
#' to silently drop them and report how many were dropped.
#'
#' @param fit             output of `fit_full_cox()`
#' @param newdata         data.frame containing the predictors and (optionally)
#'                        StudyID/time/status columns
#' @param drop_incomplete drop rows missing any predictor (default FALSE)
#' @param verbose         print summary (default TRUE)
#' @param scaling         which scaling parameters to apply to `vars_to_scale`
#'                        before predicting:
#'                        \itemize{
#'                          \item "train" (default) - use the training-time
#'                            mean/sd stored in `fit$scale_params`. This is
#'                            the standard external-validation behavior.
#'                          \item "self" - recompute mean/sd from `newdata`
#'                            itself (after complete-case filtering on the
#'                            predictor block) and apply those instead. Useful
#'                            as a diagnostic when the new cohort has a very
#'                            different predictor distribution; the resulting
#'                            LP / C-index is no longer a clean external
#'                            validation but a "moment-matched" transport.
#'                          \item "self_nmr" - use training-time mean/sd for all
#'                            scaled variables except `nmr_vars`, whose mean/sd
#'                            are recomputed from `newdata`.
#'                          \item "none" - skip scaling entirely. Only sensible
#'                            if `vars_to_scale` was NULL at fit time.
#'                        }
#' @param nmr_vars         NMR predictor names to self-scale when
#'                         `scaling = "self_nmr"`. Variables absent from the
#'                         fitted model are ignored.
#'
#' @return data.frame with at least `lp` and (if available)
#'   `StudyID`, `time`, `status`. The data frame carries an attribute
#'   `"scaling_used"` documenting which mode produced the LPs.
predict_lp_full <- function(fit, newdata,
                            drop_incomplete = FALSE,
                            verbose = TRUE,
                            scaling = c("train", "self", "self_nmr", "none"),
                            nmr_vars = NULL) {
  if (!is.data.frame(newdata)) stop("newdata must be a data.frame")

  scaling <- match.arg(scaling)

  predictors    <- fit$predictors
  studyid_col   <- fit$studyid_col
  time_col      <- fit$time_col
  status_col    <- fit$status_col
  vars_to_scale <- fit$vars_to_scale

  missing_predictors <- setdiff(predictors, names(newdata))
  if (length(missing_predictors) > 0) {
    stop("newdata is missing required predictors: ",
         paste(missing_predictors, collapse = ", "))
  }

  pred_complete <- stats::complete.cases(newdata[, predictors, drop = FALSE])
  n_drop <- sum(!pred_complete)

  if (n_drop > 0) {
    if (!drop_incomplete) {
      stop(sprintf(
        "%d rows in newdata have NAs in predictors. Set drop_incomplete=TRUE to drop them.",
        n_drop
      ))
    }
    if (verbose) {
      cat(sprintf("Dropping %d of %d rows with NA predictors\n",
                  n_drop, nrow(newdata)))
    }
    newdata <- newdata[pred_complete, , drop = FALSE]
  }

  # Decide which scale_params to use ------------------------------------
  scale_params_used <- NULL

  if (scaling == "train") {
    if (!is.null(fit$scale_params)) {
      scale_params_used <- fit$scale_params
    }
  } else if (scaling == "self") {
    if (is.null(vars_to_scale) || length(vars_to_scale) == 0L) {
      if (verbose) {
        cat("scaling='self' requested but fit has no vars_to_scale; ",
            "no scaling will be applied.\n", sep = "")
      }
    } else {
      scale_params_used <- compute_scale_params(newdata, vars_to_scale)
      if (verbose) {
        cat(sprintf(
          "scaling='self': using %s mean/sd computed on newdata (%d rows, %d vars).\n",
          "newdata-derived", nrow(newdata), length(vars_to_scale)
        ))
        cat("  NOTE: results are no longer a clean external validation;\n",
            "        treat as a moment-matched diagnostic.\n", sep = "")
      }
    }
  } else if (scaling == "self_nmr") {
    if (is.null(nmr_vars) || length(nmr_vars) == 0L) {
      stop("scaling='self_nmr' requires a non-empty `nmr_vars` vector.")
    }
    if (is.null(fit$scale_params)) {
      stop("scaling='self_nmr' requires training-time scale parameters in the fit.")
    }

    scale_params_used <- fit$scale_params
    nmr_to_self_scale <- intersect(unique(nmr_vars), vars_to_scale)
    if (length(nmr_to_self_scale) > 0L) {
      self_nmr_params <- compute_scale_params(newdata, nmr_to_self_scale)
      scale_params_used[nmr_to_self_scale] <- self_nmr_params
    }

    if (verbose) {
      cat(sprintf(
        paste0(
          "scaling='self_nmr': using training mean/sd for %d variable(s) ",
          "and newdata mean/sd for %d NMR variable(s) (%d rows).\n"
        ),
        length(scale_params_used) - length(nmr_to_self_scale),
        length(nmr_to_self_scale),
        nrow(newdata)
      ))
    }
  }
  # scaling == "none" -> leave scale_params_used = NULL

  newdata_scaled <- newdata
  if (!is.null(scale_params_used)) {
    newdata_scaled <- apply_scaling(newdata_scaled, scale_params_used)
  }

  if (fit$model_engine == "glmnet") {
    x_new <- make_x_matrix(newdata_scaled, predictors)
    lp <- as.numeric(
      predict(fit$model, newx = x_new, s = fit$lambda_used, type = "link")
    )
  } else {
    newdata_scaled <- coerce_predictors_numeric(newdata_scaled, predictors)
    lp <- as.numeric(predict(fit$model, newdata = newdata_scaled, type = "lp"))
  }

  out <- data.frame(lp = lp, stringsAsFactors = FALSE)
  if (studyid_col %in% names(newdata)) {
    out[[studyid_col]] <- newdata[[studyid_col]]
  }
  if (time_col %in% names(newdata)) {
    out[[time_col]] <- newdata[[time_col]]
  }
  if (status_col %in% names(newdata)) {
    out[[status_col]] <- newdata[[status_col]]
  }

  id_cols <- intersect(c(studyid_col, time_col, status_col), names(out))
  out <- out[, c(id_cols, "lp"), drop = FALSE]

  attr(out, "scaling_used")   <- scaling
  attr(out, "scale_params")   <- scale_params_used
  attr(out, "self_scaled_vars") <- if (scaling == "self_nmr") {
    nmr_to_self_scale
  } else if (scaling == "self") {
    vars_to_scale
  } else {
    character()
  }
  out
}

# ----------------------------------------------------------------------------
# Convenience: external-cohort C-index from a saved fit
# ----------------------------------------------------------------------------
#' Compute concordance for a saved fit on a new cohort
#'
#' `newdata` must contain `time_col`, `status_col`, and all predictors.
#'
#' @param scaling Forwarded to `predict_lp_full()`. Use `"train"` (default) for
#'   a clean external validation. Use `"self"` as a diagnostic when the new
#'   cohort's predictor distribution differs substantially from the training
#'   cohort and you want to see how much of the C-index drop is attributable
#'   to scale/location shift versus genuine signal loss.
external_cindex <- function(fit, newdata,
                            drop_incomplete = TRUE,
                            verbose = TRUE,
                            scaling = c("train", "self", "self_nmr", "none"),
                            nmr_vars = NULL) {
  scaling <- match.arg(scaling)

  scored <- predict_lp_full(fit, newdata,
                            drop_incomplete = drop_incomplete,
                            verbose = verbose,
                            scaling = scaling,
                            nmr_vars = nmr_vars)

  if (!fit$time_col %in% names(scored) || !fit$status_col %in% names(scored)) {
    stop("newdata must contain '", fit$time_col, "' and '",
         fit$status_col, "' columns to compute concordance.")
  }

  surv_obj <- survival::Surv(scored[[fit$time_col]], scored[[fit$status_col]])
  c_result <- survival::concordance(surv_obj ~ I(-scored$lp))
  se <- sqrt(c_result$var)

  list(
    c_index      = c_result$concordance,
    se           = se,
    ci_lower     = c_result$concordance - 1.96 * se,
    ci_upper     = c_result$concordance + 1.96 * se,
    n            = nrow(scored),
    n_events     = sum(scored[[fit$status_col]]),
    scaling_used = scaling,
    scored       = scored
  )
}
