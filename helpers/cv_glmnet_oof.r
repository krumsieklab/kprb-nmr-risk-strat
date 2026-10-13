# ============================================================================
# Out-of-Fold (OOF) Cross-Validation for Cox Models
#
# The coxph engine is intended for low-dimensional SOC / SOC+PRS models where
# regularisation is unnecessary.  It still performs the same outer-fold OOF
# cross-validation (scaling per-fold, predicting held-out LPs, computing an
# overall C-index), so every downstream consumer that reads `oof_predictions`,
# `c_index`, `ci_lower`, `ci_upper`, `n`, or `n_events` works without changes.
#
# Differences when model_engine = "coxph":
#   - Parameters specific to glmnet (inner_nfolds, alpha, type.measure,
#     standardize, penalty.factor, nlambda, n_cores, parallel_outer) are
#     IGNORED.  A message is printed when verbose = TRUE.
#   - fold_models[[i]] contains a coxph object (not cv.glmnet).
#     plot(fold_model) will produce Schoenfeld-residual diagnostics, NOT
#     a lambda-vs-deviance curve.  Use plot_cv_folds() only for glmnet.
#   - fold_info keeps the same column names for RDS compatibility, but
#     lambda_min, lambda_used, and cvm_selected are NA_real_.
#   - Top-level lambda_values / lambda_used are NULL.
#   - alpha and inner_nfolds are NA in the return list.
#   - A new element `model_engine` is added to the return list for both
#     engines so downstream code can branch if needed.
#
# ============================================================================

library(survival)
library(glmnet)
library(foreach)
library(doParallel)

# ----------------------------------------------------------------------------
# Parallel helpers
# ----------------------------------------------------------------------------
cleanup_all_clusters <- function(verbose = FALSE) {
  tryCatch({
    if (verbose) cat("Unregistering parallel backend...\n")
    foreach::registerDoSEQ()

    # Aggressively close leftover PSOCK connections that may linger after
    # stopCluster().  On macOS, socket FDs from dead clusters can accumulate
    # across many create/destroy cycles within a single R session, eventually
    # exhausting the 128-connection limit and causing makeCluster() to hang.
    all_cons <- showConnections(all = FALSE)
    if (nrow(all_cons) > 0) {
      sock_rows <- grep("sockconn", all_cons[, "class"], ignore.case = TRUE)
      for (i in sock_rows) {
        close(getConnection(as.integer(rownames(all_cons)[i])))
      }
    }

    invisible(gc(verbose = FALSE))
  }, error = function(e) {
    if (verbose) warning("Error during cleanup: ", e$message)
  })
}

setup_parallel_cluster <- function(n_cores = NULL, max_cores = NULL) {
  if (is.null(n_cores)) {
    n_cores <- max(1L, parallel::detectCores() - 1L)
  } else {
    n_cores <- max(1L, as.integer(n_cores))
  }

  # Never spawn more PSOCK workers than there are tasks. Extra idle workers
  # still receive a serialized copy of `dat` via foreach, which is the main
  # RAM spike for NMR-sized design matrices.
  if (!is.null(max_cores)) {
    n_cores <- min(n_cores, max(1L, as.integer(max_cores)))
  }

  cl <- parallel::makeCluster(n_cores)
  doParallel::registerDoParallel(cl)

  list(cluster = cl, n_cores = n_cores)
}

stop_parallel_cluster <- function(cl, verbose = FALSE) {
  if (!is.null(cl)) {
    tryCatch({
      parallel::stopCluster(cl)
    }, error = function(e) {
      if (verbose) warning("Error stopping cluster: ", e$message)
    })
  }
  cleanup_all_clusters(verbose = FALSE)
  invisible(NULL)
}

# ----------------------------------------------------------------------------
# Scaling helpers
# ----------------------------------------------------------------------------
compute_scale_params <- function(data, vars_to_scale) {
  scale_params <- vector("list", length(vars_to_scale))
  names(scale_params) <- vars_to_scale

  for (var in vars_to_scale) {
    scale_params[[var]] <- list(
      center = mean(data[[var]], na.rm = TRUE),
      scale  = stats::sd(data[[var]], na.rm = TRUE)
    )
  }

  scale_params
}

apply_scaling <- function(data, scale_params) {
  data_scaled <- data

  for (var in names(scale_params)) {
    s <- scale_params[[var]]$scale
    m <- scale_params[[var]]$center

    if (is.finite(s) && !is.na(s) && s > 0) {
      data_scaled[[var]] <- (data_scaled[[var]] - m) / s
    }
  }

  data_scaled
}

# ----------------------------------------------------------------------------
# Matrix helpers (glmnet path)
# ----------------------------------------------------------------------------
make_x_matrix <- function(data, predictors) {
  xdf <- data[, predictors, drop = FALSE]

  for (nm in names(xdf)) {
    if (is.factor(xdf[[nm]])) {
      xdf[[nm]] <- as.numeric(as.character(xdf[[nm]]))
    } else if (is.logical(xdf[[nm]])) {
      xdf[[nm]] <- as.integer(xdf[[nm]])
    }
  }

  x <- data.matrix(xdf)
  storage.mode(x) <- "double"

  if (anyNA(x)) {
    stop("Predictor matrix contains NA after conversion.")
  }

  colnames(x) <- predictors
  x
}

build_penalty_factor <- function(predictors, penalty.factor = NULL) {
  if (is.null(penalty.factor)) {
    return(NULL)
  }

  pf <- rep(1, length(predictors))
  names(pf) <- predictors
  pf[names(penalty.factor)] <- penalty.factor[names(penalty.factor)]
  unname(pf)
}

safe_lambda_for_prediction <- function(fold_model, verbose = FALSE) {
  lam <- fold_model$lambda.min

  if (length(lam) != 1L || is.na(lam) || !is.finite(lam)) {
    finite_lambda <- fold_model$lambda[is.finite(fold_model$lambda)]

    if (length(finite_lambda) == 0L) {
      stop("No finite lambda available for prediction.")
    }

    lam <- max(finite_lambda)

    if (verbose) {
      warning(
        sprintf(
          "Non-finite lambda.min encountered; using largest finite lambda instead: %.6f",
          lam
        )
      )
    }
  }

  lam
}

get_cvm_at_lambda <- function(fold_model, lambda_value) {
  idx <- which.min(abs(fold_model$lambda - lambda_value))
  fold_model$cvm[idx]
}

# cv.glmnet via do.call() stores the training design matrix in $call (and
# again in $glmnet.fit$call). For NMR-sized specs that is hundreds of MB per
# fold; returning five of them through PSOCK and saveRDS is what produces
# "error writing to connection" after C-index already printed. plot/predict/
# coef still work without the captured call.
strip_glmnet_call <- function(model) {
  if (inherits(model, "cv.glmnet")) {
    model$call <- NULL
    if (!is.null(model$glmnet.fit)) model$glmnet.fit$call <- NULL
  }
  model
}

# ----------------------------------------------------------------------------
# Data-frame helpers (coxph path)
# ----------------------------------------------------------------------------
# Coerces factors/logicals in `predictors` to numeric, mirroring the
# conversions that make_x_matrix() performs for the glmnet path.  This keeps
# the model matrices identical across engines.
coerce_predictors_numeric <- function(data, predictors) {
  for (nm in predictors) {
    if (is.factor(data[[nm]])) {
      data[[nm]] <- as.numeric(as.character(data[[nm]]))
    } else if (is.logical(data[[nm]])) {
      data[[nm]] <- as.integer(data[[nm]])
    }
  }
  data
}

build_coxph_formula <- function(time_col, status_col, predictors) {
  rhs <- paste(predictors, collapse = " + ")
  as.formula(paste0("Surv(", time_col, ", ", status_col, ") ~ ", rhs))
}

# ----------------------------------------------------------------------------
# Single fold worker — glmnet
# ----------------------------------------------------------------------------
fit_one_fold_glmnet <- function(
  fold,
  dat,
  fold_ids,
  predictors,
  time_col,
  status_col,
  vars_to_scale,
  inner_nfolds,
  alpha,
  type.measure,
  seed,
  standardize,
  penalty.factor,
  nlambda,
  verbose = FALSE
) {
  test_idx  <- which(fold_ids == fold)
  train_idx <- which(fold_ids != fold)

  train_dat <- dat[train_idx, , drop = FALSE]
  test_dat  <- dat[test_idx, , drop = FALSE]

  if (!is.null(vars_to_scale)) {
    scale_params <- compute_scale_params(train_dat, vars_to_scale)
    train_dat <- apply_scaling(train_dat, scale_params)
    test_dat  <- apply_scaling(test_dat, scale_params)
  }

  x_train <- make_x_matrix(train_dat, predictors)
  x_test  <- make_x_matrix(test_dat, predictors)
  y_train <- survival::Surv(train_dat[[time_col]], train_dat[[status_col]])

  pf <- build_penalty_factor(predictors, penalty.factor)

  set.seed(seed + fold)

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

  fold_model <- strip_glmnet_call(do.call(glmnet::cv.glmnet, glmnet_args))

  lambda_used <- safe_lambda_for_prediction(fold_model, verbose = verbose)

  lp_test <- as.numeric(
    predict(fold_model, newx = x_test, s = lambda_used, type = "link")
  )

  list(
    test_idx      = test_idx,
    lp_test       = lp_test,
    fold_model    = fold_model,
    lambda_min    = fold_model$lambda.min,
    lambda_used   = lambda_used,
    cvm_selected  = get_cvm_at_lambda(fold_model, lambda_used),
    sd_lp_test    = stats::sd(lp_test),
    n_train       = length(train_idx),
    n_test        = length(test_idx)
  )
}

# ----------------------------------------------------------------------------
# Single fold worker — coxph
# ----------------------------------------------------------------------------
fit_one_fold_coxph <- function(
  fold,
  dat,
  fold_ids,
  predictors,
  time_col,
  status_col,
  vars_to_scale,
  seed,
  verbose = FALSE
) {
  test_idx  <- which(fold_ids == fold)
  train_idx <- which(fold_ids != fold)

  train_dat <- dat[train_idx, , drop = FALSE]
  test_dat  <- dat[test_idx, , drop = FALSE]

  if (!is.null(vars_to_scale)) {
    scale_params <- compute_scale_params(train_dat, vars_to_scale)
    train_dat <- apply_scaling(train_dat, scale_params)
    test_dat  <- apply_scaling(test_dat, scale_params)
  }

  train_dat <- coerce_predictors_numeric(train_dat, predictors)
  test_dat  <- coerce_predictors_numeric(test_dat, predictors)

  fml <- build_coxph_formula(time_col, status_col, predictors)

  set.seed(seed + fold)

  fold_model <- survival::coxph(fml, data = train_dat)

  lp_test <- as.numeric(
    predict(fold_model, newdata = test_dat, type = "lp")
  )

  list(
    test_idx    = test_idx,
    lp_test     = lp_test,
    fold_model  = fold_model,
    sd_lp_test  = stats::sd(lp_test),
    n_train     = length(train_idx),
    n_test      = length(test_idx)
  )
}

# ----------------------------------------------------------------------------
# Main OOF function
# ----------------------------------------------------------------------------
#' @param model_engine  Character, either "glmnet" (default) or "coxph".
#'   When "coxph", fits a standard Cox PH model per fold instead of a
#'   regularised glmnet model.  Parameters that only apply to glmnet
#'   (inner_nfolds, alpha, type.measure, standardize, penalty.factor,
#'   nlambda, n_cores, parallel_outer) are silently ignored (a verbose
#'   message is printed).
cv_glmnet_cox_oof <- function(
  df,
  predictors,
  time_col = "time",
  status_col = "status",
  studyid_col = "StudyID",
  vars_to_scale = NULL,
  outer_nfolds = 5,
  inner_nfolds = 10,
  alpha = 0,
  type.measure = "C",
  seed = 123,
  n_cores = NULL,
  standardize = FALSE,
  penalty.factor = NULL,
  nlambda = 50,
  parallel_outer = TRUE,
  verbose = TRUE,
  model_engine = "glmnet"
) {
  # ------------------------------------------------------------------
  # Validate model_engine
  # ------------------------------------------------------------------
  model_engine <- match.arg(model_engine, choices = c("glmnet", "coxph"))

  cleanup_all_clusters()

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

  if (outer_nfolds < 2L) {
    stop("outer_nfolds must be >= 2")
  }

  # ------------------------------------------------------------------
  # Engine-specific validation
  # ------------------------------------------------------------------
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

  if (model_engine == "coxph" && verbose) {
    cat("\n[model_engine = 'coxph'] Ignoring glmnet-specific parameters:\n")
    cat("  inner_nfolds, alpha, type.measure, standardize,\n")
    cat("  penalty.factor, nlambda, n_cores, parallel_outer\n")
  }

  # ------------------------------------------------------------------
  # Banner
  # ------------------------------------------------------------------
  if (verbose) {
    cat("\n============================================\n")
    cat("Out-of-Fold Cross-Validation\n")
    cat("============================================\n")
    cat("Model engine:", model_engine, "\n")
    cat("Predictors:", length(predictors), "\n")
    if (!is.null(vars_to_scale)) {
      cat("Variables to scale:", length(vars_to_scale), "\n")
      cat("  (", paste(vars_to_scale, collapse = ", "), ")\n", sep = "")
    }
    cat("Outer folds:", outer_nfolds, "\n")
    if (model_engine == "glmnet") {
      cat("Inner folds (lambda tuning):", inner_nfolds, "\n")
      cat("Alpha:", alpha, "\n")
      cat("Type measure:", type.measure, "\n")
      cat("Parallelization:",
          ifelse(parallel_outer, "OUTER folds", "SEQUENTIAL outer folds"), "\n")
    }
    cat("Complete cases:", n_complete, "of", n_total,
        sprintf("(%.1f%%)\n", 100 * n_complete / n_total))
    cat("Seed:", seed, "\n\n")
  }

  # ------------------------------------------------------------------
  # Outer fold assignments (identical across engines)
  # ------------------------------------------------------------------
  set.seed(seed)
  fold_ids <- sample(rep(seq_len(outer_nfolds), length.out = n_complete))

  oof_predictions <- data.frame(
    StudyID = dat[[studyid_col]],
    fold    = fold_ids,
    lp      = NA_real_,
    time    = dat[[time_col]],
    status  = dat[[status_col]],
    stringsAsFactors = FALSE
  )

  # ------------------------------------------------------------------
  # Dispatch to engine
  # ------------------------------------------------------------------
  if (model_engine == "glmnet") {
    result <- run_oof_glmnet(
      dat = dat,
      fold_ids = fold_ids,
      oof_predictions = oof_predictions,
      predictors = predictors,
      time_col = time_col,
      status_col = status_col,
      vars_to_scale = vars_to_scale,
      outer_nfolds = outer_nfolds,
      inner_nfolds = inner_nfolds,
      alpha = alpha,
      type.measure = type.measure,
      seed = seed,
      n_cores = n_cores,
      standardize = standardize,
      penalty.factor = penalty.factor,
      nlambda = nlambda,
      parallel_outer = parallel_outer,
      verbose = verbose
    )
  } else {
    result <- run_oof_coxph(
      dat = dat,
      fold_ids = fold_ids,
      oof_predictions = oof_predictions,
      predictors = predictors,
      time_col = time_col,
      status_col = status_col,
      vars_to_scale = vars_to_scale,
      outer_nfolds = outer_nfolds,
      seed = seed,
      verbose = verbose
    )
  }

  # ------------------------------------------------------------------
  # Overall C-index (identical across engines)
  # ------------------------------------------------------------------
  oof_predictions <- result$oof_predictions

  if (anyNA(oof_predictions$lp)) {
    stop("OOF predictions contain NA values after fold aggregation.")
  }

  if (verbose) {
    cat("\n--- Computing Overall C-Index ---\n")
  }

  surv_obj <- with(oof_predictions, survival::Surv(time, status))
  c_result <- survival::concordance(surv_obj ~ I(-lp), data = oof_predictions)

  c_index  <- c_result$concordance
  se       <- sqrt(c_result$var)
  ci_lower <- c_index - 1.96 * se
  ci_upper <- c_index + 1.96 * se
  n_events <- sum(oof_predictions$status)

  if (verbose) {
    cat("  C-index:", round(c_index, 4), "\n")
    cat("  95% CI: [", round(ci_lower, 4), ",", round(ci_upper, 4), "]\n")
    cat("  n:", nrow(oof_predictions), "\n")
    cat("  Events:", n_events, "\n")
    cat("============================================\n")
    cat("OOF Cross-Validation Complete!\n")
    cat("============================================\n\n")
  }

  cleanup_all_clusters()

  # ------------------------------------------------------------------
  # Assemble return list
  # ------------------------------------------------------------------
  out <- list(
    oof_predictions = oof_predictions,
    fold_models     = result$fold_models,
    fold_info       = result$fold_info,
    c_index         = c_index,
    se              = se,
    ci_lower        = ci_lower,
    ci_upper        = ci_upper,
    n               = nrow(oof_predictions),
    n_events        = n_events,
    predictors      = predictors,
    outer_nfolds    = outer_nfolds,
    model_engine    = model_engine
  )

  if (model_engine == "glmnet") {
    out$lambda_values <- result$lambda_values
    out$lambda_used   <- result$lambda_used
    out$alpha         <- alpha
    out$inner_nfolds  <- inner_nfolds
  } else {
    out$lambda_values <- NULL
    out$lambda_used   <- NULL
    out$alpha         <- NA_real_
    out$inner_nfolds  <- NA_integer_
  }

  out
}

# ============================================================================
# Engine: glmnet (unchanged logic from v2)
# ============================================================================
run_oof_glmnet <- function(
  dat,
  fold_ids,
  oof_predictions,
  predictors,
  time_col,
  status_col,
  vars_to_scale,
  outer_nfolds,
  inner_nfolds,
  alpha,
  type.measure,
  seed,
  n_cores,
  standardize,
  penalty.factor,
  nlambda,
  parallel_outer,
  verbose
) {
  fold_models   <- vector("list", outer_nfolds)
  lambda_values <- rep(NA_real_, outer_nfolds)
  lambda_used   <- rep(NA_real_, outer_nfolds)
  cvm_selected  <- rep(NA_real_, outer_nfolds)
  sd_lp_test    <- rep(NA_real_, outer_nfolds)
  n_train_vec   <- rep(NA_integer_, outer_nfolds)
  n_test_vec    <- rep(NA_integer_, outer_nfolds)

  if (parallel_outer) {
    cl <- NULL
    on.exit(stop_parallel_cluster(cl), add = TRUE)
    cl_info <- setup_parallel_cluster(n_cores, max_cores = outer_nfolds)
    cl <- cl_info$cluster
    n_cores_used <- cl_info$n_cores

    if (verbose) {
      cat("Running outer folds in parallel on", n_cores_used, "cores\n")
    }

    fold_results <- foreach(
      fold = seq_len(outer_nfolds),
      .packages = c("survival", "glmnet"),
      .export = c(
        "compute_scale_params",
        "apply_scaling",
        "make_x_matrix",
        "build_penalty_factor",
        "safe_lambda_for_prediction",
        "get_cvm_at_lambda",
        "strip_glmnet_call",
        "fit_one_fold_glmnet"
      )
    ) %dopar% {
      fit_one_fold_glmnet(
        fold = fold,
        dat = dat,
        fold_ids = fold_ids,
        predictors = predictors,
        time_col = time_col,
        status_col = status_col,
        vars_to_scale = vars_to_scale,
        inner_nfolds = inner_nfolds,
        alpha = alpha,
        type.measure = type.measure,
        seed = seed,
        standardize = standardize,
        penalty.factor = penalty.factor,
        nlambda = nlambda,
        verbose = FALSE
      )
    }

    # Free workers before the master copies fold_models out of fold_results.
    # on.exit remains as a safety net if foreach errors.
    stop_parallel_cluster(cl)
    cl <- NULL

  } else {
    fold_results <- vector("list", outer_nfolds)

    for (fold in seq_len(outer_nfolds)) {
      if (verbose) {
        cat("\n--- Processing Fold", fold, "of", outer_nfolds, "---\n")
      }

      fold_results[[fold]] <- fit_one_fold_glmnet(
        fold = fold,
        dat = dat,
        fold_ids = fold_ids,
        predictors = predictors,
        time_col = time_col,
        status_col = status_col,
        vars_to_scale = vars_to_scale,
        inner_nfolds = inner_nfolds,
        alpha = alpha,
        type.measure = type.measure,
        seed = seed,
        standardize = standardize,
        penalty.factor = penalty.factor,
        nlambda = nlambda,
        verbose = verbose
      )
    }
  }

  for (fold in seq_len(outer_nfolds)) {
    fr <- fold_results[[fold]]

    oof_predictions$lp[fr$test_idx] <- fr$lp_test
    fold_models[[fold]] <- fr$fold_model
    names(fold_models)[fold] <- paste0("fold_", fold)

    lambda_values[fold] <- fr$lambda_min
    lambda_used[fold]   <- fr$lambda_used
    cvm_selected[fold]  <- fr$cvm_selected
    sd_lp_test[fold]    <- fr$sd_lp_test
    n_train_vec[fold]   <- fr$n_train
    n_test_vec[fold]    <- fr$n_test

    if (verbose) {
      cat(sprintf(
        "Fold %d: lambda.min=%s | lambda.used=%.6f | CV score=%.4f | sd(lp_test)=%.4g\n",
        fold,
        ifelse(is.finite(fr$lambda_min), sprintf("%.6f", fr$lambda_min), "Inf"),
        fr$lambda_used,
        fr$cvm_selected,
        fr$sd_lp_test
      ))
    }
  }

  list(
    oof_predictions = oof_predictions,
    fold_models = fold_models,
    lambda_values = lambda_values,
    lambda_used = lambda_used,
    fold_info = data.frame(
      fold = seq_len(outer_nfolds),
      n_train = n_train_vec,
      n_test = n_test_vec,
      lambda_min = lambda_values,
      lambda_used = lambda_used,
      cvm_selected = cvm_selected,
      sd_lp_test = sd_lp_test,
      stringsAsFactors = FALSE
    )
  )
}

# ============================================================================
# Engine: coxph
# ============================================================================
# Runs sequential outer folds (coxph is fast; parallel overhead not warranted).
run_oof_coxph <- function(
  dat,
  fold_ids,
  oof_predictions,
  predictors,
  time_col,
  status_col,
  vars_to_scale,
  outer_nfolds,
  seed,
  verbose
) {
  fold_models <- vector("list", outer_nfolds)
  sd_lp_test  <- rep(NA_real_, outer_nfolds)
  n_train_vec <- rep(NA_integer_, outer_nfolds)
  n_test_vec  <- rep(NA_integer_, outer_nfolds)

  for (fold in seq_len(outer_nfolds)) {
    if (verbose) {
      cat("\n--- Processing Fold", fold, "of", outer_nfolds, "(coxph) ---\n")
    }

    fr <- fit_one_fold_coxph(
      fold = fold,
      dat = dat,
      fold_ids = fold_ids,
      predictors = predictors,
      time_col = time_col,
      status_col = status_col,
      vars_to_scale = vars_to_scale,
      seed = seed,
      verbose = verbose
    )

    oof_predictions$lp[fr$test_idx] <- fr$lp_test
    fold_models[[fold]] <- fr$fold_model
    names(fold_models)[fold] <- paste0("fold_", fold)

    sd_lp_test[fold]  <- fr$sd_lp_test
    n_train_vec[fold] <- fr$n_train
    n_test_vec[fold]  <- fr$n_test

    if (verbose) {
      cat(sprintf(
        "Fold %d: sd(lp_test)=%.4g | n_train=%d | n_test=%d\n",
        fold, fr$sd_lp_test, fr$n_train, fr$n_test
      ))
    }
  }

  list(
    oof_predictions = oof_predictions,
    fold_models = fold_models,
    fold_info = data.frame(
      fold = seq_len(outer_nfolds),
      n_train = n_train_vec,
      n_test = n_test_vec,
      lambda_min = NA_real_,
      lambda_used = NA_real_,
      cvm_selected = NA_real_,
      sd_lp_test = sd_lp_test,
      stringsAsFactors = FALSE
    )
  )
}
