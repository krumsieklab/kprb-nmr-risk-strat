# ============================================================================
# Transport a frozen (UKBB-trained) Cox model to an external cohort (KPRB)
#
# Purpose
#   Take a model saved by fit_full_cox() (see helpers/fit_full_cox.r) and apply
#   its FROZEN coefficients to a new cohort to produce linear predictors (LPs),
#   without any refitting.  This is the UKBB -> KPRB external-validation step
#   that mirrors the inputs used to compute the KPRB internal-OOF C-index
#   (StudyID, time, status, lp).
#
# Why a wrapper around predict_lp_full()/external_cindex()?
#   The two cohorts encode the SAME concepts under DIFFERENT column names, e.g.
#     UKBB predictor   <-  KPRB source column
#     ------------------------------------------------
#     bmi              <-  BMI_interp
#     sbp              <-  SBP_interp
#     hba1c            <-  HGBA1C_new        (raw observed value; NA if missing)
#     Alanine          <-  Ala               (only NMR name that differs)
#   Additionally, KPRB carries missingness-aware blocks (M_a1c, A1C_obs, T_a1c)
#   that DO NOT exist in the UKBB models.  The UKBB models use plain observed
#   values, so here we simply rename the KPRB observed columns to the UKBB
#   predictor names and let predict_lp_full() complete-case filter on them.
#
# Standardization (project summary sections 3 & 6.2)
#   NMR / markers are standardized WITHIN cohort before transport (serum/plasma
#   + abundance shifts), so the default is scaling = "self" (recompute mean/sd
#   on the external cohort).  Use scaling = "train" for a strict frozen-scaling
#   external validation.
#
# DEPENDENCIES
#   Source helpers/cv_glmnet_oof.r and helpers/fit_full_cox.r first; this
#   file reuses predict_lp_full() and (optionally) external_cindex() from them.
# ============================================================================

library(survival)

.transport_required_helpers <- c("predict_lp_full")
.transport_missing <- .transport_required_helpers[
  !vapply(.transport_required_helpers, exists, logical(1))
]
if (length(.transport_missing) > 0) {
  stop(
    "transport_lp.r requires helpers from fit_full_cox.r. Missing: ",
    paste(.transport_missing, collapse = ", "),
    ". Please source('helpers/fit_full_cox.r') (and cv_glmnet_oof.r) first."
  )
}
rm(.transport_required_helpers, .transport_missing)

# ----------------------------------------------------------------------------
# Rename external-cohort columns to the fit's predictor names.
# ----------------------------------------------------------------------------
#' @param newdata     external-cohort data.frame (e.g. KPRB)
#' @param rename_map  named character vector mapping
#'                    c(<fit_predictor_name> = <newdata_column_name>), e.g.
#'                    c(bmi = "BMI_interp", sbp = "SBP_interp",
#'                      hba1c = "HGBA1C_new", Alanine = "Ala").
#'                    Only include predictors whose external name differs; any
#'                    predictor already present under its fit name is left as-is.
#' @param overwrite   if FALSE (default) refuse to clobber an existing column
#'                    with a different meaning; set TRUE to allow overwrite.
align_external_predictors <- function(newdata, rename_map = NULL,
                                      overwrite = FALSE) {
  if (!is.data.frame(newdata)) stop("newdata must be a data.frame")
  if (is.null(rename_map)) return(newdata)
  if (is.null(names(rename_map)) || any(!nzchar(names(rename_map)))) {
    stop("rename_map must be a NAMED vector: c(fit_predictor = newdata_column).")
  }

  df <- newdata
  for (target in names(rename_map)) {
    src <- rename_map[[target]]
    if (!src %in% names(df)) {
      stop("rename_map source column '", src, "' not found in newdata.")
    }
    if (target %in% names(df) && !identical(target, src) && !overwrite) {
      stop("Target column '", target, "' already exists in newdata; ",
           "set overwrite = TRUE to replace it via rename_map.")
    }
    df[[target]] <- df[[src]]
  }
  df
}

# ----------------------------------------------------------------------------
# Score ONE frozen model on an external cohort -> tidy LP data.frame (+ C-index)
# ----------------------------------------------------------------------------
#' @param fit             a fit_full_cox() object (frozen model)
#' @param newdata         external cohort with StudyID/time/status + predictors
#'                        (predictors may be under external names; see rename_map)
#' @param rename_map      see align_external_predictors()
#' @param scaling         "self" (default, within-cohort), "self_nmr"
#'                        (training scaling except self-scaled NMR), "train",
#'                        or "none"
#' @param nmr_vars        NMR predictor names used by `scaling = "self_nmr"`
#' @param drop_incomplete drop rows missing any predictor (default TRUE)
#' @param compute_cindex  also compute the transported C-index (default TRUE)
#' @param verbose         print progress (default TRUE)
#'
#' @return list with:
#'   $scored    data.frame(StudyID, time, status, lp)   [names per fit$*_col]
#'   $c_index, $se, $ci_lower, $ci_upper, $n, $n_events  (if compute_cindex)
#'   $scaling_used, $n_predictors
transport_frozen_lp <- function(fit, newdata,
                                rename_map = NULL,
                                scaling = c("self", "self_nmr", "train", "none"),
                                nmr_vars = NULL,
                                drop_incomplete = TRUE,
                                compute_cindex = TRUE,
                                verbose = TRUE) {
  scaling <- match.arg(scaling)

  df <- align_external_predictors(newdata, rename_map)

  missing_preds <- setdiff(fit$predictors, names(df))
  if (length(missing_preds) > 0) {
    stop("newdata is missing ", length(missing_preds),
         " predictor(s) after rename: ",
         paste(utils::head(missing_preds, 10), collapse = ", "),
         if (length(missing_preds) > 10) " ..." else "")
  }

  scored <- predict_lp_full(
    fit, df,
    drop_incomplete = drop_incomplete,
    verbose = verbose,
    scaling = scaling,
    nmr_vars = nmr_vars
  )

  out <- list(
    scored        = scored,
    scaling_used  = scaling,
    n_predictors  = length(fit$predictors)
  )

  time_col   <- fit$time_col
  status_col <- fit$status_col
  have_outcome <- all(c(time_col, status_col) %in% names(scored))

  if (compute_cindex && have_outcome) {
    surv_obj <- survival::Surv(scored[[time_col]], scored[[status_col]])
    cc <- survival::concordance(surv_obj ~ I(-scored$lp))
    se <- sqrt(cc$var)
    out$c_index  <- cc$concordance
    out$se       <- se
    out$ci_lower <- cc$concordance - 1.96 * se
    out$ci_upper <- cc$concordance + 1.96 * se
    out$n        <- nrow(scored)
    out$n_events <- sum(scored[[status_col]])

    if (verbose) {
      cat(sprintf(
        "  Transported C-index: %.4f [%.4f, %.4f] | n=%d, events=%d | scaling=%s\n",
        out$c_index, out$ci_lower, out$ci_upper, out$n, out$n_events, scaling
      ))
    }
  } else if (compute_cindex && verbose) {
    cat("  (time/status not both present in scored output; skipping C-index)\n")
  }

  out
}

# ----------------------------------------------------------------------------
# Shared complete-case subset across multiple frozen models (fair comparisons)
# ----------------------------------------------------------------------------
#' Keep rows complete on the UNION of predictors required by every spec in
#' `models[specs]`, plus outcome/id columns.  Mirrors
#' `filter_complete_cases_for_specs()` in the KPRB/UKBB training notebooks so
#' BASE / CLIN / MARKER / *_NMR transports are evaluated on the same people.
#'
#' @return list with $data (filtered frame), $n_total, $n_complete, $n_dropped,
#'   $predictors_union, $study_ids
subset_shared_complete_cases <- function(newdata, models,
                                         rename_map = NULL,
                                         specs = names(models),
                                         extra_vars = NULL,
                                         verbose = TRUE) {
  if (is.null(names(models))) stop("`models` must be a named list of fits.")
  missing_specs <- setdiff(specs, names(models))
  if (length(missing_specs) > 0) {
    stop("Unknown spec(s): ", paste(missing_specs, collapse = ", "))
  }

  df <- align_external_predictors(newdata, rename_map)

  all_preds <- unique(unlist(
    lapply(models[specs], function(fit) fit$predictors),
    use.names = FALSE
  ))

  fit_ref <- models[[specs[[1]]]]
  outcome_cols <- unique(c(
    fit_ref$studyid_col, fit_ref$time_col, fit_ref$status_col, extra_vars
  ))
  cols_needed <- unique(c(outcome_cols, all_preds))

  missing_cols <- setdiff(cols_needed, names(df))
  if (length(missing_cols) > 0) {
    stop("newdata is missing ", length(missing_cols),
         " column(s) needed for the shared complete-case filter: ",
         paste(utils::head(missing_cols, 10), collapse = ", "),
         if (length(missing_cols) > 10) " ..." else "")
  }

  cc <- stats::complete.cases(df[, cols_needed, drop = FALSE])
  n_total <- nrow(df)
  n_complete <- sum(cc)
  n_dropped <- n_total - n_complete

  if (n_complete == 0L) {
    stop("No shared complete cases remain across specs ",
         paste(specs, collapse = ", "), ".")
  }

  if (verbose) {
    cat(sprintf(
      "Shared complete-case filter: %d / %d rows kept (%d dropped)\n",
      n_complete, n_total, n_dropped
    ))
    cat(sprintf(
      "  Union of %d predictor(s) across %d spec(s)\n",
      length(all_preds), length(specs)
    ))
  }

  list(
    data             = df[cc, , drop = FALSE],
    n_total          = n_total,
    n_complete       = n_complete,
    n_dropped        = n_dropped,
    predictors_union = all_preds,
    study_ids        = df[[fit_ref$studyid_col]][cc]
  )
}

# ----------------------------------------------------------------------------
# Score ALL specs of one outcome (a named list of frozen models) on a cohort.
# ----------------------------------------------------------------------------
#' Applies transport_frozen_lp() to each frozen model in `models`.
#'
#' By default (`shared_complete_cases = TRUE`) all specs are scored on the SAME
#' participants: rows complete on the union of predictors across every spec
#' (e.g. marker + NMR covariates for T2D).  Set `shared_complete_cases = FALSE`
#' to complete-case filter per spec independently (N may differ).
#'
#' @param models      named list of fit_full_cox() objects, e.g. output of
#'                    load_ukbb_frozen_models("t2d_a1c").  List names become the
#'                    keys of the returned per-spec results.
#' @param newdata     external cohort (KPRB) frame
#' @param rename_map  see align_external_predictors() (shared across specs; the
#'                    superset of all cohort->fit renames is fine)
#' @param specs       which model list entries to score (default: all)
#' @param shared_complete_cases  if TRUE (default), subset once to rows complete
#'                    on the union of all spec predictors before scoring
#' @param scaling,nmr_vars,drop_incomplete,verbose forwarded to
#'                    transport_frozen_lp()
#'                    (`drop_incomplete` applies only when
#'                    `shared_complete_cases = FALSE`)
#'
#' @return list with:
#'   $per_spec  named list; each element is the transport_frozen_lp() result
#'   $lp_tables named list of tidy data.frames (StudyID, time, status, lp)
#'   $summary   data.frame(spec, n, n_events, c_index, ci_lower, ci_upper)
#'   $complete_case_info  output of subset_shared_complete_cases() when shared
transport_all_specs <- function(models, newdata,
                                rename_map = NULL,
                                specs = names(models),
                                shared_complete_cases = TRUE,
                                scaling = c("self", "self_nmr", "train", "none"),
                                nmr_vars = NULL,
                                drop_incomplete = TRUE,
                                verbose = TRUE) {
  scaling <- match.arg(scaling)
  if (is.null(names(models))) stop("`models` must be a named list of fits.")

  complete_case_info <- NULL
  scoring_data <- newdata
  scoring_rename_map <- rename_map

  if (shared_complete_cases) {
    complete_case_info <- subset_shared_complete_cases(
      newdata     = newdata,
      models      = models,
      rename_map  = rename_map,
      specs       = specs,
      verbose     = verbose
    )
    scoring_data <- complete_case_info$data
    scoring_rename_map <- NULL
  }

  per_spec  <- vector("list", length(specs)); names(per_spec)  <- specs
  lp_tables <- vector("list", length(specs)); names(lp_tables) <- specs
  summ <- vector("list", length(specs))

  for (sp in specs) {
    if (verbose) cat("\n== spec:", sp, "==\n")
    res <- transport_frozen_lp(
      fit             = models[[sp]],
      newdata         = scoring_data,
      rename_map      = scoring_rename_map,
      scaling         = scaling,
      nmr_vars        = nmr_vars,
      drop_incomplete = if (shared_complete_cases) FALSE else drop_incomplete,
      compute_cindex  = TRUE,
      verbose         = verbose
    )
    per_spec[[sp]]  <- res
    lp_tables[[sp]] <- res$scored
    summ[[sp]] <- data.frame(
      spec     = sp,
      n        = if (!is.null(res$n)) res$n else nrow(res$scored),
      n_events = if (!is.null(res$n_events)) res$n_events else NA_integer_,
      c_index  = if (!is.null(res$c_index)) res$c_index else NA_real_,
      ci_lower = if (!is.null(res$ci_lower)) res$ci_lower else NA_real_,
      ci_upper = if (!is.null(res$ci_upper)) res$ci_upper else NA_real_,
      stringsAsFactors = FALSE
    )
  }

  summary_df <- do.call(rbind, summ)
  if (shared_complete_cases) {
    n_vals <- unique(summary_df$n)
    if (length(n_vals) != 1L) {
      stop("Shared complete-case scoring failed: specs have different N (",
           paste(n_vals, collapse = ", "), ").")
    }
    if (verbose) {
      cat(sprintf(
        "\nAll %d spec(s) scored on the same %d participant(s) (%d events).\n",
        length(specs), n_vals, unique(summary_df$n_events)
      ))
    }
  }

  list(
    per_spec           = per_spec,
    lp_tables          = lp_tables,
    summary            = summary_df,
    complete_case_info = complete_case_info
  )
}
