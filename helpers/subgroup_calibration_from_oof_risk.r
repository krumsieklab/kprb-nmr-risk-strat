# ============================================================================
# subgroup_calibration_from_oof_risk()
# ============================================================================
# Stratified calibration from PRECOMPUTED fully OOF absolute risks.
#
# Unlike subgroup_calibration() (helpers/subgroup_calibration.r), this helper
# does not refit coxph() or estimate a baseline hazard from linear predictors.
# It wraps compute_calibration_df_from_oof_risk() and runs it separately for
# each (stratum, model) using a supplied risk column (e.g. pred_risk_5yr).
#
# Per row we report:
#   calibration_score    : 1 - grouped ICI (higher is better)
#   ici                  : integrated calibration index (n-weighted MAE per group)
#   mean_abs_err         : unweighted group-level mean absolute calibration error
#   max_abs_err          : maximum absolute calibration error across groups
#
# Returns a LONG-format data.frame with columns:
#   strat_name, stratum, model, t_star, risk_col,
#   observed_risk_tstar, observed_risk_lower, observed_risk_upper,
#   calibration_score, ici, mean_abs_err, max_abs_err,
#   n, n_events, stratum_n_label
# ============================================================================

if (!exists(".subgroup_resolve_models", mode = "function")) {
  source(here::here("helpers", "subgroup_utils.r"))
}
if (!exists("compute_calibration_df_from_oof_risk", mode = "function")) {
  source(here::here("helpers", "calibration_from_oof_risk.r"))
}


# Extract an OOF frame that already contains the precomputed risk column.
# Accepts the same model-entry types as .subgroup_get_oof(), but does not
# require an `lp` column.
.subgroup_get_risk_oof <- function(model_arg, input_list, risk_col,
                                   arg_name = "model") {

  if (is.character(model_arg) && length(model_arg) == 1) {
    if (is.null(input_list)) {
      stop(sprintf(
        "%s is a string key ('%s') but disease_config$input_list is NULL.",
        arg_name, model_arg
      ))
    }
    if (!(model_arg %in% names(input_list))) {
      stop(sprintf(
        "%s key '%s' not found in input_list. Available: %s",
        arg_name, model_arg, paste(names(input_list), collapse = ", ")
      ))
    }
    model_arg <- input_list[[model_arg]]
  }

  if (is.data.frame(model_arg)) {
    oof <- model_arg
  } else if (is.list(model_arg) && is.data.frame(model_arg$oof_predictions)) {
    oof <- model_arg$oof_predictions
  } else {
    stop(sprintf(
      "%s must be a string key into input_list, a data.frame, or a list with $oof_predictions.",
      arg_name
    ))
  }

  required <- c("StudyID", "time", "status", risk_col)
  missing  <- setdiff(required, names(oof))
  if (length(missing) > 0) {
    stop(sprintf(
      "%s OOF data.frame is missing required columns: %s",
      arg_name, paste(missing, collapse = ", ")
    ))
  }

  oof
}


# Build a merged data.frame with:
#   StudyID, time, status, risk_<model1>, risk_<model2>, ..., <metadata>
#
# Subjects are the INTERSECTION across all models + metadata, so per-stratum
# subsets are identical across models.
.subgroup_build_df_from_risk <- function(disease_config, risk_col,
                                         verbose = TRUE) {

  if (is.null(disease_config$metadata_df) ||
      !is.data.frame(disease_config$metadata_df)) {
    stop("disease_config$metadata_df must be a data.frame.")
  }

  meta <- disease_config$metadata_df
  if (!("StudyID" %in% names(meta))) {
    stop("disease_config$metadata_df must contain a 'StudyID' column.")
  }

  meta   <- .subgroup_add_sex_binary_if_missing(meta)
  models <- .subgroup_resolve_models(disease_config)

  model_oofs <- lapply(seq_along(models), function(i) {
    label <- names(models)[i]
    oof   <- .subgroup_get_risk_oof(
      models[[i]], disease_config$input_list, risk_col, arg_name = label
    )
    risk_name <- paste0("risk_", label)
    oof %>%
      dplyr::filter(!is.na(.data[[risk_col]])) %>%
      dplyr::transmute(
        StudyID,
        time,
        status,
        !!risk_name := .data[[risk_col]]
      )
  })
  names(model_oofs) <- names(models)

  common_ids <- Reduce(intersect, c(
    lapply(model_oofs, function(x) x$StudyID),
    list(meta$StudyID)
  ))

  if (length(common_ids) == 0) {
    stop("No common StudyIDs across all models and metadata_df.")
  }

  first_label <- names(model_oofs)[1]
  merged <- model_oofs[[first_label]] %>%
    dplyr::filter(StudyID %in% common_ids)

  if (length(model_oofs) >= 2) {
    for (i in 2:length(model_oofs)) {
      label     <- names(model_oofs)[i]
      risk_name <- paste0("risk_", label)
      other     <- model_oofs[[i]] %>%
        dplyr::filter(StudyID %in% common_ids) %>%
        dplyr::select(StudyID, dplyr::all_of(risk_name))
      merged <- merged %>% dplyr::inner_join(other, by = "StudyID")
    }
  }

  merged <- merged %>% dplyr::inner_join(meta, by = "StudyID")

  if (verbose) {
    cat(sprintf(
      "  Built merged df from precomputed %s: n=%d, models=[%s]\n",
      risk_col, nrow(merged), paste(names(model_oofs), collapse = ", ")
    ))
  }

  if (nrow(merged) == 0) {
    stop("Merge produced 0 rows. Check StudyID compatibility between OOF risks and metadata_df.")
  }

  attr(merged, "model_names") <- names(model_oofs)
  attr(merged, "risk_col")    <- risk_col
  merged
}


subgroup_calibration_from_oof_risk <- function(
  disease_config,
  strata_configs = NULL,
  t_star         = 5,
  risk_col       = NULL,
  n_groups       = 10,
  min_n          = 50,
  min_events     = 10,
  verbose        = TRUE
) {

  if (is.null(risk_col)) {
    risk_col <- sprintf("pred_risk_%dyr", as.integer(round(t_star)))
  }

  strata_configs <- .subgroup_resolve_strata_arg(strata_configs, disease_config)
  models         <- .subgroup_resolve_models(disease_config)
  model_labels   <- names(models)

  if (verbose) {
    cat("\n============================================\n")
    cat("Subgroup Calibration from Precomputed OOF Risk\n")
    cat("============================================\n")
    if (!is.null(disease_config$label)) {
      cat(sprintf("Disease     : %s\n", disease_config$label))
    }
    cat(sprintf("Models      : %s\n", paste(model_labels, collapse = ", ")))
    cat(sprintf("t_star      : %.2f\n", t_star))
    cat(sprintf("risk_col    : %s\n", risk_col))
    cat(sprintf("n_groups    : %d\n", n_groups))
    cat(sprintf("Strata      : %s\n",
                paste(vapply(strata_configs, `[[`, character(1), "name"),
                      collapse = ", ")))
  }

  df <- .subgroup_build_df_from_risk(disease_config, risk_col, verbose = verbose)
  .subgroup_validate_strata(strata_configs, df)
  strata_configs <- .subgroup_resolve_strata(strata_configs, df)

  iter <- .subgroup_iter_strata(df, strata_configs,
                                min_n = min_n, min_events = min_events)

  if (length(iter) == 0) {
    warning("No strata met min_n / min_events thresholds; returning empty df.")
    return(data.frame())
  }

  pluck_attr <- function(x, a) {
    if (is.null(x)) return(NA_real_)
    val <- attr(x, a)
    if (is.null(val)) NA_real_ else val
  }

  stratum_km_risk <- function(slice) {
    km_formula <- survival::Surv(time, status) ~ 1
    km  <- survival::survfit(km_formula, data = slice)
    kms <- summary(km, times = t_star, extend = TRUE)

    if (length(kms$surv) == 0 || all(is.na(kms$surv))) {
      return(list(
        observed_risk_tstar = NA_real_,
        observed_risk_lower = NA_real_,
        observed_risk_upper = NA_real_
      ))
    }

    list(
      observed_risk_tstar = 1 - kms$surv[1],
      observed_risk_lower = 1 - kms$upper[1],
      observed_risk_upper = 1 - kms$lower[1]
    )
  }

  rows <- list()
  row_i <- 1

  for (item in iter) {
    df_sub  <- item$df_sub
    n_sub   <- nrow(df_sub)
    n_ev    <- sum(df_sub$status, na.rm = TRUE)
    km_risk <- stratum_km_risk(df_sub)

    for (label in model_labels) {
      risk_name <- paste0("risk_", label)
      slice <- df_sub %>%
        dplyr::transmute(
          StudyID,
          time,
          status,
          !!risk_col := .data[[risk_name]]
        )

      cal <- tryCatch(
        compute_calibration_df_from_oof_risk(
          slice,
          t_star   = t_star,
          risk_col = risk_col,
          n_groups = n_groups,
          verbose  = FALSE
        ),
        error = function(e) {
          if (verbose) {
            cat(sprintf("    [calib failed | %s | %s | %s]: %s\n",
                        item$strat_name, item$stratum, label,
                        conditionMessage(e)))
          }
          NULL
        }
      )

      rows[[row_i]] <- data.frame(
        strat_name          = item$strat_name,
        stratum             = item$stratum,
        model               = label,
        t_star              = t_star,
        risk_col            = risk_col,
        observed_risk_tstar = km_risk$observed_risk_tstar,
        observed_risk_lower = km_risk$observed_risk_lower,
        observed_risk_upper = km_risk$observed_risk_upper,
        calibration_score   = pluck_attr(cal, "calibration_score"),
        ici                 = pluck_attr(cal, "ici"),
        mean_abs_err        = pluck_attr(cal, "mean_abs_calib_error"),
        max_abs_err         = pluck_attr(cal, "max_abs_calib_error"),
        n                   = n_sub,
        n_events            = n_ev,
        stringsAsFactors    = FALSE
      )
      row_i <- row_i + 1
    }
  }

  out <- dplyr::bind_rows(rows)
  out <- .subgroup_finalize_df(out, strata_configs, model_levels = model_labels)

  if (verbose) {
    cat("\nResults:\n")
    print(out %>% dplyr::select(strat_name, stratum, model,
                                observed_risk_tstar,
                                calibration_score, ici,
                                mean_abs_err, n, n_events))
    cat("============================================\n\n")
  }

  out
}
