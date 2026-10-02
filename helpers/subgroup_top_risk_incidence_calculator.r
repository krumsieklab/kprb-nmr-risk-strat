# ============================================================================
# subgroup_top_risk_incidence_calculator()
# ============================================================================
# Compute stratified top-X% predicted-risk observed incidence for an arbitrary
# number of OOF Cox models across user-defined strata. NO delta computation.
#
# For each (stratum, model, top_prop):
#   - LP -> absolute risk at t_star (baseline hazard within the stratum)
#   - Rank descending by predicted risk
#   - Select top n_top = ceiling(n_stratum * top_prop)
#   - KM observed incidence at t_star within the top group
#   - Raw + IPCW event capture (top events / total events in stratum)
#
# Returns a LONG-format data.frame with one row per (stratum, model, top_prop):
#   strat_name, stratum, model, top_prop, top_label,
#   n, n_events, n_top, n_events_top_by_tstar,
#   obs_incidence, obs_lower, obs_upper,
#   event_capture_raw, event_capture_ipcw,
#   pred_risk_mean_top, pred_risk_median_top,
#   calibration_slope, baseline_risk_tstar,
#   stratum_n_label
# ============================================================================

if (!exists(".subgroup_build_df", mode = "function")) {
  source(here::here("helpers", "subgroup_utils.r"))
}
if (!exists("top_risk_incidence_calculator", mode = "function")) {
  source(here::here("helpers", "top_risk_incidence_calculator.r"))
}

subgroup_top_risk_incidence_calculator <- function(
  disease_config,
  strata_configs = NULL,
  top_props      = c(0.01, 0.05, 0.10),
  t_star         = 5,
  recalibrate    = FALSE,
  min_n          = 50,
  min_events     = 5,
  verbose        = TRUE
) {

  strata_configs <- .subgroup_resolve_strata_arg(strata_configs, disease_config)
  models         <- .subgroup_resolve_models(disease_config)
  model_labels   <- names(models)

  if (verbose) {
    cat("\n============================================\n")
    cat("Subgroup Top-Risk Incidence Analysis\n")
    cat("============================================\n")
    if (!is.null(disease_config$label)) {
      cat(sprintf("Disease     : %s\n", disease_config$label))
    }
    cat(sprintf("Models      : %s\n", paste(model_labels, collapse = ", ")))
    cat(sprintf("t_star      : %.2f\n", t_star))
    cat(sprintf("top_props   : %s\n", paste(top_props, collapse = ", ")))
    cat(sprintf("Recalibrate : %s\n", ifelse(recalibrate, "TRUE", "FALSE")))
    cat(sprintf("Strata      : %s\n",
                paste(vapply(strata_configs, `[[`, character(1), "name"),
                      collapse = ", ")))
  }

  df <- .subgroup_build_df(disease_config, verbose = verbose)
  .subgroup_validate_strata(strata_configs, df)
  strata_configs <- .subgroup_resolve_strata(strata_configs, df)

  iter <- .subgroup_iter_strata(df, strata_configs,
                                min_n = min_n, min_events = min_events)

  if (length(iter) == 0) {
    warning("No strata met min_n / min_events thresholds; returning empty df.")
    return(data.frame())
  }

  keep_cols <- c(
    "model_key", "top_prop", "top_label",
    "n_top", "n_events_top_by_tstar",
    "obs_incidence", "obs_lower", "obs_upper",
    "event_capture_raw", "event_capture_ipcw",
    "pred_risk_mean_top", "pred_risk_median_top",
    "calibration_slope", "baseline_risk_tstar"
  )

  rows <- list()
  row_i <- 1

  for (item in iter) {
    df_sub <- item$df_sub
    n_sub  <- nrow(df_sub)
    n_ev   <- sum(df_sub$status, na.rm = TRUE)

    # Build a per-model slice list (StudyID, time, status, lp) keyed by label.
    slices <- lapply(model_labels, function(label) {
      lp_col <- paste0("lp_", label)
      df_sub %>%
        dplyr::transmute(StudyID, time, status, lp = .data[[lp_col]])
    })
    names(slices) <- model_labels

    res <- tryCatch(
      top_risk_incidence_calculator(
        input                = slices,
        disease_name         = item$stratum,
        model_keys           = model_labels,
        model_labels         = model_labels,
        reference_model_key  = model_labels[1],  # ignored downstream; we drop deltas
        top_props            = top_props,
        t_star               = t_star,
        recalibrate          = recalibrate,
        return_wide          = FALSE,
        verbose              = FALSE
      ),
      error = function(e) {
        if (verbose) {
          cat(sprintf("    [top-risk failed | %s | %s]: %s\n",
                      item$strat_name, item$stratum, conditionMessage(e)))
        }
        NULL
      }
    )

    if (is.null(res)) next

    res_use <- as.data.frame(res)[, keep_cols, drop = FALSE]
    res_use$strat_name <- item$strat_name
    res_use$stratum    <- item$stratum
    res_use$n          <- n_sub
    res_use$n_events   <- n_ev
    res_use <- res_use %>% dplyr::rename(model = model_key)

    rows[[row_i]] <- res_use
    row_i <- row_i + 1
  }

  out <- dplyr::bind_rows(rows)

  if (nrow(out) == 0) {
    warning("No strata produced top-risk results; returning empty df.")
    return(out)
  }

  out <- out %>%
    dplyr::select(
      strat_name, stratum, model, top_prop, top_label,
      n, n_events, n_top, n_events_top_by_tstar,
      obs_incidence, obs_lower, obs_upper,
      event_capture_raw, event_capture_ipcw,
      pred_risk_mean_top, pred_risk_median_top,
      calibration_slope, baseline_risk_tstar
    )

  out <- .subgroup_finalize_df(out, strata_configs, model_levels = model_labels)

  if (nrow(out) > 0) {
    out$top_label <- factor(
      out$top_label,
      levels = unique(out$top_label[order(out$top_prop)])
    )
  }

  if (verbose) {
    cat("\nResults:\n")
    print(out %>% dplyr::select(strat_name, stratum, model, top_label,
                                obs_incidence, obs_lower, obs_upper,
                                n, n_events))
    cat("============================================\n\n")
  }

  out
}
