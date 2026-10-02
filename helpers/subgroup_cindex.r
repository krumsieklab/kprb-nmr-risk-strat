# ============================================================================
# subgroup_cindex()
# ============================================================================
# Compute stratified discrimination (Harrell's C-index) for an arbitrary
# number of OOF Cox models across user-defined strata. NO delta computation.
#
# For each (stratum, model):
#   - Per-model C-index with 95% CI from survival::concordance variance
#
# Returns a LONG-format data.frame with columns:
#   strat_name, stratum, model,
#   c_index, se, ci_lower, ci_upper,
#   n, n_events, stratum_n_label
# ============================================================================

if (!exists(".subgroup_build_df", mode = "function")) {
  source(here::here("helpers", "subgroup_utils.r"))
}

subgroup_cindex <- function(
  disease_config,
  strata_configs = NULL,
  min_n          = 20,
  min_events     = 5,
  verbose        = TRUE
) {

  if (!requireNamespace("survival", quietly = TRUE)) {
    stop("Package 'survival' is required.")
  }

  strata_configs <- .subgroup_resolve_strata_arg(strata_configs, disease_config)
  models         <- .subgroup_resolve_models(disease_config)
  model_labels   <- names(models)

  if (verbose) {
    cat("\n============================================\n")
    cat("Subgroup C-index Analysis\n")
    cat("============================================\n")
    if (!is.null(disease_config$label)) {
      cat(sprintf("Disease : %s\n", disease_config$label))
    }
    cat(sprintf("Models  : %s\n", paste(model_labels, collapse = ", ")))
    cat(sprintf("Strata  : %s\n",
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

  rows <- list()
  row_i <- 1

  for (item in iter) {
    df_sub   <- item$df_sub
    surv_obj <- survival::Surv(df_sub$time, df_sub$status)
    n_sub    <- nrow(df_sub)
    n_ev     <- sum(df_sub$status, na.rm = TRUE)

    for (label in model_labels) {
      lp_col <- paste0("lp_", label)
      lp_vec <- df_sub[[lp_col]]

      conc <- survival::concordance(surv_obj ~ I(-lp_vec))
      cidx <- conc$concordance
      se   <- sqrt(conc$var)

      rows[[row_i]] <- data.frame(
        strat_name = item$strat_name,
        stratum    = item$stratum,
        model      = label,
        c_index    = cidx,
        se         = se,
        ci_lower   = cidx - 1.96 * se,
        ci_upper   = cidx + 1.96 * se,
        n          = n_sub,
        n_events   = n_ev,
        stringsAsFactors = FALSE
      )
      row_i <- row_i + 1
    }
  }

  out <- dplyr::bind_rows(rows)
  out <- .subgroup_finalize_df(out, strata_configs, model_levels = model_labels)

  if (verbose) {
    cat("\nResults:\n")
    print(out %>% dplyr::select(strat_name, stratum, model,
                                c_index, ci_lower, ci_upper,
                                n, n_events))
    cat("============================================\n\n")
  }

  out
}
