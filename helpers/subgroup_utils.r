# ============================================================================
# subgroup_utils.r
# ============================================================================
# Shared utilities used by:
#   - subgroup_cindex.r
#   - subgroup_calibration.r
#   - subgroup_calibration_from_oof_risk.r
#   - subgroup_top_risk_incidence_calculator.r
#
# A "disease_config" is a named list, e.g.:
#   t2d <- list(
#     label       = 'T2D',
#     metadata_df = meta,                  # data.frame with StudyID + strata vars
#     input_list  = t2d_list,              # collect_prs_oof_predictions(...) output
#
#     # ---- Option A: explicit named models list (preferred) ----
#     models = list(
#       "marker_SOC"      = "marker_SOC",      # label = key into input_list
#       "marker_SOC+NMR"  = "marker_SOC+NMR"
#     ),
#
#     # ---- Option B: legacy 2-model shorthand (still supported) ----
#     soc_model = 'marker_SOC',
#     nmr_model = 'marker_SOC+NMR',
#
#     pgs_id         = 'PGS002308',           # optional
#     strata_configs = strata_configs         # optional default strata config
#   )
#
# Each entry in `models` can be:
#   - a character string: looked up in `input_list[[key]]`
#   - a data.frame: used directly
#   - a list with $oof_predictions: $oof_predictions is used
#
# All three subgroup helpers return LONG-format tables with one row per
# (strat_name, stratum, model[, top_prop]) combination. No delta columns.
# ============================================================================


`%||%` <- function(a, b) if (is.null(a)) b else a


# ----------------------------------------------------------------------------
# Extract the OOF data.frame from a flexible model argument.
# ----------------------------------------------------------------------------
.subgroup_get_oof <- function(model_arg, input_list, arg_name = "model") {

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

  required <- c("time", "status", "lp")
  missing  <- setdiff(required, names(oof))
  if (length(missing) > 0) {
    stop(sprintf(
      "%s OOF data.frame is missing required columns: %s",
      arg_name, paste(missing, collapse = ", ")
    ))
  }

  oof
}


# ----------------------------------------------------------------------------
# Resolve disease_config -> named list of model entries (label -> key/df/list).
# Accepts:
#   - disease_config$models : named list/vector  (preferred)
#   - disease_config$soc_model + $nmr_model      (legacy 2-model shorthand)
# ----------------------------------------------------------------------------
.subgroup_resolve_models <- function(disease_config) {

  if (!is.null(disease_config$models)) {
    models <- disease_config$models

    if (is.null(names(models)) || any(names(models) == "")) {
      # If a bare character vector is supplied, use values as labels.
      if (is.character(models)) {
        models <- as.list(models)
        names(models) <- unlist(models)
      } else {
        stop("`disease_config$models` must be a NAMED list/vector (label -> key/df).")
      }
    }
    return(models)
  }

  if (!is.null(disease_config$soc_model) && !is.null(disease_config$nmr_model)) {
    label_soc <- if (is.character(disease_config$soc_model)) disease_config$soc_model else "SOC"
    label_nmr <- if (is.character(disease_config$nmr_model)) disease_config$nmr_model else "NMR"

    out <- list()
    out[[label_soc]] <- disease_config$soc_model
    out[[label_nmr]] <- disease_config$nmr_model
    return(out)
  }

  stop("disease_config must contain either `models` (named list) or `soc_model` + `nmr_model`.")
}


# ----------------------------------------------------------------------------
# Auto-derive sex_binary from `sex`/`gender` if missing.
# ----------------------------------------------------------------------------
.subgroup_add_sex_binary_if_missing <- function(meta) {
  if ("sex_binary" %in% names(meta)) return(meta)

  if ("sex" %in% names(meta)) {
    meta$sex_binary <- ifelse(
      as.character(meta$sex) %in% c("M", "1", "Male", "male", "m"), 1L, 0L
    )
  } else if ("gender" %in% names(meta)) {
    meta$sex_binary <- ifelse(
      as.character(meta$gender) %in% c("M", "1", "Male", "male", "m"), 1L, 0L
    )
  }

  meta
}


# ----------------------------------------------------------------------------
# Build a single merged data.frame with:
#   StudyID, time, status, lp_<model1>, lp_<model2>, ..., <metadata columns>
#
# Subjects are the INTERSECTION across all models + metadata, so per-stratum
# subsets are guaranteed identical across models.
#
# Returns the merged data.frame with attr `model_names` storing model labels.
# ----------------------------------------------------------------------------
.subgroup_build_df <- function(disease_config, verbose = TRUE) {

  if (is.null(disease_config$metadata_df) || !is.data.frame(disease_config$metadata_df)) {
    stop("disease_config$metadata_df must be a data.frame.")
  }

  meta <- disease_config$metadata_df
  if (!("StudyID" %in% names(meta))) {
    stop("disease_config$metadata_df must contain a 'StudyID' column.")
  }

  meta   <- .subgroup_add_sex_binary_if_missing(meta)
  models <- .subgroup_resolve_models(disease_config)

  # Resolve each model entry to a (StudyID, time, status, lp) data.frame.
  model_oofs <- lapply(seq_along(models), function(i) {
    label <- names(models)[i]
    oof   <- .subgroup_get_oof(models[[i]], disease_config$input_list, arg_name = label)
    if (!("StudyID" %in% names(oof))) {
      stop(sprintf("Model '%s' OOF data.frame must contain a 'StudyID' column.", label))
    }
    oof %>% dplyr::select(StudyID, time, status, lp)
  })
  names(model_oofs) <- names(models)

  # Take intersection of StudyIDs across all models + metadata.
  common_ids <- Reduce(intersect, c(
    lapply(model_oofs, function(x) x$StudyID),
    list(meta$StudyID)
  ))

  if (length(common_ids) == 0) {
    stop("No common StudyIDs across all models and metadata_df.")
  }

  # First model: keep time/status + rename lp -> lp_<label>.
  first_label <- names(model_oofs)[1]
  merged <- model_oofs[[first_label]] %>%
    dplyr::filter(StudyID %in% common_ids)
  names(merged)[names(merged) == "lp"] <- paste0("lp_", first_label)

  # Subsequent models: join lp_<label> only (assume time/status identical).
  if (length(model_oofs) >= 2) {
    for (i in 2:length(model_oofs)) {
      label <- names(model_oofs)[i]
      lp_col <- paste0("lp_", label)
      other  <- model_oofs[[i]] %>%
        dplyr::filter(StudyID %in% common_ids) %>%
        dplyr::transmute(StudyID, !!lp_col := lp)
      merged <- merged %>% dplyr::inner_join(other, by = "StudyID")
    }
  }

  # Merge with metadata.
  merged <- merged %>% dplyr::inner_join(meta, by = "StudyID")

  if (verbose) {
    cat(sprintf(
      "  Built merged df: n=%d, models=[%s]\n",
      nrow(merged), paste(names(model_oofs), collapse = ", ")
    ))
  }

  if (nrow(merged) == 0) {
    stop("Merge produced 0 rows. Check StudyID compatibility between OOF predictions and metadata_df.")
  }

  attr(merged, "model_names") <- names(model_oofs)
  merged
}


# ----------------------------------------------------------------------------
# Resolve PRS sentinel cutoffs ('median', 'p80', 'p20', 'p90', 'p10').
# ----------------------------------------------------------------------------
.subgroup_resolve_strata <- function(strata_configs, df) {

  lapply(strata_configs, function(s) {
    if (!is.null(s$cutoff) && is.character(s$cutoff)) {
      x <- df[[s$var]]
      if (is.null(x)) {
        stop(sprintf(
          "Strata config '%s' references variable '%s' which is not in the merged df.",
          s$name %||% s$var, s$var
        ))
      }
      cut_label <- s$cutoff
      s$cutoff <- switch(
        tolower(cut_label),
        "median" = stats::median(x, na.rm = TRUE),
        "p10"    = as.numeric(stats::quantile(x, 0.10, na.rm = TRUE)),
        "p20"    = as.numeric(stats::quantile(x, 0.20, na.rm = TRUE)),
        "p80"    = as.numeric(stats::quantile(x, 0.80, na.rm = TRUE)),
        "p90"    = as.numeric(stats::quantile(x, 0.90, na.rm = TRUE)),
        stop(sprintf(
          "Unrecognized sentinel cutoff '%s' for strata '%s'. Use 'median','p10','p20','p80','p90', or a numeric value.",
          cut_label, s$name %||% s$var
        ))
      )
    }
    s
  })
}


# ----------------------------------------------------------------------------
# Validate that all stratum variables exist in df.
# ----------------------------------------------------------------------------
.subgroup_validate_strata <- function(strata_configs, df) {

  missing_vars <- unique(unlist(lapply(strata_configs, function(s) {
    if (!(s$var %in% names(df))) s$var else NULL
  })))

  if (length(missing_vars) > 0) {
    stop(sprintf(
      "Strata config references variable(s) not in metadata_df: %s",
      paste(missing_vars, collapse = ", ")
    ))
  }

  invisible(TRUE)
}


# ----------------------------------------------------------------------------
# Iterate strata: for each stratum, return a data.frame subset + labels.
#
# Returns a list of lists with elements:
#   $strat_name : facet-level name (e.g. 'Age', 'Ethnicity')
#   $stratum    : within-stratum label (e.g. 'Age < 55', 'Asian')
#   $df_sub     : the subset data.frame
# ----------------------------------------------------------------------------
.subgroup_iter_strata <- function(df, strata_configs, min_n = 20, min_events = 5) {

  out <- list()
  idx <- 1

  for (strat in strata_configs) {

    if (!is.null(strat$cutoff)) {
      df_strat <- df %>%
        dplyr::filter(!is.na(.data[[strat$var]])) %>%
        dplyr::mutate(stratum = ifelse(
          .data[[strat$var]] <= strat$cutoff,
          strat$low, strat$high
        ))
      groups <- c(strat$low, strat$high)

    } else if (!is.null(strat$groups)) {
      df_strat <- df %>%
        dplyr::filter(
          !is.na(.data[[strat$var]]),
          as.character(.data[[strat$var]]) %in% names(strat$groups)
        ) %>%
        dplyr::mutate(stratum = unname(strat$groups[as.character(.data[[strat$var]])]))
      groups <- unname(strat$groups)

    } else {
      stop(sprintf(
        "Stratum '%s' must have either a 'cutoff' (binary) or 'groups' (categorical) field.",
        strat$name
      ))
    }

    for (g in groups) {
      df_sub <- df_strat %>% dplyr::filter(stratum == g)
      if (nrow(df_sub) < min_n || sum(df_sub$status, na.rm = TRUE) < min_events) {
        next
      }
      out[[idx]] <- list(
        strat_name = strat$name,
        stratum    = g,
        df_sub     = df_sub
      )
      idx <- idx + 1
    }
  }

  out
}


# ----------------------------------------------------------------------------
# Resolve `strata_configs` argument: explicit > disease_config$strata_configs.
# ----------------------------------------------------------------------------
.subgroup_resolve_strata_arg <- function(strata_configs, disease_config) {
  if (!is.null(strata_configs)) return(strata_configs)
  if (!is.null(disease_config$strata_configs)) return(disease_config$strata_configs)
  stop("No strata_configs provided. Pass `strata_configs` directly or set `disease_config$strata_configs`.")
}


# ----------------------------------------------------------------------------
# Sort the result df by strat_name (order in strata_configs) and model
# (order in models list), and add a convenience `stratum_n_label` column.
# ----------------------------------------------------------------------------
.subgroup_finalize_df <- function(df, strata_configs, model_levels = NULL) {
  if (nrow(df) == 0) return(df)

  strat_levels <- vapply(strata_configs, function(s) s$name, character(1))

  df <- df %>%
    dplyr::mutate(
      strat_name      = factor(strat_name, levels = strat_levels),
      stratum_n_label = paste0(stratum, " (n=", n, ")")
    )

  if (!is.null(model_levels) && "model" %in% names(df)) {
    df <- df %>%
      dplyr::mutate(model = factor(model, levels = model_levels)) %>%
      dplyr::arrange(strat_name, stratum, model)
  } else {
    df <- df %>% dplyr::arrange(strat_name, stratum)
  }

  df
}
