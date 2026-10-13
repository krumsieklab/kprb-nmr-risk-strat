## ============================================================
## Pre-enrollment A1C variability (SD / CV), leakage-safe

## Variability metrics (per patient, over the contributing points):
##   - mean
##   - sd
##   - cv = sd / mean      (level-normalized; preferred for A1C)
##   - n  = # of contributing measurements (after same-day collapse)
##   - recency_days = days-before-enrollment of the MOST RECENT
##                    contributing measurement (freshness of the
##                    trajectory's endpoint; the spread can be stale
##                    even when it is computable)
##
## Output uses the M / value / T (Missingness / Observed / Time)
## zero-gated encoding:
##   - M_<tag>_var   = 1 if variability is computable (n >= min_points
##                     and sd/mean are finite), else 0
##   - <OBS>_sd_obs  = sd if M==1 else 0
##   - <OBS>_cv_obs  = cv if M==1 else 0
##   - T_<tag>_var   = recency_days if M==1 else 0
## plus raw (un-gated) QC columns and the kept mean / N.
##
## Leakage safety: only measurements with
##   0 (or 1) <= days_before <= lookback_days
## are used, so nothing on/after enrollment leaks in. Multiple draws
## on the same day are collapsed to their median BEFORE counting, so
## a busy lab day does not inflate N (mirrors locf_exact_impute.r).
##
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(rlang)
})

# ------------------------------------------------------------
# Main: compute pre-enrollment variability for one analyte.
#
# Args
#   labs_df         long lab table: id_col, value_col, date_col
#   meta_df         enrollment snapshot: id_col, index_date_col
#   id_col          patient id column (default "StudyID")
#   value_col       measured value column in labs_df (default "Result_C")
#   date_col        measurement date column in labs_df (default "RESULT_DATE")
#   index_date_col  enrollment/index date column in meta_df
#                     (default "KPRB_SAMPLE_DATE")
#   lookback_days   max days before index to include (default 730 = 2y)
#   min_points      minimum contributing points required to call the
#                     metric computable (default 2; 3+ is more stable)
#   include_index_day  if TRUE, a same-day (days_before == 0) measurement
#                     counts as pre-enrollment (default TRUE)
#   value_range     c(min, max) plausibility filter applied to values;
#                     guards against sentinels (the raw A1C feed carries
#                     -1e100 codes). Set NULL to disable. Default c(3, 20).
#   tag             lowercase token for M_/T_ columns (default "a1c")
#   obs_label       UPPER token for *_sd_obs / *_cv_obs (default "A1C")
#
# Returns one row per StudyID present in meta_df:
#   <id>, n_<tag>_var, <OBS>_var_mean,
#   <OBS>_sd_raw, <OBS>_cv_raw, <tag>_var_recency_days_raw,
#   M_<tag>_var, <OBS>_sd_obs, <OBS>_cv_obs, T_<tag>_var
# left_join() this onto meta by id_col.
# ------------------------------------------------------------
compute_a1c_variability <- function(labs_df,
                                    meta_df,
                                    id_col            = "StudyID",
                                    value_col         = "Result_C",
                                    date_col          = "RESULT_DATE",
                                    index_date_col    = "KPRB_SAMPLE_DATE",
                                    lookback_days     = 730,
                                    min_points        = 2,
                                    include_index_day = TRUE,
                                    value_range       = c(3, 20),
                                    tag               = "a1c",
                                    obs_label         = "A1C") {

  # ---- validate ----
  if (!is.numeric(lookback_days) || length(lookback_days) != 1 ||
      is.na(lookback_days) || lookback_days < 0) {
    stop("`lookback_days` must be a single non-negative number.")
  }
  if (!is.numeric(min_points) || length(min_points) != 1 ||
      is.na(min_points) || min_points < 2) {
    stop("`min_points` must be a single number >= 2 (sd needs >= 2 points).")
  }
  for (col in c(id_col, value_col, date_col)) {
    if (!col %in% names(labs_df)) stop(sprintf("labs_df is missing column '%s'.", col))
  }
  for (col in c(id_col, index_date_col)) {
    if (!col %in% names(meta_df)) stop(sprintf("meta_df is missing column '%s'.", col))
  }

  value_sym <- rlang::sym(value_col)
  date_sym  <- rlang::sym(date_col)
  idx_sym   <- rlang::sym(index_date_col)
  id_sym    <- rlang::sym(id_col)

  lower_bound <- if (isTRUE(include_index_day)) 0 else 1

  # ---- output column names ----
  n_col       <- paste0("n_", tag, "_var")
  mean_col    <- paste0(obs_label, "_var_mean")
  sd_raw_col  <- paste0(obs_label, "_sd_raw")
  cv_raw_col  <- paste0(obs_label, "_cv_raw")
  rec_raw_col <- paste0(tag, "_var_recency_days_raw")
  m_col       <- paste0("M_", tag, "_var")
  sd_obs_col  <- paste0(obs_label, "_sd_obs")
  cv_obs_col  <- paste0(obs_label, "_cv_obs")
  t_col       <- paste0("T_", tag, "_var")

  # ---- 1. clean labs: numeric value, real dates, plausibility filter ----
  labs_clean <- labs_df %>%
    dplyr::transmute(
      !!id_col   := .data[[id_col]],
      value      = suppressWarnings(as.numeric(.data[[value_col]])),
      meas_date  = as.Date(.data[[date_col]])
    ) %>%
    dplyr::filter(!is.na(.data[[id_col]]), !is.na(value), !is.na(meas_date))

  if (!is.null(value_range)) {
    if (length(value_range) != 2 || any(is.na(value_range)) ||
        value_range[1] >= value_range[2]) {
      stop("`value_range` must be c(min, max) with min < max, or NULL.")
    }
    labs_clean <- labs_clean %>%
      dplyr::filter(value >= value_range[1], value <= value_range[2])
  }

  # ---- 2. collapse same-day measurements (median) so N is not inflated ----
  labs_daily <- labs_clean %>%
    dplyr::group_by(.data[[id_col]], meas_date) %>%
    dplyr::summarise(value = stats::median(value, na.rm = TRUE), .groups = "drop")

  # ---- 3. attach enrollment date, keep backward-only window ----
  meta_idx <- meta_df %>%
    dplyr::transmute(
      !!id_col   := .data[[id_col]],
      index_date = as.Date(.data[[index_date_col]])
    ) %>%
    dplyr::filter(!is.na(.data[[id_col]]), !is.na(index_date)) %>%
    dplyr::distinct(.data[[id_col]], .keep_all = TRUE)

  windowed <- labs_daily %>%
    dplyr::inner_join(meta_idx, by = id_col) %>%
    dplyr::mutate(days_before = as.numeric(index_date - meas_date)) %>%
    dplyr::filter(days_before >= lower_bound, days_before <= lookback_days)

  # ---- 4. per-patient variability over contributing points ----
  per_patient <- windowed %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::summarise(
      !!n_col      := dplyr::n(),
      !!mean_col   := mean(value, na.rm = TRUE),
      !!sd_raw_col := stats::sd(value, na.rm = TRUE),       # NA if n < 2
      # recency = days-before of the MOST RECENT contributing point
      !!rec_raw_col := min(days_before, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      # CV = SD / mean (level-normalized); guard against mean ~ 0
      !!cv_raw_col := dplyr::if_else(
        is.finite(.data[[mean_col]]) & .data[[mean_col]] != 0,
        .data[[sd_raw_col]] / .data[[mean_col]],
        NA_real_
      )
    )

  # ---- 5. assemble one row per meta patient + M/obs/T zero-gating ----
  out <- meta_idx %>%
    dplyr::select(dplyr::all_of(id_col)) %>%
    dplyr::left_join(per_patient, by = id_col) %>%
    dplyr::mutate(
      !!n_col := dplyr::coalesce(.data[[n_col]], 0L),
      # computable iff enough points AND finite spread/level
      !!m_col := dplyr::if_else(
        .data[[n_col]] >= min_points &
          is.finite(.data[[sd_raw_col]]) &
          is.finite(.data[[mean_col]]),
        1L, 0L
      )
    ) %>%
    dplyr::mutate(
      !!sd_obs_col := dplyr::if_else(.data[[m_col]] == 1L,
                                     as.numeric(.data[[sd_raw_col]]), 0),
      !!cv_obs_col := dplyr::if_else(.data[[m_col]] == 1L,
                                     as.numeric(.data[[cv_raw_col]]), 0),
      !!t_col      := dplyr::if_else(.data[[m_col]] == 1L,
                                     as.numeric(.data[[rec_raw_col]]), 0)
    )

  # tidy column order: id, N, mean, raw QC, then gated block
  out %>%
    dplyr::select(
      dplyr::all_of(c(
        id_col, n_col, mean_col,
        sd_raw_col, cv_raw_col, rec_raw_col,
        m_col, sd_obs_col, cv_obs_col, t_col
      ))
    )
}

# ------------------------------------------------------------
# Lightweight QC printer for the variability output.
# ------------------------------------------------------------
summarize_a1c_variability <- function(var_df,
                                      tag       = "a1c",
                                      obs_label = "A1C",
                                      label     = "A1C variability") {
  m_col      <- paste0("M_", tag, "_var")
  n_col      <- paste0("n_", tag, "_var")
  sd_obs_col <- paste0(obs_label, "_sd_obs")
  cv_obs_col <- paste0(obs_label, "_cv_obs")
  t_col      <- paste0("T_", tag, "_var")

  cat("\n=============================\n")
  cat("Summary:", label, "\n")
  cat("=============================\n")
  cat(sprintf("Patients total            : %d\n", nrow(var_df)))
  cat(sprintf("Computable (M==1)         : %d (%.1f%%)\n",
              sum(var_df[[m_col]] == 1L, na.rm = TRUE),
              100 * mean(var_df[[m_col]] == 1L, na.rm = TRUE)))

  cat("\nN contributing measurements (computable only):\n")
  print(summary(var_df[[n_col]][var_df[[m_col]] == 1L]))

  cat("\nSD (computable only):\n")
  print(summary(var_df[[sd_obs_col]][var_df[[m_col]] == 1L]))

  cat("\nCV (computable only):\n")
  print(summary(var_df[[cv_obs_col]][var_df[[m_col]] == 1L]))

  cat("\nRecency days of endpoint (computable only):\n")
  print(summary(var_df[[t_col]][var_df[[m_col]] == 1L]))

  invisible(NULL)
}
