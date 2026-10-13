## ============================================================
## Pre-enrollment eGFR slope (renal function decline), leakage-safe
## ------------------------------------------------------------
## Estimates each patient's eGFR trajectory slope from creatinine
## measurements taken STRICTLY BEFORE enrollment, inside a bounded
## lookback window. NO imputation: only observed values are used.
##
## Pipeline (per patient):
##   1. creatinine -> eGFR FIRST (2021 CKD-EPI, race-free), using
##      age_at_result (age at the RESULT_DATE, NOT at enrollment).
##      Slope must be computed on eGFR, not creatinine: creatinine
##      slope has the wrong sign and nonlinear scaling.
##   2. regress eGFR on time (in years) over pre-enrollment points;
##      slope = coefficient on time, units mL/min/1.73m^2 per YEAR.
##   3. M / value / T zero-gated encoding for downstream models.
##
## SIGN CONVENTION (important):
##   time axis = calendar-forward years relative to enrollment
##   (yrs = (meas_date - enroll_date) / 365.25, so pre-enrollment
##   points are negative and enrollment = 0). Therefore:
##       NEGATIVE slope = DECLINING renal function (eGFR falling
##       as time approaches enrollment) -- the clinically intuitive
##       direction. A positive slope means improving eGFR.
##   (Regressing on "years_before_enroll" instead would flip the
##   sign; we use forward time so decline reads as negative.)
##
## Leakage safety: STRICT pre-enrollment only (days_before > 0, i.e.
## meas_date < enroll_date). Same-day and post-enrollment values are
## dropped. Multiple draws on the same day are collapsed to their
## median before counting, so a busy lab day does not inflate N.
##
## Lookback window: trajectories need a longer window than point
## estimates to accumulate enough draws for a stable slope. Default
## 5 years (1826 days). State/adjust via `lookback_days`.
##
## Slope is fit with the closed-form OLS estimator
##   slope = sum((t-tbar)(y-ybar)) / sum((t-tbar)^2)
## inside dplyr -- much faster than per-group lm() on a 55K x ~1.2M
## row feed, and numerically identical to coef(lm(egfr ~ yrs))[2].
##
## Save as: helpers/compute_egfr_slope.r
## Use via: source(here('helpers', 'compute_egfr_slope.r'))
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(rlang)
})

# ------------------------------------------------------------
# 2021 CKD-EPI creatinine eGFR (race-free), vectorized.
#   scr      serum creatinine (mg/dL)
#   age      age in years (use age AT THE MEASUREMENT for longitudinal)
#   sex_num  1 = male, 0 = female
# Matches the single-value equation already used in the notebook.
# ------------------------------------------------------------
ckd_epi_2021_egfr <- function(scr, age, sex_num) {
  kappa <- ifelse(sex_num == 0, 0.7, 0.9)
  alpha <- ifelse(sex_num == 0, -0.241, -0.302)
  142 *
    pmin(scr / kappa, 1)^alpha *
    pmax(scr / kappa, 1)^(-1.200) *
    0.9938^age *
    ifelse(sex_num == 0, 1.012, 1)
}

# ------------------------------------------------------------
# Resolve a sex/gender column to sex_num (1 = male, 0 = female).
# Accepts factor/character/numeric codings ('M'/'F', 1/0, etc.).
# ------------------------------------------------------------
.resolve_sex_num <- function(x) {
  xc <- toupper(trimws(as.character(x)))
  ifelse(xc %in% c("M", "1", "MALE", "TRUE"), 1L,
         ifelse(xc %in% c("F", "0", "FEMALE", "FALSE"), 0L, NA_integer_))
}

# ------------------------------------------------------------
# Closed-form OLS slope of y on x; NA if < 2 points or no x-variation.
# ------------------------------------------------------------
.ols_slope <- function(x, y) {
  ok <- is.finite(x) & is.finite(y)
  x <- x[ok]; y <- y[ok]
  n <- length(x)
  if (n < 2) return(NA_real_)
  xc <- x - mean(x)
  denom <- sum(xc * xc)
  if (!is.finite(denom) || denom <= 0) return(NA_real_)  # no time spread
  sum(xc * (y - mean(y))) / denom
}

# ------------------------------------------------------------
# Build the per-measurement longitudinal eGFR table (pre-enrollment).
# Returns one row per (patient, measurement day) with eGFR + timing.
# Useful on its own for QC / plotting trajectories.
#
# Args
#   labs_df         long creatinine table: id_col, value_col, date_col
#   meta_df         enrollment snapshot: id_col, index_date_col,
#                     age_col (age at enrollment), sex_col
#   value_col       creatinine column (mg/dL), default "Result_C"
#   date_col        measurement date, default "RESULT_DATE"
#   index_date_col  enrollment date, default "KPRB_SAMPLE_DATE"
#   age_col         age at enrollment, default "age_at_sample"
#   sex_col         sex/gender column in meta_df, default "sex"
#   lookback_days   max days before enrollment, default 1826 (5y)
#   creat_range     plausibility filter for creatinine; guards against
#                     sentinels/errors (raw feed has 0 and 333). Default
#                     c(0.1, 20). Set NULL to disable.
#
# Columns: <id>, meas_date, days_before, yrs_to_enroll (<=... < 0),
#          creatinine, age_at_result, sex_num, egfr
# ------------------------------------------------------------
compute_longitudinal_egfr <- function(labs_df,
                                      meta_df,
                                      id_col         = "StudyID",
                                      value_col      = "Result_C",
                                      date_col       = "RESULT_DATE",
                                      index_date_col = "KPRB_SAMPLE_DATE",
                                      age_col        = "age_at_sample",
                                      sex_col        = "sex",
                                      lookback_days  = 1826,
                                      creat_range    = c(0.1, 20)) {

  for (col in c(id_col, value_col, date_col)) {
    if (!col %in% names(labs_df)) stop(sprintf("labs_df is missing column '%s'.", col))
  }
  for (col in c(id_col, index_date_col, age_col, sex_col)) {
    if (!col %in% names(meta_df)) stop(sprintf("meta_df is missing column '%s'.", col))
  }
  if (!is.numeric(lookback_days) || length(lookback_days) != 1 ||
      is.na(lookback_days) || lookback_days <= 0) {
    stop("`lookback_days` must be a single positive number.")
  }

  # ---- clean creatinine: numeric, real dates, plausibility filter ----
  labs_clean <- labs_df %>%
    dplyr::transmute(
      !!id_col  := .data[[id_col]],
      creat     = suppressWarnings(as.numeric(.data[[value_col]])),
      meas_date = as.Date(.data[[date_col]])
    ) %>%
    dplyr::filter(!is.na(.data[[id_col]]), !is.na(creat), !is.na(meas_date))

  if (!is.null(creat_range)) {
    if (length(creat_range) != 2 || any(is.na(creat_range)) ||
        creat_range[1] >= creat_range[2]) {
      stop("`creat_range` must be c(min, max) with min < max, or NULL.")
    }
    labs_clean <- labs_clean %>%
      dplyr::filter(creat >= creat_range[1], creat <= creat_range[2])
  }

  # ---- collapse same-day measurements (median) so N is not inflated ----
  labs_daily <- labs_clean %>%
    dplyr::group_by(.data[[id_col]], meas_date) %>%
    dplyr::summarise(creat = stats::median(creat, na.rm = TRUE), .groups = "drop")

  # ---- enrollment anchors (date, age, sex) ----
  meta_idx <- meta_df %>%
    dplyr::transmute(
      !!id_col       := .data[[id_col]],
      index_date     = as.Date(.data[[index_date_col]]),
      age_at_sample  = as.numeric(.data[[age_col]]),
      sex_num        = .resolve_sex_num(.data[[sex_col]])
    ) %>%
    dplyr::filter(!is.na(.data[[id_col]]), !is.na(index_date),
                  !is.na(age_at_sample), !is.na(sex_num)) %>%
    dplyr::distinct(.data[[id_col]], .keep_all = TRUE)

  # ---- strict pre-enrollment window + per-measurement eGFR ----
  labs_daily %>%
    dplyr::inner_join(meta_idx, by = id_col) %>%
    dplyr::mutate(days_before = as.numeric(index_date - meas_date)) %>%
    dplyr::filter(days_before > 0, days_before <= lookback_days) %>%  # STRICT <
    dplyr::mutate(
      age_at_result = age_at_sample - days_before / 365.25,
      yrs_to_enroll = as.numeric(meas_date - index_date) / 365.25,    # < 0, forward time
      egfr          = ckd_epi_2021_egfr(creat, age_at_result, sex_num)
    ) %>%
    dplyr::transmute(
      !!id_col := .data[[id_col]],
      meas_date, days_before, yrs_to_enroll,
      creatinine = creat, age_at_result, sex_num, egfr
    ) %>%
    dplyr::arrange(.data[[id_col]], meas_date)
}

# ------------------------------------------------------------
# Main: per-patient pre-enrollment eGFR slope with M/value/T encoding.
#
# Extra args (beyond compute_longitudinal_egfr):
#   min_points   minimum contributing points to call slope computable.
#                Default 3 (>= 2 is the floor; 3+ strongly preferred --
#                2 points is just a line through noise).
#   tag          lowercase token for M_/T_ columns (default "egfr")
#   obs_label    UPPER token for the slope value column (default "EGFR")
#
# Returns one row per StudyID present in meta_df:
#   <id>, n_egfr_slope, egfr_slope_raw, egfr_slope_recency_days_raw,
#   M_egfr_slope, EGFR_slope_obs, T_egfr_slope
# left_join() this onto meta by id_col.
# ------------------------------------------------------------
compute_egfr_slope <- function(labs_df,
                               meta_df,
                               id_col         = "StudyID",
                               value_col      = "Result_C",
                               date_col       = "RESULT_DATE",
                               index_date_col = "KPRB_SAMPLE_DATE",
                               age_col        = "age_at_sample",
                               sex_col        = "sex",
                               lookback_days  = 1826,   # 5 years
                               min_points     = 3,
                               creat_range    = c(0.1, 20),
                               tag            = "egfr",
                               obs_label      = "EGFR") {

  if (!is.numeric(min_points) || length(min_points) != 1 ||
      is.na(min_points) || min_points < 2) {
    stop("`min_points` must be a single number >= 2 (slope needs >= 2 points).")
  }

  # output column names
  n_col       <- paste0("n_", tag, "_slope")
  slope_raw   <- paste0(tag, "_slope_raw")
  rec_raw_col <- paste0(tag, "_slope_recency_days_raw")
  m_col       <- paste0("M_", tag, "_slope")
  slope_obs   <- paste0(obs_label, "_slope_obs")
  t_col       <- paste0("T_", tag, "_slope")

  # per-measurement eGFR (strict pre-enrollment, windowed)
  meas <- compute_longitudinal_egfr(
    labs_df        = labs_df,
    meta_df        = meta_df,
    id_col         = id_col,
    value_col      = value_col,
    date_col       = date_col,
    index_date_col = index_date_col,
    age_col        = age_col,
    sex_col        = sex_col,
    lookback_days  = lookback_days,
    creat_range    = creat_range
  )

  # per-patient slope (closed-form OLS on forward-time years)
  per_patient <- meas %>%
    dplyr::group_by(.data[[id_col]]) %>%
    dplyr::summarise(
      !!n_col       := dplyr::n(),
      !!slope_raw   := .ols_slope(yrs_to_enroll, egfr),
      # recency = days-before of the MOST RECENT contributing point
      !!rec_raw_col := min(days_before, na.rm = TRUE),
      .groups = "drop"
    )

  # skeleton over ALL meta patients + M/value/T zero-gating
  meta_ids <- meta_df %>%
    dplyr::filter(!is.na(.data[[id_col]])) %>%
    dplyr::distinct(.data[[id_col]])

  meta_ids %>%
    dplyr::left_join(per_patient, by = id_col) %>%
    dplyr::mutate(
      !!n_col := dplyr::coalesce(.data[[n_col]], 0L),
      !!m_col := dplyr::if_else(
        .data[[n_col]] >= min_points & is.finite(.data[[slope_raw]]),
        1L, 0L
      )
    ) %>%
    dplyr::mutate(
      !!slope_obs := dplyr::if_else(.data[[m_col]] == 1L,
                                    as.numeric(.data[[slope_raw]]), 0),
      !!t_col     := dplyr::if_else(.data[[m_col]] == 1L,
                                    as.numeric(.data[[rec_raw_col]]), 0)
    ) %>%
    dplyr::select(
      dplyr::all_of(c(
        id_col, n_col, slope_raw, rec_raw_col,
        m_col, slope_obs, t_col
      ))
    )
}

# ------------------------------------------------------------
# Lightweight QC printer for the slope output.
# ------------------------------------------------------------
summarize_egfr_slope <- function(slope_df,
                                 tag       = "egfr",
                                 obs_label = "EGFR",
                                 label     = "eGFR slope") {
  n_col     <- paste0("n_", tag, "_slope")
  m_col     <- paste0("M_", tag, "_slope")
  slope_obs <- paste0(obs_label, "_slope_obs")
  t_col     <- paste0("T_", tag, "_slope")

  cat("\n=============================\n")
  cat("Summary:", label, "\n")
  cat("=============================\n")
  cat(sprintf("Patients total            : %d\n", nrow(slope_df)))
  cat(sprintf("Computable (M==1)         : %d (%.1f%%)\n",
              sum(slope_df[[m_col]] == 1L, na.rm = TRUE),
              100 * mean(slope_df[[m_col]] == 1L, na.rm = TRUE)))

  cat("\nN contributing measurements (computable only):\n")
  print(summary(slope_df[[n_col]][slope_df[[m_col]] == 1L]))

  cat("\nSlope mL/min/1.73m^2/yr (computable only; negative = decline):\n")
  print(summary(slope_df[[slope_obs]][slope_df[[m_col]] == 1L]))

  cat("\nRecency days of endpoint (computable only):\n")
  print(summary(slope_df[[t_col]][slope_df[[m_col]] == 1L]))

  invisible(NULL)
}
