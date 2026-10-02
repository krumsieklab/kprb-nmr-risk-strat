## ============================================================
## Longitudinal eGFR history from creatinine labs
## ------------------------------------------------------------
## Converts every creatinine measurement into an eGFR value using
## the 2021 CKD-EPI creatinine equation (race-free).
##
## This is intentionally different from compute_egfr_slope.r:
##   - no pre-enrollment restriction
##   - no lookback window
##   - no slope/modeling gate
##
## It returns the full longitudinal lab history with eGFR attached,
## useful for patient-centric plots and downstream longitudinal scores.
##
## Save as: helpers/compute_egfr_history.r
## Use via: source(here('helpers', 'compute_egfr_history.r'))
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
})

# ------------------------------------------------------------
# Resolve sex/gender to CKD-EPI coding: 1 = male, 0 = female.
# Accepts factor/character/numeric codings ('M'/'F', 1/0, etc.).
# ------------------------------------------------------------
.egfr_history_resolve_sex_num <- function(x) {
  xc <- toupper(trimws(as.character(x)))
  ifelse(
    xc %in% c("M", "1", "MALE", "TRUE"), 1L,
    ifelse(xc %in% c("F", "0", "FEMALE", "FALSE"), 0L, NA_integer_)
  )
}

# ------------------------------------------------------------
# 2021 CKD-EPI creatinine eGFR (race-free), vectorized.
#   scr      serum creatinine (mg/dL)
#   age      age in years at the measurement date
#   sex_num  1 = male, 0 = female
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
# Main: full longitudinal eGFR history.
#
# Args
#   labs_df              long creatinine lab table
#   meta_df              enrollment snapshot table
#   id_col               patient id column
#   value_col            creatinine result column (mg/dL)
#   date_col             creatinine result date column
#   index_date_col       sample/enrollment date in meta_df
#   age_col              age at sample/enrollment in meta_df
#   sex_col              sex/gender column in meta_df
#   test_type_col        optional lab type column to keep if present
#   creat_range          plausibility filter for creatinine; set NULL to disable
#   collapse_same_day    if TRUE, collapse same-day creatinine values to median
#
# Returns
#   StudyID, Test_Type (if present), RESULT_DATE, KPRB_SAMPLE_DATE,
#   creatinine, age_at_sample, sex, sex_num, year_diff_result,
#   age_at_result, days_from_sample, days_before_sample, egfr
# ------------------------------------------------------------
compute_egfr_history <- function(labs_df,
                                 meta_df,
                                 id_col            = "StudyID",
                                 value_col         = "Result_C",
                                 date_col          = "RESULT_DATE",
                                 index_date_col    = "KPRB_SAMPLE_DATE",
                                 age_col           = "age_at_sample",
                                 sex_col           = "sex",
                                 test_type_col     = "Test_Type",
                                 creat_range       = c(0.1, 20),
                                 collapse_same_day = FALSE) {

  for (col in c(id_col, value_col, date_col)) {
    if (!col %in% names(labs_df)) {
      stop(sprintf("labs_df is missing column '%s'.", col))
    }
  }
  for (col in c(id_col, index_date_col, age_col, sex_col)) {
    if (!col %in% names(meta_df)) {
      stop(sprintf("meta_df is missing column '%s'.", col))
    }
  }

  keep_test_type <- !is.null(test_type_col) && test_type_col %in% names(labs_df)

  labs_clean <- labs_df %>%
    dplyr::transmute(
      !!id_col := .data[[id_col]],
      Test_Type = if (keep_test_type) as.character(.data[[test_type_col]]) else NA_character_,
      RESULT_DATE = as.Date(.data[[date_col]]),
      creatinine = suppressWarnings(as.numeric(.data[[value_col]]))
    ) %>%
    dplyr::filter(!is.na(.data[[id_col]]), !is.na(RESULT_DATE), !is.na(creatinine))

  if (!is.null(creat_range)) {
    if (length(creat_range) != 2 || any(is.na(creat_range)) || creat_range[1] >= creat_range[2]) {
      stop("`creat_range` must be c(min, max) with min < max, or NULL.")
    }
    labs_clean <- labs_clean %>%
      dplyr::filter(creatinine >= creat_range[1], creatinine <= creat_range[2])
  }

  if (isTRUE(collapse_same_day)) {
    labs_clean <- labs_clean %>%
      dplyr::group_by(.data[[id_col]], RESULT_DATE) %>%
      dplyr::summarise(
        Test_Type = dplyr::first(Test_Type),
        creatinine = stats::median(creatinine, na.rm = TRUE),
        .groups = "drop"
      )
  }

  meta_idx <- meta_df %>%
    dplyr::transmute(
      !!id_col := .data[[id_col]],
      KPRB_SAMPLE_DATE = as.Date(.data[[index_date_col]]),
      age_at_sample = as.numeric(.data[[age_col]]),
      sex = .data[[sex_col]],
      sex_num = .egfr_history_resolve_sex_num(.data[[sex_col]])
    ) %>%
    dplyr::filter(
      !is.na(.data[[id_col]]),
      !is.na(KPRB_SAMPLE_DATE),
      !is.na(age_at_sample),
      !is.na(sex_num)
    ) %>%
    dplyr::distinct(.data[[id_col]], .keep_all = TRUE)

  labs_clean %>%
    dplyr::inner_join(meta_idx, by = id_col) %>%
    dplyr::mutate(
      days_from_sample = as.numeric(RESULT_DATE - KPRB_SAMPLE_DATE),
      days_before_sample = as.numeric(KPRB_SAMPLE_DATE - RESULT_DATE),
      year_diff_result = days_from_sample / 365.25,
      age_at_result = age_at_sample + year_diff_result,
      egfr = ckd_epi_2021_egfr(creatinine, age_at_result, sex_num)
    ) %>%
    dplyr::select(
      dplyr::all_of(id_col),
      Test_Type,
      RESULT_DATE,
      KPRB_SAMPLE_DATE,
      creatinine,
      age_at_sample,
      age_at_result,
      sex,
      sex_num,
      year_diff_result,
      days_from_sample,
      days_before_sample,
      egfr
    ) %>%
    dplyr::arrange(.data[[id_col]], RESULT_DATE)
}
