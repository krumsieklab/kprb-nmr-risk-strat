locf_exact_impute <- function(df,
                          value_col      = "Result_C",
                          index_date_col = "KPRB_SAMPLE_DATE",
                          date_col       = "MEASUREMENT_DATE",
                          exact_window   = 0,        # Days for "Exactish" window (+/- days)
                          locf_window    = 730) {    # Days for LOCF fallback (default 2 years)

  value_sym <- rlang::sym(value_col)
  idx_sym   <- rlang::sym(index_date_col)
  date_sym  <- rlang::sym(date_col)

  imputed_col <- paste0("imputed_", value_col)

  # -------------------------------------------------------------------------
  # Validate parameters
  # -------------------------------------------------------------------------
  if (!is.numeric(exact_window) || length(exact_window) != 1 || is.na(exact_window) || exact_window < 0) {
    stop("`exact_window` must be a single non-negative number.")
  }
  if (!is.numeric(locf_window) || length(locf_window) != 1 || is.na(locf_window) || locf_window < 0) {
    stop("`locf_window` must be a single non-negative number.")
  }

  # -------------------------------------------------------------------------
  # 0. Make sure both date columns are real Date objects
  # -------------------------------------------------------------------------
  df <- df %>%
    dplyr::mutate(
      dplyr::across(
        .cols = dplyr::all_of(c(index_date_col, date_col)),
        .fns  = ~ as.Date(.x)
      )
    )

  # -------------------------------------------------------------------------
  # 1. Keep every StudyID that arrived with an index date.
  #    People with no usable measurement are joined back after imputation
  #    so a later missingness block can zero-gate them instead of dropping them.
  # -------------------------------------------------------------------------
  cohort <- df %>%
    dplyr::filter(!is.na(!!idx_sym)) %>%
    dplyr::group_by(StudyID) %>%
    dplyr::summarise(
      !!index_date_col := dplyr::first(!!idx_sym),
      .groups = "drop"
    )

  # -------------------------------------------------------------------------
  # 2. Collapse multiple measurements on the same day (median)
  # -------------------------------------------------------------------------
  df <- df %>%
    dplyr::filter(!is.na(!!value_sym),
                  !is.na(!!idx_sym),
                  !is.na(!!date_sym)) %>%
    dplyr::group_by(StudyID, !!date_sym) %>%
    dplyr::summarise(
      !!value_col      := median(!!value_sym, na.rm = TRUE),
      !!index_date_col := dplyr::first(!!idx_sym),
      .groups = "drop"
    )

  # -------------------------------------------------------------------------
  # 3. Imputation logic with NO future influence
  # Priority: 1) Exact/Exactish  2) LOCF  3) NA
  # -------------------------------------------------------------------------
  imputed <- df %>%
    dplyr::group_by(StudyID) %>%
    dplyr::summarise(
      INDEX_DATE = dplyr::first(!!idx_sym),

      # Only look backward or same-day
      before_date = {
        idx <- dplyr::first(!!idx_sym)
        bd  <- (!!date_sym)[!!date_sym <= idx]
        if (length(bd)) max(bd) else as.Date(NA)
      },

      before_val = if (!is.na(before_date)) {
        median((!!value_sym)[!!date_sym == before_date], na.rm = TRUE)
      } else {
        NA_real_
      },

      .groups = "drop"
    ) %>%
    dplyr::mutate(
      before_dist = abs(as.numeric(INDEX_DATE - before_date)),

      !!imputed_col := dplyr::case_when(
        # 1. EXACT: same-day measurement
        !is.na(before_date) & before_date == INDEX_DATE ~ before_val,

        # 2. EXACTISH: prior measurement within exact_window
        !is.na(before_date) & before_dist <= exact_window ~ before_val,

        # 3. LOCF: last observation within locf_window
        !is.na(before_val) & before_dist <= locf_window ~ before_val,

        # 4. NA: no valid backward imputation
        TRUE ~ NA_real_
      ),

      method = dplyr::case_when(
        !is.na(before_date) & before_date == INDEX_DATE ~ "Exact",
        !is.na(before_date) & before_dist <= exact_window ~ "Exactish",
        !is.na(before_val)  & before_dist <= locf_window ~ "LOCF",
        TRUE ~ "None"
      ),

      days_before = dplyr::case_when(
        method %in% c("Exact", "Exactish", "LOCF") & !is.na(before_date) ~ before_dist,
        TRUE ~ NA_real_
      ),

      days_after = NA_real_,

      before_date_used = dplyr::case_when(
        method %in% c("Exact", "Exactish", "LOCF") ~ before_date,
        TRUE ~ as.Date(NA)
      ),

      after_date_used = as.Date(NA),

      before_val_used = dplyr::case_when(
        method %in% c("Exact", "Exactish", "LOCF") ~ before_val,
        TRUE ~ NA_real_
      ),

      after_val_used = NA_real_
    ) %>%
    dplyr::select(-before_dist, -before_date, -before_val) %>%
    dplyr::rename(!!index_date_col := INDEX_DATE) %>%
    dplyr::select(
      StudyID,
      !!idx_sym,
      !!imputed_col,
      method,
      days_before,
      days_after,
      before_date_used,
      after_date_used,
      before_val_used,
      after_val_used
    )

  cohort %>%
    dplyr::left_join(
      imputed %>% dplyr::select(-dplyr::all_of(index_date_col)),
      by = "StudyID"
    ) %>%
    dplyr::mutate(
      method = dplyr::if_else(is.na(method), "None", method)
    ) %>%
    dplyr::select(
      StudyID,
      dplyr::all_of(index_date_col),
      dplyr::all_of(imputed_col),
      method,
      days_before,
      days_after,
      before_date_used,
      after_date_used,
      before_val_used,
      after_val_used
    )
}