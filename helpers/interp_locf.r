
# Estimate BMI / SBP at the KPRB serum-collection index date
library(dplyr)

estimate_index_value <- function(
    df,
    id_col                  = "StudyID",
    value_col               = "Result_C",
    index_date_col          = "KPRB_SAMPLE_DATE",
    measurement_date_col    = "MEASUREMENT_DATE",
    proximal_days           = 90,
    interp_before_days      = 365,
    interp_after_days       = 365,
    locf_days               = 730,
    output_col              = NULL
) {

  required_cols <- c(
    id_col, value_col, index_date_col, measurement_date_col
  )

  missing_cols <- setdiff(required_cols, names(df))

  if (length(missing_cols) > 0) {
    stop(
      "Missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  }

  check_nonnegative_scalar <- function(x, name) {
    if (
      length(x) != 1L ||
      !is.numeric(x) ||
      is.na(x) ||
      !is.finite(x) ||
      x < 0
    ) {
      stop("`", name, "` must be one finite non-negative number.")
    }
  }

  check_nonnegative_scalar(proximal_days,      "proximal_days")
  check_nonnegative_scalar(interp_before_days, "interp_before_days")
  check_nonnegative_scalar(interp_after_days,  "interp_after_days")
  check_nonnegative_scalar(locf_days,          "locf_days")

  if (!is.numeric(df[[value_col]])) {
    stop("`", value_col, "` must be numeric.")
  }

  if (is.null(output_col)) {
    output_col <- paste0("imputed_", value_col)
  }

  as_date_safe <- function(x, name) {
    if (inherits(x, "Date")) {
      return(x)
    }

    if (inherits(x, c("POSIXct", "POSIXlt"))) {
      return(as.Date(x))
    }

    if (is.numeric(x)) {
      stop(
        "`", name, "` is numeric. Convert it explicitly to Date before ",
        "calling estimate_index_value()."
      )
    }

    out <- suppressWarnings(as.Date(as.character(x)))
    bad <- !is.na(x) & is.na(out)

    if (any(bad)) {
      stop(
        "Some non-missing values in `", name,
        "` could not be converted to Date."
      )
    }

    out
  }

  dat <- df %>%
    transmute(
      .id = .data[[id_col]],
      .value = .data[[value_col]],
      .index = as_date_safe(.data[[index_date_col]], index_date_col),
      .date = as_date_safe(
        .data[[measurement_date_col]],
        measurement_date_col
      )
    )

  if (any(is.na(dat$.id))) {
    stop("Missing participant IDs detected in `", id_col, "`.")
  }

  index_dates <- dat %>%
    filter(!is.na(.index)) %>%
    distinct(.id, .index)

  multiple_index_dates <- index_dates %>%
    count(.id, name = "n_index_dates") %>%
    filter(n_index_dates > 1L)

  if (nrow(multiple_index_dates) > 0L) {
    stop(
      nrow(multiple_index_dates),
      " participant(s) have more than one distinct index date."
    )
  }

  ids_in_source <- dat %>%
    distinct(.id)

  index_tbl <- ids_in_source %>%
    left_join(index_dates, by = ".id")

  n_missing_index <- sum(is.na(index_tbl$.index))

  if (n_missing_index > 0L) {
    stop(
      n_missing_index,
      " participant(s) do not have a non-missing index date."
    )
  }

  # Collapse multiple measurements on the same date.
  daily <- dat %>%
    filter(
      !is.na(.date),
      !is.na(.value)
    ) %>%
    group_by(.id, .date) %>%
    summarise(
      .value = median(.value, na.rm = TRUE),
      .groups = "drop"
    )

  long <- index_tbl %>%
    left_join(daily, by = ".id")

  # Find exact, nearest-before, and nearest-after measurements.
  anchors <- long %>%
    group_by(.id, .index) %>%
    summarise(

      exact_value = {
        ok <- !is.na(.date) & .date == .index
        if (any(ok)) .value[ok][1L] else NA_real_
      },

      before_date = {
        ok <- !is.na(.date) & .date < .index
        if (any(ok)) max(.date[ok]) else as.Date(NA)
      },

      before_value = {
        ok <- !is.na(.date) & .date < .index

        if (any(ok)) {
          d <- .date[ok]
          v <- .value[ok]
          nearest_date <- max(d)
          v[d == nearest_date][1L]
        } else {
          NA_real_
        }
      },

      after_date = {
        ok <- !is.na(.date) & .date > .index
        if (any(ok)) min(.date[ok]) else as.Date(NA)
      },

      after_value = {
        ok <- !is.na(.date) & .date > .index

        if (any(ok)) {
          d <- .date[ok]
          v <- .value[ok]
          nearest_date <- min(d)
          v[d == nearest_date][1L]
        } else {
          NA_real_
        }
      },

      .groups = "drop"
    )

  # Apply hierarachy: 'Exact', followed by 'Linear_interp', then 'Index_proximal_pre', 'Index_proximal_post', 'LOCF', and 'None'
  result <- anchors %>%
    mutate(
      days_before = as.numeric(.index - before_date),
      days_after  = as.numeric(after_date - .index),

      can_interpolate =
        !is.na(before_value) &
        !is.na(after_value) &
        !is.na(days_before) &
        !is.na(days_after) &
        days_before <= interp_before_days &
        days_after  <= interp_after_days,

      interpolation_weight_after = if_else(
        can_interpolate,
        days_before / (days_before + days_after),
        NA_real_
      ),

      interpolated_value = if_else(
        can_interpolate,
        before_value +
          interpolation_weight_after * (after_value - before_value),
        NA_real_
      ),

      proximal_before =
        !is.na(before_value) &
        !is.na(days_before) &
        days_before <= proximal_days,

      proximal_after =
        !is.na(after_value) &
        !is.na(days_after) &
        days_after <= proximal_days,

      can_locf =
        !is.na(before_value) &
        !is.na(days_before) &
        days_before <= locf_days,

      estimated_value = case_when(
        # 1. Exact
        !is.na(exact_value) ~ exact_value,

        # 2. Linear interpolation
        can_interpolate ~ interpolated_value,

        # 3. Proximal PRE
        proximal_before ~ before_value,

        # 4. Proximal POST
        proximal_after ~ after_value,

        # 5. LOCF
        can_locf ~ before_value,

        # 6. Missing
        TRUE ~ NA_real_
      ),

      method = case_when(
        !is.na(exact_value) ~ "Exact",
        can_interpolate ~ "Linear_interp",
        proximal_before ~ "Index_proximal_pre",
        proximal_after ~ "Index_proximal_post",
        can_locf ~ "LOCF",
        TRUE ~ "None"
      ),

      # Record only measurements actually used.
      before_date_used = case_when(
        method == "Exact"              ~ .index,
        method == "Linear_interp"      ~ before_date,
        method == "Index_proximal_pre" ~ before_date,
        method == "LOCF"               ~ before_date,
        TRUE                           ~ as.Date(NA)
      ),

      before_value_used = case_when(
        method == "Exact"              ~ exact_value,
        method == "Linear_interp"      ~ before_value,
        method == "Index_proximal_pre" ~ before_value,
        method == "LOCF"               ~ before_value,
        TRUE                           ~ NA_real_
      ),

      after_date_used = case_when(
        method == "Linear_interp"       ~ after_date,
        method == "Index_proximal_post" ~ after_date,
        TRUE                            ~ as.Date(NA)
      ),

      after_value_used = case_when(
        method == "Linear_interp"       ~ after_value,
        method == "Index_proximal_post" ~ after_value,
        TRUE                            ~ NA_real_
      ),

      days_before_used = case_when(
        method == "Exact"              ~ 0,
        method == "Linear_interp"      ~ days_before,
        method == "Index_proximal_pre" ~ days_before,
        method == "LOCF"               ~ days_before,
        TRUE                           ~ NA_real_
      ),

      days_after_used = case_when(
        method == "Exact"               ~ 0,
        method == "Linear_interp"       ~ days_after,
        method == "Index_proximal_post" ~ days_after,
        TRUE                            ~ NA_real_
      )
    )

  result %>%
    transmute(
      !!id_col := .id,
      !!index_date_col := .index,
      !!output_col := estimated_value,
      method,
      days_before_used,
      days_after_used,
      before_date_used,
      after_date_used,
      before_value_used,
      after_value_used,
      interpolation_weight_after = if_else(
        method == "Linear_interp",
        interpolation_weight_after,
        NA_real_
      )
    )
}


# =============================================================================
# QC helpers
# =============================================================================

summarize_index_value_methods <- function(x, method_col = "method") {
  x %>%
    count(.data[[method_col]], name = "n") %>%
    mutate(
      pct = 100 * n / sum(n)
    ) %>%
    rename(
      method = all_of(method_col)
    ) %>%
    mutate(
      method = factor(
        method,
        levels = c(
          "Exact",
          "Linear_interp",
          "Index_proximal_pre",
          "Index_proximal_post",
          "LOCF",
          "None"
        )
      )
    ) %>%
    arrange(method) %>%
    mutate(method = as.character(method))
}


summarize_index_value_timing <- function(x) {
  safe_median <- function(z) {
    if (all(is.na(z))) NA_real_ else median(z, na.rm = TRUE)
  }

  safe_quantile <- function(z, p) {
    if (all(is.na(z))) {
      NA_real_
    } else {
      unname(quantile(z, p, na.rm = TRUE))
    }
  }

  x %>%
    group_by(method) %>%
    summarise(
      n = n(),
      pct = 100 * n / nrow(x),
      median_days_before = safe_median(days_before_used),
      q1_days_before = safe_quantile(days_before_used, 0.25),
      q3_days_before = safe_quantile(days_before_used, 0.75),
      median_days_after = safe_median(days_after_used),
      q1_days_after = safe_quantile(days_after_used, 0.25),
      q3_days_after = safe_quantile(days_after_used, 0.75),
      .groups = "drop"
    )
}

