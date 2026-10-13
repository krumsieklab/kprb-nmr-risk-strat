## ============================================================
## Longitudinal FIB-4 history from AST / ALT / platelet labs
## ------------------------------------------------------------
## FIB-4 = (age * AST) / (platelets * sqrt(ALT))
##
## The three component labs are drawn on different days and rarely
## coincide exactly, so a naive same-day join throws away most of the
## data. To capture as many timepoints as possible for trend
## visualization, this helper:
##
##   1. ANCHORS on AST measurements. AST is the binding constraint
##      (rarest of the three components), so anchoring on it maximizes
##      the number of computable FIB-4 points without inventing data.
##   2. For each AST draw, pulls the NEAREST ALT and NEAREST platelet
##      within +/- `window_days`. A FIB-4 point is only formed when
##      both partners exist inside the window (no imputation).
##   3. Uses the MEDIAN of the three contributing measurement dates as
##      the FIB-4 timepoint, and computes age at that date as
##      age_at_sample + (fib4_date - enrollment_date) / 365.25.
##
## This is intentionally different from the modeling FIB-4: it is a
## full longitudinal history (no pre-enrollment restriction, no
## lookback cap) meant for patient-centric trend line plots.
##
## Save as: helpers/compute_fib4_history.r
## Use via: source(here::here('helpers', 'compute_fib4_history.r'))
## ============================================================

suppressPackageStartupMessages({
  library(dplyr)
})

# ------------------------------------------------------------
# Standardize a component lab table to: <id>, comp_date, value.
# Coerces types, drops NA, and applies an optional plausibility filter.
# ------------------------------------------------------------
.fib4_prep_component <- function(df, id_col, value_col, date_col,
                                 value_range, collapse_same_day) {
  for (col in c(id_col, value_col, date_col)) {
    if (!col %in% names(df)) {
      stop(sprintf("component lab table is missing column '%s'.", col))
    }
  }

  out <- df %>%
    dplyr::transmute(
      !!id_col  := .data[[id_col]],
      comp_date = as.Date(.data[[date_col]]),
      value     = suppressWarnings(as.numeric(.data[[value_col]]))
    ) %>%
    dplyr::filter(!is.na(.data[[id_col]]), !is.na(comp_date), !is.na(value))

  if (!is.null(value_range)) {
    if (length(value_range) != 2 || any(is.na(value_range)) ||
        value_range[1] >= value_range[2]) {
      stop("value range must be c(min, max) with min < max, or NULL.")
    }
    out <- out %>%
      dplyr::filter(value >= value_range[1], value <= value_range[2])
  }

  if (isTRUE(collapse_same_day)) {
    out <- out %>%
      dplyr::group_by(.data[[id_col]], comp_date) %>%
      dplyr::summarise(value = stats::median(value, na.rm = TRUE), .groups = "drop")
  }

  out
}

# ------------------------------------------------------------
# Row-wise median of three Date vectors (returns Date).
# Median of 3 values = sum - min - max (the middle value).
# ------------------------------------------------------------
.fib4_median_date <- function(d1, d2, d3) {
  n1 <- as.numeric(d1); n2 <- as.numeric(d2); n3 <- as.numeric(d3)
  mid <- (n1 + n2 + n3) - pmin(n1, n2, n3) - pmax(n1, n2, n3)
  as.Date(mid, origin = "1970-01-01")
}

# ------------------------------------------------------------
# For each anchor (id, anchor_date), find the single nearest partner
# measurement within +/- window_days. Returns one row per anchor that
# has a qualifying partner: <id>, anchor_date, <out_value>, <out_date>.
# ------------------------------------------------------------
.fib4_nearest_within_window <- function(anchors, partner, id_col,
                                        window_days, out_value, out_date) {
  anchors %>%
    dplyr::inner_join(partner, by = id_col, relationship = "many-to-many") %>%
    dplyr::mutate(gap = abs(as.numeric(comp_date - anchor_date))) %>%
    dplyr::filter(gap <= window_days) %>%
    dplyr::group_by(.data[[id_col]], anchor_date) %>%
    dplyr::slice_min(gap, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::transmute(
      !!id_col   := .data[[id_col]],
      anchor_date,
      !!out_value := value,
      !!out_date  := comp_date
    )
}

# ------------------------------------------------------------
# Main: full longitudinal FIB-4 history.
#
# Args
#   ast_df, alt_df, platelet_df  long component lab tables
#   meta_df                      enrollment snapshot (age + index date)
#   id_col                       patient id column
#   value_col                    result column in each lab table
#   date_col                     measurement date column in each lab table
#   index_date_col               enrollment date in meta_df
#   age_col                      age at enrollment in meta_df
#   window_days                  +/- match window around each AST anchor
#                                  (default 90; widen to capture more points)
#   ast_range / alt_range /
#     platelet_range             plausibility filters; set NULL to disable
#   collapse_same_day            median-collapse same-day draws per component
#
# Returns one row per computable FIB-4 timepoint:
#   <id>, fib4_date, days_from_sample, age_at_result,
#   AST_new, ALT_new, PLATELETS_new,
#   ast_date, alt_date, platelet_date, max_gap_days, fib4
# ------------------------------------------------------------
compute_fib4_history <- function(ast_df,
                                 alt_df,
                                 platelet_df,
                                 meta_df,
                                 id_col            = "StudyID",
                                 value_col         = "Result_C",
                                 date_col          = "LAB_DATE",
                                 index_date_col    = "KPRB_SAMPLE_DATE",
                                 age_col           = "age_at_sample",
                                 window_days       = 90,
                                 ast_range         = c(1, 2000),
                                 alt_range         = c(1, 2000),
                                 platelet_range    = c(1, 2000),
                                 collapse_same_day = TRUE) {

  for (col in c(id_col, index_date_col, age_col)) {
    if (!col %in% names(meta_df)) stop(sprintf("meta_df is missing column '%s'.", col))
  }
  if (!is.numeric(window_days) || length(window_days) != 1 ||
      is.na(window_days) || window_days < 0) {
    stop("`window_days` must be a single non-negative number.")
  }

  # ---- standardize components ----
  ast <- .fib4_prep_component(ast_df, id_col, value_col, date_col,
                              ast_range, collapse_same_day)
  alt <- .fib4_prep_component(alt_df, id_col, value_col, date_col,
                              alt_range, collapse_same_day)
  plt <- .fib4_prep_component(platelet_df, id_col, value_col, date_col,
                              platelet_range, collapse_same_day)

  # ---- anchor on AST (binding constraint) ----
  anchors <- ast %>%
    dplyr::transmute(!!id_col := .data[[id_col]],
                     anchor_date = comp_date,
                     AST_new = value)

  # ---- nearest ALT and platelet within the window ----
  alt_match <- .fib4_nearest_within_window(
    anchors %>% dplyr::select(dplyr::all_of(id_col), anchor_date),
    alt, id_col, window_days, out_value = "ALT_new", out_date = "alt_date"
  )
  plt_match <- .fib4_nearest_within_window(
    anchors %>% dplyr::select(dplyr::all_of(id_col), anchor_date),
    plt, id_col, window_days, out_value = "PLATELETS_new", out_date = "platelet_date"
  )

  # ---- enrollment anchors (age + index date) ----
  meta_idx <- meta_df %>%
    dplyr::transmute(
      !!id_col      := .data[[id_col]],
      index_date    = as.Date(.data[[index_date_col]]),
      age_at_sample = as.numeric(.data[[age_col]])
    ) %>%
    dplyr::filter(!is.na(.data[[id_col]]), !is.na(index_date), !is.na(age_at_sample)) %>%
    dplyr::distinct(.data[[id_col]], .keep_all = TRUE)

  # ---- assemble complete AST/ALT/platelet triples + FIB-4 ----
  anchors %>%
    dplyr::rename(ast_date = anchor_date) %>%
    dplyr::inner_join(alt_match %>% dplyr::rename(ast_date = anchor_date),
                      by = c(id_col, "ast_date")) %>%
    dplyr::inner_join(plt_match %>% dplyr::rename(ast_date = anchor_date),
                      by = c(id_col, "ast_date")) %>%
    dplyr::inner_join(meta_idx, by = id_col) %>%
    dplyr::mutate(
      # median of the three contributing dates = the FIB-4 timepoint
      fib4_date = .fib4_median_date(ast_date, alt_date, platelet_date),
      max_gap_days = pmax(
        abs(as.numeric(alt_date - ast_date)),
        abs(as.numeric(platelet_date - ast_date))
      ),
      days_from_sample = as.numeric(fib4_date - index_date),
      age_at_result    = age_at_sample + days_from_sample / 365.25,
      fib4 = (age_at_result * AST_new) / (PLATELETS_new * sqrt(ALT_new))
    ) %>%
    dplyr::transmute(
      !!id_col := .data[[id_col]],
      fib4_date, days_from_sample, age_at_result,
      AST_new, ALT_new, PLATELETS_new,
      ast_date, alt_date, platelet_date, max_gap_days, fib4
    ) %>%
    dplyr::arrange(.data[[id_col]], fib4_date)
}
