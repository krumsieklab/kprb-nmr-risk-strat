# Baseline statin and blood-pressure treatment indicators.
#
# Standalone copy of the 120-day rules in R/derive_baseline_medication_status.R.
# Pass meta and rx in. Sourcing this file does not search for those objects.
#
# meta: StudyID, KPRB_SAMPLE_DATE
# rx:   StudyID, GENERIC_NM, RX_DATE
#
# statin_treatment and bp_treatment are 1 if a qualifying dispensing falls on
# the index date or within 120 days before it, 0 if not, and NA if the index
# date is missing. People with an index date and no dispensing are 0.

statin_names <- c(
  "ATORVASTATIN",
  "FLUVASTATIN",
  "LOVASTATIN",
  "PITAVASTATIN",
  "PRAVASTATIN",
  "ROSUVASTATIN",
  "SIMVASTATIN",
  "SIMVASTIN"
)

bp_med_names <- c(
  # ACE inhibitors
  "BENAZEPRIL",
  "CAPTOPRIL",
  "ENALAPRIL",
  "LISINOPRIL",
  "RAMIPRIL",
  # ARBs
  "CANDESARTAN",
  "IRBESARTAN",
  "LOSARTAN",
  "OLMESARTAN",
  "TELMISARTAN",
  "VALSARTAN",
  # Calcium-channel blockers
  "AMLODIPINE",
  "DILTIAZEM",
  "FELODIPINE",
  "ISRADIPINE",
  "NIFEDIPINE",
  "VERAPAMIL",
  # Thiazide and thiazide-like diuretics
  "CHLORTHALIDONE",
  "HYDROCHLOROTH",
  "HYDROCHLOROTHIAZIDE",
  "INDAPAMIDE",
  "METOLAZONE",
  # Loop diuretics
  "BUMETANIDE",
  "ETHACRYNIC",
  "FUROSEMIDE",
  "TORSEMIDE",
  # Potassium-sparing diuretics and mineralocorticoid receptor antagonists
  "AMILORIDE",
  "SPIRONOLACT",
  "SPIRONOLACTONE",
  "TRIAMTERENE",
  # Central antihypertensives
  "CLONIDINE",
  "GUANFACINE",
  "METHYLDOPA",
  # Direct vasodilators
  "HYDRALAZINE",
  "MINOXIDIL",
  # Alpha blockers and other adrenergic blockers
  "DOXAZOSIN",
  "PHENOXYBENZAMINE",
  "PHENTOLAMINE",
  "PRAZOSIN",
  "TERAZOSIN",
  # Beta blockers. BETAXOLOL and TIMOLOL are omitted because the observed
  # records may represent ophthalmic formulations rather than systemic BP treatment.
  "ACEBUTOLOL",
  "ATENOLOL",
  "BISOPROLOL",
  "CARVEDILOL",
  "ESMOLOL",
  "LABETALOL",
  "METOPROLOL",
  "NADOLOL",
  "NEBIVOLOL",
  "PINDOLOL",
  "PROPRANOLOL"
)

normalize_generic_name <- function(x) {
  x <- as.character(x)
  x <- toupper(x)
  x <- trimws(x)
  gsub("\\s+", " ", x)
}

parse_date_for_rx_status <- function(x) {
  if (inherits(x, "Date")) {
    return(x)
  }

  if (inherits(x, c("POSIXct", "POSIXlt"))) {
    return(as.Date(x))
  }

  if (is.numeric(x)) {
    return(as.Date(x, origin = "1970-01-01"))
  }

  x_chr <- trimws(as.character(x))
  x_chr[x_chr == ""] <- NA_character_
  parsed <- suppressWarnings(as.Date(x_chr))
  date_formats <- c("%m/%d/%Y", "%m/%d/%y", "%Y/%m/%d", "%d-%b-%Y", "%d%b%Y")

  for (fmt in date_formats) {
    needs_parse <- !is.na(x_chr) & is.na(parsed)
    if (!any(needs_parse)) {
      break
    }
    parsed[needs_parse] <- suppressWarnings(as.Date(x_chr[needs_parse], format = fmt))
  }

  parsed
}

integer_status <- function(x, index_missing) {
  out <- as.integer(dplyr::coalesce(x, FALSE))
  out[index_missing] <- NA_integer_
  out
}

derive_baseline_bp_statin_status <- function(meta, rx) {
  if (missing(meta) || missing(rx)) {
    stop(
      "`derive_baseline_bp_statin_status()` requires both `meta` and `rx`.",
      call. = FALSE
    )
  }

  message("Medication mapping review:")
  print(list(
    statin_names = statin_names,
    bp_med_names = bp_med_names
  ))

  missing_meta <- setdiff(c("StudyID", "KPRB_SAMPLE_DATE"), names(meta))
  missing_rx <- setdiff(c("StudyID", "GENERIC_NM", "RX_DATE"), names(rx))
  if (length(missing_meta) > 0) {
    stop("meta is missing required column(s): ", paste(missing_meta, collapse = ", "), call. = FALSE)
  }
  if (length(missing_rx) > 0) {
    stop("rx is missing required column(s): ", paste(missing_rx, collapse = ", "), call. = FALSE)
  }
  if (!is.atomic(meta$StudyID) || !is.atomic(rx$StudyID) || is.list(meta$StudyID) || is.list(rx$StudyID)) {
    stop("StudyID columns are not compatible atomic vectors in meta and rx.", call. = FALSE)
  }

  meta_index <- tibble::as_tibble(meta) %>%
    dplyr::mutate(
      StudyID = as.character(.data$StudyID),
      KPRB_SAMPLE_DATE = parse_date_for_rx_status(.data$KPRB_SAMPLE_DATE)
    ) %>%
    dplyr::group_by(.data$StudyID) %>%
    dplyr::summarise(
      KPRB_SAMPLE_DATE = dplyr::first(.data$KPRB_SAMPLE_DATE),
      n_distinct_sample_dates = dplyr::n_distinct(.data$KPRB_SAMPLE_DATE, na.rm = TRUE),
      .groups = "drop"
    )

  if (any(meta_index$n_distinct_sample_dates > 1L)) {
    warning(
      "Some duplicated StudyID values in meta have different KPRB_SAMPLE_DATE values; using the first date per StudyID.",
      call. = FALSE
    )
  }

  rx_aug <- tibble::as_tibble(rx) %>%
    dplyr::mutate(
      StudyID = as.character(.data$StudyID),
      RX_DATE = parse_date_for_rx_status(.data$RX_DATE),
      generic_nm_std = normalize_generic_name(.data$GENERIC_NM),
      is_statin = .data$generic_nm_std %in% statin_names,
      is_bp_med = .data$generic_nm_std %in% bp_med_names
    ) %>%
    dplyr::left_join(
      meta_index %>% dplyr::select(StudyID, KPRB_SAMPLE_DATE),
      by = "StudyID"
    ) %>%
    dplyr::mutate(
      days_before_index = as.integer(.data$KPRB_SAMPLE_DATE - .data$RX_DATE)
    )

  by_participant <- rx_aug %>%
    dplyr::group_by(.data$StudyID) %>%
    dplyr::summarise(
      statin_treatment = any(
        .data$is_statin & .data$days_before_index >= 0L & .data$days_before_index <= 120L,
        na.rm = TRUE
      ),
      bp_treatment = any(
        .data$is_bp_med & .data$days_before_index >= 0L & .data$days_before_index <= 120L,
        na.rm = TRUE
      ),
      .groups = "drop"
    )

  # Participants without rx rows are coded 0 when the index date is observed.
  meta_index %>%
    dplyr::select(StudyID, KPRB_SAMPLE_DATE) %>%
    dplyr::left_join(by_participant, by = "StudyID") %>%
    dplyr::mutate(
      statin_treatment = integer_status(.data$statin_treatment, is.na(.data$KPRB_SAMPLE_DATE)),
      bp_treatment = integer_status(.data$bp_treatment, is.na(.data$KPRB_SAMPLE_DATE))
    ) %>%
    dplyr::select(StudyID, statin_treatment, bp_treatment)
}
