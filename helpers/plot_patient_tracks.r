plot_patient_tracks <- function(
  risk_name,
  patient_id,
  meta,
  bmi = NULL, bp = NULL, a1c = NULL, ldl = NULL, hdl = NULL,
  gluran = NULL, glu_f = NULL, creatinine = NULL, egfr = NULL,
  alt = NULL, fib4 = NULL, rx = NULL,
  df_soc = NULL,
  df_nmr = NULL,
  df_base = NULL,  # legacy alias for df_soc
  # df_nm REMOVED: legacy alias for df_nmr
  risk_horizons = c(1L, 3L, 5L),
  date_start = NULL, date_end = NULL,
  use_thresholds = TRUE,
  drop_empty_tracks = TRUE,
  clinical_specs = NULL,
  show_placeholder_risk_panels = TRUE
) {
  if (!requireNamespace("dplyr", quietly = TRUE)) stop("Please install.packages('dplyr')")
  if (!requireNamespace("ggplot2", quietly = TRUE)) stop("Please install.packages('ggplot2')")
  if (!requireNamespace("patchwork", quietly = TRUE)) stop("Please install.packages('patchwork')")

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  EGFR_SPEC <- "eGFR (mL/min/1.73m2)"
  EGFR_YLAB <- "eGFR\n(mL/min/1.73m2)"

  # MARKER-only = SOC; MARKER + NMR = NMR.
  # Legacy aliases:
  #   df_base -> df_soc
  #   (df_nm   -> df_nmr) [REMOVED]
  df_soc <- df_soc %||% df_base

  # Remove all references to df_nm. No longer used.

  hz <- suppressWarnings(as.integer(unique(risk_horizons[!is.na(risk_horizons) & risk_horizons > 0])))
  hz <- hz[hz > 0]
  if (!length(hz)) stop("`risk_horizons` must contain at least one positive horizon (years).")

  # ---------- helpers for optional tracks ----------
  empty_track <- function() {
    data.frame(
      StudyID = character(),
      CLIN_DATE = as.Date(character()),
      CLIN_VAL = numeric(),
      stringsAsFactors = FALSE
    )
  }

  filter_track <- function(df) {
    if (is.null(df)) return(empty_track())
    df %>%
      dplyr::filter(.data$StudyID == patient_id) %>%
      dplyr::mutate(CLIN_DATE = as.Date(.data$CLIN_DATE))
  }

  empty_rx <- function() {
    data.frame(
      StudyID = character(),
      CLIN_DATE = as.Date(character()),
      GENERIC_NM = character(),
      stringsAsFactors = FALSE
    )
  }

  # ---------- default clinical specs ----------
  if (is.null(clinical_specs)) {
    clinical_specs <- list(
      "HbA1c (%)" = list(
        y_limits = c(5.0, 12.0),
        bands = data.frame(
          ymin  = c(-Inf, 5.7, 6.5),
          ymax  = c(5.7, 6.5, Inf),
          label = c("Normal", "Prediabetes", "Diabetes"),
          fill  = c("#e8f5e9", "#fff8e1", "#ffebee"),
          stringsAsFactors = FALSE
        ),
        lines = c(5.7, 6.5),
        label_bands = TRUE,
        ylabel_pos = c(5.1, 6.25, 6.9)
      ),
      "LDL-C (mg/dL)" = list(
        y_limits = c(50, 220),
        bands = data.frame(
          ymin  = c(-Inf, 100, 130, 160, 190),
          ymax  = c(100, 130, 160, 190, Inf),
          label = c("Optimal", "Near/Above Opt", "Borderline High", "High", "Very High"),
          fill  = c("#e8f5e9", "#e3f2fd", "#fff8e1", "#ffe0b2", "#ffebee"),
          stringsAsFactors = FALSE
        ),
        lines = c(100, 130, 160, 190),
        label_bands = TRUE
      ),
      "Creatinine (mg/dL)" = list(
        y_limits = c(0.4, 3.5),
        bands = NULL,
        lines = c(1.3),
        label_bands = FALSE
      ),
      "eGFR (mL/min/1.73m2)" = list(
        y_limits = c(0, 120),
        bands = data.frame(
          ymin  = c(-Inf, 15, 30, 45, 60, 90),
          ymax  = c(15, 30, 45, 60, 90, Inf),
          label = c("G5 (failure)", "G4 (severe)", "G3b", "G3a", "G2 (mild)", "G1 (normal)"),
          fill  = c("#b71c1c", "#ef9a9a", "#ffe0b2", "#fffde7", "#e8f5e9", "#e3f2fd"),
          stringsAsFactors = FALSE
        ),
        lines = c(15, 30, 45, 60, 90),
        label_bands = TRUE,
        ylabel_pos = c(7, 22, 37, 52, 75, 105)
      ),
      "FIB-4" = list(
        y_limits = c(0, 5),
        bands = data.frame(
          ymin  = c(-Inf, 1.3, 2.67),
          ymax  = c(1.3, 2.67, Inf),
          label = c("Low risk", "Indeterminate", "High risk"),
          fill  = c("#e8f5e9", "#fffde7", "#ef9a9a"),
          stringsAsFactors = FALSE
        ),
        lines = c(1.3, 2.67),
        label_bands = TRUE,
        ylabel_pos = c(0.65, 2.0, 3.8)
      ),
      "Fasting Glucose (mg/dL)" = list(
        y_limits = c(50, 200),
        bands = data.frame(
          ymin  = c(-Inf, 70, 100, 126),
          ymax  = c(70, 100, 126, Inf),
          label = c("Low", "Normal", "Prediabetes", "Diabetes"),
          fill  = c("#bbdefb", "#e8f5e9", "#fffde7", "#ef9a9a"),
          stringsAsFactors = FALSE
        ),
        lines = c(70, 100, 126),
        label_bands = TRUE,
        ylabel_pos = c(60, 90, 112, 160)
      )
    )
  }

  # ---------- patient meta ----------
  patient_meta <- meta %>% dplyr::filter(.data$StudyID == patient_id)
  if (nrow(patient_meta) == 0) stop("No meta row for patient_id: ", patient_id)

  # Display / event-line name is MASLD_DATE; older metadata still uses NAFLD_DATE.
  if ("NAFLD_DATE" %in% names(patient_meta) && !"MASLD_DATE" %in% names(patient_meta)) {
    names(patient_meta)[names(patient_meta) == "NAFLD_DATE"] <- "MASLD_DATE"
  }

  date_cols   <- c("KPRB_SAMPLE_DATE", "HYPERTENS_DATE", "DEATH_DATE", "DIAB_DATE", "CKD_DATE", "MASLD_DATE")
  date_labels <- c("Enrollment", "Hypertension", "Death", "Diabetes", "CKD", "MASLD")
  date_colors <- c("red", "orange", "black", "purple", "brown", "green")

  vertical_lines <- data.frame()
  enrollment_date <- NULL

  for (i in seq_along(date_cols)) {
    if (!date_cols[i] %in% names(patient_meta)) next

    val <- patient_meta[[date_cols[i]]][1]
    if (is.null(val) || is.na(val) || identical(val, "NA")) next

    dt <- tryCatch(as.Date(val), error = function(e) NA)
    if (is.na(dt)) next

    vertical_lines <- rbind(
      vertical_lines,
      data.frame(
        date = dt,
        label = date_labels[i],
        color = date_colors[i],
        stringsAsFactors = FALSE
      )
    )

    if (date_cols[i] == "KPRB_SAMPLE_DATE") {
      enrollment_date <- dt
    }
  }

  date_cols_all <- names(patient_meta)[
    grepl("DATE", names(patient_meta), ignore.case = TRUE) &
      names(patient_meta) != "KPRB_SAMPLE_DATE"
  ]

  subtitle_pieces <- vapply(date_cols_all, function(col) {
    val <- patient_meta[[col]][1]
    if (is.null(val) || is.na(val) || identical(val, "NA")) return(NA_character_)

    d <- tryCatch(as.Date(val), error = function(e) NA)
    if (is.na(d)) return(NA_character_)

    paste0(col, " = ", format(d, "%Y-%m-%d"))
  }, character(1))

  subtitle_text <- paste(stats::na.omit(subtitle_pieces), collapse = "\n")
  if (identical(subtitle_text, "")) subtitle_text <- "No clinical DATE columns available"

  # ---------- clinical layer helper ----------
  make_clin_layers <- function(var_label, x_limits, specs, use_thresholds_flag = TRUE) {
    if (!isTRUE(use_thresholds_flag)) return(list())

    spec <- specs[[var_label]]
    if (is.null(spec)) return(list())

    layers <- list()

    if (!is.null(spec$bands) && nrow(spec$bands)) {
      band_df <- transform(spec$bands, xmin = x_limits[1], xmax = x_limits[2])

      layers <- c(
        layers,
        list(
          ggplot2::geom_rect(
            data = band_df,
            ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = fill),
            inherit.aes = FALSE,
            alpha = 0.12
          ),
          ggplot2::scale_fill_identity()
        )
      )

      if (isTRUE(spec$label_bands)) {
        label_df <- band_df

        if (!is.null(spec$ylabel_pos)) {
          if (length(spec$ylabel_pos) != nrow(band_df)) {
            warning(sprintf("ylabel_pos mismatch for '%s'; using midpoints.", var_label))
            spec$ylabel_pos <- NULL
          }
        }

        if (is.null(spec$ylabel_pos)) {
          label_df$y <- mapply(function(ymin, ymax) {
            if (is.infinite(ymin) && is.infinite(ymax)) return(NA_real_)
            if (is.infinite(ymin)) return(ymax - 0.1 * abs(ymax))
            if (is.infinite(ymax)) return(ymin + 0.1 * abs(ymin))
            (ymin + ymax) / 2
          }, label_df$ymin, label_df$ymax)
        } else {
          label_df$y <- spec$ylabel_pos
        }

        label_df$x <- band_df$xmax - 1
        label_df <- label_df[is.finite(label_df$y), , drop = FALSE]

        if (nrow(label_df)) {
          layers <- c(
            layers,
            list(
              ggplot2::geom_text(
                data = label_df,
                ggplot2::aes(x = x, y = y, label = label),
                hjust = 1.02,
                vjust = 0.5,
                size = 3,
                inherit.aes = FALSE
              )
            )
          )
        }
      }
    }

    if (!is.null(spec$lines) && length(spec$lines)) {
      layers <- c(
        layers,
        list(
          ggplot2::geom_hline(
            yintercept = spec$lines,
            linetype = 3,
            linewidth = 0.8,
            color = "grey30"
          )
        )
      )
    }

    layers
  }

  coord_xy <- function(xm, ym = NULL, expand_xy = TRUE) {
    exp <- if (isFALSE(expand_xy)) FALSE else TRUE

    if (!is.null(xm) && length(xm) == 2 && !is.null(ym) && length(ym) == 2) {
      return(ggplot2::coord_cartesian(xlim = xm, ylim = ym, expand = exp))
    }

    if (!is.null(xm) && length(xm) == 2) {
      return(ggplot2::coord_cartesian(xlim = xm, expand = exp))
    }

    if (!is.null(ym) && length(ym) == 2) {
      return(ggplot2::coord_cartesian(ylim = ym, expand = exp))
    }

    NULL
  }

  # ---------- filter to patient ----------
  bmi_p   <- filter_track(bmi)
  bp_p    <- filter_track(bp)
  a1c_p   <- filter_track(a1c)
  ldl_p   <- filter_track(ldl)
  hdl_p   <- filter_track(hdl)
  glu_p   <- filter_track(gluran)
  glu_f_p <- filter_track(glu_f)
  crt_p   <- filter_track(creatinine)
  egfr_p  <- filter_track(egfr)
  alt_p   <- filter_track(alt)
  fib4_p  <- filter_track(fib4)

  rx_p <- if (is.null(rx)) {
    empty_rx()
  } else {
    rx %>%
      dplyr::filter(.data$StudyID == patient_id) %>%
      dplyr::mutate(CLIN_DATE = as.Date(.data$CLIN_DATE))
  }

  # ---------- x-limits & breaks ----------
  if (!is.null(date_start) && !is.null(date_end)) {
    x_limits <- c(as.Date(date_start), as.Date(date_end))
  } else {
    all_dates <- c(
      bmi_p$CLIN_DATE,
      bp_p$CLIN_DATE,
      a1c_p$CLIN_DATE,
      ldl_p$CLIN_DATE,
      hdl_p$CLIN_DATE,
      glu_f_p$CLIN_DATE,
      glu_p$CLIN_DATE,
      crt_p$CLIN_DATE,
      egfr_p$CLIN_DATE,
      alt_p$CLIN_DATE,
      fib4_p$CLIN_DATE,
      rx_p$CLIN_DATE,
      vertical_lines$date
    )

    all_dates <- all_dates[!is.na(all_dates)]

    x_limits <- if (length(all_dates)) {
      c(min(all_dates, na.rm = TRUE) - 30, max(all_dates, na.rm = TRUE) + 30)
    } else {
      c(Sys.Date() - 365, Sys.Date() + 30)
    }
  }

  start_y <- as.integer(format(x_limits[1], "%Y"))
  end_y   <- as.integer(format(x_limits[2], "%Y"))
  x_breaks_year <- as.Date(paste0(seq(start_y, end_y, by = 1), "-01-01"))

  # ---------- clinical panel ----------
  line_panel <- function(dat, ylab, spec_key = ylab, x_breaks = NULL, vlines = NULL, x_limits = NULL,
                         specs = clinical_specs, show_x_axis = FALSE) {
    x_minor_breaks <- if (!is.null(x_limits)) {
      seq(
        from = as.Date(format(x_limits[1], "%Y-%m-01")),
        to   = as.Date(format(x_limits[2], "%Y-%m-01")),
        by   = "6 months"
      )
    } else {
      NULL
    }

    if (nrow(dat) == 0 && isTRUE(drop_empty_tracks)) return(NULL)

    spec <- specs[[spec_key]]
    ylab_ylim <- if (!is.null(spec)) spec$y_limits else NULL

    if (nrow(dat) == 0) {
      p <- ggplot2::ggplot() +
        ggplot2::labs(y = ylab, x = NULL, subtitle = paste0(spec_key, " (no data)")) +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::theme(
          panel.grid = ggplot2::element_blank(),
          axis.text.x  = if (show_x_axis) ggplot2::element_text() else ggplot2::element_blank(),
          axis.ticks.x = if (show_x_axis) ggplot2::element_line() else ggplot2::element_blank(),
          axis.line.x  = if (show_x_axis) ggplot2::element_line() else ggplot2::element_blank(),
          axis.title.y = ggplot2::element_text(angle = 0, vjust = 0.5, hjust = 1, lineheight = 0.9)
        )

      if (!is.null(x_limits) && show_x_axis) {
        p <- p +
          ggplot2::scale_x_date(
            breaks = x_breaks %||% NULL,
            date_labels = "%Y",
            minor_breaks = x_minor_breaks
          ) +
          ggplot2::theme(
            panel.grid.minor = ggplot2::element_line(color = "grey90", linewidth = 0.3)
          )
      }

      if (!is.null(vlines) && nrow(vlines) > 0) {
        for (i in seq_len(nrow(vlines))) {
          p <- p + ggplot2::geom_vline(
            xintercept = vlines$date[i],
            linetype = 2,
            color = vlines$color[i],
            alpha = 0.7
          )
        }
      }

      if (!is.null(x_limits)) {
        cxy <- coord_xy(x_limits, ylab_ylim, expand_xy = TRUE)
        if (!is.null(cxy)) p <- p + cxy
      }

      return(p)
    }

    base_layers <- make_clin_layers(
      spec_key,
      x_limits,
      specs,
      use_thresholds_flag = use_thresholds
    )

    p <- ggplot2::ggplot() +
      base_layers +
      ggplot2::geom_line(
        data = dat,
        ggplot2::aes(.data$CLIN_DATE, .data$CLIN_VAL),
        na.rm = TRUE,
        linewidth = 0.5,
        color = "steelblue"
      ) +
      ggplot2::geom_point(
        data = dat,
        ggplot2::aes(.data$CLIN_DATE, .data$CLIN_VAL),
        size = 1.2,
        na.rm = TRUE,
        alpha = 0.5
      ) +
      ggplot2::labs(y = ylab, x = NULL) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_blank(),
        plot.subtitle = ggplot2::element_blank(),
        axis.text.x  = if (show_x_axis) ggplot2::element_text() else ggplot2::element_blank(),
        axis.ticks.x = if (show_x_axis) ggplot2::element_line() else ggplot2::element_blank(),
        axis.line.x  = if (show_x_axis) ggplot2::element_line() else ggplot2::element_blank(),
        axis.title.y = ggplot2::element_text(angle = 0, vjust = 0.5, hjust = 1, lineheight = 0.9)
      )

    if (!is.null(vlines) && nrow(vlines) > 0) {
      for (i in seq_len(nrow(vlines))) {
        p <- p + ggplot2::geom_vline(
          xintercept = vlines$date[i],
          linetype = 2,
          color = vlines$color[i],
          alpha = 0.7
        )
      }
    }

    if (!is.null(x_breaks) && show_x_axis) {
      p <- p +
        ggplot2::scale_x_date(
          breaks = x_breaks,
          date_labels = "%Y",
          minor_breaks = x_minor_breaks
        ) +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_line(color = "grey90", linewidth = 0.3)
        )
    }

    if (!is.null(x_limits)) {
      cxy <- coord_xy(x_limits, ylab_ylim, expand_xy = !is.null(ylab_ylim))
      if (!is.null(cxy)) p <- p + cxy
    }

    p
  }

  # ---------- risk projection panel ----------
  risk_panel <- function(risk_df, ylab, enrollment_date,
                         x_breaks = NULL, vlines = NULL, x_limits = NULL,
                         show_x_axis = FALSE, horizons = hz, cohort_label = "") {
    x_minor_breaks <- if (!is.null(x_limits)) {
      seq(
        from = as.Date(format(x_limits[1], "%Y-%m-01")),
        to   = as.Date(format(x_limits[2], "%Y-%m-01")),
        by   = "6 months"
      )
    } else {
      NULL
    }

    pid_chr <- trimws(as.character(patient_id))

    patient_risk <- risk_df %>%
      dplyr::mutate(.study_id_chr = trimws(as.character(.data$StudyID))) %>%
      dplyr::filter(.data$.study_id_chr == pid_chr)

    missing_enrollment <- is.null(enrollment_date) ||
      (inherits(enrollment_date, "Date") && length(enrollment_date) == 1L && is.na(enrollment_date))

    missing_patient <- nrow(patient_risk) == 0

    silent_drop <- isTRUE(drop_empty_tracks) && !isTRUE(show_placeholder_risk_panels)

    if ((missing_enrollment || missing_patient) && silent_drop) {
      return(NULL)
    }

    if (missing_enrollment || missing_patient) {
      sub_msg <- if (missing_enrollment) {
        sprintf(
          paste0(
            "Enrollment missing for %s. ",
            "Set KPRB_SAMPLE_DATE / KPRB sample date to anchor risk timelines."
          ),
          pid_chr
        )
      } else {
        sprintf(
          paste0(
            "\"%s\" not in SOC/NMR OOF table (%s; n modeled = %d). ",
            "These RDS usually cover Kaiser-5k / primary-care modeling—not all KP55k IDs."
          ),
          pid_chr,
          cohort_label,
          nrow(risk_df)
        )
      }

      p <- ggplot2::ggplot() +
        ggplot2::labs(y = ylab, x = NULL, subtitle = sub_msg) +
        ggplot2::theme_bw(base_size = 11) +
        ggplot2::theme(
          plot.subtitle = ggplot2::element_text(
            hjust = 0.5,
            colour = "firebrick",
            size = 11,
            margin = ggplot2::margin(t = 24, b = 8),
            lineheight = 1.05
          ),
          panel.grid = ggplot2::element_blank(),
          axis.text.x  = if (show_x_axis) ggplot2::element_text() else ggplot2::element_blank(),
          axis.ticks.x = if (show_x_axis) ggplot2::element_line() else ggplot2::element_blank(),
          axis.line.x  = if (show_x_axis) ggplot2::element_line() else ggplot2::element_blank(),
          axis.title.y = ggplot2::element_text(angle = 0, vjust = 0.5, hjust = 1)
        )

      if (!is.null(vlines) && nrow(vlines) > 0) {
        for (i in seq_len(nrow(vlines))) {
          p <- p + ggplot2::geom_vline(
            xintercept = vlines$date[i],
            linetype = 2,
            color = vlines$color[i],
            alpha = 0.7
          )
        }
      }

      if (!is.null(x_limits) && show_x_axis) {
        p <- p +
          ggplot2::scale_x_date(
            breaks = x_breaks %||% NULL,
            date_labels = "%Y",
            minor_breaks = x_minor_breaks
          ) +
          ggplot2::theme(
            panel.grid.minor = ggplot2::element_line(color = "grey90", linewidth = 0.3)
          )
      }

      if (!is.null(x_limits)) {
        cxy <- coord_xy(x_limits, c(0, 100), expand_xy = TRUE)
        if (!is.null(cxy)) p <- p + cxy
      }

      return(p)
    }

    patient_risk <- patient_risk[1L, , drop = FALSE]

    horizons <- horizons[!is.na(horizons) & horizons > 0]

    risk_segments <- data.frame(
      risk_col = sprintf("risk_%dy", as.integer(round(horizons))),
      years = horizons,
      stringsAsFactors = FALSE
    )

    risk_to_color <- function(risk) {
      if (is.na(risk)) return("#CCCCCC")

      risk <- pmin(pmax(risk, 0), 1)

      ramp <- grDevices::colorRamp(c("#2ecc71", "#f7e134", "#e74c3c"))
      col <- ramp(risk)

      grDevices::rgb(col[1] / 255, col[2] / 255, col[3] / 255)
    }

    seg_data <- data.frame()

    for (i in seq_len(nrow(risk_segments))) {
      col_name <- risk_segments$risk_col[i]
      if (!col_name %in% names(patient_risk)) next

      risk_val <- suppressWarnings(as.numeric(patient_risk[[col_name]][1]))
      if (is.na(risk_val)) next

      # Expected input is probability scale, 0 to 1.
      # If accidentally supplied as percent scale, tolerate 0 to 100.
      if (risk_val > 1 && risk_val <= 100) {
        risk_val <- risk_val / 100
      }

      risk_val <- pmin(pmax(risk_val, 0), 1)

      end_date <- enrollment_date + 365.25 * risk_segments$years[i]

      seg_data <- rbind(
        seg_data,
        data.frame(
          xmin = enrollment_date,
          xmax = end_date,
          y = risk_val * 100,
          label = sprintf("%dyr: %.1f%%", risk_segments$years[i], risk_val * 100),
          color = risk_to_color(risk_val),
          years = risk_segments$years[i],
          stringsAsFactors = FALSE
        )
      )
    }

    if (nrow(seg_data) == 0) {
      if (isTRUE(drop_empty_tracks)) return(NULL)

      return(
        ggplot2::ggplot() +
          ggplot2::labs(y = ylab, x = NULL, subtitle = paste0(ylab, " (no valid risk data)")) +
          ggplot2::theme_bw(base_size = 12) +
          ggplot2::theme(
            panel.grid = ggplot2::element_blank(),
            axis.text.x  = if (show_x_axis) ggplot2::element_text() else ggplot2::element_blank(),
            axis.ticks.x = if (show_x_axis) ggplot2::element_line() else ggplot2::element_blank(),
            axis.line.x  = if (show_x_axis) ggplot2::element_line() else ggplot2::element_blank(),
            axis.title.y = ggplot2::element_text(angle = 0, vjust = 0.5, hjust = 1)
          )
      )
    }

    p <- ggplot2::ggplot() +
      ggplot2::geom_segment(
        data = seg_data,
        ggplot2::aes(x = xmin, xend = xmax, y = y, yend = y, color = color),
        linewidth = 1.5,
        alpha = 0.8
      ) +
      ggplot2::scale_color_identity() +
      ggplot2::geom_text(
        data = seg_data,
        ggplot2::aes(x = xmax, y = y, label = label),
        hjust = 0,
        vjust = 0.5,
        size = 3,
        inherit.aes = FALSE,
        nudge_x = 20
      ) +
      ggplot2::labs(y = ylab, x = NULL) +
      ggplot2::scale_y_continuous(
        breaks = seq(0, 100, 25),
        labels = function(x) paste0(x, "%"),
        expand = ggplot2::expansion(mult = c(0, 0.08))
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        panel.grid.minor = ggplot2::element_blank(),
        plot.title = ggplot2::element_blank(),
        plot.subtitle = ggplot2::element_blank(),
        axis.text.x  = if (show_x_axis) ggplot2::element_text() else ggplot2::element_blank(),
        axis.ticks.x = if (show_x_axis) ggplot2::element_line() else ggplot2::element_blank(),
        axis.line.x  = if (show_x_axis) ggplot2::element_line() else ggplot2::element_blank(),
        axis.title.y = ggplot2::element_text(angle = 0, vjust = 0.5, hjust = 1)
      )

    if (!is.null(vlines) && nrow(vlines) > 0) {
      for (i in seq_len(nrow(vlines))) {
        p <- p + ggplot2::geom_vline(
          xintercept = vlines$date[i],
          linetype = 2,
          color = vlines$color[i],
          alpha = 0.7
        )
      }
    }

    if (!is.null(x_breaks) && show_x_axis) {
      p <- p +
        ggplot2::scale_x_date(
          breaks = x_breaks,
          date_labels = "%Y",
          minor_breaks = x_minor_breaks
        ) +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_line(color = "grey90", linewidth = 0.3)
        )
    }

    if (!is.null(x_limits)) {
      cxy <- coord_xy(x_limits, c(0, 100), expand_xy = TRUE)
      if (!is.null(cxy)) p <- p + cxy
    }

    p
  }

  # ---------- medication panel ----------
  rx_panel <- function(rx_dat, x_breaks = NULL, x_limits = NULL, vlines = NULL, show_x_axis = TRUE) {
    x_minor_breaks <- if (!is.null(x_breaks) && length(x_breaks) > 0) {
      seq(
        from = as.Date(format(min(x_breaks), "%Y-%m-01")),
        to   = as.Date(format(max(x_breaks), "%Y-%m-01")),
        by   = "6 months"
      )
    } else {
      NULL
    }

    if (nrow(rx_dat) == 0) {
      p <- ggplot2::ggplot() +
        ggplot2::theme_bw(base_size = 12) +
        ggplot2::labs(subtitle = "RX (no data)") +
        ggplot2::theme(
          panel.grid = ggplot2::element_blank(),
          axis.text.x  = if (show_x_axis) ggplot2::element_text() else ggplot2::element_blank(),
          axis.ticks.x = if (show_x_axis) ggplot2::element_line() else ggplot2::element_blank(),
          axis.line.x  = if (show_x_axis) ggplot2::element_line() else ggplot2::element_blank()
        )

      if (!is.null(x_breaks) && show_x_axis) {
        p <- p +
          ggplot2::scale_x_date(
            breaks = x_breaks,
            date_labels = "%Y",
            minor_breaks = x_minor_breaks
          ) +
          ggplot2::theme(
            panel.grid.minor = ggplot2::element_line(color = "grey90", linewidth = 0.3)
          )
      }

      if (!is.null(vlines) && nrow(vlines) > 0) {
        for (i in seq_len(nrow(vlines))) {
          p <- p + ggplot2::geom_vline(
            xintercept = vlines$date[i],
            linetype = 2,
            color = vlines$color[i],
            alpha = 0.7
          )
        }
      }

      if (!is.null(x_limits)) {
        cxy <- coord_xy(x_limits, NULL, expand_xy = TRUE)
        if (!is.null(cxy)) p <- p + cxy
      }

      return(p)
    }

    rx_processed <- rx_dat %>%
      dplyr::mutate(CLIN_DATE = as.Date(.data$CLIN_DATE)) %>%
      dplyr::group_by(.data$GENERIC_NM, .data$CLIN_DATE) %>%
      dplyr::distinct(.data$CLIN_DATE, .keep_all = TRUE) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(
        med_rank = as.numeric(
          factor(.data$GENERIC_NM, levels = sort(unique(.data$GENERIC_NM)))
        )
      ) %>%
      dplyr::group_by(.data$CLIN_DATE, .data$GENERIC_NM) %>%
      dplyr::mutate(med_rank = .data$med_rank + (dplyr::row_number() - 1) * 0.3) %>%
      dplyr::ungroup()

    unique_meds <- sort(unique(rx_processed$GENERIC_NM))

    p <- ggplot2::ggplot(
      rx_processed,
      ggplot2::aes(.data$CLIN_DATE, .data$med_rank, color = .data$GENERIC_NM)
    ) +
      ggplot2::geom_hline(
        data = data.frame(med_rank = seq_along(unique_meds)),
        ggplot2::aes(yintercept = .data$med_rank),
        color = "gray90",
        linewidth = 0.3,
        inherit.aes = FALSE
      ) +
      ggplot2::geom_point(size = 2, alpha = 0.7) +
      ggplot2::labs(y = "Rx", x = "Date") +  # was "Medications"; REVERT that string to restore the old label
      ggplot2::scale_y_continuous(
        breaks = seq_along(unique_meds),
        labels = unique_meds
      ) +
      ggplot2::theme_bw(base_size = 12) +
      ggplot2::theme(
        panel.grid = ggplot2::element_blank(),
        plot.title = ggplot2::element_blank(),
        plot.subtitle = ggplot2::element_blank(),
        axis.title.y = ggplot2::element_text(angle = 0, vjust = 0.5, hjust = 1),
        axis.text.y = ggplot2::element_text(size = 10),
        axis.ticks.y = ggplot2::element_line(linewidth = 0.3),
        legend.position = "none"
      )

    if (!is.null(vlines) && nrow(vlines) > 0) {
      for (i in seq_len(nrow(vlines))) {
        p <- p + ggplot2::geom_vline(
          xintercept = vlines$date[i],
          linetype = 2,
          color = vlines$color[i],
          alpha = 0.7
        )
      }
    }

    if (!is.null(x_breaks) && show_x_axis) {
      p <- p +
        ggplot2::scale_x_date(
          breaks = x_breaks,
          date_labels = "%Y",
          minor_breaks = x_minor_breaks
        ) +
        ggplot2::theme(
          panel.grid.minor = ggplot2::element_line(color = "grey90", linewidth = 0.3)
        )
    }

    y_rx_lim <- c(0.5, length(unique_meds) + 0.5)
    # cxy <- coord_xy(if (!is.null(x_limits)) x_limits else NULL, y_rx_lim, expand_xy = FALSE)
    cxy <- coord_xy(if (!is.null(x_limits)) x_limits else NULL, y_rx_lim, expand_xy = TRUE)
    if (!is.null(cxy)) p <- p + cxy

    p
  }

  # ---------- decide which panel shows the x-axis ----------
  presence <- c(
    BMI      = !is.null(bmi) && (nrow(bmi_p) > 0 || !drop_empty_tracks),
    SBP      = !is.null(bp) && (nrow(bp_p) > 0 || !drop_empty_tracks),
    A1C      = !is.null(a1c) && (nrow(a1c_p) > 0 || !drop_empty_tracks),
    LDL      = !is.null(ldl) && (nrow(ldl_p) > 0 || !drop_empty_tracks),
    HDL      = !is.null(hdl) && (nrow(hdl_p) > 0 || !drop_empty_tracks),
    FGLU     = !is.null(glu_f) && (nrow(glu_f_p) > 0 || !drop_empty_tracks),
    RGLU     = !is.null(gluran) && (nrow(glu_p) > 0 || !drop_empty_tracks),
    CRT      = !is.null(creatinine) && (nrow(crt_p) > 0 || !drop_empty_tracks),
    EGFR     = !is.null(egfr) && (nrow(egfr_p) > 0 || !drop_empty_tracks),
    ALT      = !is.null(alt) && (nrow(alt_p) > 0 || !drop_empty_tracks),
    FIB4     = !is.null(fib4) && (nrow(fib4_p) > 0 || !drop_empty_tracks),
    RISK_SOC = !is.null(df_soc),
    RISK_NMR = !is.null(df_nmr),
    RX       = !is.null(rx) && (nrow(rx_p) > 0 || !drop_empty_tracks)
  )

  show_axis <- stats::setNames(rep(FALSE, length(presence)), names(presence))
  if (any(presence)) show_axis[max(which(presence))] <- TRUE

  # ---------- build panels ----------
  p_bmi <- if (!is.null(bmi)) {
    line_panel(
      bmi_p, "BMI",
      vlines = vertical_lines,
      x_limits = x_limits,
      x_breaks = x_breaks_year,
      show_x_axis = show_axis[["BMI"]]
    )
  } else NULL

  p_bp <- if (!is.null(bp)) {
    line_panel(
      bp_p, "SBP",
      vlines = vertical_lines,
      x_limits = x_limits,
      x_breaks = x_breaks_year,
      show_x_axis = show_axis[["SBP"]]
    )
  } else NULL

  p_a1c <- if (!is.null(a1c)) {
    line_panel(
      a1c_p, "HbA1c (%)",
      vlines = vertical_lines,
      x_limits = x_limits,
      x_breaks = x_breaks_year,
      show_x_axis = show_axis[["A1C"]]
    )
  } else NULL

  p_ldl <- if (!is.null(ldl)) {
    line_panel(
      ldl_p, "LDL-C (mg/dL)",
      vlines = vertical_lines,
      x_limits = x_limits,
      x_breaks = x_breaks_year,
      show_x_axis = show_axis[["LDL"]]
    )
  } else NULL

  p_hdl <- if (!is.null(hdl)) {
    line_panel(
      hdl_p, "HDL-C (mg/dL)",
      vlines = vertical_lines,
      x_limits = x_limits,
      x_breaks = x_breaks_year,
      show_x_axis = show_axis[["HDL"]]
    )
  } else NULL

  p_glu_f <- if (!is.null(glu_f)) {
    line_panel(
      glu_f_p, "Fasting Glucose (mg/dL)",
      vlines = vertical_lines,
      x_limits = x_limits,
      x_breaks = x_breaks_year,
      show_x_axis = show_axis[["FGLU"]]
    )
  } else NULL

  p_glu <- if (!is.null(gluran)) {
    line_panel(
      glu_p, "Rand Gluc (mg/dL)",
      vlines = vertical_lines,
      x_limits = x_limits,
      x_breaks = x_breaks_year,
      show_x_axis = show_axis[["RGLU"]]
    )
  } else NULL

  p_crt <- if (!is.null(creatinine)) {
    line_panel(
      crt_p, "Creatinine (mg/dL)",
      vlines = vertical_lines,
      x_limits = x_limits,
      x_breaks = x_breaks_year,
      show_x_axis = show_axis[["CRT"]]
    )
  } else NULL

  p_egfr <- if (!is.null(egfr)) {
    line_panel(
      egfr_p, EGFR_YLAB, spec_key = EGFR_SPEC,
      vlines = vertical_lines,
      x_limits = x_limits,
      x_breaks = x_breaks_year,
      show_x_axis = show_axis[["EGFR"]]
    )
  } else NULL

  p_alt <- if (!is.null(alt)) {
    line_panel(
      alt_p, "ALT (U/L)",
      vlines = vertical_lines,
      x_limits = x_limits,
      x_breaks = x_breaks_year,
      show_x_axis = show_axis[["ALT"]]
    )
  } else NULL

  p_fib4 <- if (!is.null(fib4)) {
    line_panel(
      fib4_p, "FIB-4",
      vlines = vertical_lines,
      x_limits = x_limits,
      x_breaks = x_breaks_year,
      show_x_axis = show_axis[["FIB4"]]
    )
  } else NULL

  p_risk_soc <- if (!is.null(df_soc)) {
    risk_panel(
      df_soc,
      paste0(risk_name, " Risk\n(SOC)\n%"),
      enrollment_date,
      x_breaks = x_breaks_year,
      vlines = vertical_lines,
      x_limits = x_limits,
      show_x_axis = show_axis[["RISK_SOC"]],
      cohort_label = "SOC (marker_oof)"
    )
  } else NULL

  p_risk_nmr <- if (!is.null(df_nmr)) {
    risk_panel(
      df_nmr,
      paste0(risk_name, " Risk\n(NMR)\n%"),
      enrollment_date,
      x_breaks = x_breaks_year,
      vlines = vertical_lines,
      x_limits = x_limits,
      show_x_axis = show_axis[["RISK_NMR"]],
      cohort_label = "NMR (marker_nmr_oof)"
    )
  } else NULL

  p_rx <- if (!is.null(rx)) {
    rx_panel(
      rx_p,
      x_breaks = x_breaks_year,
      x_limits = x_limits,
      vlines = vertical_lines,
      show_x_axis = show_axis[["RX"]]
    )
  } else NULL

  # 2026-09-17: Patchwork left-aligns every y-axis title to the longest
  # label ("Medications"), which leaves a gap between short labels (BMI,
  # SBP, HbA1c, risk, ...) and the plot. Free the left axis *title* on
  # non-medication panels so those labels sit next to their own axes.
  # Panels and tick labels stay aligned; only the title position changes.
  # REVERT: delete `tight_left_ylab()` and pass the plots through as-is
  # (e.g. `panels <- list(p_bmi, p_bp, ..., p_rx)`).
  tight_left_ylab <- function(p) {
    if (is.null(p)) return(NULL)
    patchwork::free(p, type = "label", side = "l")
  }

  panels <- list(
    tight_left_ylab(p_bmi),
    tight_left_ylab(p_bp),
    tight_left_ylab(p_a1c),
    tight_left_ylab(p_ldl),
    tight_left_ylab(p_hdl),
    tight_left_ylab(p_glu_f),
    tight_left_ylab(p_glu),
    tight_left_ylab(p_crt),
    tight_left_ylab(p_egfr),
    tight_left_ylab(p_alt),
    tight_left_ylab(p_fib4),
    tight_left_ylab(p_risk_soc),
    tight_left_ylab(p_risk_nmr),
    p_rx  # keep "Medications" aligned with drug-name ticks
  )

  panels <- Filter(Negate(is.null), panels)

  if (!length(panels)) {
    stop("No panels to plot. Check patient_id, input tracks, and drop_empty_tracks.")
  }

  layout_heights <- rep(1.2, length(panels))

  combined_plot <- (
    Reduce(`/`, panels) +
      patchwork::plot_layout(heights = layout_heights, guides = "collect")
  ) &
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        face = "bold",
        size = 14,
        hjust = 0.5,
        margin = ggplot2::margin(b = 20)
      )
    )

  enroll_str <- {
    idx <- which(vertical_lines$label == "Enrollment")
    if (length(idx) == 0) "N/A" else format(vertical_lines$date[idx][1], "%Y-%m-%d")
  }

  age_str <- {
    v <- if ("age_at_sample" %in% names(patient_meta)) patient_meta$age_at_sample[1] else NA

    if (is.null(v) || is.na(v) || identical(v, "NA")) {
      "N/A"
    } else if (is.numeric(v)) {
      sprintf("%.1f", v)
    } else {
      as.character(v)
    }
  }

  sex_str <- {
    v <- if ("gender" %in% names(patient_meta)) patient_meta$gender[1] else NA

    if (is.null(v) || is.na(v) || identical(v, "NA") || v == "") {
      "N/A"
    } else {
      as.character(v)
    }
  }

  missing_tracks <- c(
    if (!is.null(bmi) && nrow(bmi_p) == 0) "BMI",
    if (!is.null(bp) && nrow(bp_p) == 0) "SBP",
    if (!is.null(a1c) && nrow(a1c_p) == 0) "HbA1c (%)",
    if (!is.null(ldl) && nrow(ldl_p) == 0) "LDL-C (mg/dL)",
    if (!is.null(hdl) && nrow(hdl_p) == 0) "HDL-C (mg/dL)",
    if (!is.null(glu_f) && nrow(glu_f_p) == 0) "Fasting Glucose (mg/dL)",
    if (!is.null(gluran) && nrow(glu_p) == 0) "Rand Gluc (mg/dL)",
    if (!is.null(creatinine) && nrow(crt_p) == 0) "Creatinine (mg/dL)",
    if (!is.null(egfr) && nrow(egfr_p) == 0) EGFR_YLAB,
    if (!is.null(alt) && nrow(alt_p) == 0) "ALT (U/L)",
    if (!is.null(fib4) && nrow(fib4_p) == 0) "FIB-4"
  )

  if (length(missing_tracks)) {
    subtitle_text <- paste0(
      subtitle_text,
      "\nMissing: ",
      paste(missing_tracks, collapse = ", ")
    )
  }

  combined_plot +
    patchwork::plot_annotation(
      title = paste0(
        patient_id,
        "  |  Enrollment: ", enroll_str,
        "  |  Age: ", age_str,
        "  |  Sex: ", sex_str
      ),
      subtitle = subtitle_text,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(
          hjust = 0.5,
          margin = ggplot2::margin(b = 4)
        ),
        plot.subtitle = ggplot2::element_text(
          hjust = 0.5,
          margin = ggplot2::margin(t = 2)
        )
      )
    )
}