# C-index comparison of SOC, PRS, and NMR out-of-fold models.
# Reads {prefix}{base|clin|marker}[_prs][_nmr]_oof.rds and plots SOC-anchored dumbbells.

collect_prs_cindex <- function(dir, prefix) {
  groups   <- c("base", "clin", "marker")
  suffixes <- c("", "_prs", "_nmr", "_prs_nmr")

  manifest <- expand.grid(group = groups, suffix = suffixes, stringsAsFactors = FALSE)
  manifest$filename <- paste0(prefix, manifest$group, manifest$suffix, "_oof.rds")
  manifest$type <- dplyr::case_when(
    manifest$suffix == ""         ~ "SOC",
    manifest$suffix == "_prs"     ~ "SOC+PRS",
    manifest$suffix == "_nmr"     ~ "SOC+NMR",
    manifest$suffix == "_prs_nmr" ~ "SOC+PRS+NMR"
  )

  rows <- vector("list", nrow(manifest))
  for (i in seq_len(nrow(manifest))) {
    path <- file.path(dir, manifest$filename[i])
    if (!file.exists(path)) {
      warning("File not found, skipping: ", path)
      next
    }
    res <- readRDS(path)
    rows[[i]] <- data.frame(
      model    = manifest$group[i],
      type     = manifest$type[i],
      c_index  = res$c_index,
      ci_lower = res$ci_lower,
      ci_upper = res$ci_upper,
      stringsAsFactors = FALSE
    )
    rm(res)
  }

  df <- do.call(rbind, rows)
  df$model <- factor(df$model, levels = c("base", "clin", "marker"))
  df$type  <- factor(df$type, levels = c("SOC", "SOC+PRS", "SOC+NMR", "SOC+PRS+NMR"))
  df
}

plot_prs_nmr_cindex_dumbbell <- function(
  df,
  title        = NULL,
  subtitle     = NULL,
  xlim         = NULL,
  model_labels = c("base" = "BASE", "clin" = "CLIN", "marker" = "MARKER"),
  type_order   = c("SOC+PRS+NMR", "SOC+NMR", "SOC+PRS"),
  use_ggsci    = TRUE
) {
  required_cols <- c("model", "type", "c_index", "ci_lower", "ci_upper")
  missing <- setdiff(required_cols, colnames(df))
  if (length(missing) > 0) stop("Missing columns: ", paste(missing, collapse = ", "))

  matched_keys <- intersect(names(model_labels), unique(as.character(df$model)))
  model_levels <- unname(model_labels[matched_keys])

  df <- df %>%
    dplyr::mutate(model = dplyr::recode(as.character(model), !!!model_labels))
  df$model <- factor(df$model, levels = model_levels)

  df_soc <- df %>%
    dplyr::filter(type == "SOC") %>%
    dplyr::select(model, soc_c = c_index, soc_lo = ci_lower, soc_hi = ci_upper)

  df_cmp <- df %>%
    dplyr::filter(type != "SOC") %>%
    dplyr::left_join(df_soc, by = "model")

  present_types <- intersect(type_order, unique(as.character(df_cmp$type)))
  df_cmp$type <- factor(df_cmp$type, levels = rev(present_types))

  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(
      data = df_cmp,
      ggplot2::aes(y = type, yend = type, x = soc_c, xend = c_index),
      color = "grey72", linewidth = 0.65, linetype = "dotted"
    ) +
    ggplot2::geom_errorbarh(
      data = df_cmp,
      ggplot2::aes(y = type, x = soc_c, xmin = soc_lo, xmax = soc_hi),
      color = "black", height = 0.4, linewidth = 0.6
    ) +
    ggplot2::geom_errorbarh(
      data = df_cmp,
      ggplot2::aes(y = type, x = c_index, xmin = ci_lower, xmax = ci_upper),
      color = "black", height = 0.4, linewidth = 0.6
    ) +
    ggplot2::geom_point(
      data = df_cmp,
      ggplot2::aes(y = type, x = soc_c),
      color = "grey45", size = 1, shape = 1.5
    ) +
    ggplot2::geom_point(
      data = df_cmp,
      ggplot2::aes(y = type, x = c_index, color = type),
      size = 1, shape = 1.5
    ) +
    ggplot2::facet_grid(. ~ model, scales = "free_x", space = "free_x", switch = "y") +
    ggplot2::labs(x = "C-index (95% CI)", y = NULL, color = NULL, title = title, subtitle = subtitle) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::theme(
      panel.grid.major.y = ggplot2::element_blank(),
      panel.grid.minor   = ggplot2::element_blank(),
      panel.spacing.y    = ggplot2::unit(0.3, "lines"),
      strip.text.y       = ggplot2::element_text(angle = 0, face = "bold", size = 10),
      strip.background   = ggplot2::element_rect(fill = "grey93", colour = NA),
      plot.title         = ggplot2::element_text(hjust = 0.5, face = "bold", size = 12),
      plot.subtitle      = ggplot2::element_text(hjust = 0.5, size = 10),
      legend.position    = "bottom",
      legend.key.size    = ggplot2::unit(0.5, "lines")
    )

  if (use_ggsci && requireNamespace("ggsci", quietly = TRUE)) {
    p <- p + ggsci::scale_color_npg()
  } else {
    fallback <- c("#E64B35", "#4DBBD5", "#00A087")[seq_along(present_types)]
    p <- p + ggplot2::scale_color_manual(values = stats::setNames(fallback, rev(present_types)))
  }

  if (!is.null(xlim)) p <- p + ggplot2::coord_cartesian(xlim = xlim)
  p
}

plot_oof_prs_nmr_cindex <- function(df, disease = NULL, title = NULL, xlim = NULL, ...) {
  if (is.null(title) && !is.null(disease)) title <- disease
  plot_prs_nmr_cindex_dumbbell(df, title = title, xlim = xlim, ...)
}
