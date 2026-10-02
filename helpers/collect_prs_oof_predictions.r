# ============================================================================
# collect_prs_oof_predictions: Load RDS results one-at-a-time and build a
#   named list of tidy oof-prediction data frames
# ============================================================================
# Loads each of the 12 model RDS files for a given disease sequentially,
# extracts $oof_predictions (keeping StudyID and lp), annotates each with
# model / type / outcome, then immediately discards the full object before
# loading the next one.  The function returns a named list so the caller can
# inspect individual slices or rbind them all at once.
#
# Expected filename convention (inside `dir`):
#   {prefix}base_oof.rds          {prefix}clin_oof.rds
#   {prefix}base_prs_oof.rds      {prefix}clin_prs_oof.rds
#   {prefix}base_nmr_oof.rds      {prefix}clin_nmr_oof.rds
#   {prefix}base_prs_nmr_oof.rds  {prefix}clin_prs_nmr_oof.rds
#   ... (repeat for marker)
#
# Usage:
#   source(here("helpers/collect_prs_oof_predictions.r"))
#
#   t2d_list   <- collect_prs_oof_predictions(
#                   here("outputs/outputs_nofuture/prs_comparison"), "t2d_")
#   ckd_list   <- collect_prs_oof_predictions(
#                   here("outputs/outputs_nofuture/prs_comparison"), "ckd_")
#
#   # Combine all outcomes into one long data frame manually:
#   all_df <- do.call(rbind, c(t2d_list, ckd_list, nafld_list, ...))
#   rownames(all_df) <- NULL
# ============================================================================

#' Sequentially load 12 OOF result RDS files and return a named list of tidy
#' oof-prediction data frames.
#'
#' @param dir    Path to the directory containing the RDS files.
#' @param prefix Filename prefix that identifies the disease, e.g. "t2d_".
#'               The trailing underscore is stripped to derive the outcome label.
#'
#' @return A named list (length ≤ 12) of data frames.  Each element is named
#'   "{model}_{type}" (e.g. "base_SOC", "clin_SOC+PRS") and contains columns:
#'   \describe{
#'     \item{model}{Factor: base < clin < marker}
#'     \item{type}{Factor: SOC < SOC+PRS < SOC+NMR < SOC+PRS+NMR}
#'     \item{outcome}{Character: disease label derived from \code{prefix}}
#'     \item{StudyID}{Character: participant identifier from oof_predictions}
#'     \item{lp}{Numeric: linear predictor from oof_predictions}
#'     \item{time}{Numeric: follow-up time from oof_predictions}
#'     \item{status}{Integer: event indicator from oof_predictions}
#'   }
#'   The \code{time} and \code{status} columns are carried along so that
#'   downstream helpers (e.g. \code{compute_oof_risk_probability}) can fit
#'   the Cox baseline hazard at a chosen horizon.
collect_prs_oof_predictions <- function(dir, prefix) {

  # ---- derive the outcome label (strip trailing underscore) ----------------
  outcome_label <- sub("_$", "", prefix)

  # ---- build the 12-entry manifest -----------------------------------------
  groups   <- c("base", "clin", "marker")
  suffixes <- c("", "_prs", "_nmr", "_prs_nmr")   # SOC / +PRS / +NMR / +PRS+NMR

  manifest <- expand.grid(
    group  = groups,
    suffix = suffixes,
    stringsAsFactors = FALSE
  )

  manifest$filename <- paste0(prefix, manifest$group, manifest$suffix, "_oof.rds")

  manifest$type <- dplyr::case_when(
    manifest$suffix == ""         ~ "SOC",
    manifest$suffix == "_prs"     ~ "SOC+PRS",
    manifest$suffix == "_nmr"     ~ "SOC+NMR",
    manifest$suffix == "_prs_nmr" ~ "SOC+PRS+NMR"
  )

  manifest$list_key <- paste0(manifest$group, "_", manifest$type)

  # ---- iterate, loading one RDS at a time ----------------------------------
  result <- vector("list", nrow(manifest))
  names(result) <- manifest$list_key

  for (i in seq_len(nrow(manifest))) {
    path <- file.path(dir, manifest$filename[i])

    if (!file.exists(path)) {
      warning("File not found, skipping: ", path)
      next
    }

    res <- readRDS(path)

    oof <- res$oof_predictions
    rm(res)   # free memory before loading next object

    result[[i]] <- data.frame(
      model   = factor(manifest$group[i],  levels = c("base", "clin", "marker")),
      type    = factor(manifest$type[i],   levels = c("SOC", "SOC+PRS", "SOC+NMR", "SOC+PRS+NMR")),
      outcome = outcome_label,
      StudyID = oof$StudyID,
      lp      = oof$lp,
      time    = oof$time,
      status  = oof$status,
      stringsAsFactors = FALSE
    )
  }

  # Drop any NULLs left by skipped files while preserving names
  result[!vapply(result, is.null, logical(1))]
}
