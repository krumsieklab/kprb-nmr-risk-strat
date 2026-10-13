# ============================================================================
# Helper: Batch Delta C-Index Comparison via compareC
# ============================================================================
# Provides:
#   1. run_delta_cindex()
#      - single paired compareC test
#
#   2. strip_oof_result()
#      - keeps only StudyID, lp, time, status from large model result objects
#
#   3. run_delta_cindex_batch()
#      - runs a set of pairwise delta C-index comparisons
#      - optional caching
#      - optional lower-memory parallelization
#      - by default restricts all models to the shared complete-case
#        StudyID intersection before any comparisons
#
# Dependencies:
#   compareC
#
# Optional dependencies:
#   digest       - needed only if cache_dir is used
#   future       - needed only if parallel = TRUE
#   future.apply - needed only if parallel = TRUE
# ============================================================================

library(compareC)


# ----------------------------------------------------------------------------
# Restrict OOF predictions to a shared StudyID set
# ----------------------------------------------------------------------------

restrict_oof_to_ids <- function(res, ids, object_name = "result") {
  if (is.null(res$oof_predictions)) {
    stop(object_name, "$oof_predictions is missing.", call. = FALSE)
  }

  oof <- as.data.frame(res$oof_predictions, stringsAsFactors = FALSE)
  if (!"StudyID" %in% names(oof)) {
    stop(object_name, "$oof_predictions is missing StudyID.", call. = FALSE)
  }

  oof$StudyID <- as.character(oof$StudyID)
  ids <- as.character(ids)

  idx <- match(ids, oof$StudyID)
  if (anyNA(idx)) {
    stop(
      object_name, " is missing ", sum(is.na(idx)),
      " StudyID values required by the complete-case set.",
      call. = FALSE
    )
  }

  res$oof_predictions <- oof[idx, , drop = FALSE]
  rownames(res$oof_predictions) <- NULL
  res
}


get_common_oof_ids <- function(results, model_names = names(results)) {
  if (length(model_names) == 0L) {
    stop("No models supplied for complete-case StudyID intersection.", call. = FALSE)
  }

  id_lists <- lapply(model_names, function(nm) {
    res <- results[[nm]]
    if (is.null(res$oof_predictions)) {
      stop("Model '", nm, "' is missing $oof_predictions.", call. = FALSE)
    }
    as.character(res$oof_predictions$StudyID)
  })
  names(id_lists) <- model_names

  common_ids <- Reduce(intersect, id_lists)
  common_ids <- sort(unique(common_ids))

  if (length(common_ids) == 0L) {
    stop(
      "No overlapping StudyIDs across models: ",
      paste(model_names, collapse = ", "),
      call. = FALSE
    )
  }

  n_by_model <- vapply(id_lists, length, integer(1L))
  dropped_by_model <- n_by_model - length(common_ids)

  list(
    ids = common_ids,
    n = length(common_ids),
    n_by_model = n_by_model,
    dropped_by_model = dropped_by_model
  )
}

# ----------------------------------------------------------------------------
# Single comparison helper
# ----------------------------------------------------------------------------

run_delta_cindex <- function(
  res1,
  res2,
  label_model1,
  label_model2,
  id_col = "StudyID",
  lp_col = "lp",
  time_col = "time",
  status_col = "status",
  lp_higher_risk = TRUE,
  conf_level = 0.95
) {
  required_cols <- c(id_col, lp_col, time_col, status_col)

  prep_oof <- function(res, object_name) {
    if (is.null(res$oof_predictions)) {
      stop(object_name, "$oof_predictions is missing.", call. = FALSE)
    }

    oof <- as.data.frame(res$oof_predictions, stringsAsFactors = FALSE)

    missing_cols <- setdiff(required_cols, names(oof))
    if (length(missing_cols) > 0) {
      stop(
        object_name, "$oof_predictions is missing columns: ",
        paste(missing_cols, collapse = ", "),
        call. = FALSE
      )
    }

    oof <- oof[, required_cols, drop = FALSE]
    names(oof) <- c("StudyID", "lp", "time", "status")

    oof$StudyID <- as.character(oof$StudyID)

    if (anyDuplicated(oof$StudyID)) {
      dupes <- unique(oof$StudyID[duplicated(oof$StudyID)])
      stop(
        object_name, "$oof_predictions has duplicate StudyID values, e.g. ",
        paste(head(dupes, 5), collapse = ", "),
        call. = FALSE
      )
    }

    oof$lp <- as.numeric(oof$lp)
    oof$time <- as.numeric(oof$time)
    oof$status <- as.integer(as.character(oof$status))

    if (anyNA(oof[, c("StudyID", "lp", "time", "status")])) {
      stop(object_name, "$oof_predictions contains missing values.", call. = FALSE)
    }

    if (any(oof$time <= 0)) {
      stop(object_name, "$oof_predictions contains non-positive survival times.", call. = FALSE)
    }

    if (!all(oof$status %in% c(0L, 1L))) {
      stop(object_name, "$status must be coded 0 = censored, 1 = event.", call. = FALSE)
    }

    oof
  }

  oof1 <- prep_oof(res1, "res1")
  oof2 <- prep_oof(res2, "res2")

  idx <- match(oof1$StudyID, oof2$StudyID)

  if (nrow(oof1) != nrow(oof2) || anyNA(idx)) {
    missing_in_2 <- setdiff(oof1$StudyID, oof2$StudyID)
    missing_in_1 <- setdiff(oof2$StudyID, oof1$StudyID)

    stop(
      "OOF predictions are not for the exact same StudyID set.\n",
      "Missing in res2: ", paste(head(missing_in_2, 10), collapse = ", "), "\n",
      "Missing in res1: ", paste(head(missing_in_1, 10), collapse = ", "),
      call. = FALSE
    )
  }

  if (!isTRUE(all.equal(oof1$time, oof2$time[idx]))) {
    stop("Survival times differ between res1 and res2 after matching StudyID.", call. = FALSE)
  }

  if (!identical(oof1$status, oof2$status[idx])) {
    stop("Event statuses differ between res1 and res2 after matching StudyID.", call. = FALSE)
  }

  score1 <- oof1$lp
  score2 <- oof2$lp[idx]

  if (lp_higher_risk) {
    score1 <- -score1
    score2 <- -score2
  }

  # scoreY = model 2
  # scoreZ = model 1
  # est.diff_c = C(model 2) - C(model 1)
  cc <- compareC::compareC(
    timeX   = oof1$time,
    statusX = oof1$status,
    scoreY  = score2,
    scoreZ  = score1
  )

  var_delta <- cc$est.vardiff_c

  if (is.na(var_delta) || var_delta < 0) {
    stop("compareC returned an invalid variance for delta C.", call. = FALSE)
  }

  se_delta <- sqrt(var_delta)
  alpha <- 1 - conf_level
  zcrit <- qnorm(1 - alpha / 2)

  z_score <- unname(cc$zscore)

  # Useful when cc$pval underflows to 0
  log_p_value <- log(2) + pnorm(-abs(z_score), log.p = TRUE)

  data.frame(
    model_1     = label_model1,
    model_2     = label_model2,
    n           = nrow(oof1),
    n_events    = sum(oof1$status == 1L),
    c_index_1   = unname(cc$est.c["Cxz"]),
    c_index_2   = unname(cc$est.c["Cxy"]),
    delta_c     = unname(cc$est.diff_c),
    se_delta    = se_delta,
    ci_lower    = unname(cc$est.diff_c - zcrit * se_delta),
    ci_upper    = unname(cc$est.diff_c + zcrit * se_delta),
    z_score     = z_score,
    p_value     = unname(cc$pval),
    log_p_value = log_p_value,
    stringsAsFactors = FALSE
  )
}


# ----------------------------------------------------------------------------
# Strip large result objects down to the fields needed for compareC
# ----------------------------------------------------------------------------

strip_oof_result <- function(
  res,
  id_col = "StudyID",
  lp_col = "lp",
  time_col = "time",
  status_col = "status"
) {
  if (is.null(res$oof_predictions)) {
    stop("Result object is missing $oof_predictions.", call. = FALSE)
  }

  oof <- as.data.frame(res$oof_predictions, stringsAsFactors = FALSE)
  required_cols <- c(id_col, lp_col, time_col, status_col)

  missing_cols <- setdiff(required_cols, names(oof))
  if (length(missing_cols) > 0) {
    stop(
      "$oof_predictions is missing columns: ",
      paste(missing_cols, collapse = ", "),
      call. = FALSE
    )
  }

  out <- data.frame(
    StudyID = as.character(oof[[id_col]]),
    lp      = oof[[lp_col]],
    time    = oof[[time_col]],
    status  = oof[[status_col]],
    stringsAsFactors = FALSE
  )

  list(oof_predictions = out)
}


# ----------------------------------------------------------------------------
# Internal helpers for caching
# ----------------------------------------------------------------------------

.safe_cache_label <- function(x, max_chars = 80L) {
  x <- gsub("[^A-Za-z0-9._-]+", "_", as.character(x))
  substr(x, 1L, max_chars)
}


.make_delta_cache_file <- function(
  label_model1,
  label_model2,
  model_hash1,
  model_hash2,
  cache_dir,
  lp_higher_risk,
  conf_level
) {
  if (!requireNamespace("digest", quietly = TRUE)) {
    stop(
      "Package 'digest' is required when cache_dir is used. ",
      "Install it with install.packages('digest').",
      call. = FALSE
    )
  }

  cache_key <- digest::digest(
    list(
      label_model1 = label_model1,
      label_model2 = label_model2,
      model_hash1 = model_hash1,
      model_hash2 = model_hash2,
      lp_higher_risk = lp_higher_risk,
      conf_level = conf_level
    ),
    algo = "xxhash64"
  )

  file_name <- paste0(
    .safe_cache_label(label_model1),
    "__vs__",
    .safe_cache_label(label_model2),
    "__",
    cache_key,
    ".rds"
  )

  file.path(cache_dir, file_name)
}


.atomic_save_rds <- function(object, file) {
  tmp <- paste0(
    file,
    ".tmp_",
    Sys.getpid(),
    "_",
    as.integer(runif(1L, 1L, .Machine$integer.max))
  )

  saveRDS(object, tmp)

  ok <- file.rename(tmp, file)

  if (!ok) {
    unlink(tmp)
    saveRDS(object, file)
  }

  invisible(file)
}


.run_delta_cindex_cached <- function(
  res1,
  res2,
  label_model1,
  label_model2,
  cache_dir = NULL,
  force = FALSE,
  lp_higher_risk = TRUE,
  conf_level = 0.95,
  model_hash1 = NULL,
  model_hash2 = NULL
) {
  if (!is.null(cache_dir)) {
    if (!requireNamespace("digest", quietly = TRUE)) {
      stop(
        "Package 'digest' is required when cache_dir is used. ",
        "Install it with install.packages('digest').",
        call. = FALSE
      )
    }

    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)

    if (is.null(model_hash1)) {
      model_hash1 <- digest::digest(res1$oof_predictions, algo = "xxhash64")
    }

    if (is.null(model_hash2)) {
      model_hash2 <- digest::digest(res2$oof_predictions, algo = "xxhash64")
    }

    cache_file <- .make_delta_cache_file(
      label_model1 = label_model1,
      label_model2 = label_model2,
      model_hash1 = model_hash1,
      model_hash2 = model_hash2,
      cache_dir = cache_dir,
      lp_higher_risk = lp_higher_risk,
      conf_level = conf_level
    )

    if (!force && file.exists(cache_file)) {
      return(readRDS(cache_file))
    }
  }

  ans <- run_delta_cindex(
    res1 = res1,
    res2 = res2,
    label_model1 = label_model1,
    label_model2 = label_model2,
    lp_higher_risk = lp_higher_risk,
    conf_level = conf_level
  )

  if (!is.null(cache_dir)) {
    .atomic_save_rds(ans, cache_file)
  }

  ans
}


# ----------------------------------------------------------------------------
# Batch runner
# ----------------------------------------------------------------------------

run_delta_cindex_batch <- function(
  results,
  comparisons,
  strip_first = TRUE,
  complete_case = TRUE,
  cache_dir = NULL,
  force = FALSE,
  parallel = FALSE,
  workers = NULL,
  parallel_transport = c("rds", "globals"),
  worker_tmp_dir = NULL,
  cleanup_worker_tmp = TRUE,
  close_workers_on_exit = TRUE,
  lp_higher_risk = TRUE,
  conf_level = 0.95
) {
  parallel_transport <- match.arg(parallel_transport)

  if (!is.list(results) || is.null(names(results))) {
    stop("'results' must be a named list of result objects.", call. = FALSE)
  }

  if (anyNA(names(results)) || any(names(results) == "")) {
    stop("'results' must have non-empty names for all models.", call. = FALSE)
  }

  if (anyDuplicated(names(results))) {
    stop("'results' names must be unique.", call. = FALSE)
  }

  if (!all(c("model_1", "model_2") %in% names(comparisons))) {
    stop(
      "'comparisons' must contain columns named model_1 and model_2.",
      call. = FALSE
    )
  }

  comparisons <- as.data.frame(comparisons, stringsAsFactors = FALSE)

  if (nrow(comparisons) == 0L) {
    return(data.frame())
  }

  comparisons$model_1 <- as.character(comparisons$model_1)
  comparisons$model_2 <- as.character(comparisons$model_2)

  if (!"label_model1" %in% names(comparisons)) {
    comparisons$label_model1 <- comparisons$model_1
  }

  if (!"label_model2" %in% names(comparisons)) {
    comparisons$label_model2 <- comparisons$model_2
  }

  comparisons$label_model1 <- as.character(comparisons$label_model1)
  comparisons$label_model2 <- as.character(comparisons$label_model2)

  models_used <- unique(c(comparisons$model_1, comparisons$model_2))

  missing_models <- setdiff(models_used, names(results))

  if (length(missing_models) > 0) {
    stop(
      "These comparison models are missing from 'results': ",
      paste(missing_models, collapse = ", "),
      call. = FALSE
    )
  }

  if (strip_first) {
    results <- lapply(results, strip_oof_result)
    gc()
  }

  # Restrict every model used in comparisons to the shared complete-case
  # StudyID intersection so all delta-C tests use the same patient set.
  if (isTRUE(complete_case)) {
    common_info <- get_common_oof_ids(results, models_used)

    if (any(common_info$dropped_by_model > 0L)) {
      drop_msg <- paste(
        sprintf(
          "%s: %d -> %d (dropped %d)",
          names(common_info$dropped_by_model),
          common_info$n_by_model,
          common_info$n,
          common_info$dropped_by_model
        ),
        collapse = "; "
      )
      message(
        "complete_case=TRUE: restricting all comparisons to n=",
        common_info$n, " shared StudyIDs. ",
        drop_msg
      )
    } else {
      message(
        "complete_case=TRUE: all models already share n=",
        common_info$n, " StudyIDs."
      )
    }

    results[models_used] <- lapply(models_used, function(nm) {
      restrict_oof_to_ids(results[[nm]], common_info$ids, object_name = nm)
    })
  }

  if (!is.null(cache_dir)) {
    if (!requireNamespace("digest", quietly = TRUE)) {
      stop(
        "Package 'digest' is required when cache_dir is used. ",
        "Install it with install.packages('digest').",
        call. = FALSE
      )
    }

    dir.create(cache_dir, recursive = TRUE, showWarnings = FALSE)
    cache_dir <- normalizePath(cache_dir, mustWork = TRUE)

    # Important: hash each model once, not once per comparison.
    model_hashes <- vapply(
      results,
      function(x) digest::digest(x$oof_predictions, algo = "xxhash64"),
      character(1L)
    )
  } else {
    model_hashes <- setNames(rep(NA_character_, length(results)), names(results))
  }

  run_one_direct <- function(i) {
    m1 <- comparisons$model_1[i]
    m2 <- comparisons$model_2[i]

    .run_delta_cindex_cached(
      res1 = results[[m1]],
      res2 = results[[m2]],
      label_model1 = comparisons$label_model1[i],
      label_model2 = comparisons$label_model2[i],
      cache_dir = cache_dir,
      force = force,
      lp_higher_risk = lp_higher_risk,
      conf_level = conf_level,
      model_hash1 = model_hashes[[m1]],
      model_hash2 = model_hashes[[m2]]
    )
  }

  if (!parallel) {
    out_list <- lapply(seq_len(nrow(comparisons)), run_one_direct)
  } else {
    if (!requireNamespace("future", quietly = TRUE)) {
      stop(
        "Package 'future' is required when parallel = TRUE. ",
        "Install it with install.packages('future').",
        call. = FALSE
      )
    }

    if (!requireNamespace("future.apply", quietly = TRUE)) {
      stop(
        "Package 'future.apply' is required when parallel = TRUE. ",
        "Install it with install.packages('future.apply').",
        call. = FALSE
      )
    }

    if (is.null(workers)) {
      workers <- min(
        2L,
        max(1L, parallel::detectCores() - 1L),
        nrow(comparisons)
      )
    }

    workers <- as.integer(workers)

    if (is.na(workers) || workers < 1L) {
      stop("'workers' must be a positive integer.", call. = FALSE)
    }

    workers <- min(workers, nrow(comparisons))

    old_plan <- future::plan()

    on.exit({
      if (isTRUE(close_workers_on_exit)) {
        future::plan(future::sequential)
      } else {
        future::plan(old_plan)
      }
      gc()
    }, add = TRUE)

    future::plan(future::multisession, workers = workers)

    if (parallel_transport == "rds") {
      # Lower-memory path:
      # Save each stripped model once to disk, then workers read only the two
      # models needed for each comparison. This avoids exporting the full
      # results list to every worker process.
      if (is.null(worker_tmp_dir)) {
        worker_tmp_dir <- tempfile("delta_cindex_worker_")
      }

      dir.create(worker_tmp_dir, recursive = TRUE, showWarnings = FALSE)
      worker_tmp_dir <- normalizePath(worker_tmp_dir, mustWork = TRUE)

      if (isTRUE(cleanup_worker_tmp)) {
        on.exit(unlink(worker_tmp_dir, recursive = TRUE, force = TRUE), add = TRUE)
      }

      model_files <- setNames(
        file.path(
          worker_tmp_dir,
          paste0(
            seq_along(results),
            "__",
            .safe_cache_label(names(results)),
            ".rds"
          )
        ),
        names(results)
      )

      for (nm in names(results)) {
        saveRDS(results[[nm]], model_files[[nm]])
      }

      jobs <- lapply(seq_len(nrow(comparisons)), function(i) {
        m1 <- comparisons$model_1[i]
        m2 <- comparisons$model_2[i]

        list(
          file1 = model_files[[m1]],
          file2 = model_files[[m2]],
          label_model1 = comparisons$label_model1[i],
          label_model2 = comparisons$label_model2[i],
          cache_dir = cache_dir,
          force = force,
          lp_higher_risk = lp_higher_risk,
          conf_level = conf_level,
          model_hash1 = model_hashes[[m1]],
          model_hash2 = model_hashes[[m2]]
        )
      })

      rm(model_files)
      gc()

      out_list <- future.apply::future_lapply(
        jobs,
        function(job) {
          res1 <- readRDS(job$file1)
          res2 <- readRDS(job$file2)

          ans <- .run_delta_cindex_cached(
            res1 = res1,
            res2 = res2,
            label_model1 = job$label_model1,
            label_model2 = job$label_model2,
            cache_dir = job$cache_dir,
            force = job$force,
            lp_higher_risk = job$lp_higher_risk,
            conf_level = job$conf_level,
            model_hash1 = job$model_hash1,
            model_hash2 = job$model_hash2
          )

          rm(res1, res2)
          gc()

          ans
        },
        future.seed = 123
      )
    } else {
      # Faster but more memory-hungry path:
      # Sends model objects through future globals. Use only when results are small.
      jobs <- lapply(seq_len(nrow(comparisons)), function(i) {
        m1 <- comparisons$model_1[i]
        m2 <- comparisons$model_2[i]

        list(
          res1 = results[[m1]],
          res2 = results[[m2]],
          label_model1 = comparisons$label_model1[i],
          label_model2 = comparisons$label_model2[i],
          cache_dir = cache_dir,
          force = force,
          lp_higher_risk = lp_higher_risk,
          conf_level = conf_level,
          model_hash1 = model_hashes[[m1]],
          model_hash2 = model_hashes[[m2]]
        )
      })

      out_list <- future.apply::future_lapply(
        jobs,
        function(job) {
          .run_delta_cindex_cached(
            res1 = job$res1,
            res2 = job$res2,
            label_model1 = job$label_model1,
            label_model2 = job$label_model2,
            cache_dir = job$cache_dir,
            force = job$force,
            lp_higher_risk = job$lp_higher_risk,
            conf_level = job$conf_level,
            model_hash1 = job$model_hash1,
            model_hash2 = job$model_hash2
          )
        },
        future.seed = FALSE
      )
    }
  }

  out <- do.call(rbind, out_list)
  rownames(out) <- NULL

  out
}