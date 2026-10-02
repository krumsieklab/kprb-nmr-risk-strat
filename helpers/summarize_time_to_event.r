summarize_time_to_event <- function(df, status_col = "status", time_col = "time") {
  df %>%
    group_by(.data[[status_col]]) %>%
    summarise(
      n = n(),
      mean_time = mean(.data[[time_col]], na.rm = TRUE),
      median_time = median(.data[[time_col]], na.rm = TRUE),
      sd_time = sd(.data[[time_col]], na.rm = TRUE),
      min_time = min(.data[[time_col]], na.rm = TRUE),
      max_time = max(.data[[time_col]], na.rm = TRUE)
    )
}