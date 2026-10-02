prep_incident_time_to_event_df <- function(meta_df, status_col, time_col, time_filter = NULL) {
  stopifnot(status_col %in% colnames(meta_df))
  stopifnot(time_col %in% colnames(meta_df))

  # Remove 'Prevalent' from status
  df <- meta_df %>%
    dplyr::filter(.data[[status_col]] != 'Prevalent')
  
  # Apply optional time filter
  if (!is.null(time_filter) && nzchar(time_filter)) {
    # Use standard evaluation to interpret time_filter
    df <- df %>% dplyr::filter(rlang::parse_expr(time_filter) %>% rlang::eval_tidy(df))
    # ^ fallback if you want robust NSE: df <- df %>% dplyr::filter(!!rlang::parse_expr(time_filter))
  }
  
  # Add 'status' (incident=1, other=0); add 'time'
  df <- df %>%
    dplyr::mutate(
      status = ifelse(.data[[status_col]] == 'Incident', 1, 0),
      time = .data[[time_col]]
    )
  return(df)
}
