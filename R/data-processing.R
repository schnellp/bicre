#' @title Collapse redundant intervals in dataset
#'
#' @description Combines adjacent or overlapping intervals with identical
#'              data into one interval.
#'
#' @param data A data frame.
#' @param t_start_var A string indicating the name of the column in \code{data}
#'                    storing the start of each interval.
#' @param t_end_var A string indicating the name of the column in \code{data}
#'                  storing the start of each interval.
#'
#' @return A data frame in which adjacent or overlapping intervals with
#'         identical data combined into one interval.
#'
#' @examples
#' data_overlap_ints <- data.frame(
#' id = c(1, 1, 1, 1,
#'        2, 2, 2, 2),
#' t_start = c(0, 4, 7, 21,
#'             7, 14, 21, 28),
#' t_end = c(7, 11, 14, 28,
#'           14, 21, 28, 35),
#' med = c("A", "A", "B", "B",
#'         "A", "A", "A", "A")
#' )
#'
#' data_overlap_ints %>% collapse_interval_data("t_start", "t_end")
#'
#' @export
collapse_interval_data <- function(data, t_start_var, t_end_var) {
  data %>%
    arrange(!!sym(t_start_var)) %>%
    group_by(across(-all_of(c(t_start_var, t_end_var)))) %>%
    mutate(.index = c(0,
                      cumsum(as.numeric(lead(!!sym(t_start_var))) >
                               cummax(as.numeric(!!sym(t_end_var))))[-n()])) %>%
    group_by(across(-all_of(c(t_start_var, t_end_var)))) %>%
    summarize({{ t_start_var }} := min(!!sym(t_start_var)),
              {{ t_end_var }} := max(!!sym(t_end_var), na.rm = TRUE),
              .groups = "drop") %>%
    select(-.index)
}

#' @title Fill in gaps between intervals
#'
#' @description Creates additional data rows representing
#'              gaps between existing intervals.
#'
#' @param data A data frame.
#' @param id Group identifier in \code{data}.
#' @param t_start_var The column in \code{data}
#'                    storing the start of each interval.
#' @param t_end_var The column in \code{data}
#'                  storing the end of each interval.
#' @param fill A named list specifying the value of each remaining variable
#'             that should be used for newly added rows.
#'
#' @return A data frame containing the original intervals plus
#'         the intervals implied by the gaps between them.
#'
#' @export
complete_interval_data <- function(data, id, t_start_var, t_end_var,
                                   fill = NA) {

  var_names <- data %>%
    select(-c({{ id }}, {{ t_start_var }}, {{ t_end_var }})) %>%
    colnames()

  data_list <- data %>% split(f = data %>% select({{ id }}))
  data_filled_list <- list()

  ids <- data %>% pull({{ id }}) %>% unique() %>% sort()

  for (i in ids) {

    nodes <- data_list[[i]] %>%
      select({{ t_start_var }}, {{ t_end_var }}) %>%
      unlist() %>%
      unique() %>%
      sort()

    data_filled_list[[i]] <- data.frame(
      t_start = head(nodes, -1),
      t_end = tail(nodes, -1)
    ) %>%
      mutate({{ id }} := i, .before = t_start)

    for (var_name in var_names) {
      if (!is.list(fill)) {
        fill_val <- fill
      } else if (is.null(fill[[var_name]])) {
        fill_val <- NA
      } else {
        fill_val <- fill[[var_name]]
      }

      data_filled_list[[i]] <-
        data_filled_list[[i]] %>%
        mutate(!!sym(var_name) := fill_val)


      for (j in 1 : (length(nodes) - 1)) {
        node_start <- nodes[j]
        node_end <- nodes[j + 1]

        rows <- data_list[[i]] %>%
          filter({{ t_start_var }} <= node_start, {{ t_end_var }} >= node_end)

        if (nrow(rows) == 1) {
          data_filled_list[[i]][j, var_name] <- rows[1, var_name]
        } else if (nrow(rows) > 1) {
          print(rows)
          print(nodes)
          print(node_start)
          print(node_end)
          stop("Multiple rows covering interval")
        }
      }
    }

    rownames(data_filled_list[[i]]) <- NULL
  }



  do.call(rbind.data.frame, data_filled_list)
}

combine_variables <- function() {

}


co_events <- function(data_covariates, data_events,
                      id, t_start, t_end, e_min, e_max, e_type = NULL) {

}
