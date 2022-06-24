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
#' @param t_start_var The column in \code{data}
#'                    storing the start of each interval.
#' @param t_end_var The column in \code{data}
#'                  storing the end of each interval.
#' @param fill A named list specifying the value of each remaining variable
#'             that should be used for newly added rows.
#' @param new_nodes A vector of new interval endpoints that should
#'                  be added regardless of whether they appear in the data.
#'
#' @return A data frame containing the original intervals plus
#'         the intervals implied by the gaps between them.
#'
#' @export
#'
#'
complete_interval_data_single <- function(data, t_start_var, t_end_var,
                                          fill = NA, new_nodes = c()) {

  var_names <- data %>%
    select(-c({{ t_start_var }}, {{ t_end_var }})) %>%
    colnames()


  nodes <- data %>%
    select({{ t_start_var }}, {{ t_end_var }}) %>%
    unlist() %>%
    c(new_nodes) %>%
    unique() %>%
    sort()

  data_filled <- data.frame(
    t_start = head(nodes, -1),
    t_end = tail(nodes, -1)
  )

  for (var_name in var_names) {
    if (!is.list(fill)) {
      fill_val <- fill
    } else if (is.null(fill[[var_name]])) {
      fill_val <- NA
    } else {
      fill_val <- fill[[var_name]]
    }

    data_filled <-
      data_filled %>%
      mutate(!!sym(var_name) := fill_val)


    for (j in 1 : (length(nodes) - 1)) {
      node_start <- nodes[j]
      node_end <- nodes[j + 1]

      rows <- data %>%
        filter({{ t_start_var }} <= node_start, {{ t_end_var }} >= node_end)

      if (nrow(rows) == 1) {
        data_filled[j, var_name] <- rows[1, var_name]
      } else if (nrow(rows) > 1) {
        print(rows)
        print(nodes)
        print(node_start)
        print(node_end)
        stop("Multiple rows covering interval")
      }
    }
  }

  rownames(data_filled) <- NULL

  data_filled
}

complete_interval_data <- function(data, id, t_start_var, t_end_var,
                                   fill = NA, new_nodes = c()) {

  data_list <- data %>% split(f = data %>% select({{ id }}))
  data_filled_list <- list()

  ids <- data %>% pull({{ id }}) %>% unique() %>% sort() %>% as.character()

  for (i in ids) {
    id_fill_list <- list()
    id_fill_list[[deparse(substitute(id))]] <- i

    data_filled_list[[i]] <-
      data_list[[i]] %>%
      complete_interval_data_single(
        t_start_var = {{ t_start_var }},
        t_end_var = {{ t_end_var }},
        fill = c(fill, id_fill_list),
        new_nodes = new_nodes
      )
  }

  do.call(rbind.data.frame, data_filled_list)
}

merge_interval_data_single <- function(data, new_data,
                                       t_start_var, t_end_var, new_var,
                                       fill = NA) {

  nodes <- union(
    data %>%
      select(c({{ t_start_var }}, {{ t_end_var }})) %>%
      unlist(),
    new_data %>%
      select(c({{ t_start_var }}, {{ t_end_var }})) %>%
      unlist()
  ) %>%
    unique() %>%
    sort()

  temp_data <- new_data %>%
    select(c({{ t_start_var }}, {{ t_end_var }}, {{ new_var }})) %>%
    complete_interval_data_single(t_start_var = {{ t_start_var }},
                                  t_end_var = {{ t_end_var }},
                                  fill = fill,
                                  new_nodes = nodes)

  data <- data %>%
    complete_interval_data_single(t_start_var = {{ t_start_var }},
                                  t_end_var = {{ t_end_var }},
                                  fill = fill,
                                  new_nodes = nodes) %>%
    mutate({{ new_var }} := temp_data %>% pull({{ new_var }}))

  data
}

merge_interval_data <- function(data, new_data,
                                id,
                                t_start_var, t_end_var, new_var,
                                fill = NA) {

  data_list <- data %>% split(f = data %>% select({{ id }}))
  data_merged_list <- list()

  ids <- data %>% pull({{ id }}) %>% unique() %>% sort() %>% as.character()

  for (i in ids) {
    id_fill_list <- list()
    id_fill_list[[deparse(substitute(id))]] <- i

    data_merged_list[[i]] <-
      data_list[[i]] %>%
      merge_interval_data_single(
        new_data = new_data,
        t_start_var = {{ t_start_var }},
        t_end_var = {{ t_end_var }},
        new_var = {{ new_var }}
        fill = c(fill, id_fill_list)
      )
  }

  do.call(rbind.data.frame, data_merged_list)
}

combine_variables <- function() {

}


co_events <- function(data_covariates, data_events,
                      id, t_start, t_end, e_min, e_max, e_type = NULL) {

}
