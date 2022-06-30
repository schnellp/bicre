#' @title Collapse redundant intervals in dataset
#'
#' @description Combines adjacent or overlapping intervals with identical
#'              data into one interval.
#'
#' @param data A data frame.
#' @param t_start_var Column in \code{data} storing the start of each interval.
#' @param t_end_var Column in \code{data} storing the start of each interval.
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
#' data_overlap_ints %>% collapse_interval_data(t_start, t_end)
#'
#' @export
collapse_interval_data <- function(data, t_start_var, t_end_var) {
  data_collapsed <- data %>%
    arrange({{ t_start_var }}) %>%
    group_by(across(-c({{ t_start_var }}, {{ t_end_var }}))) %>%
    mutate(.index = c(0,
                      cumsum(as.numeric(lead({{ t_start_var }})) >
                               cummax(as.numeric({{ t_end_var }})))[-n()])) %>%
    group_by(across(-c({{ t_start_var }}, {{ t_end_var }}))) %>%
    summarize({{ t_start_var }} := min({{ t_start_var }}),
              {{ t_end_var }} := max({{ t_end_var }}, na.rm = TRUE),
              .groups = "drop") %>%
    select(-.index)

  data_collapsed[, colnames(data_collapsed) %in% colnames(data)] <-
    data_collapsed[, colnames(data)]
}

#' @title Single-id helper for \code{complete_interval_data}
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

  data_filled[, colnames(data)]
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

#' @title Single-id helper for \code{merge_interval_data}
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

  merged_data <- data %>%
    complete_interval_data_single(t_start_var = {{ t_start_var }},
                                  t_end_var = {{ t_end_var }},
                                  fill = fill,
                                  new_nodes = nodes) %>%
    mutate({{ new_var }} := temp_data %>% pull({{ new_var }}))

  merged_data
}

#' @title Merge new time-dependent variable into interval dataset
#'
#' @description Merges an additional column from a new dataset
#'              into an existing interval dataset,
#'              subdividing intervals as necessary.
#'
#' @param data A data frame.
#' @param new_data A data frame containing the new variable to be merged.
#' @param t_start_var The column in \code{data} and \code{new_data}
#'                    storing the start of each interval.
#' @param t_end_var The column in \code{data} and \code{new_data}
#'                  storing the end of each interval.
#' @param new_var The column in \code{new_data} to be merged into \code{data}.
#' @param fill A named list specifying the value of each remaining variable
#'             that should be used for newly added rows.
#'
#' @return A data frame containing the original columns of \code{data} plus
#'         the new column from \code{new_data} with subdivided intervals
#'         to accommodate all variables.
#'
#' @export
merge_interval_data <- function(data, new_data,
                                id,
                                t_start_var, t_end_var, new_var,
                                fill = NA) {

  data_list <- data %>% split(f = data %>% select({{ id }}))
  new_data_list <- new_data %>% split(f = new_data %>% select({{ id }}))
  data_merged_list <- list()

  ids <- data %>% pull({{ id }}) %>% unique() %>% sort() %>% as.character()
  new_ids <- new_data %>% pull({{ id }}) %>% unique() %>% sort() %>% as.character()

  if (!all(new_ids %in% ids)) {
    warning(paste("Ids in new_data missing from data and will be ignored:",
                  paste0(setdiff(new_ids, ids), collapse = ", ")))
  }

  if (!all(ids %in% new_ids)) {
    warning(paste("Ids in data missing from new_data and will be filled according to `fill` argument:",
                  paste0(setdiff(ids, new_ids), collapse = ", ")))
  }

  for (i in ids) {
    id_fill_list <- list()
    id_fill_list[[deparse(substitute(id))]] <- i

    if (is.null(new_data_list[[i]])) {
      var_name <- as.character(substitute(new_var))

      if (!is.list(fill)) {
        fill_val <- fill
      } else if (is.null(fill[[var_name]])) {
        fill_val <- NA
      } else {
        fill_val <- fill[[var_name]]
      }

      new_data_list[[i]] <- data_list[[i]] %>%
        select({{ id }}, {{ t_start_var }}, {{ t_end_var }}) %>%
        mutate({{ new_var }} := fill[[var_name]])
    }

    data_merged_list[[i]] <-
      data_list[[i]] %>%
      merge_interval_data_single(
        new_data = new_data_list[[i]],
        t_start_var = {{ t_start_var }},
        t_end_var = {{ t_end_var }},
        new_var = {{ new_var }},
        fill = c(fill, id_fill_list)
      )
  }

  data_merged <- do.call(rbind.data.frame, data_merged_list)
  rownames(data_merged) <- NULL

  data_merged[, c(colnames(data), as.character(substitute(new_var)))]
}

#' @title Create \code{co_events} object.
#'
#' @description Combine time-varying covariate data and censored event data
#'              into a single computation-friendly object.
#'
#' @param data_covariates A data frame containing time-varying covariates
#'                        over disjoint intervals within units,
#'                        one row per interval, including columns for
#'                        unit id, interval start, and interval end.
#' @param data_events A data frame containing censored event counts
#'                    in intervals possibly overlapping within units,
#'                    one row per interval, including columns for
#'                    unit id and minimum and maximum event counts.
#' @param id The column in \code{data_covariates} and \code{data_events}
#'           containing IDs linking units between the datasets.
#' @param t_start The column in \code{data_covariates} and \code{data_events}
#'                containing the (open) start time of each interval.
#' @param t_end The column in \code{data_covariates} and \code{data_events}
#'              containing the (closed) end time of each interval.
#' @param e_min The column in \code{data_events}
#'              containing the minimum number of events in each interval.
#' @param e_max The column in \code{data_events}
#'              containing the maximum number of events in each interval.
#' @param fill A named list specifying the value of each time-varying covariate
#'             that should be used for newly added rows.
#'
#' @return An object of class \code{co_events}.
#'
#' @export
co_events <- function(data_covariates, data_events,
                      id, t_start, t_end, e_min, e_max,
                      fill = NA) {

  ids_covariates <- data_covariates %>%
    pull({{ id }}) %>%
    unique() %>%
    sort()

  ids_events <- data_events %>%
    pull({{ id }}) %>%
    unique() %>%
    sort()

  if (!all(ids_events %in% ids_covariates)) {
    ids_missing <- setdiff(ids_events, ids_covariates)

    warning(paste("Dropping IDs in data_events not in data_covariates:",
                  paste0(ids_missing, collapse = ", ")))

  }

  data_covariates_list <- data_covariates %>%
    split(f = data_covariates %>% select({{ id }}))

  data_events_list <- data_events %>%
    split(f = data_events %>% select({{ id }}))

  co_events <- list()

  # ignores IDs in data_events missing from data_covariates
  for (i in ids_covariates) {

    # make sure covariate intervals cover event intervals
    covariate_nodes <- data_covariates_list[[i]] %>%
      select(c({{ t_start }}, {{ t_end }})) %>%
      unlist() %>%
      unique()

    event_nodes <- data_events_list[[i]] %>%
      select(c({{ t_start }}, {{ t_end }})) %>%
      unlist() %>%
      unique()

    data_covariates_list[[i]] <- data_covariates_list[[i]] %>%
      complete_interval_data_single({{ t_start }}, {{ t_end }},
                                    fill = fill,
                                    new_nodes = c(min(event_nodes, covariate_nodes),
                                                  max(event_nodes, covariate_nodes))) %>%
      mutate({{ id }} := i)

    # pair covariate and event datasets
    co_events[[i]] <- list(
      covariates = data_covariates_list[[i]],
      events = data_events_list[[i]]
    )
  }

  class(co_events) <- c(class(co_events), "co_events")

  co_events
}
