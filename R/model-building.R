#' @title Build model frames from \code{co_events} object.
#'
#' @description For each element in a \code{co_events} object,
#'              builds a model frame from \code{covariates}, and constructs
#'              data frames for intervals with standardized column names
#'              to be used internally in model fitting functions.
#'
#' @param co_events A \code{co_events} object.
#' @param formula An object of class \code{formula} symbolically describing
#'                the model to be fit.
#' @param contrasts An optional list to be passed to \code{model.matrix}.
#' @param ... Additional arguments to be passed to \code{model.frame}.
#'
#' @return An object of class \code{co_events_frame}.
#'         In addition to the components of the given \code{co_events} object,
#'         each list element has the following sub-elements:
#'         \describe{
#'           \item{\code{X}}{a model matrix built from \code{covariates}
#'           and the \code{formula} argument;}
#'           \item{\code{cov_t}}{the \code{t_start} and \code{t_end}
#'           columns of \code{covariates} but with standardized variable
#'           names \code{"t_start"} and \code{"t_end"};}
#'           \item{\code{ev}}{same as the \code{events} component
#'           but with standardized variable names \code{"t_start"},
#'           \code{"t_end"}, \code{"e_min"}, and \code{"e_max"}.}
#'         }
#' @export

co_events_frame <- function(co_events, formula, contrasts = NULL, ...) {

  # get the function call with full argument names,
  # ignoring arguments passed to ...
  mf <- match.call(expand.dots = FALSE)

  # get indices of arguments matching the following names.
  # if a listed name isn't found in the arguments, return a 0 for that index.
  m <- match(c("formula", "co_events"),
             table = names(mf), nomatch = 0L)

  # trim the function call to only the function name
  # and the arguments matched in m
  mf <- mf[c(1L, m)]

  # replace the original function name with "model.frame"
  mf[[1L]] <- quote(stats::model.frame)

  # rename co_events entry "data"
  names(mf)[names(mf) == "co_events"] <- "data"
  ce_frame <- list()
  for (i in 1 : length(co_events)) {
    ce_frame[[i]] <- list()
    ce_frame[[i]]$covariates <- co_events[[i]]$covariates %>% select(-any_of(attr(co_events, "special_cols")))
    ce_frame[[i]]$events <- co_events[[i]]$events %>% select(-any_of(attr(co_events, "special_cols")["id"]))
    mf[["data"]] <- co_events[[i]]$covariates
    mfi <- eval(mf, envir = parent.frame())
    mt <- attr(mfi, "terms")
    ce_frame[[i]]$X <- model.matrix(mt, mfi, contrasts)

    cov_t <- data.frame(
      t_start = co_events[[i]]$covariates[, attr(co_events, "special_cols")["t_start"]],
      t_end = co_events[[i]]$covariates[, attr(co_events, "special_cols")["t_end"]]
    )
    ce_frame[[i]]$cov_t <- cov_t

    ev <- co_events[[i]]$events[, c(attr(co_events, "special_cols")["t_start"],
                                    attr(co_events, "special_cols")["t_end"],
                                    attr(co_events, "special_cols")["e_min"],
                                    attr(co_events, "special_cols")["e_max"])]
    colnames(ev) <- c("t_start", "t_end", "e_min", "e_max")
    ce_frame[[i]]$ev <- ev

  }

  class(ce_frame) <- c(class(co_events), "co_events_frame")

  ce_frame
}



simplify_ev <- function(ev){

  ### 1. make a union of all F intervals
  ev_zero <- ev %>% filter(e_max == 0)
  ev_impute <- ev %>% filter(e_max > 0)

  if(nrow(ev_zero) > 1){
    ev_zero <- ev_zero %>% collapse_interval_data(t_start, t_end) %>% as.data.frame()
  }


  ### 2. shorten >0-event intervals with one-sided overlap with
  ### zero-event intervals.
  if (nrow(ev_zero) > 0 && nrow(ev_impute) > 0) {
    for (i in 1 : nrow(ev_zero)) {
      e <- ev_zero %>% slice(i)
      ev_impute <- ev_impute %>%
        mutate(
          t_end = if_else( # t_end but not t_start in current interval
            t_start < e$t_start &
              t_end > e$t_start &
              t_end <= e$t_end,
            as.double(e$t_start),
            as.double(t_end)),
          t_start = if_else( # t_start but not t_end in current interval
            t_start >= e$t_start &
              t_start < e$t_end &
              t_end > e$t_end,
            as.double(e$t_end),
            as.double(t_start)
          )
        )
    }
  }

  ### 3. delete >0-event intervals covered by zero-event intervals.
  if (nrow(ev_zero) > 0 && nrow(ev_impute) > 0){
    ev_impute_delete_row <- c()
    for(row_num in 1:nrow(ev_impute)){
      e <- ev_impute %>% slice(row_num)
      if(nrow(ev_zero %>% filter(t_start <= e$t_start & t_end >= e$t_end)) > 0){
        if(e$e_min != 0){
          stop("contradiction between T interval and F interval covering it")
        }
        ev_impute_delete_row <- c(ev_impute_delete_row, row_num)
      }
    }
    if(length(ev_impute_delete_row) > 0){
      ev_impute <- ev_impute[-ev_impute_delete_row,]
    }

  }

  ev_impute <- ev_impute[!duplicated(ev_impute), ]
  ### 4. delete upper or lower bound only T intervals
  ### whose information has already represented by other intervals
  if (nrow(ev_impute) > 1){
    ev_upper_only_row <- with(ev_impute, e_max != Inf & e_min ==0) %>% which
    ev_lower_only_row <- with(ev_impute, e_max == Inf & e_min !=0) %>% which
    ev_both_row <- (1:nrow(ev_impute))[-c(ev_upper_only_row, ev_lower_only_row)]

    delete_row_upper <- c()
    if (length(ev_upper_only_row) > 0){
      for(row_num in ev_upper_only_row){
        e <- ev_impute %>% slice(row_num)
        if(sum(ev_impute %>% slice(c(ev_both_row, ev_upper_only_row)) %>%
               apply(1, function(x) x["t_start"] <= e$t_start &
                     x["t_end"] >= e$t_end &
                     x["e_max"] <= e$e_max) ) > 1){
          delete_row_upper <- c(delete_row_upper, row_num)
        }
      }
    }

    delete_row_lower <- c()
    if (length(ev_lower_only_row) > 0){
      for(row_num in ev_lower_only_row){
        e <- ev_impute %>% slice(row_num)
        if(sum(ev_impute %>% slice(c(ev_both_row, ev_lower_only_row)) %>%
               apply(1, function(x) x["t_start"] >= e$t_start &
                     x["t_end"] <= e$t_end &
                     x["e_min"] >= e$e_min) ) > 1){
          delete_row_lower <- c(delete_row_lower, row_num)
        }
      }

    }

    if (length(c(delete_row_upper, delete_row_lower)) > 0){
      ev_impute <- ev_impute[-c(delete_row_upper, delete_row_lower),]
    }
  }

  ### 5. A T interval's upper bound can be changed by the smaller upper bound of
  ### any T interval strictly covers the former T interval.
  ### A T interval's lower bound can be changed by the bigger lower bound of
  ### any T interval strictly covered by the former T interval.

  if (nrow(ev_impute) > 1){
    for(i in 1 : nrow(ev_impute)){
      e <- ev_impute %>% slice(i)
      ev_impute$e_max[i] <- ev_impute %>% filter(t_start <= e$t_start &
                             t_end >= e$t_end &
                             e_max <= e$e_max) %>% pull(e_max) %>% min
      ev_impute$e_min[i] <- ev_impute %>% filter(t_start >= e$t_start &
                             t_end <= e$t_end &
                             e_min >= e$e_min) %>% pull(e_min) %>% max
    }
  }

  # if (nrow(ev_zero) > 0){
  #   ev_zero_outer_row <- c()
  #   for (i in 1:nrow(ev_zero)){
  #     e <- ev_zero %>% slice(i)
  #     with(ev_impute, (t_start < e$t_start & t_end > e$t_end))
  #   }
  # }

  ev <- rbind(ev_zero, ev_impute)  %>% arrange(t_start, t_end, e_min, e_max)
  ev[!duplicated(ev), ]
}



partition_ev <- function(ev) {

  ### partition into self-overlapping units using column unit

  ev <- ev %>%
    mutate(
      unit = as.numeric(NA)
    )

  u <- 1
  ev$unit[1] <- 1
  t_end <- ev$t_end[1]

  if (nrow(ev) > 1) {
    for (i in 2 : nrow(ev)) {
      if (ev$t_start[i] >= t_end) {
        u <- u + 1
      }
      ev$unit[i] <- u
      t_end <- max(t_end, ev$t_end[i])
    }
  }

  ev
}

compound_impute_unit_prepare <- function(compound_impute_unit){
  # for compound imputation unit first distinguish between
  # Unit disjoint sub intervals vs Unit disjoint check intervals
  impute_unit <- compound_impute_unit
  impute_unit$sim <- FALSE
  t_end <- impute_unit$t_end[1]
  impute_unit$sim[1] <- TRUE
  for (j in 2 : nrow(impute_unit)) {
    if (impute_unit$t_start[j] >= t_end) {
      t_end <- impute_unit$t_end[j]
      impute_unit$sim[j] <- TRUE
    }
  }

  sim_split <- impute_unit %>%
    split(f = impute_unit %>% select(sim))


  # Prepare for 1st round simulation
  round_disjoint <- vector("list", length = nrow(sim_split$`TRUE`))
  inner_check_row_all <- c()
  # check_interv_index <- 1:nrow(sim_split$`FALSE`)
  for(j in 1 : nrow(sim_split$`TRUE`)){
    # For each Unit disjoint sub F-interval
    # Nothing needs preparing since no simulation there
    # For each Unit disjoint sub T-interval, find
    # a.	Inner Unit check sub interval covered by it
    # b.  Outer Unit check sub interval overlapping with it

    e_sim <- sim_split$`TRUE`[j,]
    round_disjoint[[j]]$sim <- e_sim

    # initialized inner_check_row for later inner_check_row_all updates
    inner_check_row <- c()


    if(e_sim$e_max != 0){
      inner_check_row <- with(sim_split$`FALSE`, t_start >= e_sim$t_start &  t_end <= e_sim$t_end)
      check_row_all <- with(sim_split$`FALSE`, t_end >  e_sim$t_start & t_start < e_sim$t_end)


      outer_check_row <- check_row_all & (!inner_check_row)

      round_disjoint[[j]]$inner_check <- sim_split$`FALSE`[inner_check_row,]
      # for outer check, only upper bound is used hence delete those without upper bounds
      # since check interval still covers other area, so lower bound can not contribute
      round_disjoint[[j]]$outer_check <- sim_split$`FALSE`[outer_check_row,] %>% filter(e_max < e_sim$e_max)
    }else{
      round_disjoint[[j]]$inner_check <- round_disjoint[[j]]$outer_check <- data.frame(matrix(nrow = 0, ncol = 6))
    }


    if(any(inner_check_row)){
      inner_check_row_all <- c(inner_check_row_all, which(inner_check_row))
    }
  }


  # Prepare for 2st round simulation
  left_over <- data.frame(t_start = sim_split$`TRUE`$t_end[-nrow(sim_split$`TRUE`)],
            t_end = sim_split$`TRUE`$t_start[-1])
  unit_end <- max(impute_unit$t_end)
  left_over <- rbind(left_over, data.frame(t_start = last(sim_split$`TRUE`$t_end),
                                  t_end = unit_end) )
  left_over <- left_over[with(left_over, t_start != t_end) ,]


  # Pick out outer unit check intervals, these have to be checked again in the final check
  if(length(inner_check_row_all) != 0){
    outer_check_all <- sim_split$`FALSE`[-inner_check_row_all,]
  }else{
    outer_check_all <- sim_split$`FALSE`
  }


  if(nrow(left_over) != 0){
    round_leftover <- vector("list", length = nrow(left_over))
    for(j in 1:nrow(left_over)){
      # For each left over interval, find
      #a.	Outer Unit check sub interval covering this left over interval
      #b.	Outer Unit check sub interval overlapping with it
      #These two must be T- intervals since all inner F intervals must be contained
      #in a T-interval or disjoint with a T-interval for each T-interval after simplification.

      e_leftover <- left_over[j,]
      round_leftover[[j]]$leftover <- e_leftover

      covered_row <- with(outer_check_all, t_start <= e_leftover$t_start &  t_end >= e_leftover$t_end)
      overlap_row <- with(outer_check_all, t_end >  e_leftover$t_start & t_start < e_leftover$t_end)
      # for check sub interval covering the left over, only upper bound is useful
      # since check interval still covers other area
      round_leftover[[j]]$covered <- outer_check_all[covered_row,] %>% filter(e_max != Inf)
      # for check sub interval strictly overlapping with the left over, only upper bound is useful
      # since check interval still covers other area
      round_leftover[[j]]$overlap_strict <- outer_check_all[overlap_row & (!covered_row), ] %>% filter(e_max != Inf)

    }
  }else{
    round_leftover <- NULL
  }

  return(list(round_disjoint = round_disjoint,
              round_leftover = round_leftover,
              outer_check_all = outer_check_all,
              unit_range = data.frame(t_start = impute_unit$t_start[1],
                                      t_end = unit_end)))


}

seq_samp_prepare_ev <- function(ev){

  # partition by imputation unit
  ev <- ev %>%
    split(f = ev %>% select(unit))


  impute_unit_ez <- list()
  ez_index <- 1
  impute_unit_compound <- list()
  compound_index <- 1
  for (i in 1:length(ev)){
    if(nrow(ev[[i]]) == 1){
      # find out easy imputation unit that only contains one T interval
      # ev[[i]]$sim <- T
      impute_unit_ez[[ez_index]] <- ev[[i]]
      ez_index <- ez_index + 1
    }else{
      impute_unit_compound[[compound_index]] <- compound_impute_unit_prepare(ev[[i]])
      compound_index <- compound_index + 1
    }

  }


  # impute_unit_ez <- impute_unit_ez %>% do.call(rbind,.)


  # covered_ints <- rbind(impute_unit_ez[,c(1,2)],
  #      impute_unit_compound %>% lapply("[[","unit_range") %>% do.call(rbind,.)) %>%
  #   collapse_interval_data(t_start, t_end)
  ### record covered intervals

  impute_unit_single <- list(impute_unit_ez = impute_unit_ez,
                             impute_unit_compound = impute_unit_compound)

  covered_ints <- list()

  for (i in 1 : length(ev)) {
    covered_ints[[i]] <- data.frame(
      t_start = min(ev[[i]]$t_start),
      t_end = max(ev[[i]]$t_end)
    )
  }

  attributes(impute_unit_single)$covered_ints <- do.call("rbind", covered_ints) %>%
    collapse_interval_data(t_start, t_end) %>% as.data.frame()


  return(impute_unit_single)

}


#' build imputation units for a single id without tiny intervals
#' @export
build_imputation_units_single <- function(ev) {
  ev %>% simplify_ev %>% partition_ev %>% seq_samp_prepare_ev
}



#' @title Build imputation units
#'
#' @description An imputation unit is a set of one or more pairwise overlapping
#'              intervals on which the event history may be imputed
#'              independently of all intervals outside of their union.
#'
#' @param ev The \code{ev} member of a \code{co_events_frame} object
#' @export


build_imputation_units <- function(co_events_frame, tiny_diff = NULL, prepare_style = "normal", verbose = FALSE, length_rate_perc_prod = 0.1) {
  if(is.null(tiny_diff)){
    if(verbose){
      print("Preprocessing WITHOUT tiny intervals")
    }
    for (i in 1 : length(co_events_frame)) {
      if(verbose){
        if(i %% 100 == 0)
        print(paste("Now data preprocessing for individual ", i, " out of ",  length(co_events_frame), " total",sep = ""))
      }
      co_events_frame[[i]]$ev <- co_events_frame[[i]]$ev %>%
        build_imputation_units_single()
    }
  }else{
    if(verbose){
      print("Preprocessing WITH tiny intervals")
    }
    for (i in 1 : length(co_events_frame)) {
      if(verbose){
        if(i %% 100 == 0)
        print(paste("Now data preprocessing for individual ", i, " out of ",  length(co_events_frame), " total",sep = ""))
      }
      co_events_frame[[i]]$ev <- co_events_frame[[i]]$ev %>%
        build_imputation_units_single_w_tiny_interval(tiny_diff = tiny_diff, prepare_style = prepare_style, length_rate_perc_prod = length_rate_perc_prod)
    }
  }


  co_events_frame
}
