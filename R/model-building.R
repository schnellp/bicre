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

  for (i in 1 : length(co_events)) {
    mf[["data"]] <- co_events[[i]]$covariates
    mfi <- eval(mf, envir = parent.frame())
    mt <- attr(mfi, "terms")
    co_events[[i]]$X <- model.matrix(mt, mfi, contrasts)

    cov_t <- data.frame(
      t_start = co_events[[i]]$covariates[, attr(co_events, "special_cols")["t_start"]],
      t_end = co_events[[i]]$covariates[, attr(co_events, "special_cols")["t_end"]]
    )
    co_events[[i]]$cov_t <- cov_t

    ev <- co_events[[i]]$events[, c(attr(co_events, "special_cols")["t_start"],
                                    attr(co_events, "special_cols")["t_end"],
                                    attr(co_events, "special_cols")["e_min"],
                                    attr(co_events, "special_cols")["e_max"])]
    colnames(ev) <- c("t_start", "t_end", "e_min", "e_max")
    co_events[[i]]$ev <- ev

  }

  class(co_events) <- c(class(co_events), "co_events_frame")

  co_events
}

#' @title Build imputation units
#'
#' @description An imputation unit is a set of one or more pairwise overlapping
#'              intervals on which the event history may be imputed
#'              independently of all intervals outside of their union.
#'
#' @param ev The \code{ev} member of a \code{co_events_frame} object
build_imputation_units_single <- function(ev) {

  ### shorten >0-event intervals with one-sided overlap with
  ### zero-event intervals.
  ev_zero <- ev %>% filter(e_max == 0)
  ev_impute <- ev %>% filter(e_max > 0)

  if (nrow(ev_zero) > 0 && nrow(ev_impute) > 0) {
    for (i in 1 : nrow(ev_zero)) {
      e <- ev_zero %>% slice(i)
      ev_impute <- ev_impute %>%
        mutate(
          t_end = if_else( # t_end but not t_start in current interval
            t_start < e$t_start &
              t_end > e$t_start &
              t_end <= e$t_end,
            e$t_start,
            t_end),
          t_start = if_else( # t_start but not t_end in current interval
            t_start >= e$t_start &
              t_start < e$t_end &
              t_end > e$t_end,
            e$t_end,
            t_start
          )
        )
    }
  }



  ev <- rbind(ev_zero, ev_impute)

  ### partition into self-overlapping units

  ev <- ev %>%
    arrange(t_start, t_end, e_min, e_max) %>%
    mutate(
      unit = as.numeric(NA),
      noncontig = FALSE,
      t_start_trunc = t_start,
      t_end_trunc = t_end
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

  ev <- ev %>%
    split(f = ev %>% select(unit))

  ### find maximal non-overlapping set within each unit

  for (i in 1 : length(ev)) {
    t_end <- ev[[i]]$t_end[1]
    ev[[i]]$noncontig[1] <- TRUE

    if (nrow(ev[[i]]) > 1) {
      for (j in 2 : nrow(ev[[i]])) {
        if (ev[[i]]$t_start[j] >= t_end) {
          t_end <- ev[[i]]$t_end[j]
          ev[[i]]$noncontig[j] <- TRUE
        } else {
          ev[[i]]$t_start_trunc[j] <- t_end
        }
      }

      t_start <- Inf
      for (j in nrow(ev[[i]]) : 2) { # noncontig[1] == TRUE
        if (ev[[i]]$noncontig[j] == FALSE && ev[[i]]$t_end[j] > t_start) {
          ev[[i]]$t_end_trunc[j] <- t_start
        } else if (ev[[i]]$noncontig[j] == TRUE) {
          t_start <- ev[[i]]$t_start[j]
        }
      }
    }

    ev[[i]] <- ev[[i]] %>%
      split(f = ev[[i]] %>% select(noncontig))
  }

  ev
}

build_imputation_units <- function(co_events_frame) {
  for (i in 1 : length(co_events_frame)) {
    co_events_frame[[i]]$ev <- co_events_frame[[i]]$ev %>%
      build_imputation_units_single()
  }

  co_events_frame
}
