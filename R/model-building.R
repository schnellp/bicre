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

  for (i in 1 : length(ce)) {
    mf[["data"]] <- ce[[i]]$covariates
    mfi <- eval(mf, envir = parent.frame())
    mt <- attr(mfi, "terms")
    ce[[i]]$X <- model.matrix(mt, mfi, contrasts)

    cov_t <- data.frame(
      t_start = ce[[i]]$covariates[, attr(ce, "special_cols")$t_start],
      t_end = ce[[i]]$covariates[, attr(ce, "special_cols")$t_end]
    )
    ce[[i]]$cov_t <- cov_t

    ev <- ce[[i]]$events[, c(attr(ce, "special_cols")$t_start,
                             attr(ce, "special_cols")$t_end,
                             attr(ce, "special_cols")$e_min,
                             attr(ce, "special_cols")$e_max)]
    colnames(ev) <- c("t_start", "t_end", "e_min", "e_max")
    ce[[i]]$ev <- ev

  }

  class(ce) <- c(class(ce), "co_events_frame")

  ce
}

