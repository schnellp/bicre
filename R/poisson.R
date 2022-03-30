#' @title Weibull expected cumulative event count
#'
#' @description Expected cumulative event count of Poisson process
#'              with Weibull hazard.
#'              Default parameterization matches that of \code{stats::Weibull}.
#'
#' @param t Time at which to evaluate expected cumulative count.
#' @param shape Shape of Weibull distribution, as in \code{stats::Weibull}.
#' @param scale Scale of Weibull distribution, as in \code{stats::Weibull}.
#' @param k Alternative parameterization (shape).
#' @param b Alternative parameterization (scale).
#'
#' @return Expected cumulative event count at time \code{t}.
#'
#' @examples
#' expect_cum_weibull(4, 2, 3)
#' expect_cum_weibull(4, k = 2, b = 1 / 3^2)
expect_cum_weibull <- function(t,
                               shape = 1, scale = 1,
                               k = shape, b = scale^(-shape)) {
  b * t^k
}

#' @title Inverse of Weibull expected cumulative event count
#'
#' @description Inverse of expected cumulative event count of Poisson process
#'              with Weibull hazard.
#'              Default parameterization matches that of \code{stats::Weibull}.
#'
#' @param m Expected (mean) Weibull cumulative event count.
#' @param shape Shape of Weibull distribution, as in \code{stats::Weibull}.
#' @param scale Scale of Weibull distribution, as in \code{stats::Weibull}.
#' @param k Alternative parameterization (shape).
#' @param b Alternative parameterization (scale).
#'
#' @return Time at which Weibull expected cumulative event count equals \code{m}.
#'
#' @examples
#' expect_cum_weibull_inverse(expect_cum_weibull(4, 2, 3))
expect_cum_weibull_inverse <- function(m,
                                       shape = 1, scale = 1,
                                       k = shape, b = scale^(-shape)) {
  (m / b)^(1 / k)
}

#' @title Weibull expected cumulative event count with time-varying covariates
#'
#' @description
#' Expected cumulative event count
#' given Weibull baseline expected cumulative event count and
#' piecewise constant multiplicative factor exp(lin_pred).
#'
#' @param t Time at which to evaluate expected cumulative count.
#' @param lin_pred Vector of log scaling factors on baseline expected
#'                 cumulative event count function over intervals.
#' @param t_breaks Vector of right (closed) endpoints of intervals
#'                 over which lin_pred applies.
#'                 Elements must be strictly increasing and at least one
#'                 must be \code{>= t}.
#' @param shape Shape of Weibull distribution, as in \code{stats::Weibull}.
#' @param scale Scale of Weibull distribution, as in \code{stats::Weibull}.
#' @param k Alternative parameterization (shape).
#' @param b Alternative parameterization (scale).
#'
#' @return Expected cumulative event count at time \code{t}.
#'
#' @examples
#' expect_cum_weibull_tvc(3, lin_pred = 1 : 5, t_breaks = 1 : 5, k = 1, b = 1)
expect_cum_weibull_tvc <- function(t,
                                   lin_pred = 0,
                                   t_breaks = t,
                                   shape = 1, scale = 1,
                                   k = shape, b = scale^(-shape)) {

  stopifnot(t >= 0)
  stopifnot(length(t) == 1)
  stopifnot(t <= tail(t_breaks, 1))
  stopifnot(length(lin_pred) == length(t_breaks))

  t_breaks <- t_breaks[1 : (sum(t_breaks < t) + 1)]
  lin_pred <- lin_pred[1 : length(t_breaks)]

  t_breaks[length(t_breaks)] <- t

  sum(diff(c(0, expect_cum_weibull(t_breaks, k = k, b = b))) * exp(lin_pred))
}

#' @title Inverse of Weibull expected cumulative event count with time-varying covariates
#'
#' @description
#' Inverse of expected cumulative event count
#' given Weibull baseline expected cumulative event count and
#' piecewise constant multiplicative factor exp(lin_pred).
#'
#' @param m Expected (mean) Weibull cumulative event count.
#' @param lin_pred Vector of log scaling factors on baseline expected
#'                 cumulative event count function over intervals.
#' @param t_breaks Vector of right (closed) endpoints of intervals
#'                 over which lin_pred applies.
#'                 Elements must be strictly increasing and at least one
#'                 must be \code{>= t}.
#' @param shape Shape of Weibull distribution, as in \code{stats::Weibull}.
#' @param scale Scale of Weibull distribution, as in \code{stats::Weibull}.
#' @param k Alternative parameterization (shape).
#' @param b Alternative parameterization (scale).
#'
#' @return Time at which Weibull expected cumulative event count equals \code{m}.
#'
#' @examples
#' expect_cum_weibull_tvc_inverse(
#'   expect_cum_weibull_tvc(3, lin_pred = 1 : 5, t_breaks = 1 : 5, k = 1, b = 1),
#'   lin_pred = 1 : 5, t_breaks = 1 : 5, k = 1, b = 1)
expect_cum_weibull_tvc_inverse <- function(m,
                                           lin_pred = 0,
                                           t_breaks = Inf,
                                           shape = 1, scale = 1,
                                           k = shape, b = scale^(-shape)) {

  stopifnot(m >= 0)
  stopifnot(length(m) == 1)
  stopifnot(length(lin_pred) == length(t_breaks))

  m_breaks <- cumsum(diff(c(0, expect_cum_weibull(t_breaks, k = k, b = b))) * exp(lin_pred))

  stopifnot(m <= tail(m_breaks, 1))

  left_index <- which.max(m < m_breaks)
  t_left <- c(0, t_breaks)[left_index]

  ((m - c(0, m_breaks)[left_index]) / (b * exp(lin_pred[left_index])) + t_left^k)^(1 / k)
}


