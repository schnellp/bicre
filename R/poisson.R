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

  stopifnot(all(t >= 0))
  stopifnot(all(t <= tail(t_breaks, 1)))
  stopifnot(length(lin_pred) == length(t_breaks))

  sapply(t, function(t) {
    t_breaks <- t_breaks[1 : (sum(t_breaks < t) + 1)]
    lin_pred <- lin_pred[1 : length(t_breaks)]

    t_breaks <- c(t_breaks[-length(t_breaks)], t)

    sum(diff(c(0, expect_cum_weibull(t_breaks, k = k, b = b))) * exp(lin_pred))
  })
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
  stopifnot(length(lin_pred) == length(t_breaks))

  m_breaks <- cumsum(diff(c(0, expect_cum_weibull(t_breaks, k = k, b = b))) * exp(lin_pred))

  stopifnot(m <= tail(m_breaks, 1))

  sapply(m, function(m) {
    left_index <- which.max(m < m_breaks)
    t_left <- c(0, t_breaks)[left_index]

    ((m - c(0, m_breaks)[left_index]) / (b * exp(lin_pred[left_index])) + t_left^k)^(1 / k)
  })
}

#' @title Truncated Poisson generation
#'
#' @description
#' Generates random draws from un-, left-, right-, or double-truncated
#' Poisson distribution with parameter \code{lambda} truncated to \code{[min, max]}
#' (inclusive). The Poisson distribution is parameterized identically to
#' \code{stats::rpois}, but \code{lambda} cannot be a vector.
#'
#' @param n Number of random values to return.
#' @param lambda Non-negative (scalar) mean of un-truncated distribution.
#' @param min The minimum allowed return value.
#' @param max The maximum allowed return value.
#' @param parallel_draws For left-only-truncated draws, the number of
#'                       simultaneous draws by the rejection sampler.
#' @param max_tries For left-only-truncated draws, the maximum number of
#'                  iterations of the rejection sampler.
#'
#' @return A vector of draws.
#'
#' @examples
#' rpois_trunc(100, lambda = 1)
#' rpois_trunc(100, lambda = 1, min = 1)
#' rpois_trunc(100, lambda = 1, max = 10)
#' rpois_trunc(100, lambda = 1, min = 1, max = 10)
rpois_trunc <- function(n, lambda, min = 0, max = Inf,
                        parallel_draws = n * 2, max_tries = 1e3) {

  stopifnot(length(lambda) == 1)
  stopifnot(lambda >= 0)

  if (lambda == 0) {
    if (min > 0) {
      stop("Improper distribution (lambda = 0, min > 0).")
    } else {
      return(rep(0, n))
    }
  }

  if (max < Inf) {
    return(
      sample(min : max, size = n, replace = TRUE,
             prob = dpois(min : max, lambda = lambda))
      )
  } else if (min > 0) {
    x <- numeric(0)
    for (i in 1 : max_tries) {
      m <- max(ceiling(min - 1 - lambda), 0)
      y <- rpois(parallel_draws, lambda)
      x_pro <- m + y
      accept_prob <- exp(lfactorial(y) - lfactorial(x_pro) + lfactorial(min) -
                           lfactorial(min - m)) * as.numeric(x_pro >= min)

      x <- c(x, x_pro[as.logical(rbinom(parallel_draws, 1, accept_prob))])

      if (length(x) >= n) {
        return(x[1 : n])
      }
    }

    stop("Maximum number of rejection sampler iterations reached.")
  } else {
    return(rpois(n, lambda = lambda))
  }
}

#' @title Simulate non-homogeneous Poisson process via inversion
#'
#' @param t_start Start time of interval on which to simulate.
#' @param t_end End time of interval on which to simulate.
#' @param expect_cum_FUN Cumulative expected event count function.
#'                       Expected to be vectorized.
#' @param expect_cum_inverse_FUN Inverse cumulative expected event count function.
#'                               Expected to be vectorized.
#' @param ... Additional arguments passed to \code{expect_cum_FUN} and
#'            \code{expect_cum_inverse_FUN}
#'
#' @return A sorted vector of event times simulated from the specified process.
#'
#' @examples
#' simulate_nonhomog_inversion(
#'   t_start = 0, t_end = 100,
#'   expect_cum_FUN = expect_cum_weibull_tvc,
#'   expect_cum_inverse_FUN = expect_cum_weibull_tvc_inverse,
#'   k = 2, b = 1,
#'   t_breaks = c(50, 100),
#'   lin_pred = c(log(10), log(2)))
simulate_nonhomog_inversion <- function(t_start, t_end,
                                        expect_cum_FUN, expect_cum_inverse_FUN,
                                        ...) {

  m_start <- expect_cum_FUN(t_start, ...)
  m_end <- expect_cum_FUN(t_end, ...)

  N <- rpois(1, lambda = m_end - m_start)
  z <- sort(runif(N, min = m_start, max = m_end))

  expect_cum_inverse_FUN(z, ...)
}
