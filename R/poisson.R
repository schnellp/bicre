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
#' @param ... Accepts other arguments; not currently used except
#'            to prevent errors from e.g., simulate_nonhomog_inversion
#'
#' @return Expected cumulative event count at time \code{t}.
#'
#' @examples
#' expect_cum_weibull(4, 2, 3)
#' expect_cum_weibull(4, k = 2, b = 1 / 3^2)
#'
#' @export
expect_cum_weibull <- function(t,
                               shape = 1, scale = 1,
                               k = shape, b = scale^(-shape),
                               ...) {
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
#' @param ... Accepts other arguments; not currently used except
#'            to prevent errors from e.g., simulate_nonhomog_inversion
#'
#' @return Time at which Weibull expected cumulative event count equals \code{m}.
#'
#' @examples
#' expect_cum_weibull_inverse(expect_cum_weibull(4, 2, 3))
#'
#' @export
expect_cum_weibull_inverse <- function(m,
                                       shape = 1, scale = 1,
                                       k = shape, b = scale^(-shape),
                                       ...) {
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
#' @param ... Accepts other arguments; not currently used except
#'            to prevent errors from e.g., simulate_nonhomog_inversion
#'
#' @return Expected cumulative event count at time \code{t}.
#'
#' @examples
#' expect_cum_weibull_tvc(3, lin_pred = 1 : 5, t_breaks = 1 : 5, k = 1, b = 1)
#'
#' @export
expect_cum_weibull_tvc <- function(t,
                                   lin_pred = 0,
                                   t_breaks = Inf,
                                   shape = 1, scale = 1,
                                   k = shape, b = scale^(-shape),
                                   ...) {

  if (length(t) == 0) {
    return(numeric(0))
  }

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



#' @title expect_cum_weibull_tvc written in Rcpp
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
#' @param ... Accepts other arguments; not currently used except
#'            to prevent errors from e.g., simulate_nonhomog_inversion
#'
#' @return Expected cumulative event count at time \code{t}.
#'
#' @export
expect_cum_weibull_tvc_Rcpp <- function(t,
                                         lin_pred = 0,
                                         t_breaks = Inf,
                                         shape = 1, scale = 1,
                                         k = shape, b = scale^(-shape),
                                         ...){
  expect_cum_weibull_tvc_cpp(t, lin_pred,t_breaks, k, b)
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
#' @param ... Accepts other arguments; not currently used except
#'            to prevent errors from e.g., simulate_nonhomog_inversion
#'
#' @return Time at which Weibull expected cumulative event count equals \code{m}.
#'
#' @examples
#' expect_cum_weibull_tvc_inverse(
#'   expect_cum_weibull_tvc(3, lin_pred = 1 : 5, t_breaks = 1 : 5, k = 1, b = 1),
#'   lin_pred = 1 : 5, t_breaks = 1 : 5, k = 1, b = 1)
#'
#' @export
expect_cum_weibull_tvc_inverse <- function(m,
                                           lin_pred = 0,
                                           t_breaks = Inf,
                                           shape = 1, scale = 1,
                                           k = shape, b = scale^(-shape),
                                           ...) {

  if (length(m) == 0) {
    return(numeric(0))
  }

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

#' @title expect_cum_weibull_tvc_inverse written in Rcpp
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
#' @param ... Accepts other arguments; not currently used except
#'            to prevent errors from e.g., simulate_nonhomog_inversion
#'
#' @return Time at which Weibull expected cumulative event count equals \code{m}.
#'
#' @export
expect_cum_weibull_tvc_inverse_Rcpp <- function(m,
                                                lin_pred = 0,
                                                t_breaks = Inf,
                                                shape = 1, scale = 1,
                                                k = shape, b = scale^(-shape),
                                                ...){
  expect_cum_weibull_tvc_inverse_cpp(m, lin_pred,t_breaks, k, b)
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
#' @param fail_mode Logical. If \code{max_tries} is exceeded, should the
#'                  function return the mode (\code{floor(lambda)}) of the
#'                  distribution (with a warning) rather than throwing an error?
#' @param rejection_sampler Logical (default `FALSE`). Avoid specialized algorithm
#'                          and use a simple rejection sampler instead. Available
#'                          only for `n = 1`.
#' @param ... Accepts other arguments; not currently used except
#'            to prevent errors from e.g., simulate_nonhomog_inversion
#'
#' @details
#' If \code{max < Inf} then \code{sample} is used to draw from the finite sample space.
#' If \code{min > 0} and \code{max == Inf} then the improved rejection sampler
#' of Geyer (2021) is used.
#'
#' @return A vector of draws.
#'
#' @examples
#' rpois_trunc(100, lambda = 1)
#' rpois_trunc(100, lambda = 1, min = 1)
#' rpois_trunc(100, lambda = 1, max = 10)
#' rpois_trunc(100, lambda = 1, min = 1, max = 10)
#'
#' @references
#' Geyer, CJ. "Lower-Truncated Poisson and Negative Binomial Distributions"
#' 2021.
#' \url{https://cran.r-project.org/web/packages/aster/vignettes/trunc.pdf}
#'
#' @export
rpois_trunc <- function(n, lambda, min = 0, max = Inf,
                        parallel_draws = n * 2, max_tries = 1e3,
                        fail_mode = FALSE,
                        rejection_sampler = FALSE,
                        ...) {

  # stopifnot(length(lambda) == 1)
  # stopifnot(lambda >= 0)

  if(length(lambda) != 1){
    stop("0 or more than 1 lambda has been given")
  }

  if(lambda < 0){
    stop("a negative lambda has been given")
  }

  if (lambda == 0) {
    if (min > 0) {
      stop("Improper distribution (lambda = 0, min > 0).")
    } else {
      return(rep(0, n))
    }
  }

  if (min == max) {
    return(min)
  } else if (rejection_sampler == TRUE & n == 1) {
    for (try in 1 : max_tries) {
      x <- rpois(n = 1, lambda = lambda)
      if (min <= x & x <= max) {
        return(n)
      }
    }

    if (fail_mode) {
      warning(paste("Maximum number of rejection sampler iterations reached.",
                    "Returning mode."))
      if (lambda > max) {
        return(rep(floor(max), n))
      } else {
        return(rep(floor(max(min, lambda)), n))
      }
    } else {
      stop("Maximum number of rejection sampler iterations reached.")
    }
  } else if (max < Inf) {
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
      accept_prob <- pmin(1, exp(lfactorial(y) - lfactorial(x_pro) + lfactorial(min) -
                           lfactorial(min - m)) * as.numeric(x_pro >= min))
      x <- c(x, x_pro[as.logical(rbinom(parallel_draws, 1, accept_prob))])

      if (length(x) >= n) {
        return(x[1 : n])
      }
    }

    if (fail_mode) {
      warning(paste("Maximum number of rejection sampler iterations reached.",
                    "Returning mode."))
      if (lambda > max) {
        return(rep(floor(max), n))
      } else {
        return(rep(floor(max(min, lambda)), n))
      }
    } else {
      stop("Maximum number of rejection sampler iterations reached.")
    }

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
#' @param count_min Minimum number of events.
#' @param count_max Maximum number of events.
#' @param sorted_times (Logical) Should times be sorted before returning?
#' @param ... Additional arguments passed to \code{expect_cum_FUN} and
#'            \code{expect_cum_inverse_FUN}
#'
#' @details
#' Arguments \code{count_min} and \code{count_max} defines constraints on
#' realizations of the Poisson process. Probability law of the constrained
#' process is the same as would be derived via rejection sampling.
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
#'
#' @export
simulate_nonhomog_inversion <- function(t_start, t_end,
                                        expect_cum_FUN, expect_cum_inverse_FUN,
                                        count_min = 0, count_max = Inf,
                                        sorted_times = FALSE,
                                        ...) {

  m_start <- expect_cum_FUN(t_start, ...)
  m_end <- expect_cum_FUN(t_end, ...)

  if(count_min == 1 & count_max == Inf){
    z <- generate_events_no0_poisson1(m_start, m_end)
  }else{
    N <- rpois_trunc(1, lambda = m_end - m_start,
                     min = count_min, max = count_max,
                     ...)
    z <- runif(N, min = m_start, max = m_end)

    if(sorted_times){
      z <- sort(z)
    }
  }


  expect_cum_inverse_FUN(z, ...)
}

# Simulate non-homogeneous Poisson process via inversion with condition on
# multiple non-contiguous intervals
simulate_nonhomog_inversion_multi_interval <- function(interval_df,
                                                       expect_cum_FUN, expect_cum_inverse_FUN,
                                                       count_min = 0, count_max = Inf,
                                                       sorted_times = FALSE,
                                                       ...){
  m_starts <- expect_cum_FUN(interval_df$t_start, ...)
  m_ends <- expect_cum_FUN(interval_df$t_end, ...)
  int_length <-  m_ends - m_starts
  int_length_cumsum <- cumsum(int_length)
  total_length <- last(int_length_cumsum)

  if(count_min == 1 & count_max == Inf){
    z_raw <- generate_events_no0_poisson1(0, total_length)
  }else{
    N <-  rpois_trunc(1, lambda = total_length,
                      min = count_min, max = count_max,
                      ...)
    z_raw <- runif(N, min = 0, max = total_length)
    if(sorted_times){
      z_raw <- sort(z_raw)
    }
  }

  z_raw_int_index <- findInterval(z_raw,int_length_cumsum, left.open = TRUE) + 1
  z_raw_add <- z_raw - c(0, int_length_cumsum)[z_raw_int_index]
  z <- m_starts[z_raw_int_index] + z_raw_add


  expect_cum_inverse_FUN(z, ...)
}

