% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/poisson.R
\name{simulate_nonhomog_inversion}
\alias{simulate_nonhomog_inversion}
\title{Simulate non-homogeneous Poisson process via inversion}
\usage{
simulate_nonhomog_inversion(
  t_start,
  t_end,
  expect_cum_FUN,
  expect_cum_inverse_FUN,
  count_min = 0,
  count_max = Inf,
  ...
)
}
\arguments{
\item{t_start}{Start time of interval on which to simulate.}

\item{t_end}{End time of interval on which to simulate.}

\item{expect_cum_FUN}{Cumulative expected event count function.
Expected to be vectorized.}

\item{expect_cum_inverse_FUN}{Inverse cumulative expected event count function.
Expected to be vectorized.}

\item{count_min}{Minimum number of events.}

\item{count_max}{Maximum number of events.}

\item{...}{Additional arguments passed to \code{expect_cum_FUN} and
\code{expect_cum_inverse_FUN}}
}
\value{
A sorted vector of event times simulated from the specified process.
}
\description{
Simulate non-homogeneous Poisson process via inversion
}
\details{
Arguments \code{count_min} and \code{count_max} defines constraints on
realizations of the Poisson process. Probability law of the constrained
process is the same as would be derived via rejection sampling.
}
\examples{
simulate_nonhomog_inversion(
  t_start = 0, t_end = 100,
  expect_cum_FUN = expect_cum_weibull_tvc,
  expect_cum_inverse_FUN = expect_cum_weibull_tvc_inverse,
  k = 2, b = 1,
  t_breaks = c(50, 100),
  lin_pred = c(log(10), log(2)))

}
