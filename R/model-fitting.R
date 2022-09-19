# exponential model only (not weibull!)
log_likelihood_pois <- function(event_times,
                                observed_intervals,
                                log_rates,
                                t_breaks) {

  if (length(event_times) == 0) {
    log_rates_at_events <- 0
  } else {
    log_rates_at_events <- sapply(event_times, function(z) {
      log_rates[z <= t_breaks][1]
    })
  }



  interval_neg_expect_cum <-
    expect_cum_weibull_tvc(observed_intervals$t_start,
                           lin_pred = log_rates,
                           t_breaks = t_breaks) -
    expect_cum_weibull_tvc(observed_intervals$t_end,
                           lin_pred = log_rates,
                           t_breaks = t_breaks)

  sum(log_rates_at_events, interval_neg_expect_cum)
}

log_prior_coef <- function(coef) {
  0
}

#' @title Bayesian imputation for censored recurrent events
#'
#' @description Fits a non-homogeneous Poisson process model
#'              to interval count data with complex structure and censoring
#'              using Bayesian data augmentation.
#'
#' @param formula An object of class \code{formula}:
#'                a symbolic description of the model to be fit.
#'                Details of the model specification are given
#'                in the "Details" section.
#' @param co_events An object of class \code{co_events} containing
#'                  the pre-structured outcome and covariate data
#'                  for the model.
#' @param n_burn The number of burn-in iterations for the MCMC sampler.
#' @param n_keep The number of retained iterations from the MCMC sampler.
#'
#' @details
#' The \code{formula} argument does not require a left-hand side,
#' and any LHS supplied is ignored.
#' The variables in the right-hand side must be available in the
#' \code{co_events} argument.
#'
#' @return An object of the class \code{bicre}.
#'
#' @export
bicre <- function(formula, co_events,
                  n_burn = 100, n_keep = 1000,
                  z_force = NULL,
                  verbose = FALSE) {

  n_iter <- n_burn + n_keep

  iu <- co_events %>%
    co_events_frame(formula = formula) %>%
    build_imputation_units()

  coef_dim <- ncol(iu[[1]]$X)
  coef <- numeric(coef_dim)

  trace <- matrix(NA, nrow = n_iter, ncol = coef_dim)
  colnames(trace) <- colnames(iu[[1]]$X)
  trace[1, ] <- coef

  for (iter in 2 : n_iter) {
    if (verbose) {
      print(iter)
    }

    # imputation
    if (is.null(z_force)) {
      z <- sapply(iu,
                  impute_single_id,
                  coef = coef,
                  expect_cum_FUN = expect_cum_weibull_tvc,
                  expect_cum_inverse_FUN = expect_cum_weibull_tvc_inverse,
                  fail_mode = TRUE,
                  simplify = FALSE)
    } else {
      z <- z_force
    }


    # regression coefficients

    coef_cur <- coef
    coef_pro <- coef + rnorm(coef_dim, mean = 0, sd = 0.04)

    log_lik_cur <- sum(sapply(
      1 : length(z),
      function(i) {
        log_likelihood_pois(
          event_times = z[[i]],
          observed_intervals = attributes(iu[[i]]$ev)$covered_ints,
          log_rates = iu[[i]]$X %*% coef_cur,
          t_breaks = iu[[i]]$cov_t$t_end
        )
      }))
    log_lik_pro <- sum(sapply(
      1 : length(z),
      function(i) {
        log_likelihood_pois(
          event_times = z[[i]],
          observed_intervals = attributes(iu[[i]]$ev)$covered_ints,
          log_rates = iu[[i]]$X %*% coef_pro,
          t_breaks = iu[[i]]$cov_t$t_end
        )
      }))

    log_prior_cur <- log_prior_coef(coef_cur)
    log_prior_pro <- log_prior_coef(coef_pro)

    if (rbinom(1, 1, prob = min(1,
                                exp(log_lik_pro + log_prior_pro -
                                    log_lik_cur - log_prior_cur)))) {
      coef <- coef_pro
    }

    # record
    trace[iter, ] <- coef
  }

  class(trace) <- "bicre"
  trace

}

