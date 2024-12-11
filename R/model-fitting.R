# likelihood and prior functions

log_lik_uis_var <- function(uis, uis_var){
  sum(dgamma(uis, shape =  1/uis_var, scale = uis_var, log = T))
}

normal_coef <- function(coef, coef_mean = 0, coef_sd = 5){
  sum(dnorm(coef, mean = coef_mean, sd = coef_sd, log = T ))
}

# not lognormal densities (missing jacobian)
# because random walk moves are on log scale
# without change-of-variables
log_normal_k <- function(k, logk_mean = 0, logk_sd = 3){
  dnorm(log(k), mean = logk_mean, sd = logk_sd, log = T)
}

log_normal_b <- function(b, logb_mean = 0, logb_sd = 3){
  dnorm(log(b), mean = logb_mean, sd = logb_sd, log = T)
}

log_normal_uis_var <- function(uis_var, log_uis_var_mean = 0, log_uis_var_sd = 3){
  dnorm(log(uis_var), mean = log_uis_var_mean, sd = log_uis_var_sd, log = T)
}

Xmatrix_lik <- function(event_time_list, iu, mc_cores = 1L){

  Xmatrix_list <-  mclapply(1:length(iu), function(ii) iu[[ii]]$covariates[
    1 + findInterval(event_time_list[[ii]],
                     iu[[ii]]$cov_t$t_end,left.open = T),,drop=FALSE],
    mc.cores = mc_cores)

  (rbindlist(Xmatrix_list)  |>
      mutate("(Intercept)" = 1, .before = 1) |> as.matrix())[,colnames(iu[[1]]$X)]
}



event.times.loglik <- function(k, coef, uis, z_counts.sum, log.z.sum, XAllMatrix,
                               iu, mc_cores = 1L){
  log(k) * z_counts.sum + (k - 1) * log.z.sum +  sum(XAllMatrix %*% coef) +
    (1:length(iu) |> mclapply(function(i){
      expect_cum_weibull_tvc_Rcpp(attributes(iu[[i]]$ev)$covered_ints$t_start,
                                  lin_pred = iu[[i]]$X %*% coef,
                                  t_breaks = iu[[i]]$cov_t$t_end,
                                  k = k,
                                  b = uis[i]) -
        expect_cum_weibull_tvc_Rcpp(attributes(iu[[i]]$ev)$covered_ints$t_end,
                                    lin_pred = iu[[i]]$X %*% coef,
                                    t_breaks = iu[[i]]$cov_t$t_end,
                                    k = k,
                                    b = uis[i])
    }, mc.cores = mc_cores) |> unlist() |> sum())



}


#' @title Bayesian Data augmentation for censored recurrent events with Weibull baseline hazard
#' @description This function applies the Bayesian data augmentation method on the
#'  input censored recurrent events data assuming
#' the recurrent events are from a Poisson process with the following intensity:
#'
#' \deqn{\lambda(t|u_i, \boldsymbol{x(t)}) = u_i \cdot bkt^{k-1} \cdot e^{\boldsymbol{x(t)\beta}}}
#' Where \eqn{u_i} is the random effect, \eqn{bkt^{k-1}} is the Weibull baeline hazard,
#' \eqn{\boldsymbol{x(t)}} is the covariates matrix, and \eqn{\boldsymbol{\beta}} is the covariates coefficient vector.
#'
# from \code{data_covariates} and \code{data_events} or the prerpocessed data \code{iu})
# Although \link{bicre} already contain steps for preprocessing data,
# it is recommended to do this outside it so you don't need to rerun the preprocessing
# every time you run  \link{bicre}.
#'
#' @param formula A formula with a format of ~ \eqn{A_1 + A_2 + \cdots + A_n} where \eqn{A_1, A_2, ..., A_n}
#' are chosen column names of \code{data_covariates} representing user chosen predictors for the analysis.
#'  E.g. if \code{data_covariates} contain predictor column names "trt", "age", "gender"
#' a formula of ~ trt + age means using "trt" and "age" column in the \code{data_covariates} dataframe
#' @param data_covariates A dataframe containing information for covariates values in different time intervals.
#' The first column and the last two columns are named "id", "t_start", and "t_end".
#' The other columns are covariate values.
#' Each row represents individual "id"'s covariate values between time "t_start" and "t_end".
#' E.g. if \code{data_covariates} contains predictor columns "trt", "age" except for "id", "t_start", and "t_end" columns.
#' The first row of \code{data_covariates} is 1, TRUE, 30, 0, 10. Then it means the individual with id 1 are have "trt" value 1
#' and "age" of 30 between time 0 and time 10.
#' @param data_events A dataframe containing observations for the study outcome recurrent events.
#' \code{data_events} has 5 columns with names respectively "id", "t_start", "t_end", "e_min", and "e_max".
#' Each row of \code{data_events} is an observation representing that the recurrent event counts is between
#'  "e_min" and "e_max" for the time interval ("e_min", "e_max"] for individual "id".
#' @param id The name of the special `id` column in `data_covariates` and `data_events`.
#'           See \link{co_events}.
#' @param t_start The name of the special `t_start` column in `data_covariates` and `data_events`.
#'                See \link{co_events}.
#' @param t_end The name of the special `t_end` column in `data_covariates` and `data_events`.
#'              See \link{co_events}.
#' @param e_min The name of the special `e_min` column in `data_events`.
#'                See \link{co_events}.
#' @param e_max The name of the special `e_max` column in `data_events`.
#'                See \link{co_events}.
#' @param fill A named list specifying the value of each time-varying covariate
#'             that should be used for newly added rows. See \link{co_events}
#' @param check_cov_cover_ev Logical value whether to check covariate time range cover event time range.  See \link{co_events}
#' @param n_burn_in Number of MCMC iterations used as burn-in period.
#' @param n_keep Number of MCMC iterations used for inference.
#' @param prior_dist_k A R function for the probability density function of k's prior distribution.
#' This function takes in the current value of k and output k's prior probability density at current value of k.
#' Default prior is a log normal distribution with mean 0 and standard deviation 3.
#' @param prior_dist_b A R function for the probability density function of b's prior distribution.
#' This function takes in the current value of b and output b's prior probability density at current value of b.
#' Default prior is a log normal distribution with mean 0 and standard deviation 3.
#' @param prior_dist_coef A R function for the probability density function of the covariate coefficients vector \eqn{\boldsymbol{\beta}}.
#' This function takes in the current vector value of covariate coefficients \eqn{\boldsymbol{\beta}} and output
#' its prior probability density.
#' Default prior assumes covariate coefficients are indepdent and all follow the same normal distribution with mean 0 and standard deviation 5.
#' @param prior_dist_uis_var A R function for the probability density function of \eqn{\sigma^2}'s prior distribution.
#'  Here \eqn{\sigma^2} is the random effects' variance.
#' This function takes in the current value of random effects variance and output \eqn{\sigma^2}'s prior probability density at its current value.
#' Default prior is a log normal distribution with mean 0 and standard deviation 3.
#' @param truncate_uis If FALSE (default) or 0, do not truncate the random effect distribution.
#'                     If >0, left-truncate at given value.
#' @param seed The seed to set for running the MCMC iteration.
#' @param tiny_diff The threshold for differentiate between T-zones of tiny length (tiny T-zones) and other T-zones.
#' If a T-zone has a length smaller than tiny_diff, it is treated as a tiny T-zone for the preprocessing otherwise it is treated as other T-zones.
#' @param prepare_style The method used to rank non-tiny T-zones when determining the disjoint set of compound imputation-units.
#' Available options are:
#' \enumerate{
#' \item "target_emin": order non-tiny T-zones by decreasing minimum count when they are chosen into the disjoint set.
#' This is the default option since it is the fastest by our simulation study.
#' \item "target_min": order T-zones by decreasing minimum count rate when they are chosen into the disjoint set.
#' \item "target_max": order T-zones by increasing maximum count rate when they are chosen into the disjoint set.
#' \item "normal": order T-zones by "t_start" when they are chosen into the disjoint set.
#' }
#' @param iu An object of class \code{co_events_frame} containing special intervals and sets for optimized rejection sampling.
#' This object is derived using \link{build_imputation_units} function.
#' When this argument has been given a value, the following arguments are all ignored:
#'  \code{data_covariates}, \code{data_events}, \code{fill}, \code{check_cov_cover_ev}, \code{formula}, \code{tiny_diff}, \code{prepare_style }
#' Otherwise \code{iu}'s default value is calculated from these arguments.
#' @param trace_start A list containing the starting values of parameters.
#' This list has 3 elements named \code{k}, \code{uis_var}, and \code{coef}, they are seperately
#' \enumerate{
#' \item The starting value of k. Default to 1
#' \item The starting value of random effects variance \eqn{\sigma^2}. Default to 0.5
#' \item The starting value of covariate coefficients \eqn{\boldsymbol{\beta}}. This is a vector. Default to a vector with all values being 0.
#' }
#'
#'
#' @param uis_start A vector with a length equal to the number of individuals.
#' This vector specifies the starting values of random effects.
#' Default to a vector with all values being 1.
#' @param fail_mode Logical. See \link{rpois_trunc}.
#' @param run_and_save Logical. Whether to save the results every \code{iter_per_save} number of iterations. Default TRUE
#' @param run_and_print Logical. Whether to print iteration numbers every \code{iter_per_print} number of iterations. Default FALSE
#' @param iter_per_save Integer.
#' When \code{run_and_save} is TRUE, save the results every \code{iter_per_save} number of iterations.
#' Otherwise this argument is ignored.
#' Default 1000
#' @param iter_per_print Integer.
#' When \code{run_and_print} is TRUE,
#' print the current iteration number every \code{iter_per_print} number of iterations.
#' Otherwise this argument is ignored.
#' Default to be equal to \code{iter_per_save}.
#' @param save_folder The path to save the results when \code{iter_per_save} is TRUE. Otherwise this argument is ignored.
#' @param mc_cores The number of cores to run on parallel for the algorithm. Default to 1. mc_cores > 1 is only avaliable for R on Linux or macOS.
#' @param keep_going Logical. When there is errors, whether repeat the current iteration and continue use \link{bicre_continue}.
#'
#'
#' @return An list with two elements.
#'
#' \describe{
#' \item{\code{trace}}{A dataframe that lists the posterior samples of all parameters;}
#' \item{\code{accept_rate}}{A named vector lists the acceptance rates of k,b and \eqn{\boldsymbol{\beta}} and \eqn{\sigma^2};}
#' }
#'
#' @export
#'
#'
bicre <- function(formula, data_covariates, data_events,
                          id = id,
                          t_start = t_start,
                          t_end = t_end,
                          e_min = e_min,
                          e_max = e_max,
                          fill = NA, check_cov_cover_ev = TRUE,
                          n_burn_in = 5000, n_keep = 10000,
                          prior_dist_k = log_normal_k,
                          prior_dist_b = log_normal_b,
                          prior_dist_coef = normal_coef,
                          prior_dist_uis_var =  log_normal_uis_var,
                          truncate_uis = FALSE,
                          seed = 11374,
                          tiny_diff = NULL,
                          prepare_style = "normal",
                          iu = NULL,
                          trace_start = list(k = 1,
                                             uis_var = 0.5,
                                             coef = rep(0,times = ncol(iu[[1]]$X))),
                          uis_start = rep(1, iu |> length()),
                          fail_mode = FALSE,
                          run_and_save = TRUE,
                          run_and_print = FALSE,
                          iter_per_save = 1000,
                          iter_per_print = iter_per_save,
                          save_folder = NULL,
                          mc_cores = 1L,
                          keep_going = FALSE){

  # decisions on print and save iterations
  if(run_and_save){
    if(!is.null(save_folder)){
      dir.create(file.path(save_folder), showWarnings = FALSE, recursive = TRUE)
      save_prefix <- file.path(save_folder, "trace")
    }else{
      save_prefix <- "trace"
    }

    if(iter_per_print == iter_per_save){
      run_and_print <- FALSE
    }
  }


  # prepare imputation units if not provided
  if (is.null(iu)) {
    ce <- co_events(
      data_covariates = data_covariates,
      data_events = data_events,
      id = {{ id }},
      t_start = {{ t_start }},
      t_end = {{ t_end }},
      e_min = {{ e_min }},
      e_max = {{ e_max }},
      fill = fill,
      check_cov_cover_ev = check_cov_cover_ev,
      special_cols = c(
        "id" = deparse(substitute(id)),
        "t_start" = deparse(substitute(t_start)),
        "t_end" = deparse(substitute(t_end)),
        "e_min" = deparse(substitute(e_min)),
        "e_max" = deparse(substitute(e_max))
      ))
    cef <- co_events_frame(ce, formula = formula)
    iu <- build_imputation_units(cef, tiny_diff = tiny_diff, prepare_style = prepare_style)
  }


  # preparation for MCMC
  n_iter <- n_burn_in + n_keep
  trace <- matrix(NA, nrow = n_iter, ncol = unlist(trace_start) |> length())
  colnames(trace) <- c("log(k)", "log(uis_var)", "log(b)" ,colnames(iu[[1]]$X)[-1])
  trace[1, ] <- unlist(trace_start)
  trace[1, c("log(k)", "log(uis_var)")]  <- log(trace[1, c("log(k)", "log(uis_var)")])

  accept <- rep(0, 2)
  names(accept) <- c("uis_var", "kcoef")

  list2env(trace_start, envir = environment())
  uis <- uis_start

  est.cov <- diag(1 + length(coef)) * 0.25
  S.logkcoef <- t(chol(est.cov))
  S_uis_var <-  0.25

  #### MCMC run
  set.seed(seed)



  we_catch <- tryCatch(expr = {
    not_finish <- TRUE
    for (iter in 2 : n_iter) {

      if(run_and_print){
        if(iter %% iter_per_print == 0){
          print(iter)
        }
      }




      # imputation


      if(iter == 2){
        z_list <- list()

        for (i in 1:length(iu)) {
          t <- 0
          repeat{
            t <- t + 1
            event_times_i <-impute_single_id_initialize(co_events_frame_single = iu[[i]],
                                                        b = uis[i], coef = coef, expect_cum_FUN = expect_cum_weibull_tvc_Rcpp,
                                                        expect_cum_inverse_FUN = expect_cum_weibull_tvc_inverse_Rcpp,
                                                        k = k)
            if(class(event_times_i) ==  "numeric"){
              z_list[[i]] <- event_times_i
              break
            }else if(event_times_i == "Sequential sampling: larger ui needed"){
              uis[i] <- uis[i] * 2
            }else if(event_times_i == "Sequential sampling: smaller ui needed"){
              uis[i] <- uis[i] / 1.5
            }else{
              stop("imputation initialize has unexpected errors")
            }

            if(t > 15){
              stop("imputation initialize failed to find a proper starting value")
            }

          }

        }

      }else{

        z_list <- mcmapply(impute_single_id, co_events_frame_single = iu, b = uis,
                         MoreArgs = list(coef = coef, expect_cum_FUN = expect_cum_weibull_tvc_Rcpp,
                                         expect_cum_inverse_FUN = expect_cum_weibull_tvc_inverse_Rcpp,
                                         k = k, fail_mode = fail_mode),
                         mc.cores = mc_cores)
      }


      z_counts <- z_list |> sapply(length)

      # uis update
      shapes <- z_counts + 1/uis_var
      scales <- uis_var /
        (1 + uis_var * (
          mcmapply(function(i){
            sum( - expect_cum_weibull_tvc_Rcpp(attributes(iu[[i]]$ev)$covered_ints$t_start,
                                               lin_pred = iu[[i]]$X %*% coef,
                                               t_breaks = iu[[i]]$cov_t$t_end,
                                               k = k) +
                   expect_cum_weibull_tvc_Rcpp(attributes(iu[[i]]$ev)$covered_ints$t_end,
                                               lin_pred = iu[[i]]$X %*% coef,
                                               t_breaks = iu[[i]]$cov_t$t_end,
                                               k = k))
          },1:length(iu), mc.cores = mc_cores)
        )
        )

      if (truncate_uis) {
        uis <- ltrgamma(length(z_counts), shapes = shapes, scales = scales,
                        truncate = truncate_uis)
      } else {
        uis <- rgamma(length(z_counts), shape = shapes, scale = scales)
      }



      # uis_var update
      u_uis_var <- rnorm(1,mean = 0, sd = 1)
      uis_var_ori <- uis_var
      uis_var_pro <- (log(uis_var_ori) + S_uis_var * u_uis_var) |> exp()

      log_lik_uis_var_ori <- log_lik_uis_var(uis, uis_var_ori) +
        prior_dist_uis_var(uis_var_ori)
      log_lik_uis_var_pro <- log_lik_uis_var(uis, uis_var_pro) +
        prior_dist_uis_var(uis_var_pro)

      acceptance_prob_uis_var <- min(1, exp((log_lik_uis_var_pro) - (log_lik_uis_var_ori)))

      if (rbinom(1, 1, prob = acceptance_prob_uis_var)) {
        uis_var <- uis_var_pro
        accept["uis_var"] <- accept["uis_var"] + 1
      }

      if(iter <= n_burn_in){
        S_uis_var <-  ramcmc::adapt_S(S_uis_var, u_uis_var, acceptance_prob_uis_var, iter - 1) |>  as.vector()
      }


      # k and coef update
      kcoef.ori <- c(k, coef)
      logkcoef.ori <- c(log(k), coef)
      u.logkcoef <- rnorm(length(logkcoef.ori), 0, sd = 1)
      logkcoef.pro <- logkcoef.ori + S.logkcoef %*% u.logkcoef
      kcoef.pro <- c(exp(logkcoef.pro[1]), logkcoef.pro[-1])


      # likelihood calculation preparation
      z_counts.sum <- sum(z_counts)
      log.z.sum <- sum(log(z_list |> unlist()))
      XAllMatrix <- Xmatrix_lik(z_list, iu, mc_cores = mc_cores)

      loglik.kcoef.ori <- event.times.loglik(k = kcoef.ori[1], coef = kcoef.ori[-1], uis = uis,
                                             z_counts.sum = z_counts.sum, log.z.sum = log.z.sum,
                                             XAllMatrix = XAllMatrix, iu = iu, mc_cores = mc_cores) +
        prior_dist_k(kcoef.ori[1]) +
        prior_dist_b(exp(kcoef.ori[2])) +
        prior_dist_coef(kcoef.ori[-c(1,2)])



      loglik.kcoef.pro <- event.times.loglik(k = kcoef.pro[1], coef = kcoef.pro[-1], uis = uis,
                                             z_counts.sum = z_counts.sum, log.z.sum = log.z.sum,
                                             XAllMatrix = XAllMatrix, iu = iu, mc_cores = mc_cores) +
        prior_dist_k(kcoef.pro[1]) +
        prior_dist_b(exp(kcoef.pro[2])) +
        prior_dist_coef(kcoef.pro[-c(1,2)])

      acceptance.prob.kcoef <- min(1, exp((loglik.kcoef.pro ) - (loglik.kcoef.ori )))

      if(rbinom(1,1,prob = acceptance.prob.kcoef)) {
        k <- kcoef.pro[1]
        coef <- kcoef.pro[-1]
        accept["kcoef"] <- accept["kcoef"] + 1
      }

      if(iter <= n_burn_in){
        S.logkcoef <- ramcmc::adapt_S(S.logkcoef, u.logkcoef, acceptance.prob.kcoef, iter - 1)
      }

      ###  record the traces
      trace[iter,] <- c(log(k), log(uis_var), coef)
      if(run_and_save){
        if(iter %% iter_per_save == 0){
          print(iter)
          save(list = ls(all.names = TRUE),
               file = paste(save_prefix, "_target_", n_iter, "seed_", seed, ".Rda", sep = ""))
        }
      }

    }
    not_finish <- FALSE
  },
  error=function(e){print(e)},
  finally = {      if(not_finish & run_and_save){
      print(iter)
      save(list = ls(all.names = TRUE),
           file = paste(save_prefix, "_target_", n_iter, "seed_", seed, ".Rda", sep = ""))
  }})

  if(inherits(we_catch, c("error","warning"))){
    if(keep_going){
      return(bicre_continue(paste(save_prefix, "_target_", n_iter, "seed_", seed, ".Rda", sep = ""),
             seed, keep_going = keep_going))
    }else{
      stop(we_catch)
    }
  }else{
    print("target number iterations reached")
    return(list(trace = trace,
                  accept_rate = round(accept/n_iter, 2)))
  }



}



#' @title Continue Bayesian imputation for censored recurrent events with weibull baseline hazard
#'
#' @param stop_file_dir The file path for the results from last unfinished run of the Bayesian data augmentation method.
#' @param seed The seed to set for running the MCMC iteration.
#' @param mc_cores The number of cores to run on parallel for the algorithm. Default to 1. mc_cores > 1 is only avaliable for R on Linux or macOS.
#' @param keep_going Logical. When there is errors, whether repeat the current iteration and continue use this function \link{bicre_continue}.
#' @return An list with two elements.
#'
#' \describe{
#' \item{\code{trace}}{A dataframe that lists the posterior samples of all parameters;}
#' \item{\code{accept_rate}}{A named vector lists the acceptance rates of k,b and \eqn{\boldsymbol{\beta}} and \eqn{\sigma^2};}
#' }
#' @export
bicre_continue <- function(stop_file_dir, seed, mc_cores = 1L, keep_going = FALSE){
  load(stop_file_dir, envir = environment())
  iter_current <- iter
  set.seed(seed)

  we_catch <- tryCatch(expr = {
    not_finish <- TRUE
    for (iter in iter_current : n_iter) {

      if(run_and_print){
        if(iter %% iter_per_print == 0){
          print(iter)
        }
      }




      # imputation


      if(iter == 2){
        z_list <- list()

        for (i in 1:length(iu)) {
          t <- 0
          repeat{
            t <- t + 1
            event_times_i <-impute_single_id_initialize(co_events_frame_single = iu[[i]],
                                                        b = uis[i], coef = coef, expect_cum_FUN = expect_cum_weibull_tvc_Rcpp,
                                                        expect_cum_inverse_FUN = expect_cum_weibull_tvc_inverse_Rcpp,
                                                        k = k)
            if(class(event_times_i) ==  "numeric"){
              z_list[[i]] <- event_times_i
              break
            }else if(event_times_i == "Sequential sampling: larger ui needed"){
              uis[i] <- uis[i] + 1
            }else if(event_times_i == "Sequential sampling: smaller ui needed"){
              uis[i] <- uis[i] - 1
            }else{
              stop("imputation initialize has unexpected errors")
            }

            if(t > 15){
              stop("imputation initialize failed to find a proper starting value")
            }

          }

        }

      }else{

        z_list <- mcmapply(impute_single_id, co_events_frame_single = iu, b = uis,
                           MoreArgs = list(coef = coef, expect_cum_FUN = expect_cum_weibull_tvc_Rcpp,
                                           expect_cum_inverse_FUN = expect_cum_weibull_tvc_inverse_Rcpp,
                                           k = k, fail_mode = fail_mode),
                           mc.cores = mc_cores)
      }


      z_counts <- z_list |> sapply(length)

      # uis update
      shapes <- z_counts + 1/uis_var
      scales <- uis_var /
        (1 + uis_var * (
          mcmapply(function(i){
            sum( - expect_cum_weibull_tvc_Rcpp(attributes(iu[[i]]$ev)$covered_ints$t_start,
                                               lin_pred = iu[[i]]$X %*% coef,
                                               t_breaks = iu[[i]]$cov_t$t_end,
                                               k = k) +
                   expect_cum_weibull_tvc_Rcpp(attributes(iu[[i]]$ev)$covered_ints$t_end,
                                               lin_pred = iu[[i]]$X %*% coef,
                                               t_breaks = iu[[i]]$cov_t$t_end,
                                               k = k))
          },1:length(iu), mc.cores = mc_cores)
        )
        )

      if (truncate_uis) {
        uis <- ltrgamma(length(z_counts), shapes = shapes, scales = scales,
                        truncate = truncate_uis)
      } else {
        uis <- rgamma(length(z_counts), shape = shapes, scale = scales)
      }


      # uis_var update
      u_uis_var <- rnorm(1,mean = 0, sd = 1)
      uis_var_ori <- uis_var
      uis_var_pro <- (log(uis_var_ori) + S_uis_var * u_uis_var) |> exp()

      log_lik_uis_var_ori <- log_lik_uis_var(uis, uis_var_ori) +
        prior_dist_uis_var(uis_var_ori)
      log_lik_uis_var_pro <- log_lik_uis_var(uis, uis_var_pro) +
        prior_dist_uis_var(uis_var_pro)

      acceptance_prob_uis_var <- min(1, exp((log_lik_uis_var_pro) - (log_lik_uis_var_ori)))

      if (rbinom(1, 1, prob = acceptance_prob_uis_var)) {
        uis_var <- uis_var_pro
        accept["uis_var"] <- accept["uis_var"] + 1
      }

      if(iter <= n_burn_in){
        S_uis_var <-  ramcmc::adapt_S(S_uis_var, u_uis_var, acceptance_prob_uis_var, iter - 1) |>  as.vector()
      }


      # k and coef update
      kcoef.ori <- c(k, coef)
      logkcoef.ori <- c(log(k), coef)
      u.logkcoef <- rnorm(length(logkcoef.ori), 0, sd = 1)
      logkcoef.pro <- logkcoef.ori + S.logkcoef %*% u.logkcoef
      kcoef.pro <- c(exp(logkcoef.pro[1]), logkcoef.pro[-1])


      # likelihood calculation preparation
      z_counts.sum <- sum(z_counts)
      log.z.sum <- sum(log(z_list |> unlist()))
      XAllMatrix <- Xmatrix_lik(z_list, iu, mc_cores = mc_cores)

      loglik.kcoef.ori <- event.times.loglik(k = kcoef.ori[1], coef = kcoef.ori[-1], uis = uis,
                                             z_counts.sum = z_counts.sum, log.z.sum = log.z.sum,
                                             XAllMatrix = XAllMatrix, iu = iu, mc_cores = mc_cores) +
        prior_dist_k(kcoef.ori[1]) +
        prior_dist_b(exp(kcoef.ori[2])) +
        prior_dist_coef(kcoef.ori[-c(1,2)])



      loglik.kcoef.pro <- event.times.loglik(k = kcoef.pro[1], coef = kcoef.pro[-1], uis = uis,
                                             z_counts.sum = z_counts.sum, log.z.sum = log.z.sum,
                                             XAllMatrix = XAllMatrix, iu = iu, mc_cores = mc_cores) +
        prior_dist_k(kcoef.pro[1]) +
        prior_dist_b(exp(kcoef.pro[2])) +
        prior_dist_coef(kcoef.pro[-c(1,2)])

      acceptance.prob.kcoef <- min(1, exp((loglik.kcoef.pro ) - (loglik.kcoef.ori )))

      if(rbinom(1,1,prob = acceptance.prob.kcoef)) {
        k <- kcoef.pro[1]
        coef <- kcoef.pro[-1]
        accept["kcoef"] <- accept["kcoef"] + 1
      }

      if(iter <= n_burn_in){
        S.logkcoef <- ramcmc::adapt_S(S.logkcoef, u.logkcoef, acceptance.prob.kcoef, iter - 1)
      }

      ###  record the traces
      trace[iter,] <- c(log(k), log(uis_var), coef)
      if(run_and_save){
        if(iter %% iter_per_save == 0){
          print(iter)
          save(list = ls(all.names = TRUE),
               file = paste(save_prefix, "_target_", n_iter, "seed_", seed, ".Rda", sep = ""))
        }
      }

    }
    not_finish <- FALSE
  },
  # warning=function(w){print(w)},
  error=function(e){print(e)},
  finally = {if(not_finish & run_and_save){
    print(iter)
    save(list = ls(all.names = TRUE),
         file = paste(save_prefix, "_target_", n_iter, "seed_", seed, ".Rda", sep = ""))
  }})

  if(inherits(we_catch, c("error","warning"))){
    if(keep_going){
      return(bicre_continue(paste(save_prefix, "_target_", n_iter, "seed_", seed, ".Rda", sep = ""),
                                    seed, keep_going = keep_going))
    }else{
      stop(we_catch)
    }
  }else{
    print("target number iterations reached")
    return(list(trace = trace,
                accept_rate = round(accept/n_iter, 2)))
  }
}
