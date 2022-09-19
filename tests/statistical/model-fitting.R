library(dplyr)

M <- 100
est <- upper <- lower <- matrix(NA, nrow = M, ncol = 3)
for (m in 5 : M) {
  set.seed(m)
  print(paste(Sys.time(), m))

  beta <- c(-2, 1, 0)

  N <- 100
  max_time <- 10
  lag <- 1
  log_observation_rate <- 2.4
  treatment_length <- 5


  treatment_start <- runif(N, 0, max_time)
  data_treatment <- data.frame(
    id = 1 : N,
    trt = TRUE,
    t_start = treatment_start,
    t_end = pmin(treatment_start + treatment_length, max_time)
  )

  data_subjects <- data.frame(
    id = 1 : N,
    sex = sample(c("female", "male"), N, replace = TRUE)
  )

  data_treatment_complete <- data_treatment %>%
    complete_interval_data(id = id, t_start = t_start, t_end = t_end,
                           fill = list(trt = FALSE), new_nodes = c(0, max_time)) %>%
    mutate(id = as.numeric(id))

  data_cov <- data_subjects %>%
    left_join(data_treatment_complete, by = "id")

  mf <- model.frame( ~ trt + sex, data = data_cov)
  X <- model.matrix(mf, data = data_cov)

  lin_pred <- drop(X %*% beta)

  data_gen <- cbind(
    data_treatment_complete,
    lin_pred
  )

  data_events <- list()
  data_events_full <- list()

  for (i in 1 : N) {
    sel <- which(data_gen$id == i)

    z <- bicre::simulate_nonhomog_inversion(
      t_start = 0,
      t_end = max_time,
      expect_cum_FUN = expect_cum_weibull_tvc,
      expect_cum_inverse_FUN = expect_cum_weibull_tvc_inverse,
      lin_pred = data_gen$lin_pred[sel],
      t_breaks = data_gen$t_end[sel])

    n_obs <- max(rpois(1, exp(log_observation_rate)), 1)
    t_obs <- sort(runif(n_obs, 0, max_time))
    # n_obs <- 10
    # t_obs <- 1 : 10
    y <- logical(n_obs)
    for (j in 1 : n_obs) {
      y[j] <- any(t_obs[j] - lag < z & z < t_obs[j])
    }

    data_events[[i]] <- data.frame(id = i,
                                   t_start = pmax(t_obs - lag, 0),
                                   t_end = t_obs,
                                   e_min = if_else(y, 1, 0),
                                   e_max = if_else(y, Inf, 0))

    if (length(z) > 0) {
      v <- logical(length(z))
      for (k in 1 : length(z)) {
        v[k] <- any(t_obs - lag < z[k] & z[k] < t_obs)
      }

      data_events_full[[i]] <- z[v]
    } else {
      data_events_full[[i]] <- numeric(0)
    }

  }

  data_events <- do.call("rbind", data_events)

  ce <- co_events(data_cov, data_events,
                  id, t_start, t_end, e_min, e_max,
                  fill = list("trt" = FALSE))

  n_burn <- 1000
  n_keep <- 10000

  # trace <- bicre( ~ trt + sex, ce, n_burn = n_burn, n_keep = n_keep,
  #                 z_force = data_events_full)

  trace <- bicre( ~ trt + sex, ce, n_burn = n_burn, n_keep = n_keep,
                  z_force = NULL)

  est[m, ] <- colMeans(trace[n_burn + (1 : n_keep), ])
  lower[m, ] <- apply(trace[n_burn + (1 : n_keep), ], 2, quantile, prob = c(0.025))
  upper[m, ] <- apply(trace[n_burn + (1 : n_keep), ], 2, quantile, prob = c(0.975))

  plot(trace[, "(Intercept)"])
}

colnames(est) <- colnames(lower) <- colnames(upper) <- colnames(trace)

hist(est[, "(Intercept)"])
hist(est[, "trtTRUE"])
hist(est[, "sex"])

for (p in 1 : ncol(est)) {
  plot(est[, p], ylim = range(lower[, p], upper[, p], na.rm = TRUE),
       main = colnames(est)[p])
  segments(1 : M, lower[, p], 1 : M, upper[, p],
           col = ifelse(lower[, p] < beta[p] & upper[, p] > beta[p], "black", "red"))
  abline(h = beta[p], lty = 2)
}

