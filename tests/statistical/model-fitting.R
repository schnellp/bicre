set.seed(0)

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

for (i in 1 : N) {
  sel <- which(data_gen$id == i)

  z <- simulate_nonhomog_inversion(
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
}

data_events <- do.call("rbind", data_events)

ce <- co_events(data_cov, data_events,
                id, t_start, t_end, e_min, e_max,
                fill = list("trt" = FALSE))

trace <- bicre( ~ trt + sex, ce, n_burn = 1000, n_keep = 10000)

