test_that("expect_cum_tvc calls expect_cum_weibull correctly", {
  expect_equal(expect_cum_weibull_tvc(4, shape = 2, scale = 3),
               expect_cum_weibull(4, 2, 3))
  expect_equal(expect_cum_weibull_tvc(4, k = 2, b = 3),
               expect_cum_weibull(4, k = 2, b = 3))
})

test_that("expect_cum_tvc telescopes correctly", {
  lin_pred <- 1 : 5
  t_breaks <- 1 : 5

  expect_equal(
    expect_cum_weibull_tvc(3, lin_pred = lin_pred, t_breaks = t_breaks, k = 1, b = 1),
    expect_cum_weibull(1, k = 1, b = 1) * exp(lin_pred[1]) +
      (expect_cum_weibull(2, k = 1, b = 1) * exp(lin_pred[2]) - expect_cum_weibull(1, k = 1, b = 1) * exp(lin_pred[2])) +
      (expect_cum_weibull(3, k = 1, b = 1) * exp(lin_pred[3]) - expect_cum_weibull(2, k = 1, b = 1) * exp(lin_pred[3]))
  )

  expect_equal(expect_cum_weibull_tvc(0), 0)
})

test_that("expect_cum_tvc validates input", {
  # no t_breaks >= t
  expect_error(expect_cum_weibull_tvc(6, lin_pred = 1 : 5, t_breaks = 1 : 5))

  # lin_pred and t_breaks not of same length
  expect_error(expect_cum_weibull_tvc(3, lin_pred = 1 : 4, t_breaks = 1 : 5))
})

test_that("expect_cum_weibull_tvc vectorized correctly", {
  lin_pred <- 1 : 5
  t_breaks <- 1 : 5

  expect_equal(
    expect_cum_weibull_tvc(c(1.5, 4), lin_pred = lin_pred, t_breaks = t_breaks, k = 2, b = 3),
    c(
      expect_cum_weibull_tvc(1.5, lin_pred = lin_pred, t_breaks = t_breaks, k = 2, b = 3),
      expect_cum_weibull_tvc(4, lin_pred = lin_pred, t_breaks = t_breaks, k = 2, b = 3)
    )
  )
})

test_that("expect_cum_weibull_tvc_inverse inverts correctly", {
  lin_pred <- 1 : 5
  t_breaks <- 1 : 5

  expect_equal(
    expect_cum_weibull_tvc_inverse(
      expect_cum_weibull_tvc(2, lin_pred = lin_pred, t_breaks = t_breaks, k = 2, b = 3),
      lin_pred = lin_pred, t_breaks = t_breaks, k = 2, b = 3),
    2)

  expect_equal(
    expect_cum_weibull_tvc_inverse(
      expect_cum_weibull_tvc(1.5, lin_pred = lin_pred, t_breaks = t_breaks, k = 2, b = 3),
      lin_pred = lin_pred, t_breaks = t_breaks, k = 2, b = 3),
    1.5)

  expect_equal(
    expect_cum_weibull_tvc_inverse(
      expect_cum_weibull_tvc(0.5, lin_pred = lin_pred, t_breaks = t_breaks, k = 2, b = 3),
      lin_pred = lin_pred, t_breaks = t_breaks, k = 2, b = 3),
    0.5)

  expect_equal(
    expect_cum_weibull_tvc_inverse(
      expect_cum_weibull_tvc(0, lin_pred = lin_pred, t_breaks = t_breaks, k = 2, b = 3),
      lin_pred = lin_pred, t_breaks = t_breaks, k = 2, b = 3),
    0)
})

test_that("expect_cum_weibull_inverse inverts correctly", {
  expect_equal(
    expect_cum_weibull_inverse(expect_cum_weibull(4, 2, 3), 2, 3),
    4
  )

  expect_equal(
    expect_cum_weibull_inverse(0),
    0
  )

  expect_equal(
    expect_cum_weibull_inverse(Inf),
    Inf
  )
})

test_that("expect_cum_weibull_inverse vectorized correctly", {
  expect_equal(
    expect_cum_weibull_inverse(c(30, 800), k = 2, b = 3),
    c(
      expect_cum_weibull_inverse(c(30), k = 2, b = 3),
      expect_cum_weibull_inverse(c(800), k = 2, b = 3)
    )
  )
})

test_that("simulate_nonhomog_inversion simulates correctly", {
  set.seed(0)

  z <- simulate_nonhomog_inversion(
    t_start = 0, t_end = 100,
    expect_cum_FUN = expect_cum_weibull_tvc,
    expect_cum_inverse_FUN = expect_cum_weibull_tvc_inverse,
    k = 2, b = 1,
    t_breaks = c(50, 100),
    lin_pred = c(log(10), log(2))
  )

  expect_equal(
    expect_cum_weibull_tvc(50,
                           lin_pred = c(log(10), log(2)),
                           t_breaks = c(50, 100),
                           k = 2, b = 1),
    sum(z < 50),
    tolerance = 0.01
  )

  expect_equal(
    expect_cum_weibull_tvc(50,
                           lin_pred = c(log(10), log(2)),
                           t_breaks = c(50, 100),
                           k = 2, b = 1),
    sum(z < 50),
    tolerance = 0.01
  )

})

test_that("simulate_nonhomog_inversion respects bounds", {
  expect_lte(max(
    simulate_nonhomog_inversion(
      t_start = 1, t_end = 10,
      expect_cum_FUN = expect_cum_weibull,
      expect_cum_inverse_FUN = expect_cum_weibull_inverse,
      k = 2, b = 1,
      count_min = 0, count_max = 10)
    ), 10)

  expect_gte(min(
    simulate_nonhomog_inversion(
      t_start = 1, t_end = 10,
      expect_cum_FUN = expect_cum_weibull,
      expect_cum_inverse_FUN = expect_cum_weibull_inverse,
      k = 2, b = 1,
      count_min = 0, count_max = 10)
  ), 1)

  expect_gte(min(
    simulate_nonhomog_inversion(
      t_start = 1, t_end = 10,
      expect_cum_FUN = expect_cum_weibull,
      expect_cum_inverse_FUN = expect_cum_weibull_inverse,
      k = 2, b = 1,
      count_min = 0, count_max = Inf)
  ), 1)
})

test_that("rpois_trunc draws from correct un-truncated distribution", {
  set.seed(0)

  x_test <- rpois_trunc(10000, lambda = 1)
  x_ref <- rpois(10000, lambda = 1)

  expect_equal(mean(x_test == 4), mean(x_ref == 4), tolerance = 0.05)
})

test_that("rpois_trunc draws from correct left-truncated distribution", {
  set.seed(0)

  x_test <- rpois_trunc(10000, lambda = 1, min = 1)
  x_ref <- rpois(10000, lambda = 1)
  x_ref <- x_ref[x_ref >= 1]

  expect_equal(sum(x_test == 0), 0)
  expect_equal(mean(x_test == 4), mean(x_ref == 4), tolerance = 0.05)
})

test_that("rpois_trunc draws from correct right-truncated distribution", {
  set.seed(0)

  x_test <- rpois_trunc(10000, lambda = 1, max = 10)
  x_ref <- rpois(10000, lambda = 1)
  x_ref <- x_ref[x_ref <= 10]

  expect_equal(sum(x_test > 10), 0)
  expect_equal(mean(x_test == 4), mean(x_ref == 4), tolerance = 0.05)
})

test_that("rpois_trunc draws from correct double-truncated distribution", {
  set.seed(0)

  x_test <- rpois_trunc(10000, lambda = 1, min = 1, max = 10)
  x_ref <- rpois(10000, lambda = 1)
  x_ref <- x_ref[x_ref >= 1 & x_ref <= 10]

  expect_equal(sum(x_test < 1 | x_test > 10), 0)
  expect_equal(mean(x_test == 4), mean(x_ref == 4), tolerance = 0.05)
})

test_that("rpois_trunc handles edge cases correctly", {
  set.seed(0)

  x <- rpois_trunc(n = 10, lambda = 0)
  expect_true(all(x == 0))
  expect_length(x, 10)

  expect_error(rpois_trunc(n = 10, lambda = 0, min = 1))
})
