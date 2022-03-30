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
