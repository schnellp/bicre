library(dplyr)

data_overlap_ints <- data.frame(
  id = c(1, 1, 1, 1,
         2, 2, 2, 2),
  t_start = c(0, 4, 7, 21,
              7, 14, 21, 28),
  t_end = c(7, 11, 14, 28,
            14, 21, 28, 35),
  med = c("A", "A", "B", "B",
          "A", "A", "A", "A")
)

data_collapsed_ints_benchmark <- data.frame(
  id = c(1, 1, 1,
         2),
  t_start = c(0, 7, 21,
              7),
  t_end = c(11, 14, 28,
            35),
  med = c("A", "B", "B",
          "A")
)

data_collapsed_ints <- data_overlap_ints %>%
  collapse_interval_data("t_start", "t_end")

test_that("collapse_interval_data produces correct intervals", {
  expect_equal(data_collapsed_ints_benchmark$id,
               data_collapsed_ints$id)

  expect_equal(data_collapsed_ints_benchmark$t_start,
               data_collapsed_ints$t_start)

  expect_equal(data_collapsed_ints_benchmark$t_end,
               data_collapsed_ints$t_end)
})

test_that("collapse_interval_data keeps data intact", {
  expect_equal(data_collapsed_ints_benchmark$med,
               data_collapsed_ints$med)
})

