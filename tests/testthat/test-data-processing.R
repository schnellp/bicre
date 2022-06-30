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
  collapse_interval_data(t_start, t_end)

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

test_that("collapse_interval_data is invariant to input data order", {
  set.seed(0)
  expect_equal(data_overlap_ints %>%
                 collapse_interval_data(t_start, t_end),
               data_overlap_ints %>%
                 slice_sample(prop = 1) %>%
                 collapse_interval_data(t_start, t_end))
})

test_that("collapse_interval_data is idempotent", {
  expect_equal(data_overlap_ints %>%
                 collapse_interval_data(t_start, t_end),
               data_overlap_ints %>%
                 collapse_interval_data(t_start, t_end) %>%
                 collapse_interval_data(t_start, t_end))
})

##########################
### complete_intervals ###
##########################

data_long <- data.frame(
  id = c(1, 1, 1,
         2),
  t_start = c(0, 7, 21,
              7),
  t_end = c(11, 14, 28,
            35),
  med = c("A", "B", "B",
          "A")
)

data_a <- data_long %>%
  filter(med == "A") %>%
  mutate(A = (med == "A")) %>%
  select(-med) %>%
  collapse_interval_data(t_start_var = t_start, t_end_var = t_end)

test_that("complete_interval_data is idempotent", {
  expect_equal(
    data_a %>%
      complete_interval_data(id = id, t_start_var = t_start, t_end_var = t_end,
                             fill = list(A = FALSE), new_nodes = 42),
    data_a %>%
      complete_interval_data(id = id, t_start_var = t_start, t_end_var = t_end,
                             fill = list(A = FALSE), new_nodes = 42) %>%
      complete_interval_data(id = id, t_start_var = t_start, t_end_var = t_end,
                             fill = list(A = FALSE), new_nodes = 42)
  )
})

test_that("complete_interval_data retains column ordering", {
  expect_equal(
    data_a %>%
      complete_interval_data(id = id, t_start_var = t_start, t_end_var = t_end,
                             fill = list(A = FALSE), new_nodes = 42) %>%
      colnames(),
    data_a %>%
      colnames()
  )
})

###########################
### merge_interval_data ###
###########################

data_merged_benchmark <- data.frame(
  id = as.character(c(1, 1, 1, 1, 1,
                      2)),
  t_start = c(0, 7, 11, 14, 21,
              7),
  t_end = c(7, 11, 14, 21, 28,
            35),
  A = c(TRUE, TRUE, FALSE, FALSE, FALSE,
        TRUE),
  B = c(FALSE, TRUE, TRUE, FALSE, TRUE,
        FALSE)
)

test_that("merge_interval_data handles an id not appearing in new_data", {
  data_a <- data_long %>%
    filter(med == "A") %>%
    mutate(A = (med == "A")) %>%
    select(-med)

  data_b <- data_long %>%
    filter(med == "B") %>%
    mutate(B = (med == "B")) %>%
    select(-med)

  expect_warning(
    data_merged <- data_a %>%
      merge_interval_data(new_data = data_b,
                          id = id, t_start_var = t_start, t_end_var = t_end,
                          new_var = B,
                          fill = list("A" = FALSE,
                                      "B" = FALSE))
  )

  expect_equal(data_merged, data_merged_benchmark)
})

test_that("merge_interval_data warns when an id id appears in new_data but not data", {
  data_a <- data_long %>%
    filter(med == "A") %>%
    mutate(A = (med == "A")) %>%
    select(-med)

  data_b <- data_long %>%
    filter(med == "B") %>%
    mutate(B = (med == "B")) %>%
    select(-med)

  expect_warning(
    data_b %>%
      merge_interval_data(new_data = data_a,
                          id = id, t_start_var = t_start, t_end_var = t_end,
                          new_var = A,
                          fill = list("A" = FALSE,
                                      "B" = FALSE))
  )
})

test_that("merge_interval_data chains correctly", {
  data_long <- data.frame(
    id = "1",
    t_start = c(8,
                3, 6,
                1, 5, 8),
    t_end = c(9,
              5, 8,
              2, 6, 10),
    med = c("A",
            "B", "B",
            "C", "C", "C")
    )

  data_merged_benchmark <- data.frame(
    id = "1",
    t_start = c(1, 2, 3, 5, 6, 8, 9),
    t_end = c(2, 3, 5, 6, 8, 9, 10),
    A = c(FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, FALSE),
    B = c(FALSE, FALSE, TRUE, FALSE, TRUE, FALSE, FALSE),
    C = c(TRUE, FALSE, FALSE, TRUE, FALSE, TRUE, TRUE)
  )

  data_a <- data_long %>%
    filter(med == "A") %>%
    mutate(A = (med == "A")) %>%
    select(-med)

  data_b <- data_long %>%
    filter(med == "B") %>%
    mutate(B = (med == "B")) %>%
    select(-med)

  data_c <- data_long %>%
    filter(med == "C") %>%
    mutate(C = (med == "C")) %>%
    select(-med)

  data_merged <- data_a %>%
    merge_interval_data(data_b, id, t_start, t_end,
                        new_var = B,
                        fill = list("A" = FALSE,
                                    "B" = FALSE)) %>%
    merge_interval_data(data_c, id, t_start, t_end,
                        new_var = C,
                        fill = list("A" = FALSE,
                                    "B" = FALSE,
                                    "C" = FALSE))

  expect_equal(data_merged, data_merged_benchmark)
})
