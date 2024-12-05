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

data_test_events <- data.frame(
  id = as.character(c(1, 1, 1,
                      1, 1, 1, 1,
                      2)),
  t_start = c(20, 22, 29,
              40, 45, 50, 55,
              0),
  t_end = c(22, 32, 30,
            46, 52, 57, 60,
            4),
  e_min = c(0, 1, 0,
            1, 1, 1, 0,
            1),
  e_max = c(0, Inf, 0,
            Inf, Inf, Inf, 0,
            Inf)
)

ce <- co_events(data_merged_benchmark, data_test_events,
                id, t_start, t_end, e_min, e_max,
                fill = list("A" = FALSE,
                            "B" = FALSE))

cef <- co_events_frame(ce, ~ A + B)

iu <- build_imputation_units(cef)






test_that("impute_single satisfies constraints", {
  z <- impute_single_id(iu[[1]],
                        coef = c(0, 1, 0),
                        expect_cum_FUN = expect_cum_weibull_tvc,
                        expect_cum_inverse_FUN = expect_cum_weibull_tvc_inverse)


  counts <- sapply(1 : nrow(iu[[1]]$events),
                   function(i) {
                     sum(iu[[1]]$events$t_start[i] < z &
                           z <= iu[[1]]$events$t_end[i])
                   })

  expect_true(all(counts >= iu[[1]]$events$e_min &
                    counts <= iu[[1]]$events$e_max))
})


############# Test 2
############# each T-F disjoint no limit interval (see the definition in our candidacy paper) has been imputed only once?

data_test_events2 <- data.frame(
  id = as.character(c(1, 1, 1, 1, 1,
                      2)),
  t_start = c(1, 2, 3, 4, 6,
              0),
  t_end = c(5, 9, 8, 7, 10,
            4),
  e_min = c(2, 5, 3, 2, 4,
            1),
  e_max = c(5, 7, 6, 4, 6,
            Inf)
)

ce2 <- co_events(data_merged_benchmark, data_test_events2,
                 id, t_start, t_end, e_min, e_max,
                 fill = list("A" = FALSE,
                             "B" = FALSE))

cef2 <- co_events_frame(ce2, ~ A + B)

iu2 <- build_imputation_units(cef2)

coef = c(0, 1, 0)
expect_cum_FUN = expect_cum_weibull_tvc
expect_cum_inverse_FUN = expect_cum_weibull_tvc_inverse
imputation_unit <- iu2[[1]]$ev[[1]]


