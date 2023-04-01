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

minmax  <-   c(2, 5, 1, 3, 0, 2, 2, 6, 3, 10, 0, 0, 0, 2, 1, 5, 2, 7, 1, 3, 2, 5) %>% matrix(ncol = 2, byrow = T)

data_test_events <- data.frame(
  id = as.character(rep(1,11)),
  t_start = (c(-10, -3, 0, 1,2,3,4,8,9,10,13) + 10),
  t_end = (c(-5, -1, 0.5, 6,14,5,7,11,15,12,16) + 10),
  e_min = minmax[,1],
  e_max = minmax[,2]
)
data_test_events <- rbind(data_test_events, c(2, 0, 3, 0, 1))
ce <- co_events(data_merged_benchmark, data_test_events,
                id, t_start, t_end, e_min, e_max,
                fill = list("A" = FALSE,
                            "B" = FALSE))

cef <- co_events_frame(ce, ~ A + B)
iu2 <- build_imputation_units(cef)

roundone_unit <- iu2[[1]]$ev$impute_unit_compund[[1]]$roundone[[1]]

test_that("round one imputation satisfy constraints", {
  z_roundone <- roundone_impute(roundone_unit,
                  expect_cum_FUN = expect_cum_weibull_tvc,
                  expect_cum_inverse_FUN = expect_cum_weibull_tvc_inverse)

  counts_inner <- sapply(1:nrow(roundone_unit$inner_check),
                         function(i) {sum(roundone_unit$inner_check$t_start[i] < z_roundone &
                               z_roundone <= roundone_unit$inner_check$t_end[i])}
                         )
  expect_true(all(counts_inner >= roundone_unit$inner_check$e_min &
                    counts_inner <= roundone_unit$inner_check$e_max))

  counts_outer <- sapply(1:nrow(roundone_unit$outer_check),
                         function(i) {sum(roundone_unit$outer_check$t_start[i] < z_roundone &
                                            z_roundone <= roundone_unit$outer_check$t_end[i])}
  )

  expect_true(all(counts_outer <= roundone_unit$outer_check$e_max))
})

roundtwo_unit <- iu2[[1]]$ev$impute_unit_compund[[1]]$roundtwo[[1]]

test_that("round two imputation satisfy constraints",{
  z_roundtwo <- roundtwo_impute(roundtwo_unit,
                                expect_cum_FUN = expect_cum_weibull_tvc,
                                expect_cum_inverse_FUN = expect_cum_weibull_tvc_inverse)

  counts_covered <- sapply(1:nrow(roundtwo_unit$covered),
                           function(i) {sum(roundtwo_unit$covered$t_start[i] < z_roundtwo &
                                              z_roundtwo <= roundtwo_unit$covered$t_end[i])}
  )
  expect_true(all(counts_covered <= roundtwo_unit$covered$e_max))

  counts_overlap <- sapply(1:nrow(roundtwo_unit$overlap_strict),
                           function(i) {sum(roundtwo_unit$overlap_strict$t_start[i] < z_roundtwo &
                                              z_roundtwo <= roundtwo_unit$overlap_strict$t_end[i])}
  )

  expect_true(all(counts_overlap <= roundtwo_unit$overlap_strict$e_max))

})


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


  z2 <-  impute_single_id(iu2[[1]],
                          coef = c(0, 1, 0),
                          expect_cum_FUN = expect_cum_weibull_tvc,
                          expect_cum_inverse_FUN = expect_cum_weibull_tvc_inverse)

  counts2 <- sapply(1 : nrow(iu2[[1]]$events),
                   function(i) {
                     sum(iu2[[1]]$events$t_start[i] < z2 &
                           z2 <= iu2[[1]]$events$t_end[i])
                   })

  expect_true(all(counts2 >= iu2[[1]]$events$e_min &
                    counts2 <= iu2[[1]]$events$e_max))

})

