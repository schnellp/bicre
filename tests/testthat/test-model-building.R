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
  id = as.character(c(1, 1, 1, 1,1, 1,
                      1, 1, 1, 1,
                      1, 1, 1,
                      1, 1,
                      2)),
  t_start = c(19, 20, 21, 29, 32, 33,
              40, 45, 50, 55,
              65, 69, 67,
              76, 78,
              0),
  t_end = c(21, 22, 35, 30,36, 34,
            46, 52, 57, 60,
            75, 72, 73,
            85, 82,
            4),
  e_min = c(0, 0, 1, 0,0, 0,
            1, 1, 1, 0,
            1, 0, 2,
            2, 3,
            1),
  e_max = c(0, 0, Inf, 0,0,5,
            Inf, Inf, Inf, 0,
            5, 6, 7,
            Inf, 5,
            Inf)
)

ce <- co_events(data_merged_benchmark, data_test_events,
                id, t_start, t_end, e_min, e_max,
                fill = list("A" = FALSE,
                            "B" = FALSE))

cef <- co_events_frame(ce, ~ A + B)

ev <- cef[[1]]$ev

ev_simp <- simplify_ev(ev)

iu <- build_imputation_units(cef)


test_that("simplify steps conducted correctly",{

# test step 1
    expect_equal(
      as.numeric(ev_simp[1,]),
      as.numeric(data.frame(t_start = 19,
               t_end = 22,
               e_min = 0,
               e_max = 0))
    )

# test step 2
  expect_equal(
    as.numeric(ev_simp[2,]),
    as.numeric(data.frame(t_start = ev[2,"t_end"],
                          t_end = ev[5,"t_start"],
                          e_min = 1,
                          e_max = Inf))

  )

# test step 3
  expect_equal(
    nrow(ev_simp %>% filter(t_start >= 33 & t_end <= 34 & e_max > 0)),
    0
  )

# test step 4
  expect_equal(
    nrow(ev_simp %>% filter(t_start ==  69 & t_end  ==  72)),
    0
  )
  expect_equal(
    nrow(ev_simp %>% filter(t_start == 76 & t_end == 85 )),
    0
  )

# test step 5
  expect_equal(
    (ev_simp %>% filter(t_start == 65 & t_end == 75 ))$e_min,
    (ev %>% filter(t_start == 67 & t_end == 73 ))$e_min
  )
  expect_equal(
    (ev_simp %>% filter(t_start == 67 & t_end == 73 ))$e_max,
    (ev %>% filter(t_start == 65 & t_end == 75 ))$e_max
  )

# test no change
  expect_equal(
    simplify_ev(ev_simp),
    ev_simp
  )
  })

test_that("partition conducted correctly",{
  expect_equal(
    partition_ev(ev_simp)$unit,
    c(1,2,2,3,4,4,4,5,6,6,7)
  )
})


minmax  <-   c(0,2, 2, 6, 3, 10, 0, 0, 0, 2, 1, 5, 2, 7, 1, 3, 2, 5) %>% matrix(ncol = 2, byrow = T)

data_test_events <- data.frame(
  id = as.character(rep(1,9)),
  t_start = c(0, 1,2,3,4,8,9,10,13),
  t_end = c(0.5, 6,14,5,7,11,15,12,16),
  e_min = minmax[,1],
  e_max = minmax[,2]
)
data_test_events <- rbind(data_test_events, c(2, 0, 3, 0, 1))
ce <- co_events(data_merged_benchmark, data_test_events,
                id, t_start, t_end, e_min, e_max,
                fill = list("A" = FALSE,
                            "B" = FALSE))

cef <- co_events_frame(ce, ~ A + B)
unit_prepared1 <- cef[[1]]$ev %>% build_imputation_units_single





minmax  <-   c(0,5,2,8,1,4,0,0, 1, 3, 2,5) %>% matrix(ncol = 2, byrow = T)

data_test_events <- data.frame(
  id = as.character(rep(1,6)),
  t_start = c(1,2,3,5,7,8),
  t_end = c(4,9,5,6,10,11),
  e_min = minmax[,1],
  e_max = minmax[,2]
)
data_test_events <- rbind(data_test_events, c(2, 0, 3, 0, 1))

ce <- co_events(data_merged_benchmark, data_test_events,
                id, t_start, t_end, e_min, e_max,
                fill = list("A" = FALSE,
                            "B" = FALSE))

cef <- co_events_frame(ce, ~ A + B)
unit_prepared2 <- cef[[1]]$ev %>% build_imputation_units_single

test_that("imputation units prepared correctly",{
  expect_equal(unit_prepared1$impute_unit_compund[[1]]$roundone[[1]]$outer_check$t_start,
               c(2,5))
  expect_equal(unit_prepared1$impute_unit_compund[[1]]$roundone[[2]]$outer_check$t_start,
               c(2,9,10))
  expect_equal(unit_prepared1$impute_unit_compund[[1]]$roundtwo[[2]]$overlap_strict[,1:4] %>% as.numeric(),
               data_test_events[8,-1] %>% as.numeric())
  expect_equal(unit_prepared1$impute_unit_compund[[1]]$outer_check_all$t_start,
               c(2,5,9,10))

  expect_equal(unit_prepared2$impute_unit_compund[[1]]$roundone[[2]] %>% length,
               1)
  expect_equal(unit_prepared2$impute_unit_compund[[1]]$roundtwo[[3]]$leftover %>% as.numeric(),
               c(10,11))


})

iu <- build_imputation_units(cef)
test_that("imputation units built correctly",{
  expect_equal(iu[[1]]$ev,
               unit_prepared2)
})
