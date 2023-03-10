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



#######following code are copied and pasted from part of the sample_sequential function #####
set.seed(1)

z2 <- numeric(0)

# independently impute first round
for (i in 1 : nrow(imputation_unit$`TRUE`)) {
  ei <- imputation_unit$`TRUE`[i, ]
  
  zi <- simulate_nonhomog_inversion(
    t_start = ei$t_start,
    t_end = ei$t_end,
    expect_cum_FUN = expect_cum_FUN,
    expect_cum_inverse_FUN = expect_cum_inverse_FUN,
    count_min = ei$e_min,
    count_max = ei$e_max,
    # ...
  )
  
  z2 <- c(z2, zi)
}

# first-round consistency check

if (is.null(imputation_unit$`FALSE`)) {
  return(z2)
}

counts <- sapply(1 : nrow(imputation_unit$`FALSE`),
                 function(i) {
                   sum(imputation_unit$`FALSE`$t_start[i] < z2 &
                         z2 <= imputation_unit$`FALSE`$t_end[i])
                 })

if (any(counts > imputation_unit$`FALSE`$e_max)) {
  if (verbose) {
    print("Failed first-round check:")
    print(z2)
  }
  next
}

# impute remainder in second round
set.seed(2)
for (i in 1 : nrow(imputation_unit$`FALSE`)) {
  ei <- imputation_unit$`FALSE`[i, ]
  
  if (ei$t_start_trunc >= ei$t_end_trunc) {
    next
  }
  
  zi <- simulate_nonhomog_inversion(
    t_start = ei$t_start_trunc,
    t_end = ei$t_end_trunc,
    expect_cum_FUN = expect_cum_FUN,
    expect_cum_inverse_FUN = expect_cum_inverse_FUN,
    count_min = ei$e_min - counts[i],
    count_max = ei$e_max - counts[i],
    # ...
  )
  
  z2 <- c(z2, zi)
}

test_that("each T-F disjoint no limit interval (see the definition in our paper) has been imputed only once",
          {
          num_repeat_sim_T_F_disjoint_no_limit <-  group_size(group_by(imputation_unit$`FALSE`,
                                t_start_trunc,t_end_trunc))
          expect_true(all(num_repeat_sim_T_F_disjoint_no_limit == 1))
          }
          
          
  
)


############### Test 3
############### Noncontig FALSE interval t start end trunc covers all
############### T-F disjoint no limit interval (see the definition in our candidacy paper)?

data_test_events3 <- data.frame(
  id = as.character(c(1, 1, 1, 1, 
                      2)),
  t_start = c(1, 2, 4, 6, 
              0),
  t_end = c(3, 7, 5, 8,
            4),
  e_min = c(2, 3, 1, 1, 
            1),
  e_max = c(4, 5, 3, 3, 
            Inf)
)

ce3 <- co_events(data_merged_benchmark, data_test_events3,
                 id, t_start, t_end, e_min, e_max,
                 fill = list("A" = FALSE,
                             "B" = FALSE))

cef3 <- co_events_frame(ce3, ~ A + B)

iu3 <- build_imputation_units(cef3)


coef = c(0, 1, 0)
expect_cum_FUN = expect_cum_weibull_tvc
expect_cum_inverse_FUN = expect_cum_weibull_tvc_inverse
imputation_unit <- iu3[[1]]$ev[[1]]

#######following code are copied and pasted from part of the sample_sequential function #####
set.seed(1)
z3 <- numeric(0)

# independently impute first round
for (i in 1 : nrow(imputation_unit$`TRUE`)) {
  ei <- imputation_unit$`TRUE`[i, ]
  
  zi <- simulate_nonhomog_inversion(
    t_start = ei$t_start,
    t_end = ei$t_end,
    expect_cum_FUN = expect_cum_FUN,
    expect_cum_inverse_FUN = expect_cum_inverse_FUN,
    count_min = ei$e_min,
    count_max = ei$e_max,
    # ...
  )
  
  z3 <- c(z3, zi)
}

# first-round consistency check

if (is.null(imputation_unit$`FALSE`)) {
  return(z3)
}

counts <- sapply(1 : nrow(imputation_unit$`FALSE`),
                 function(i) {
                   sum(imputation_unit$`FALSE`$t_start[i] < z3 &
                         z3 <= imputation_unit$`FALSE`$t_end[i])
                 })

if (any(counts > imputation_unit$`FALSE`$e_max)) {
  if (verbose) {
    print("Failed first-round check:")
    print(z3)
  }
  next
}

# impute remainder in second round
set.seed(2)
for (i in 1 : nrow(imputation_unit$`FALSE`)) {
  ei <- imputation_unit$`FALSE`[i, ]
  
  if (ei$t_start_trunc >= ei$t_end_trunc) {
    next
  }
  
  zi <- simulate_nonhomog_inversion(
    t_start = ei$t_start_trunc,
    t_end = ei$t_end_trunc,
    expect_cum_FUN = expect_cum_FUN,
    expect_cum_inverse_FUN = expect_cum_inverse_FUN,
    count_min = ei$e_min - counts[i],
    count_max = ei$e_max - counts[i],
    # ...
  )
  
  z3 <- c(z3, zi)
}


test_that("Noncontig FALSE interval t start end trunc covers all T-F disjoint no limit interval (see the definition in our candidacy paper)",
          {
            get_T_F_disjoint_nolimit <- imputation_unit$`FALSE`[,c("t_start_trunc", "t_end_trunc")]  %>% as.matrix()
            dimnames(get_T_F_disjoint_nolimit) <- NULL
            true_T_F_disjoint_nolimit <- matrix(NA, ncol = 2, nrow= 0)
            for(i in 1:(nrow(imputation_unit$`TRUE`)- 1)){
              true_T_F_disjoint_nolimit <- rbind(true_T_F_disjoint_nolimit, 
                                                 c(imputation_unit$`TRUE`[i, "t_end_trunc"], 
                                                   imputation_unit$`TRUE`[i + 1, "t_start_trunc"]))
            }
            
            expect_equal(get_T_F_disjoint_nolimit, true_T_F_disjoint_nolimit)
          }
          
          
          
)