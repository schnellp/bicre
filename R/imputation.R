roundone_impute <- function(roundone_unit,
                            expect_cum_FUN,
                            expect_cum_inverse_FUN,
                            max_tries = 1e4,
                            verbose = FALSE,
                            ...){
  t <- 0
  while(t < max_tries){
    t <- t + 1
    #events simulation
    e_sim <- roundone_unit$sim
    z <- simulate_nonhomog_inversion(
      t_start = e_sim$t_start,
      t_end = e_sim$t_end,
      expect_cum_FUN = expect_cum_FUN,
      expect_cum_inverse_FUN = expect_cum_inverse_FUN,
      count_min = e_sim$e_min,
      count_max = e_sim$e_max,
      ...
    )


    #events check
    if(nrow(roundone_unit$inner_check) > 0){
      counts_inner_check <-  sapply(1 : nrow(roundone_unit$inner_check),
                                    function(i) {
                                      sum(roundone_unit$inner_check$t_start[i] < z &
                                            z <= roundone_unit$inner_check$t_end[i])
                                    })

      if(any(counts_inner_check > roundone_unit$inner_check$e_max |
             counts_inner_check < roundone_unit$inner_check$e_min)){
        if (verbose) {
          print("Failed first-round inner check:")
          print(z)
        }
        next
      }
    }

    if(nrow(roundone_unit$outer_check) > 0){
      counts_outer_check <- sapply(1 : nrow(roundone_unit$outer_check),
                                   function(i) {
                                     sum(roundone_unit$outer_check$t_start[i] < z &
                                           z <= roundone_unit$outer_check$t_end[i])
                                   })
      if(any(counts_outer_check > roundone_unit$inner_check$e_max)){
        if (verbose) {
          print("Failed first-round outer check:")
          print(z)
        }
        next
      }
    }


    return(z)


  }

  stop("round one imputations do not simulate a set of event times satisfying conditions with max tries: increase number of max tries ")
}


roundtwo_impute <- function(roundtwo_unit,
                            expect_cum_FUN,
                            expect_cum_inverse_FUN,
                            max_tries = 1e4,
                            verbose = FALSE,
                            ...){
  t <- 0
  while(t < max_tries){
    t <- t + 1

    #events simulation
    e_leftover <-  roundtwo_unit$leftover

    #Check intervals covering the left over interval provides a upper bound
    #Left over interval cannot cover any check intervals otherwise the check interval
    # should be the disjoint sub interval for simulation in round one
    #Hence only upper bound can be used for simulation in round two.
    e_leftover_upper <- ifelse(nrow(roundtwo_unit$covered) > 0,
                               min(roundtwo_unit$covered$e_max), Inf)


    z <- simulate_nonhomog_inversion(
      t_start = e_leftover$t_start,
      t_end = e_leftover$t_end,
      expect_cum_FUN = expect_cum_FUN,
      expect_cum_inverse_FUN = expect_cum_inverse_FUN,
      count_min = 0,
      count_max = e_leftover_upper,
      ...
    )

    # events check
    if(nrow(roundtwo_unit$overlap_strict) > 0){
      counts_overlap_strict <-  sapply(1 : nrow(roundtwo_unit$overlap_strict),
                                       function(i){
                                         sum(roundtwo_unit$overlap_strict$t_start[i] < z &
                                               z <= roundtwo_unit$overlap_strict$t_end[i])
                                      })

      if(any(counts_overlap_strict > roundtwo_unit$overlap_strict$e_max)){
        if (verbose) {
          print("Failed second-round strict overlapping interval check:")
          print(z)
        }
        next
      }
    }

    return(z)


  }

  stop("round two imputations do not simulate a set of event times satisfying conditions with max tries: increase number of max tries ")
}

sample_sequential <- function(imputation_unit,
                              expect_cum_FUN,
                              expect_cum_inverse_FUN,
                              max_tries = 1e4,
                              verbose = FALSE,
                              ...) {

  t <- 0

  while (t < max_tries) {
    t <- t + 1
    z <- numeric(0)

    # independently impute first round
    # For each single unit disjoint sub interval
    # 1st round check, rejection and re-impute are independent
    # Since they are disjoint

    z_roundone <- imputation_unit$roundone %>% lapply(roundone_impute,
                                        expect_cum_FUN = expect_cum_FUN,
                                        expect_cum_inverse_FUN = expect_cum_inverse_FUN,
                                        max_tries = max_tries,
                                        verbose = verbose,
                                        ...
                                        )

    # independently impute second round for left over intervals
    z_roundtwo <- imputation_unit$roundtwo %>% lapply(roundtwo_impute,
                                                      expect_cum_FUN = expect_cum_FUN,
                                                      expect_cum_inverse_FUN = expect_cum_inverse_FUN,
                                                      max_tries = max_tries,
                                                      verbose = verbose,
                                                      ...
                                                      )

    z <- unlist(c(z_roundone, z_roundtwo)) %>% sort

    # final check of all outer check intervals
    if(nrow(imputation_unit$outer_check_all) > 0){
      counts_outer_check_all <- sapply(1:nrow(imputation_unit$outer_check_all),
                                       function(i) {
                                         sum(imputation_unit$outer_check_all$t_start[i] < z &
                                               z <= imputation_unit$outer_check_all$t_end[i] )
                                       })


      if (any(counts_outer_check_all < imputation_unit$outer_check_all$e_min |
              counts_outer_check_all > imputation_unit$outer_check_all$e_max)) {
        if (verbose) {
          print("Failed second-round check:")
          print(z)
        }
        next
      }

    }


    return(z)
  }

  stop("sequential sampling do not simulate a set of event times satisfying conditions with max tries: increase number of max tries ")
  # return(NULL)
}

impute_single_id <- function(co_events_frame_single,
                             coef,
                             expect_cum_FUN,
                             expect_cum_inverse_FUN,
                             ...) {

  lin_pred <- co_events_frame_single$X %*% coef
  t_breaks <- co_events_frame_single$cov_t$t_end

  # simulate events for easy imputation units
  if(length(co_events_frame_single$ev$impute_unit_ez) > 0){
    z_easy <- co_events_frame_single$ev$impute_unit_ez %>% lapply(function(x) simulate_nonhomog_inversion(
      t_start = x$t_start,
      t_end = x$t_end,
      expect_cum_FUN = expect_cum_FUN,
      expect_cum_inverse_FUN = expect_cum_inverse_FUN,
      count_min = x$e_min,
      count_max = x$e_max,
      lin_pred = lin_pred,
      t_breaks = t_breaks,
      ...
    ))
  }else{
    z_easy <- numeric(0)
  }


  # simulate events for compound imputation units with sequential sampling
  if(length(co_events_frame_single$ev$impute_unit_compund) > 0){
    z_compound <- lapply(co_events_frame_single$ev$impute_unit_compund,
                sample_sequential,
                expect_cum_FUN = expect_cum_FUN,
                expect_cum_inverse_FUN = expect_cum_inverse_FUN,
                lin_pred = lin_pred,
                t_breaks = t_breaks,
                ...
                )
  }else{
    z_compound <- numeric(0)
  }


  z <- sort(unlist(c(z_easy, z_compound)))

  if (class(z) == "list") {
    return(numeric(0))
  } else {
    return(drop(z))
  }
}
