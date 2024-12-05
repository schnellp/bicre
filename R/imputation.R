round_disjoint_impute <- function(round_disjoint_unit,
                            expect_cum_FUN,
                            expect_cum_inverse_FUN,
                            max_tries = 2e4,
                            verbose = FALSE,
                            ...){
  t <- 0
  while(t < max_tries){
    t <- t + 1
    #events simulation
    e_sim <- round_disjoint_unit$sim
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
    if(nrow(round_disjoint_unit$inner_check) > 0){
      counts_inner_check <-  sapply(1 : nrow(round_disjoint_unit$inner_check),
                                    function(i) {
                                      sum(round_disjoint_unit$inner_check$t_start[i] < z &
                                            z <= round_disjoint_unit$inner_check$t_end[i])
                                    })

      if(any(counts_inner_check > round_disjoint_unit$inner_check$e_max |
             counts_inner_check < round_disjoint_unit$inner_check$e_min)){
        if (verbose) {
          print("Failed first-round inner check:")
          print(z)
        }
        next
      }
    }

    if(nrow(round_disjoint_unit$outer_check) > 0){
      counts_outer_check <- sapply(1 : nrow(round_disjoint_unit$outer_check),
                                   function(i) {
                                     sum(round_disjoint_unit$outer_check$t_start[i] < z &
                                           z <= round_disjoint_unit$outer_check$t_end[i])
                                   })
      if(any(counts_outer_check > round_disjoint_unit$outer_check$e_max)){
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


round_disjoint_impute_w_tiny <- function(round_disjoint_unit,
                                         expect_cum_FUN,
                                         expect_cum_inverse_FUN,
                                         max_tries = 2e4,
                                         verbose = FALSE,
                                         ...){
  t <- 0
  while(t < max_tries){
    t <- t + 1
    #events simulation
    e_sim <- round_disjoint_unit$sim
    z <- simulate_nonhomog_inversion_multi_interval(
      interval_df = e_sim$int,
      expect_cum_FUN = expect_cum_FUN,
      expect_cum_inverse_FUN = expect_cum_inverse_FUN,
      count_min = e_sim$range$e_min,
      count_max = e_sim$range$e_max,
      ...
    )


    #events check
    if(nrow(round_disjoint_unit$inner_check) > 0){
      counts_inner_check <-  sapply(1 : nrow(round_disjoint_unit$inner_check),
                                    function(i) {
                                      sum(round_disjoint_unit$inner_check$t_start[i] < z &
                                            z <= round_disjoint_unit$inner_check$t_end[i])
                                    })

      if(any(counts_inner_check > round_disjoint_unit$inner_check$e_max |
             counts_inner_check < round_disjoint_unit$inner_check$e_min)){
        if (verbose) {
          print("Failed first-round inner check:")
          print(z)
        }
        next
      }
    }

    if(nrow(round_disjoint_unit$outer_check) > 0){
      counts_outer_check <- sapply(1 : nrow(round_disjoint_unit$outer_check),
                                   function(i) {
                                     sum(round_disjoint_unit$outer_check$t_start[i] < z &
                                           z <= round_disjoint_unit$outer_check$t_end[i])
                                   })
      if(any(counts_outer_check > round_disjoint_unit$outer_check$e_max)){
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




round_leftover_impute <- function(round_leftover_unit,
                            expect_cum_FUN,
                            expect_cum_inverse_FUN,
                            max_tries = 2e4,
                            verbose = FALSE,
                            ...){
  t <- 0
  while(t < max_tries){
    t <- t + 1

    #events simulation
    e_leftover <-  round_leftover_unit$leftover

    #Check intervals covering the left over interval provides a upper bound
    #Left over interval cannot cover any check intervals otherwise the check interval
    # should be the disjoint sub interval for simulation in round one
    #Hence only upper bound can be used for simulation in round two.
    e_leftover_upper <- ifelse(nrow(round_leftover_unit$covered) > 0,
                               min(round_leftover_unit$covered$e_max), Inf)


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
    if(nrow(round_leftover_unit$overlap_strict) > 0){
      counts_overlap_strict <-  sapply(1 : nrow(round_leftover_unit$overlap_strict),
                                       function(i){
                                         sum(round_leftover_unit$overlap_strict$t_start[i] < z &
                                               z <= round_leftover_unit$overlap_strict$t_end[i])
                                      })

      if(any(counts_overlap_strict > round_leftover_unit$overlap_strict$e_max)){
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

#' Sequential sampling for a compound inputation unit
sample_sequential <- function(imputation_unit,
                              expect_cum_FUN,
                              expect_cum_inverse_FUN,
                              max_tries = 1e5,
                              verbose = FALSE,
                              round_disjoint_impute_fun,
                              ...) {

  t <- 0

  while (t < max_tries) {
    t <- t + 1
    z <- numeric(0)

    # independently impute first round
    # For each single unit disjoint sub interval
    # 1st round check, rejection and re-impute are independent
    # Since they are disjoint

    z_round_disjoint <- imputation_unit$round_disjoint %>% lapply(round_disjoint_impute_fun,
                                        expect_cum_FUN = expect_cum_FUN,
                                        expect_cum_inverse_FUN = expect_cum_inverse_FUN,
                                        max_tries = max_tries,
                                        verbose = verbose,
                                        ...
                                        )

    # independently impute second round for left over intervals
    z_round_leftover <- imputation_unit$round_leftover %>% lapply(round_leftover_impute,
                                                      expect_cum_FUN = expect_cum_FUN,
                                                      expect_cum_inverse_FUN = expect_cum_inverse_FUN,
                                                      max_tries = max_tries,
                                                      verbose = verbose,
                                                      ...
                                                      )

    z <- unlist(c(z_round_disjoint, z_round_leftover)) %>% sort

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


sample_sequential_initialize <- function(imputation_unit,
                                         expect_cum_FUN,
                                         expect_cum_inverse_FUN,
                                         max_tries = 1e4,
                                         verbose = FALSE,
                                         round_disjoint_impute_fun,
                                         ...) {

  t <- 0
  reject_reason <- matrix(NA,ncol = 2, nrow = max_tries)
  colnames(reject_reason) <- c("less_emin", "larger_emax")
  while (t < max_tries) {
    # print(t)
    t <- t + 1
    z <- numeric(0)

    # independently impute first round
    # For each single unit disjoint sub interval
    # 1st round check, rejection and re-impute are independent
    # Since they are disjoint

    z_round_disjoint <- imputation_unit$round_disjoint %>% lapply(round_disjoint_impute_fun,
                                                                  expect_cum_FUN = expect_cum_FUN,
                                                                  expect_cum_inverse_FUN = expect_cum_inverse_FUN,
                                                                  max_tries = max_tries,
                                                                  verbose = verbose,
                                                                  ...
    )

    # independently impute second round for left over intervals
    z_round_leftover <- imputation_unit$round_leftover %>% lapply(round_leftover_impute,
                                                                  expect_cum_FUN = expect_cum_FUN,
                                                                  expect_cum_inverse_FUN = expect_cum_inverse_FUN,
                                                                  max_tries = max_tries,
                                                                  verbose = verbose,
                                                                  ...
    )

    z <- unlist(c(z_round_disjoint, z_round_leftover)) %>% sort

    # final check of all outer check intervals
    if(nrow(imputation_unit$outer_check_all) > 0){
      counts_outer_check_all <- sapply(1:nrow(imputation_unit$outer_check_all),
                                       function(i) {
                                         sum(imputation_unit$outer_check_all$t_start[i] < z &
                                               z <= imputation_unit$outer_check_all$t_end[i] )
                                       })

      reject_reason[t,"less_emin"] <- sum(counts_outer_check_all < imputation_unit$outer_check_all$e_min)
      reject_reason[t,"larger_emax"] <- sum(counts_outer_check_all > imputation_unit$outer_check_all$e_max)
      if (sum(reject_reason[t,]) > 0 ) {
        if (verbose) {
          print("Failed second-round check:")
          print(z)
        }
        next
      }

    }


    return(z)
  }
  reject_reason <-  reject_reason  %>% colMeans()
  if(reject_reason["less_emin"] > reject_reason["larger_emax"]){
    return("Sequential sampling: larger ui needed")
  }else if(reject_reason["less_emin"] < reject_reason["larger_emax"]){
    return("Sequential sampling: smaller ui needed")
  }else{
    stop("sequential sampling do not simulate a set of event times satisfying conditions with max tries: increase number of max tries ")
  }


  # return(NULL)
}


#' impute events for a single patient
#' @export
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

    if(is.null(co_events_frame_single$ev$impute_unit_compund[[1]]$round_disjoint[[1]]$sim$int)){
      round_disjoint_impute_fun <- round_disjoint_impute
    }else{
      round_disjoint_impute_fun <- round_disjoint_impute_w_tiny
    }

    z_compound <- lapply(co_events_frame_single$ev$impute_unit_compund,
                sample_sequential,
                expect_cum_FUN = expect_cum_FUN,
                expect_cum_inverse_FUN = expect_cum_inverse_FUN,
                lin_pred = lin_pred,
                t_breaks = t_breaks,
                round_disjoint_impute_fun = round_disjoint_impute_fun,
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

#' impute events for a single patient for getting proper first values
#' @export

impute_single_id_initialize <- function(co_events_frame_single,
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

    if(is.null(co_events_frame_single$ev$impute_unit_compund[[1]]$round_disjoint[[1]]$sim$int)){
      round_disjoint_impute_fun <- round_disjoint_impute
    }else{
      round_disjoint_impute_fun <- round_disjoint_impute_w_tiny
    }

    z_compound <- list()
    for(j in 1:length(co_events_frame_single$ev$impute_unit_compund)){
     events_imputed <- sample_sequential_initialize(co_events_frame_single$ev$impute_unit_compund[[j]],
                                                      expect_cum_FUN = expect_cum_FUN,
                                                      expect_cum_inverse_FUN = expect_cum_inverse_FUN,
                                                      lin_pred = lin_pred,
                                                      t_breaks = t_breaks,
                                                      round_disjoint_impute_fun = round_disjoint_impute_fun,
                                                      ...)
     if(class(events_imputed) == "character"){
       return(events_imputed)
     }else{
       z_compound[[j]] <- events_imputed
     }

    }


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
