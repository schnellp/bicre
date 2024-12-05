whether_disjoint <- function(df1, df2){
  all(apply(df2,1, function(x) all(x["t_start"] >= df1["t_end"] | x["t_end"] <= df1["t_start"])))
}


whether_cover <- function(covering_df1, covered_df2){
  # note the covering_df1 must be a clean interval without overlapping and without the case where one interval's t_start
  # is equal to another one's t_end.
  # in the following use, this is guarateened.
  all(apply(covered_df2,1, function(x) any(covering_df1["t_start"] <= x["t_start"] & x["t_end"] <= covering_df1["t_end"])))
}

# whether_cover <- function(covering_easy, covered_compound){
#   all(with(covered_compound, covering_easy["t_start"] <= t_start  & t_end <= covering_easy["t_end"]))
# }
#l
# whether_cover_multi <- function(covering_easy_multi, covered_compound){
#    apply(covering_easy_multi, 1, function(x)  whether_cover(x,covered_compound))
# }

intersect_length <- function(int1, int2){
  apply(int1, 1, function(x)
    pmax(0, pmin(x["t_end"], int2[,"t_end"]) -
           pmax(x["t_start"],int2[, "t_start"]))) %>%  sum
}

sim_FALSE_change <- function(e_sim_int, FALSE_inner_check){
  for(jj in 1:nrow(FALSE_inner_check)){
    e_F <- FALSE_inner_check[jj,]
    sim_change_i <- which(e_sim_int$t_start < e_F$t_start & e_F$t_end < e_sim_int$t_end)
    e_sim_int <- e_sim_int %>% add_row(t_start = e_F$t_end, t_end = e_sim_int$t_end[sim_change_i], .after = sim_change_i)
    e_sim_int$t_end[sim_change_i] <- e_F$t_start
  }

  e_sim_int
}

tiny_interval_prepare <- function(compound_impute_unit, tiny_diff, length_rate_perc_prod = 0.10, rate_perc = 0.25){
  #First find inty interval that is hard to calculate.
  # any interval with emin not 0 could possible be that.
  # each interval can have multiple overlapping tiny intervals

  tiny_simp <- data.frame(matrix(nrow = 0, ncol = 8))
  colnames(tiny_simp) <-  c(colnames(compound_impute_unit), "compound_index", "length", "rate")
  tiny_compound <- list()
  tiny_compound_first <-  data.frame(matrix(nrow = 0, ncol = 7))

  compound_index <- 1

  for(j in 1:nrow(compound_impute_unit)){
    if(compound_impute_unit[j,]$e_min > 0){
      e <- compound_impute_unit %>% slice(j)
      post_overlap <- compound_impute_unit[-j,] %>%
        filter( t_start >= e$t_start & t_end >= e$t_end & t_start < e$t_end &
                  e_max < e$e_min &
                  # Note only the higher emin interval - the lower emax interval
                  # forms the tiny interval and only they matter,
                  # we are not calculating the difference of the intervals here
                  # (t_start - e$t_start +t_end -e$t_end) <= tiny_diff)
                  (t_start - e$t_start) <= tiny_diff)
      if(nrow(post_overlap) >0 ){
        post_tiny <- with(post_overlap, data.frame(t_start = e$t_start, t_end = t_start,
                                                   e_min = e$e_min - e_max, e_max = e$e_max, unit = unit))
        tiny_simp <- rbind(tiny_simp,
                           post_tiny)
      }



      pre_overlap <-  compound_impute_unit[-j,] %>%
        filter( t_start <=  e$t_start & t_end <= e$t_end & t_end > e$t_start &
                  e_max < e$e_min &
                  # Same here Note only the higher emin interval - the lower emax interval
                  # forms the tiny interval and only they matter,
                  # we are not calculating the difference of the intervals here
                  # (-(t_start - e$t_start +t_end -e$t_end)) <=  tiny_diff)
                  (e$t_end - t_end) <=  tiny_diff)

      if(nrow(pre_overlap) > 0){
        pre_tiny <- with(pre_overlap,
                         data.frame(t_start = t_end, t_end = e$t_end,
                                    e_min = e$e_min - e_max, e_max = e$e_max, unit = unit))
        tiny_simp <- rbind(tiny_simp,
                           pre_tiny)
      }


      covered <- compound_impute_unit[-j,] %>%
        filter( t_start >  e$t_start & t_end < e$t_end &
                  e_max < e$e_min &
                  # only when we have covering we add both ends up
                  (t_start - e$t_start + e$t_end - t_end ) <=  tiny_diff)


      if(nrow(covered) > 0){
        for(jj in 1:nrow(covered)){
          covered_tiny <- list(list(int = data.frame(t_start = c(e$t_start, covered$t_end[[jj]]),
                                                     t_end = c(covered$t_start[jj], e$t_end)),
                                    e_min = e$e_min - covered$e_max[[jj]],
                                    e_max = e$e_max,
                                    unit = e$unit))
          tiny_compound <- c(tiny_compound,
                             covered_tiny)

          tiny_compound_first <- rbind(tiny_compound_first,
                                       c(covered_tiny[[1]]$int$t_start[1],
                                         covered_tiny[[1]]$int$t_end[1],
                                         covered_tiny[[1]]$e_min,
                                         covered_tiny[[1]]$e_max,
                                         covered_tiny[[1]]$unit,
                                         compound_index,
                                         sum(covered_tiny[[1]]$int$t_end - covered_tiny[[1]]$int$t_start)) )
          compound_index <- compound_index + 1


        }
      }
      # covered_tiny <- list()
      # for(jj in 1:nrow(covered)){
      #   covered_tiny[[jj]] <- list(t_start = c(e$t_start, covered$t_end[[jj]]),
      #                              t_end = c(covered$t_start[jj], e$t_end),
      #                              e_min = e$e_min - covered$e_max[[jj]],
      #                              e_max = e$e_max,
      #                              unit = e$unit)
      #
      # # for true tiny cases, I assume they covered tiny intervals cannot be covered
      # # by those single tiny intervals
      #
      # }



    }
  }

  colnames(tiny_compound_first) <- c(colnames(compound_impute_unit), "compound_index", "length")
  tiny_compound_first$rate <- with(tiny_compound_first, e_min/length)


  compound_impute_unit$length <- with(compound_impute_unit, t_end - t_start)
  compound_impute_unit$compound_index <- 0
  compound_impute_unit$rate <- with(compound_impute_unit, e_min/length)

  if(nrow(tiny_simp) > 0){
    # for the best results , the tiny_simp and tiny_compound should simplify together
    # but simplification of tiny_compound is complicated and won't affect the results here I am not doing it.
    tiny_simp <- tiny_simp %>% simplify_ev
    tiny_simp$length <- tiny_simp$t_end - tiny_simp$t_start
    tiny_simp$compound_index <- 0
    tiny_simp$rate <- with(tiny_simp, e_min/length)
  }






  tiny_int_cal <- rbind(tiny_simp, tiny_compound_first)


  if(nrow(tiny_int_cal) > 0){
    cal_take <- c()
    for(j in 1:nrow(tiny_int_cal)){
      e <- tiny_int_cal %>% slice(j)
      if(tiny_int_cal$compound_index[j] == 0){
        e_int <- e %>% select(t_start, t_end)
      }else{
        e_int <- tiny_compound[[e$compound_index]]$int
      }

      e_exchangable_fromob <- compound_impute_unit %>%
        filter(e_min > e$e_min &
                 ( ((rate/e$rate) > rate_perc &
                     (apply(e_int, 1, function(x)
                       pmax(0, pmin(x["t_end"], compound_impute_unit$t_end) -
                              pmax(x["t_start"],compound_impute_unit$t_start))) %>% rowSums) > 0) |
                 ((rate/e$rate) *
                  ((apply(e_int, 1, function(x)
                    pmax(0, pmin(x["t_end"], compound_impute_unit$t_end) -
                           pmax(x["t_start"],compound_impute_unit$t_start))) %>% rowSums)/e$length)) > length_rate_perc_prod ) )


      e_exchangable_fromself <- tiny_int_cal %>%
        filter(e_min > e$e_min)


      if(nrow(e_exchangable_fromself) > 0){
        exchangable_fromself_take <- c()
        for(jj in 1:nrow(e_exchangable_fromself)){
          if(e_exchangable_fromself$compound_index[jj] == 0){
            intersection_length <- intersect_length(e_int, e_exchangable_fromself[jj,c("t_start", "t_end")])
          }else{
            intersection_length <-
              intersect_length(e_int,
                               tiny_compound[[e_exchangable_fromself$compound_index[jj]]]$int)
          }

          if((e_exchangable_fromself$rate[jj]/e$rate * intersection_length/e$length) > length_rate_perc_prod){
            exchangable_fromself_take <- c(exchangable_fromself_take, jj)
          }else if(intersection_length > 0 & e_exchangable_fromself$rate[jj]/e$rate > rate_perc){
            exchangable_fromself_take <- c(exchangable_fromself_take, jj)
          }

        }
        e_exchangable_fromself <- e_exchangable_fromself[exchangable_fromself_take,]

      }


      if( (nrow(e_exchangable_fromob) + nrow(e_exchangable_fromself)) == 0 ){
        cal_take <- c(cal_take, j)
      }
    }
    tiny_int_cal <- tiny_int_cal[cal_take,]
  }





  #it does not matter whether there is covering or overlapping between tiny intervals
  #Since we only need to pick out those with the most events and shortest length.
  # all others are not considered anyway in the process.
  # tiny_simp <- tiny_simp %>% simplify_ev







  #add in tiny intervals from the original observations
  if(nrow(tiny_int_cal) > 0){
    tiny_int_cal$from_ob <- 0
  }else{
    tiny_int_cal$from_ob <- numeric(0)
  }


  tiny_ob_row <- with(compound_impute_unit, which(length <= tiny_diff & e_min > 0))
  tiny_int_ob <- compound_impute_unit[tiny_ob_row,]
  tiny_int_ob$from_ob <- tiny_ob_row

  if(nrow(tiny_int_ob) > 0){
    ob_take <- c()
    for(j in 1:nrow(tiny_int_ob)){
      e <- tiny_int_ob %>% slice(j)
      e_exchangable <- compound_impute_unit[-e$from_ob,] %>%
        filter(e_min > e$e_min &
                 ( (rate/e$rate > rate_perc &
                     pmax(0, pmin(e$t_end, t_end) - pmax(e$t_start,t_start)) > 0) |
                 ((pmax(0, pmin(e$t_end, t_end) - pmax(e$t_start,t_start))/e$length) *
                 rate/e$rate) > length_rate_perc_prod) )

      if(nrow(e_exchangable) == 0){
        ob_take <- c(ob_take, j)
      }
    }
    tiny_int_ob <- tiny_int_ob[ob_take,]
  }



  # choose disjoint tiny intervals with the highest lower bounds and the shortest length
  # since potential unqualified tiny intervals is already removed, rate would be the good for choosing
  # more important tiny intervals.
  tiny_int <- rbind(tiny_int_cal, tiny_int_ob) %>% arrange(desc(rate), desc(e_min), length)
  # tiny_int <- rbind(tiny_int, tiny_int_ob) %>% arrange(desc(e_min), length)


  if(nrow(tiny_int) > 0){
    if(tiny_int$compound_index[1] > 0){
      check.overlap <- tiny_compound[[tiny_int$compound_index[1]]]$int
    }else{
      check.overlap <- tiny_int[1,c("t_start", "t_end")]
    }

    if(nrow(tiny_int) > 1){
      tiny_take <- 1

      for(i in 2:nrow(tiny_int)){
        if(tiny_int$compound_index[i]>0){
          if(whether_disjoint(check.overlap, tiny_compound[[tiny_int$compound_index[i]]]$int)){
            tiny_take <- c(tiny_take, i)
            check.overlap <- rbind(check.overlap, tiny_compound[[tiny_int$compound_index[i]]]$int)
          }
        }else{
          if(whether_disjoint(check.overlap, tiny_int[i,c("t_start", "t_end")])){
            tiny_take <- c(tiny_take, i)
            check.overlap <-  rbind(check.overlap, tiny_int[i,c("t_start", "t_end")])
          }
        }

      }


      tiny_int <-  tiny_int[tiny_take,]
    }


    #delete tiny intervals selected in original observations from non-tiny observations.
    tiny_int_fromob_row <- tiny_int$from_ob[tiny_int$from_ob != 0]
    if(length(tiny_int_fromob_row) > 0 ){
      compound_impute_unit_omit_tiny <-  compound_impute_unit[-tiny_int_fromob_row,]
    }else{
      compound_impute_unit_omit_tiny <- compound_impute_unit
    }


    ### tiny intervals' upper bound can be changed by the smaller upper bound of
    ### any T check interval strictly covers the former tiny intervals.

    for(i in 1:nrow(tiny_int)){
      e <- tiny_int %>% slice(i)
      if(tiny_int$compound_index[i] == 0){
        # the reason using compound_impute_unit here is that it at least contains
        # itself so it won't throw an error when no other interval covers it.
        tiny_int$e_max[i] <- compound_impute_unit %>%
          filter(t_start <= e$t_start &
                   t_end >= e$t_end &
                   e_max <= e$e_max) %>% pull(e_max) %>% min
      }else{
        cover_TF <-
          # whether_cover_multi(compound_impute_unit, tiny_compound[[tiny_int$compound_index[i]]]$int)
          compound_impute_unit %>% apply(1, function(x)  whether_cover(x, tiny_compound[[tiny_int$compound_index[i]]]$int))
        tiny_int$e_max[i] <- compound_impute_unit %>% filter(cover_TF & e_max <= e$e_max) %>% pull(e_max) %>% min
      }
    }

  }else{
    check.overlap <- data.frame(t_start = numeric(0), t_end = numeric(0))
    compound_impute_unit_omit_tiny <- compound_impute_unit
  }

  tiny_int <- tiny_int %>% select(-from_ob)






  return(list(tiny_int = tiny_int,
              check.overlap = check.overlap,
              tiny_compound = tiny_compound,
              compound_impute_unit_omit_tiny = compound_impute_unit_omit_tiny))

}


# tiny_int_info <- tiny_interval_prepare(ev[[1]], tiny_diff)

compound_impute_unit_prepare_w_tiny_interval <- function(tiny_int_info, prepare_style){
  list2env(tiny_int_info, envir = environment())
  if(prepare_style == "target_min"){
    compound_impute_unit_omit_tiny$rate <- with(compound_impute_unit_omit_tiny, e_min/length)
    compound_impute_unit_cal <- compound_impute_unit_omit_tiny %>% arrange(desc(rate), length) %>% select(-any_of(c("length", "rate")))
  }else if(prepare_style == "target_max"){
    compound_impute_unit_omit_tiny$rate <- with(compound_impute_unit_omit_tiny, e_max/length)
    compound_impute_unit_cal <- compound_impute_unit_omit_tiny %>% arrange(rate, desc(length)) %>% select(-any_of(c("length", "rate")))
  }else if(prepare_style == "normal"){
    compound_impute_unit_cal <- compound_impute_unit_omit_tiny %>% select(-any_of(c("length", "rate")))
  }else if(prepare_style == "target_emin"){
    compound_impute_unit_cal <- compound_impute_unit_omit_tiny %>% arrange(desc(e_min), length, desc(e_max)) %>% select(-any_of(c("length", "rate")))
  }else{
    stop("type in a proper prepare style")
  }


  impute_disjoint_intervals <- tiny_int %>% select(-any_of(c("length", "rate")))

  check_intervals_row <- c()
  if(nrow(compound_impute_unit_cal) > 0){
    for(j in 1:nrow(compound_impute_unit_cal)){
      if(whether_disjoint(check.overlap, compound_impute_unit_cal[j, c("t_start", "t_end")])){
        impute_disjoint_intervals <- rbind(impute_disjoint_intervals, compound_impute_unit_cal[j, ])
        check.overlap <- rbind(check.overlap, compound_impute_unit_cal[j, c("t_start", "t_end")])
      }else{
        check_intervals_row <- c(check_intervals_row, j)
      }
    }
  }
  # now check.overlap contain all intervals that the  impute_disjoint_intervals covers.
  # So it can be used to make the left over intervals
  disjoint_intervals_covering <- check.overlap %>% arrange(t_start)

  check_intervals <- compound_impute_unit_cal[check_intervals_row,] %>% arrange(t_start, t_end, e_min, e_max)
  impute_disjoint_intervals <- impute_disjoint_intervals %>% arrange(t_start, t_end, e_min, e_max)


  # Prepare for 1st round simulation: disjoint simulation
  FALSE_check_row <- check_intervals$e_max == 0
  round_disjoint <- vector("list", length = nrow(impute_disjoint_intervals))
  inner_check_row_all <- c()
  for(j in 1 : nrow(impute_disjoint_intervals)){
    # For each Unit disjoint sub F-interval
    # Nothing needs preparing since no simulation there
    # For each Unit disjoint sub T-interval, find
    # a.	Inner Unit check sub interval covered by it
    # b.  Outer Unit check sub interval overlapping with it
    e_sim <- impute_disjoint_intervals[j,]
    if(e_sim$compound_index > 0){
      round_disjoint[[j]]$sim <- list(int = tiny_compound[[e_sim$compound_index]]$int,
                                      range = data.frame(e_min =  tiny_compound[[e_sim$compound_index]]$e_min,
                                                         e_max =  tiny_compound[[e_sim$compound_index]]$e_max))
    }else{
      round_disjoint[[j]]$sim <- list(int = e_sim[,c("t_start", "t_end")],
                                      range = e_sim[, c("e_min", "e_max")])
    }

    # initialized inner_check_row for later inner_check_row_all updates
    inner_check_row <- c()
    if(e_sim$e_max != 0){
      # check_intervals %>% apply(1, function(x)  whether_cover(x,round_disjoint[[j]]$sim$int))
      inner_check_row <- check_intervals %>% split(1:nrow(check_intervals)) %>%
        lapply( whether_cover,covering_df1 = round_disjoint[[j]]$sim$int) %>% unlist
      check_row_all <- check_intervals %>%  apply(1, function(x) !whether_disjoint(x,round_disjoint[[j]]$sim$int) )
      outer_check_row <- check_row_all & (!inner_check_row)


      # if there is false interval as inner check interval, the information can be integrated to change the disjoint simulation
      # interval for direct simulation without the need for check
      if(any(FALSE_check_row & inner_check_row)){
        FALSE_inner_check <- check_intervals[FALSE_check_row & inner_check_row,]
        round_disjoint[[j]]$sim$int <- sim_FALSE_change(round_disjoint[[j]]$sim$int, FALSE_inner_check)
      }

      # then the inner check only contains the true check rows being contained in the disjoint simulation interval
      # Note inner_check_row is kept as it is for later use to pick check intervals in the final check
      round_disjoint[[j]]$inner_check <- check_intervals[(!FALSE_check_row) & inner_check_row,]

      # for outer check, only upper bound is used hence delete those without upper bounds
      # delete those upper bounds that are larger than or equal to the simulation upper bound. since it will always satisfies.
      # since check interval still covers other area, so lower bound can not contribute

      round_disjoint[[j]]$outer_check <- check_intervals[outer_check_row,] %>% filter(e_max < e_sim$e_max)
    }else{
      round_disjoint[[j]]$inner_check <- round_disjoint[[j]]$outer_check <-  data.frame(matrix(nrow = 0, ncol = ncol(check_intervals)))
    }

    if(any(inner_check_row)){
      inner_check_row_all <- c(inner_check_row_all, which(inner_check_row))
    }
  }

  # Prepare for 2nd round simulation: leftover simulation
  left_over <- data.frame(t_start = disjoint_intervals_covering$t_end[-nrow(disjoint_intervals_covering)],
                          t_end = disjoint_intervals_covering$t_start[-1])
  unit_start <- min(disjoint_intervals_covering$t_start[1], check_intervals$t_start[1])
  unit_end <- max(last(disjoint_intervals_covering$t_end), check_intervals$t_end)
  left_over <- rbind(data.frame(t_start = unit_start,
                                t_end = disjoint_intervals_covering$t_start[1]),
                     left_over,
                     data.frame(t_start = last(disjoint_intervals_covering$t_end),
                                t_end = unit_end) )
  left_over <- left_over[with(left_over, t_start != t_end) ,]

  # Pick out outer unit check intervals, these have to be checked again in the final check
  if(length(inner_check_row_all) != 0){
    outer_check_all <- check_intervals[-inner_check_row_all,]
  }else{
    outer_check_all <- check_intervals
  }

  if(nrow(left_over) != 0){
    round_leftover <- vector("list", length = nrow(left_over))
    for(j in 1:nrow(left_over)){
      # For each left over interval, it cannot contain any check intervals
      # otherwise that check interval can be put into simulation disjoint intervals (either T or F)
      # so only the following two scenarios:
      #a.	Outer Unit check sub interval covering this left over interval
      #b.	Outer Unit check sub interval overlapping with it
      #These two must be T- intervals since all inner F intervals must be contained
      #in a T-interval or disjoint with a T-interval for each T-interval after simplification.
      e_leftover <- left_over[j,]
      round_leftover[[j]]$leftover <- e_leftover
      covered_row <- with(outer_check_all, t_start <= e_leftover$t_start &  t_end >= e_leftover$t_end)
      overlap_row <- with(outer_check_all, t_end >  e_leftover$t_start & t_start < e_leftover$t_end)
      # for check sub interval covering the left over, only upper bound is useful
      # since check interval still covers other area
      round_leftover[[j]]$covered <- outer_check_all[covered_row,] %>% filter(e_max != Inf)
      # for check sub interval strictly overlapping with the left over, only upper bound is useful
      # since check interval still covers other area
      round_leftover[[j]]$overlap_strict <- outer_check_all[overlap_row & (!covered_row), ] %>% filter(e_max != Inf)
    }
  }else{
    round_leftover <- NULL
  }


  return(list(round_disjoint = round_disjoint,
              round_leftover = round_leftover,
              outer_check_all = outer_check_all,
              unit_range = data.frame(t_start = unit_start,
                                      t_end = unit_end)))

}

seq_samp_prepare_ev_w_tiny_interval <- function(ev, tiny_diff, prepare_style, length_rate_perc_prod){
  # partition by imputation unit
  ev <- ev %>%
    split(f = ev %>% select(unit))


  impute_unit_ez <- list()
  ez_index <- 1
  impute_unit_compound <- list()
  compound_index <- 1
  for (i in 1:length(ev)){
    if(nrow(ev[[i]]) == 1){
      # find out easy imputation unit that only contains one T interval
      # ev[[i]]$sim <- T
      impute_unit_ez[[ez_index]] <- ev[[i]]
      ez_index <- ez_index + 1
    }else{
      impute_unit_compound[[compound_index]] <-
        compound_impute_unit_prepare_w_tiny_interval(tiny_interval_prepare(ev[[i]], tiny_diff, length_rate_perc_prod), prepare_style)
      compound_index <- compound_index + 1
    }

  }


  # impute_unit_ez <- impute_unit_ez %>% do.call(rbind,.)


  # covered_ints <- rbind(impute_unit_ez[,c(1,2)],
  #      impute_unit_compound %>% lapply("[[","unit_range") %>% do.call(rbind,.)) %>%
  #   collapse_interval_data(t_start, t_end)
  ### record covered intervals

  impute_unit_single <- list(impute_unit_ez = impute_unit_ez,
                             impute_unit_compound = impute_unit_compound)

  covered_ints <- list()

  for (i in 1 : length(ev)) {
    covered_ints[[i]] <- data.frame(
      t_start = min(ev[[i]]$t_start),
      t_end = max(ev[[i]]$t_end)
    )
  }

  attributes(impute_unit_single)$covered_ints <- do.call("rbind", covered_ints) %>%
    collapse_interval_data(t_start, t_end) %>% as.data.frame()


  return(impute_unit_single)
}



#' build imputation units for a single id with tiny intervals
#' @export
build_imputation_units_single_w_tiny_interval <- function(ev, tiny_diff, prepare_style, length_rate_perc_prod){
  ev %>% simplify_ev %>% partition_ev %>% seq_samp_prepare_ev_w_tiny_interval(tiny_diff, prepare_style, length_rate_perc_prod)
}
