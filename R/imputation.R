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
    for (i in 1 : nrow(imputation_unit$`TRUE`)) {
      ei <- imputation_unit$`TRUE`[i, ]

      zi <- simulate_nonhomog_inversion(
        t_start = ei$t_start,
        t_end = ei$t_end,
        expect_cum_FUN = expect_cum_FUN,
        expect_cum_inverse_FUN = expect_cum_inverse_FUN,
        count_min = ei$e_min,
        count_max = ei$e_max,
        ...
      )

      z <- c(z, zi)
    }

    # first-round consistency check

    if (is.null(imputation_unit$`FALSE`)) {
      return(z)
    }

    counts <- sapply(1 : nrow(imputation_unit$`FALSE`),
                     function(i) {
                       sum(imputation_unit$`FALSE`$t_start[i] < z &
                             z <= imputation_unit$`FALSE`$t_end[i])
                     })

    if (any(counts > imputation_unit$`FALSE`$e_max)) {
      if (verbose) {
        print("Failed first-round check:")
        print(z)
      }
      next
    }

    # impute remainder in second round

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
        ...
      )

      z <- c(z, zi)
    }

    # second-round consistency check

    counts <- sapply(1 : nrow(imputation_unit$`FALSE`),
                     function(i) {
                       sum(imputation_unit$`FALSE`$t_start[i] < z &
                             z <= imputation_unit$`FALSE`$t_end[i])
                     })

    if (any(counts < imputation_unit$`FALSE`$e_min |
            counts > imputation_unit$`FALSE`$e_max)) {
      if (verbose) {
        print("Failed second-round check:")
        print(z)
      }
      next
    }

    counts <- sapply(1 : nrow(imputation_unit$`TRUE`),
                     function(i) {
                       sum(imputation_unit$`TRUE`$t_start[i] < z &
                             z <= imputation_unit$`TRUE`$t_end[i])
                     })

    if (any(counts < imputation_unit$`TRUE`$e_min |
            counts > imputation_unit$`TRUE`$e_max)) {
      if (verbose) {
        print("Failed second-round check:")
        print(z)
      }
      next
    }

    return(z)
  }

  return(NULL)
}

