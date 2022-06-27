
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bicre

<!-- badges: start -->
<!-- badges: end -->

Bayesian Imputation for Censored Recurrent Events.

## Installation

You can install the development version of bicre like so:

``` r
# Package not currently publicly available.
```

## Data formatting

This package works primarily with datasets composed of two data tables:

1.  A table of (possibly time-varying) covariates;
2.  A table of (possibly censored) event counts in (possibly
    overlapping) intervals.

In both the covariate and event tables, a row corresponds to an
interval, and columns are present representing the start and end times
of an interval, as well as a unit identifier (e.g., a subject ID) that
is used to link observations between the two tables. Additionally, the
event table contains columns indicating the minimum and maximum event
counts within each interval. For example, if the minimum and maximum
event counts are both `1`, then it is known that there was exactly one
event within the interval, whereas if the minimum is `1` and the maximum
is `Inf`, then all that is known is that there was at least one event in
the interval.

The following conventions are taken:

1.  Intervals are open on the left and closed on the right;
2.  Covariates are considered constant over each interval;
3.  Event data are constraints only, and do not imply any particular
    distribution of event counts (e.g., Poisson) or times (e.g.,
    uniform) within intervals.

A `co_events` object is created from the covariate and event tables
using the `co_events()` function.

## Example

This is a basic example which shows you how to solve a common problem:

``` r
# library(bicre)
## basic example code
```
