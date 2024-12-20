
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bicre

<!-- badges: start -->

<!-- badges: end -->

Bayesian Imputation for Censored Recurrent Events.

Uses Bayesian data augmentation to fit a proportional hazards Poisson
model for recurrent events which are observed as censored event counts
within intervals which may overlap. A Weibull baseline hazard is used,
with piecewise constant time-varying covariates allowed.
Gamma-distributed frailty random effects are available.

## Installation

You can install the development version of bicre like so:

``` r
remotes::install_github("schnellp/bicre", build_vignettes = TRUE)
```

The vignette provides a simulated example of data formatting and model
fitting. It may take a moment for the vignette to build (it includes
some MCMC sampling), but it should finish within one minute on a modern
laptop.

## Documentation

The vignette is accessible via:

``` r
utils::vignette("analysis-simulate", package = "bicre")
```
