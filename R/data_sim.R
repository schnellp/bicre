#' Simulated Example Data
#'
#' A simulated dataset used as an example in the vignette,
#' illustrating the assumed data format.
#'
#' @format
#' The dataset is a list of two data frames.
#'
#' ## `data_cov`
#' A data frame with 87 rows and 6 columns:
#' \describe{
#'   \item{id}{ID linking covariates in `data_cov` and events in `data_events` for the same individual}
#'   \item{trt}{A binary time-varying treatment indicator}
#'   \item{age}{A numeric time-varying age variable, e.g., in decades}
#'   \item{t_start}{The start times (e.g., in years) of each interval with constant values of time-varying covariates}
#'   \item{t_end}{The end times of those intervals}
#'   \item{sex}{A non-time-varying covariate indicating male sex assigned at birth}
#' }
#'
#' ## `data_events`
#' A data frame with 342 rows and 5 columns:
#' \describe{
#'   \item{id}{ID linking covariates in `data_cov` and events in `data_events` for the same individual}
#'   \item{t_start}{The start times (e.g., in years) of each interval on which partial information about event counts was observed}
#'   \item{t_end}{The end times of those intervals}
#'   \item{e_min}{The minimum number of events that occurred within the interval}
#'   \item{e_max}{The maximum number of events that occurred within the interval}
#' }
"data_sim"
