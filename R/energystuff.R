#' energystuff: a package for energy modeling and monitoring
#' 
#' The energystuff package serves as a repository for various R functions that perform regression fitting and hypothesis testing.
#' 
#' @section Hybrid energy monitoring:
#' `hybridmonitor.single` function performs hypothesis testing for a single utility account.
#' `hybridmonitor.multiple` function equivalently performs hypothesis testing for multiple utility accounts by stochastically constructing a posterior predictive distribution of a test statistic, against which the test data are compared.
#' 
#' 
#' @docType package
#' @name energystuff
#' @useDynLib energystuff, .registration = TRUE
#' @md
#' @importFrom Rcpp evalCpp
NULL