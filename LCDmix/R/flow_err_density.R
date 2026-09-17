# Generated from create-LCDmix.Rmd: do not edit by hand

#' Estimated error density of one flowmix component.
#'
#' NOTE: flowmix stores a variance in `sigma`, so the standard deviation is
#' sqrt(sigma_k). The work order's `dnorm(xs, 0, sigma_k)` passes a variance
#' where a standard deviation belongs. dens_est_fun() uses sqrt(), and so does
#' this function.
#' @export
flow_err_density <- function(sigma_k) {
  return(function(v) stats::dnorm(v, 0, sqrt(sigma_k)))
}
