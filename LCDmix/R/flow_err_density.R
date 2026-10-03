# Generated from create-LCDmix.Rmd: do not edit by hand

#' Estimated error density of one flowmix component
#'
#' @description
#' Returns the normal density with mean 0 and variance \code{sigma_k}. flowmix
#' stores the variance of each component in \code{sigma}, so the standard
#' deviation is \code{sqrt(sigma_k)}, as in \code{dens_est_fun()}.
#'
#' @param sigma_k Variance of the component (one element of \code{sigma} of a
#'   flowmix fit).
#'
#' @return A function of a numeric vector \code{v} that returns the densities
#'   at \code{v}.
#' @export
flow_err_density <- function(sigma_k) {
  return(function(v) stats::dnorm(v, 0, sqrt(sigma_k)))
}
