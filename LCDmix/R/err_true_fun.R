# Generated from create-LCDmix.Rmd: do not edit by hand

#' True error density of a simulated dataset
#'
#' @description
#' Evaluates the density of the error distribution that
#' \code{gen_simul_data()} used. Both components share it, and it has mean 0.
#' \code{err_true_quantile()} is its inverse; when a noise type is added,
#' change both functions.
#'
#' @param sim A dataset from \code{gen_simul_data()}.
#' @param y_grid Numeric vector of residual values.
#'
#' @return Numeric vector of densities at \code{y_grid}.
#'
#' @export
err_true_fun <- function(
  sim,
  y_grid
) {
  if (sim$noisetype == 'skewed') {
    res <- sn::dsn(y_grid, xi = -sim$mn_shift, omega = sim$omega, alpha = sim$skew_alpha)
  } else if (sim$noisetype == 'heavytail') { 
    res <- dt(y_grid  * sqrt(sim$variance), df = sim$df) * sqrt(sim$variance)
  } else if (sim$noisetype == 'laplace') {
    res <- VGAM::dlaplace(y_grid, scale = 1) 
  } else if (sim$noisetype == 'exponential') {
    res <- dexp(y_grid + 1) 
  } else { # Gaussian
    res <- dnorm(y_grid, mean = 0)
  }
  return(res)
}
