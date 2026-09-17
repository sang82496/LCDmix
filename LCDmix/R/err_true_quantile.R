# Generated from create-LCDmix.Rmd: do not edit by hand

#' True error quantile function of a simulated dataset, on the residual scale
#'
#' @description
#' The inverse of \code{err_true_fun()}. The two are kept side by side on
#' purpose: when a noise family is added, both branches must be added together.
#' Every family is mean zero, so these are quantiles of the error, not of the
#' response.
#'
#' @param sim A simulated dataset from \code{gen_simul_data()}.
#' @param p   Numeric vector of probabilities.
#'
#' @return Numeric vector of quantiles of the true error distribution.
#'
#' @export
err_true_quantile <- function(
  sim,
  p
) {
  if (sim$noisetype == 'skewed') {
    res <- sn::qsn(p, xi = -sim$mn_shift, omega = sim$omega, alpha = sim$skew_alpha)
  } else if (sim$noisetype == 'heavytail') {
    res <- qt(p, df = sim$df) / sqrt(sim$variance)
  } else if (sim$noisetype == 'laplace') {
    res <- VGAM::qlaplace(p, scale = 1)
  } else if (sim$noisetype == 'exponential') {
    res <- qexp(p) - 1
  } else { # Gaussian
    res <- qnorm(p, mean = 0)
  }
  return(res)
}
