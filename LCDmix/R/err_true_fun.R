# Generated from create-LCDmix.Rmd: do not edit by hand

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
    res <- VGAM::dlaplace(y_grid, scale = sim$laplace_scale) 
  } else { # Gaussian
    res <- dnorm(y_grid, mean = 0)
  }
  return(res)
}
