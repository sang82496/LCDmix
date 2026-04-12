# Generated from create-LCDmix.Rmd: do not edit by hand

#' @export
dens_true_fun <- function(
  sim,
  t,
  y_grid
) {
  K = 2
  if (sim$noisetype == 'skewed') {
    res <- sapply(seq_len(K), function(k) {
        mu = sim$mnmat[t,k]
        sn::dsn(y_grid - mu, xi = -sim$mn_shift, omega = sim$omega, alpha = sim$skew_alpha)
      })
    } else if (sim$noisetype == 'heavytail') { 
    res <- sapply(seq_len(K), function(k) {
        mu = sim$mnmat[t,k]
        dt((y_grid - mu) * sqrt(sim$variance), df = sim$df) * sqrt(sim$variance)
      })
    } else if (sim$noisetype == 'laplace') {
    res <- sapply(seq_len(K), function(k) {
        mu = sim$mnmat[t,k]
        VGAM::dlaplace(y_grid - mu, scale = 1) 
      })
    } else if (sim$noisetype == 'exponential') {
    res <- sapply(seq_len(K), function(k) {
        mu = sim$mnmat[t,k]
        dexp(y_grid - mu) 
      })
    } else { # Gaussian
    res <- sapply(seq_len(K), function(k) {
        mu = sim$mnmat[t,k]
        dnorm(y_grid, mean = mu)
      })
    }
  return(res)
}
