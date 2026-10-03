# Generated from create-LCDmix.Rmd: do not edit by hand

#' True component densities of a simulated dataset at one time point
#'
#' @description
#' Evaluates the density of each of the two components at time \eqn{t},
#' \eqn{f(y - \mu_{tk})} with \eqn{\mu_{tk}} = \code{sim$mnmat[t, k]} and
#' \eqn{f} the error density (\code{err_true_fun()}). The gate probabilities
#' are not included: to get the mixture density, multiply column \eqn{k} by
#' \code{sim$prob[t, k]} and sum over \eqn{k}.
#'
#' @param sim A dataset from \code{gen_simul_data()}.
#' @param t Time point.
#' @param y_grid Numeric vector of response values.
#'
#' @return A \code{length(y_grid)} \eqn{\times 2} matrix of densities (a
#'   vector of length 2 when \code{y_grid} has length 1). The number of
#'   components is fixed at 2.
#'
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
        dexp(y_grid - mu + 1) 
      })
    } else { # Gaussian
    res <- sapply(seq_len(K), function(k) {
        mu = sim$mnmat[t,k]
        dnorm(y_grid, mean = mu)
      })
    }
  return(res)
}
