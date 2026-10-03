# Generated from create-LCDmix.Rmd: do not edit by hand

#' Estimated component densities at one time point
#'
#' @description
#' Evaluates the density of each fitted component at time \eqn{t}, for an
#' LCDmix fit or a flowmix fit. For LCDmix (\code{est_res} has
#' \code{alpha_new}), component \eqn{k} is the log-concave density at
#' \eqn{y - \theta_{0k} - X_t^\top \theta_k}, and 0 outside its support. For
#' flowmix, it is the normal density with mean \code{est_res$mn[t, 1, k]} and
#' variance \code{est_res$sigma[k]}. The gate probabilities are not included.
#'
#' @param est_res The \code{iter} element of a \code{main()} fit, or a flowmix
#'   fit (fields \code{mn} and \code{sigma}).
#' @param t Time point.
#' @param y_grid Numeric vector of response values.
#' @param X Covariate matrix; used for LCDmix only.
#'
#' @return A \code{length(y_grid)} \eqn{\times K} matrix of densities.
#'
#' @export
dens_est_fun <- function(
  est_res,
  t,
  y_grid,
  X
) {
  if (!is.null(est_res$alpha_new)) { # if LCDmix 
    K <- length(est_res$g_new)
    res_est <- sapply(seq_len(K), function(k) {
        mu <- est_res$theta0_new[[k]] + sum(X[t,] * est_res$theta_new[[k]])
        logcondens::evaluateLogConDens(y_grid - mu, est_res$g_new[[k]], which = 2)[,3]   # fixT: column 3 = density
      })
    } else { # if flowmix
    mn_arr <- est_res$mn
    sigma  <- as.numeric(est_res$sigma)
    K      <- dim(mn_arr)[3]
    res_est <- sapply(seq_len(K), function(k) {
        dnorm(y_grid, mean = mn_arr[t,1,k], sd = sqrt(sigma[k]))
      })
    }
  return(res_est)
}
