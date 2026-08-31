# Generated from create-LCDmix.Rmd: do not edit by hand

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
        logcondens::evaluateLogConDens(y_grid - mu, est_res$g_new[[k]])[,3]
      })
    } else { # if flowmix
    mn_arr <- est_res$mn
    sigma  <- as.numeric(est_res$sigma)
    K      <- dim(mn_arr)[3]
    res_est <- sapply(seq_len(K), function(k) {
        dnorm(y_grid, mean = mn_arr[t,1,k], sd = sqrt(sigma[k]))
      })
    }
  return(res)
}
