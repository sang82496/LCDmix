# Generated from create-LCDmix.Rmd: do not edit by hand

#' Compute biomass‐weighted L1 (total‐variation) distance between estimated and true mixtures
#'
#' @description
#' For each time point \code{t}, evaluates the estimated mixture density 
#' (either from an LCDmix fit or a flowmix fit) and the true mixture density 
#' at the observed bin midpoints, computes the \(L^1\) distance, and then 
#' returns both the per‐time distances and the overall weighted average 
#' distance (weighted by total biomass at each \code{t}).
#'
#' @param Y_bin List of length \code{TT}; each element is a numeric vector of bin midpoints for time \code{t}.
#' @param X Numeric \code{TT × p} matrix of covariates (rows = time points).
#' @param bin_mass List of length \code{TT}; each element is a numeric vector of biomass weights aligned with \code{Y_bin[[t]]}.
#' @param est_res A fitted model object:
#'   \describe{
#'     \item{LCDmix fits}{A list containing \code{alpha_new}, \code{theta0_new}, \code{theta_new}, \code{g_new}.}
#'     \item{flowmix fits}{A list containing \code{alpha}, \code{mn}, and \code{sigma}.}
#'   }
#' @param true_res A list with components:
#'   \describe{
#'     \item{\code{prob}}{\code{TT × K} numeric matrix of true mixing proportions \eqn{\pi_{t,k}}.}
#'     \item{\code{dens_true}}{List of length \code{TT}, each an \code{n_t × K} numeric matrix of the true component densities evaluated at \code{Y_bin[[t]]}.}
#'   }
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{per_time}}{Numeric vector of length \code{TT}, the \(L^1\) distance at each time point.}
#'   \item{\code{weighted}}{Scalar, the biomass‐weighted average of those distances.}
#' }
#'
#' @examples
#' \dontrun{
#' sim <- generate_skewed_data(seed = 42)
#' Y_bin    <- sim$ylist
#' X        <- sim$X
#' bin_mass <- sim$countslist
#'
#' # Fit an LCDmix model
#' fit <- main(Y = Y_bin, X = X, biomass = bin_mass, K = sim$numclust)
#' true_res <- list(prob = sim$prob, dens_true = sim$dens_true)
#'
#' # Compute metric
#' dist <- mixture_metric(Y_bin, X, bin_mass, fit$iter, true_res)
#' plot(dist$per_time, type = "b", xlab = "Time", ylab = "L1 distance")
#' print(dist$weighted)
#' }
#' @export
mixture_metric <- function(
  sim,
  est_res
) {
  Y_bin     = sim$Y_bin
  X         = sim$X
  bin_mass  = sim$bin_mass
  TT        = length(Y_bin)
  prob_true = pi_k(X, t(sim$alpha))
  
  if (sim$noisetype == 'skewed') {
    dens_true_fun <- function(t, y_grid) {
      sapply(seq_len(K), function(k) {
        mu = sim$mnmat[t,k]
        sn::dsn(y_grid, xi = -sim$mn_shift, omega = sim$omega, alpha = sim$skew_alpha)
      })
    }
  }
  
  if (!is.null(est_res$alpha_new)) { # if LCDmix
    pi_est <- pi_k(X, est_res$alpha_new)
    K <- length(est_res$g_new)
    dens_est_fun <- function(t, y_grid) {
      sapply(seq_len(K), function(k) {
        mu <- est_res$theta0_new[[k]] + sum(X[t,] * est_res$theta_new[[k]])
        logcondens::evaluateLogConDens(y_grid - mu, est_res$g_new[[k]])[,3]
      })
    }
  } else if (!is.null(est_res$alpha)) { # if flowmix
    pi_est <- pi_k(X, est_res$alpha)
    mn_arr <- est_res$mn
    sigma  <- as.numeric(est_res$sigma)
    K      <- dim(mn_arr)[3]
    dens_est_fun <- function(t, y_grid) {
      sapply(seq_len(K), function(k) {
        dnorm(y_grid, mean = mn_arr[t,1,k], sd = sqrt(sigma[k]))
      })
    }
  } else {
    stop("`est_res` must contain either `alpha_new` (LCDmix) or `alpha` (flowmix).")
  }
  
  per_time <- numeric(TT)
  for (t in seq_len(TT)) {
    dens_est  <- dens_est_fun(t, Y_bin[[t]])
    mix_est   <- dens_est %*% pi_est[t, ]
    dens_true <- dens_true_fun(t, Y_bin[[t]])
    mix_true  <- dens_true %*% prob_true[t, ]
    per_time[t] <- sum(abs(mix_est - mix_true))
  }
  w_t      <- vapply(bin_mass, sum, numeric(1))
  metric   <- sum(w_t * per_time) / sum(w_t)
  
  return(metric)
}
