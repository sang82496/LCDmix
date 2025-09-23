# Generated from create-LCDmix.Rmd: do not edit by hand

#' M‐step update of all component regression parameters under log‐concavity
#'
#' @description
#' For each mixture component \(k\), performs the M‐step update of intercept \(\theta_{0k}\)
#' and slope vector \(\theta_k\) by invoking the LP solver with log‐concave constraints
#' (\code{mstep_update_theta_lp()}), returning updated lists of intercepts and slopes.
#'
#' @param Y_bin A list of length \code{TT}; each element is an \eqn{M_t \times 1} matrix of binned responses at time \(t\).
#' @param X A numeric \eqn{TT \times p} matrix of covariates (rows = time points).
#' @param posterior_weights A list of length \code{TT}; each element is an \eqn{M_t \times K} matrix of posterior weights.
#' @param residuals A list of length \code{TT}; each element is an \eqn{M_t \times K} matrix of residuals.
#' @param densities A list of length \code{K}, each an object returned by \code{modified_logcondens()} for one component.
#' @param responsibility_mask A list of length \code{TT}; each element is an \eqn{M_t \times K} logical matrix indicating effective bin membership.
#' @param intercepts A list of length \code{K} of current intercept parameters \(\theta_{0k}\).
#' @param slopes A list of length \code{K} of current slope vectors \(\theta_k\).
#' @param lambda_theta Nonnegative numeric L1 penalty on slopes.
#'
#' @return A list with elements:
#' \describe{
#'   \item{theta0}{List of length \code{K} of updated intercepts \(\theta_{0k}\).}
#'   \item{theta}{List of length \code{K} of updated slope vectors \(\theta_k\).}
#' }
#'
#' @examples
#' \dontrun{
#' TT      <- 3; K <- 2; p <- 2; n_bins <- 5
#' Y_bin   <- lapply(1:TT, function(t) matrix(rnorm(n_bins), ncol=1))
#' X       <- matrix(rnorm(TT*p), nrow=TT, ncol=p)
#' weights <- lapply(Y_bin, function(m) matrix(runif(length(m)*K), ncol=K))
#' resi    <- lapply(Y_bin, function(m) matrix(rnorm(length(m)*K), ncol=K))
#' mask    <- lapply(Y_bin, function(m) matrix(TRUE, nrow(m), ncol=K))
#' densities <- replicate(K,
#'   modified_logcondens(rnorm(50), w = rep(1/50,50)),
#'   simplify = FALSE
#' )
#' intercepts <- replicate(K, 0, simplify = FALSE)
#' slopes     <- replicate(K, rep(0, p), simplify = FALSE)
#' result <- mstep_update_theta_log_concave(
#'   Y_bin, X, weights, resi,
#'   densities, mask,
#'   intercepts, slopes,
#'   lambda_theta = 1e-3
#' )
#' str(result)
#' }
#' @export
mstep_update_theta_log_concave <- function(
  Y_bin,
  X,
  posterior_weights,
  residuals,
  densities,
  responsibility_mask,
  intercepts,
  slopes,
  lambda_theta
) {
  K <- length(densities)
  theta0_new <- vector("list", K)
  theta_new  <- vector("list", K)
  
  for (k in seq_len(K)) {
    tmp <- mstep_update_theta_lp(
      Y_bin               = Y_bin,
      X                   = X,
      posterior_weights   = posterior_weights,
      residuals           = residuals,
      density_k           = densities[[k]],
      responsibility_mask = responsibility_mask,
      intercept_k         = intercepts[[k]],
      slopes_k            = slopes[[k]],
      lambda_theta        = lambda_theta,
      component           = k
    )
    theta0_new[[k]] <- tmp$theta0_k
    theta_new[[k]]  <- tmp$theta_k
  }
  
  # Explicit return
  return(list(
    theta0 = theta0_new,
    theta  = theta_new
  ))
}
