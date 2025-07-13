# Generated from create-LCDmix.Rmd: do not edit by hand

#' Compute posterior mixing probabilities π_{t,k} = P(Z^{(t)} = k | X^{(t)})
#'
#' @description
#' Given a sequence of covariate vectors over \code{TT} time points and a parameter
#' matrix \code{alpha} of size \eqn{K \times (p+1)}, this function computes the
#' \eqn{TT \times K} matrix of posterior probabilities
#' \eqn{\pi_{t,k} = P(Z^{(t)} = k \mid X^{(t)})} via a row-wise softmax of the linear predictors.
#'
#' @param X A numeric matrix with \code{TT} rows (time points) and \code{p} covariate columns.
#' @param alpha A numeric matrix of dimension \eqn{K \times (p+1)}, where each row
#'   corresponds to one mixture component \code{k}. Column 1 is the intercept,
#'   and columns 2 through \eqn{p+1} are the slopes for that component.
#'
#' @return A numeric \eqn{TT \times K} matrix \code{pi_mat}, where
#'   \code{pi_mat[t, k]} is the posterior probability of component \code{k}
#'   at time point \code{t}. Each row of \code{pi_mat} sums to 1.
#'
#' @examples
#' \dontrun{
#' # Simulate 3 time points, 2 covariates, 2 mixture components
#' X     <- matrix(rnorm(3 * 2), nrow = 3, ncol = 2)
#' alpha <- matrix(
#'   c(0.1,  1.0, -0.5,   # component 1: intercept = 0.1, slopes = (1, -0.5)
#'     0.2, -1.0,  0.3),  # component 2: intercept = 0.2, slopes = (-1, 0.3)
#'   nrow = 2, byrow = TRUE
#' )
#' pi_mat <- pi_k(X, alpha)
#' stopifnot(all.equal(rowSums(pi_mat), rep(1, nrow(pi_mat))))
#' }
#'
#' @export
pi_k <- function(
  X,
  alpha
) {
  n_time <- nrow(X)
  p_cov  <- ncol(X)
  K_comp <- nrow(alpha)

  # Split alpha into intercepts (length K) and slopes (K x p)
  intercepts <- alpha[, 1]
  slopes     <- alpha[, -1, drop = FALSE]

  # Build intercept matrix: repeat intercepts for each time point
  intercept_mat <- matrix(
    intercepts,
    nrow = n_time,
    ncol = K_comp,
    byrow = TRUE
  )

  # Linear predictor: TT x K
  lin_pred <- intercept_mat + X %*% t(slopes)

  # Softmax: exp and normalize row-wise
  exp_lp <- exp(lin_pred)
  pi_mat <- exp_lp / rowSums(exp_lp)

  return(pi_mat)
}
