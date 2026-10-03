# Generated from create-LCDmix.Rmd: do not edit by hand

#' Gate probabilities of the mixture components
#'
#' @description
#' Computes, for every time point \eqn{t} and component \eqn{k},
#' \deqn{\pi_{tk} = P(Z_t = k \mid X_t) = \frac{\exp(\alpha_{k0} + X_t^\top \alpha_k)}{\sum_{l=1}^K \exp(\alpha_{l0} + X_t^\top \alpha_l)}.}
#' These are the gate (prior) probabilities: they depend on the covariates
#' only, not on the response.
#'
#' @param X A numeric matrix with \code{TT} rows (time points) and \code{p}
#'   covariate columns.
#' @param alpha A numeric \eqn{K \times (p+1)} matrix; row \eqn{k} belongs to
#'   component \eqn{k}. Column 1 is the intercept and columns 2 to
#'   \eqn{p+1} are the slopes.
#'
#' @return A numeric \eqn{TT \times K} matrix; entry \code{[t, k]} is
#'   \eqn{\pi_{tk}}. Each row sums to 1.
#'
#' @examples
#' \dontrun{
#' # 3 time points, 2 covariates, 2 components
#' X     <- matrix(rnorm(3 * 2), nrow = 3, ncol = 2)
#' alpha <- matrix(
#'   c(0.1,  1.0, -0.5,   # component 1: intercept 0.1, slopes (1, -0.5)
#'     0.2, -1.0,  0.3),  # component 2: intercept 0.2, slopes (-1, 0.3)
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
