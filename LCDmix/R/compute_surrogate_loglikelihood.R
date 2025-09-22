# Generated from create-LCDmix.Rmd: do not edit by hand

#' Compute the surrogate log-likelihood Q for the EM algorithm
#'
#' @description
#' Calculates the expected complete-data log-likelihood (the “Q” function) given
#' current parameter estimates, combining:
#' - Log-concave density contributions of residuals  
#' - Mixture weight (α) contributions  
#' - L1 penalties on α (excluding intercepts) and on regression slopes  
#'
#' Formally,
#' \deqn{Q = \frac{1}{N}\sum_{t=1}^{TT}\sum_{k=1}^{K}\sum_{i\in\mathcal{I}_{tk}}
#'   w_{t,i,k}\bigl[\log f_k(r_{t,i,k}) + \log\pi_{t,k}\bigr]
#'   - \lambda_\alpha\|\alpha_{-0}\|_1 - \lambda_\theta\|\theta\|_1,}
#' where:
#'   - \(r_{t,i,k}\) are residuals,  
#'   - \(f_k\) is the log-concave density for component \(k\),  
#'   - \(\pi_{t,k}=P(Z^{(t)}=k\mid X^{(t)})\),  
#'   - \(w_{t,i,k}\) are responsibilities,  
#'   - \(N=\sum_{t,i,k}w_{t,i,k}\).  
#'
#' @param X A \eqn{TT \times p} numeric covariate matrix (rows = time points).
#' @param densities A list of length \eqn{K}, each object returned by \code{modified_logcondens()}.
#' @param residuals A list of length \eqn{TT}, each an \eqn{M_t \times K} matrix of residuals.
#' @param slopes A list of length \eqn{K}, each a numeric vector of length \eqn{p} of slope parameters.
#' @param alpha A \eqn{K \times (p+1)} numeric matrix of mixture parameters (column 1 = intercepts).
#' @param indices A list of length \eqn{TT}, each an \eqn{M_t \times K} logical/integer matrix indicating which residuals belong to component \eqn{k}.
#' @param responsibilities A list of length \eqn{TT}, each an \eqn{M_t \times K} matrix of responsibilities \eqn{w_{t,i,k}}.
#' @param lambda_alpha Nonnegative numeric L1 penalty on non-intercept columns of \code{alpha}.
#' @param lambda_theta Nonnegative numeric L1 penalty on slope parameters \code{slopes}.
#'
#' @return A single numeric: the surrogate log-likelihood \eqn{Q} (normalized) minus the L1 penalties.
#'
#' @examples
#' \dontrun{
#' TT <- 3; p <- 2; K <- 2
#' X  <- matrix(rnorm(TT * p), nrow = TT)
#' densities <- replicate(K,
#'   modified_logcondens(rnorm(50), w = rep(1/50,50)),
#'   simplify = FALSE
#' )
#' residuals <- lapply(1:TT, function(t) matrix(rnorm(5 * K), ncol = K))
#' indices   <- lapply(residuals,  function(m) matrix(TRUE, nrow(m), ncol(m)))
#' responsibilities <- lapply(indices,
#'   function(idx) matrix(runif(length(idx)), nrow = nrow(idx))
#' )
#' slopes <- replicate(K, rnorm(p), simplify = FALSE)
#' alpha  <- matrix(rnorm(K*(p+1)), nrow = K)
#'
#' Q_val <- compute_surrogate_loglikelihood(
#'   X, densities, residuals, slopes, alpha,
#'   indices, responsibilities,
#'   lambda_alpha = 1e-3, lambda_theta = 1e-3
#' )
#' }
#' @export
compute_surrogate_loglikelihood <- function(
  X,
  densities,
  residuals,
  slopes,
  alpha,
  indices,
  responsibilities,
  lambda_alpha,
  lambda_theta
) {
  # Number of time points
  TT     <- nrow(X)
  K_comp <- length(densities)

  # Total “effective sample size”
  N_total <- sum(unlist(responsibilities))

  # Mixture probabilities: TT x K matrix
  pi_mat <- pi_k(X, alpha)

  total_ll <- 0

  for (k in seq_len(K_comp)) {
    for (t in seq_len(TT)) {
      # Select bins for component k at time t
      idx_tk    <- indices[[t]][, k]
      w_tk      <- responsibilities[[t]][idx_tk, k]
      resi_tk   <- residuals[[t]][idx_tk, k]

      if (length(w_tk) == 0) next

      # Evaluate log-density for each residual under component k
      log_dens <- suppressWarnings(
        logcondens::evaluateLogConDens(resi_tk, densities[[k]])[, 2]
      )
      finite_mask <- is.finite(log_dens)

      # Accumulate weighted log-density and mixture log-prob terms
      total_ll <- total_ll +
        sum(w_tk[finite_mask] * log_dens[finite_mask]) +
        sum(w_tk[finite_mask]) * log(pi_mat[t, k])
    }
  }

  # Normalize and subtract L1 penalties
  l1_alpha  <- sum(abs(alpha[, -1]))
  l1_theta  <- sum(abs(unlist(slopes)))
  pen_total <- lambda_alpha * l1_alpha + lambda_theta * l1_theta

  Q_val <- total_ll / N_total - pen_total

  return(Q_val)
}
