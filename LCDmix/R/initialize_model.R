# Generated from create-LCDmix.Rmd: do not edit by hand

#' Initialize mixture-of-experts model parameters via Gaussian mixture regression
#'
#' @description
#' Performs the initialization for the mixture-of-experts model pipeline:
#' 1. Fit a Gaussian mixture regression (GMR) via \code{flowmix::flowmix}.  
#' 2. Extract initial mixture parameters (\code{alpha_init}, \code{theta0_init}, \code{theta_init}).  
#' 3. Compute initial residuals, responsibilities, and posterior weights.  
#' 4. Perform one M-step update for intercepts.  
#' 5. Estimate initial log-concave component densities.  
#' 6. Compute the initial surrogate log-likelihood.
#'
#' @param Y_bin A list of length \code{TT}; each element is an \eqn{M_t \times 1} matrix of binned responses.
#' @param X A numeric matrix of dimension \eqn{TT \times p}; covariate values for each time point.
#' @param bin_mass A list of length \code{TT}; each element is a numeric vector of length \eqn{M_t} giving sum of biomass per bin.
#' @param K Integer; number of mixture components.
#' @param lambda_alpha Numeric; L1 regularization parameter for mixture weights in GMR.
#' @param lambda_theta Numeric; L1 regularization parameter for regression slopes in GMR.
#' @param resp_threshold Numeric; threshold on responsibilities for bin assignment.
#' @param maxdev Numeric or \code{NULL}; maximum deviance threshold for \code{flowmix}.
#' @param n_restarts Integer; number of random restarts in \code{flowmix::flowmix}.
#'
#' @return A list with components:
#' \describe{
#'   \item{flow}{Result object returned by \code{flowmix::flowmix}.}
#'   \item{idx_init}{List of length \code{TT} of \eqn{M_t \times K} logical matrices indicating high-confidence bin assignments.}
#'   \item{resp_init}{List of length \code{TT} of \eqn{M_t \times K} matrices of posterior probabilities (responsibilities).}
#'   \item{weight_init}{List of length \code{TT} of \eqn{M_t \times K} matrices of posterior weights (responsibilities × \code{bin_mass}).}
#'   \item{theta0_init}{List of length \code{K} of initial intercepts.}
#'   \item{theta_init}{List of length \code{K} of initial slope vectors.}
#'   \item{alpha_init}{\eqn{K \times (p+1)} matrix of initial mixture parameters from GMR.}
#'   \item{resi_init}{List of length \code{TT} of \eqn{M_t \times K} matrices of residuals.}
#'   \item{g_init}{List of length \code{K} of log-concave density objects from \code{modified_logcondens()}.}
#'   \item{Q}{Numeric; value of the surrogate log-likelihood at initialization.}
#'   \item{Q_every}{Numeric; same as \code{Q} (for consistency with iterative output).}
#' }
#'
#' @examples
#' \dontrun{
#' TT      <- 4; K <- 2; p <- 3; n_bins <- 5
#' Y_bin   <- lapply(1:TT, function(t) matrix(rnorm(n_bins), ncol = 1))
#' bin_mass<- lapply(Y_bin, function(y) runif(nrow(y)))
#' X       <- matrix(rnorm(TT * p), nrow = TT, ncol = p)
#' init    <- initialize_model(
#'   Y_bin, X, bin_mass, K,
#'   lambda_alpha   = 1e-3,
#'   lambda_theta   = 1e-3,
#'   resp_threshold = 1e-3,
#'   maxdev         = NULL,
#'   n_restarts     = 1
#' )
#' str(init)
#' }
#' @export
initialize_model <- function(
  Y_bin,
  X,
  bin_mass,
  K,
  lambda_alpha,
  lambda_theta,
  resp_threshold,
  maxdev,
  n_restarts
) {
  # 1) Fit Gaussian mixture regression via flowmix
  flow_res <- flowmix::flowmix(
    ylist        = Y_bin,
    X            = X,
    countslist   = bin_mass,
    numclust     = K,
    prob_lambda  = lambda_alpha,
    mean_lambda  = lambda_theta,
    maxdev       = maxdev,
    nrep         = n_restarts
  )
  
  message("✔ flowmix initialization complete")
  
  # 2) Extract initial parameters
  alpha_init  <- flow_res$alpha
  theta0_init <- lapply(flow_res$beta, `[[`, 1)
  theta_init  <- lapply(flow_res$beta, function(b) b[-1])
  
  # 3) Compute initial residuals
  resi_init <- compute_residuals(Y_bin, X, theta0_init, theta_init)
  
  # 4) Compute initial responsibilities and weights
  TT           <- length(Y_bin)
  resp_init    <- vector("list", TT)
  idx_init     <- vector("list", TT)
  weight_init  <- vector("list", TT)
  pi_mat       <- pi_k(X, alpha_init)  # TT x K matrix
  
  for (t in seq_len(TT)) {
    M_t      <- nrow(Y_bin[[t]])
    likeli_t <- matrix(0, nrow = M_t, ncol = K)
    for (k in seq_len(K)) {
      likeli_t[, k] <- dnorm(
        x    = resi_init[[t]][, k],
        mean = flow_res$mn[t, 1, k],
        sd   = sqrt(flow_res$sigma[k])
      ) * pi_mat[t, k]
    }
    resp_t        <- likeli_t / rowSums(likeli_t)
    resp_init[[t]]   <- resp_t
    idx_init[[t]]    <- resp_t > resp_threshold
    weight_init[[t]] <- resp_t * bin_mass[[t]]
  }
  
  # 5) One M-step update for intercepts
  theta0_init <- mstep_update_intercepts(Y_bin, X, weight_init, idx_init, theta_init)
  resi_init   <- compute_residuals(Y_bin, X, theta0_init, theta_init)
  
  # 6) Initial log-concave density estimates
  g_init <- mstep_estimate_log_concave_densities(resi_init, weight_init, idx_init)
  
  # 7) Compute initial surrogate log-likelihood
  Q_init <- compute_surrogate_loglikelihood(
    X, g_init, resi_init, theta_init,
    alpha_init, idx_init, weight_init,
    lambda_alpha, lambda_theta
  )
  names(Q_init) = 'Q'
  
  # Return initialization results
  return(list(
    flow        = flow_res,
    idx_init    = idx_init,
    resp_init   = resp_init,
    weight_init = weight_init,
    theta0_init = theta0_init,
    theta_init  = theta_init,
    alpha_init  = alpha_init,
    resi_init   = resi_init,
    g_init      = g_init,
    Q           = Q_init,
    Q_every     = Q_init
  ))
}
