# Generated from create-LCDmix.Rmd: do not edit by hand

#' Perform the E-step using log-concave component densities
#'
#' @description
#' Given current model parameters and log-concave density estimates, computes:
#'   - Soft responsibilities \(\mathrm{resp}_{t,i,k}\)  
#'   - Posterior weights (responsibilities × biomass per bin)  
#'   - Logical indices of “effectively nonzero” responsibilities  
#' for each time point \(t\) and component \(k\), thresholding responsibilities
#' below \code{resp_threshold} to zero for numerical stability and speed.
#'
#' @param X A numeric \eqn{TT \times p} matrix of covariates (rows = time points).
#' @param bin_mass A list of length \code{TT}, where each element is a numeric
#'   vector of length \eqn{M_t} giving the biomass weight per bin at time \(t\).
#' @param residuals A list of length \code{TT}; each element is an \eqn{M_t \times K}
#'   matrix of residuals for bins at time \(t\).
#' @param alpha A numeric \eqn{K \times (p+1)} matrix of mixture parameters
#'   (column 1 = intercepts, columns 2:\(p+1\) = slopes).
#' @param densities A list of length \code{K}, each an object returned by
#'   \code{modified_logcondens()} representing the log‐concave density for one component.
#' @param resp_threshold Numeric in \([0,1]\); responsibilities below this value
#'   are set to zero. Default: \code{1e-3}.
#'
#' @return A list with components:
#' \describe{
#'   \item{resp}{List of length \code{TT} of \eqn{M_t \times K} matrices of responsibilities.}
#'   \item{weight}{List of length \code{TT} of \eqn{M_t \times K} matrices of posterior weights.}
#'   \item{idx}{List of length \code{TT} of \eqn{M_t \times K} logical matrices, where
#'       \code{idx[[t]][i,k]} = \code{TRUE} if \code{resp[[t]][i,k]} > 0.}
#' }
#'
#' @examples
#' \dontrun{
#' TT <- 4; K <- 3; p <- 2
#' # Simulate residuals (5 bins × K) and biomass
#' residuals <- lapply(1:TT, function(t) matrix(rnorm(5*K), ncol = K))
#' bin_mass  <- lapply(residuals, function(m) runif(nrow(m)))
#' # Dummy densities
#' densities <- replicate(K,
#'   modified_logcondens(rnorm(100), w = rep(1/100, 100)),
#'   simplify = FALSE
#' )
#' # Random mixture parameters
#' alpha <- matrix(rnorm(K*(p+1)), nrow = K)
#' # One E‐step
#' e_res <- estep_lcd(
#'   X              = matrix(rnorm(TT*p), nrow = TT),
#'   bin_mass       = bin_mass,
#'   residuals      = residuals,
#'   alpha          = alpha,
#'   densities      = densities,
#'   resp_threshold = 1e-3
#' )
#' str(e_res)
#' }
#' @export
estep_lcd <- function(
  X,
  bin_mass,
  residuals,
  alpha,
  densities,
  resp_threshold = 1e-3
) {
  TT      <- nrow(X)
  K_comp  <- length(densities)
  # Precompute mixing probabilities
  pi_mat  <- pi_k(X, alpha)
  
  resp_list   <- vector("list", TT)
  weight_list <- vector("list", TT)
  idx_list    <- vector("list", TT)
  
  for (t in seq_len(TT)) {
    M_t     <- nrow(residuals[[t]])
    lik_mat <- matrix(0, nrow = M_t, ncol = K_comp)
    
    # Compute (log‐density × mixing probability) for each component
    for (k in seq_len(K_comp)) {
      log_dens <- suppressWarnings(
        logcondens::evaluateLogConDens(
          residuals[[t]][, k],
          densities[[k]]
        )[, 3]
      )
      lik_mat[, k] <- log_dens * pi_mat[t, k]
    }
    
    # Normalize to get soft responsibilities
    resp_t <- lik_mat / rowSums(lik_mat)
    # Threshold small probabilities
    resp_t[resp_t < resp_threshold] <- 0
    
    # Indices of effectively nonzero responsibilities
    idx_t <- resp_t > 0
    # Posterior weights = responsibilities × biomass per bin
    weight_t <- resp_t * bin_mass[[t]]
    
    resp_list[[t]]   <- resp_t
    weight_list[[t]] <- weight_t
    idx_list[[t]]    <- idx_t
  }
  
  return(list(
    resp   = resp_list,
    weight = weight_list,
    idx    = idx_list
  ))
}
