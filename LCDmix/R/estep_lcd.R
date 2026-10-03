# Generated from create-LCDmix.Rmd: do not edit by hand

#' E-step with log-concave component densities
#'
#' @description
#' Computes, for each time point \eqn{t}, bin \eqn{i} and component \eqn{k},
#' the responsibility
#' \deqn{r_{tik} = \frac{\pi_{tk} \hat f_k(u_{tik})}{\sum_l \pi_{tl} \hat f_l(u_{til})},}
#' with each density evaluated at the residual of the bin center.
#' Responsibilities below \code{resp_threshold} are set to 0. The others are
#' not renormalized, so a row can sum to slightly less than 1. A bin outside
#' the support of every component gives 0/0, that is \code{NaN}.
#'
#' @param X A numeric \eqn{TT \times p} covariate matrix.
#' @param bin_mass A list of length \code{TT}; each element is a numeric vector
#'   of length \eqn{M_t}, the weight of each bin.
#' @param residuals A list of length \code{TT}; each element is an
#'   \eqn{M_t \times K} matrix of residuals.
#' @param alpha A numeric \eqn{K \times (p+1)} matrix of gate coefficients
#'   (column 1 holds the intercepts).
#' @param densities A list of length \code{K} of \code{modified_logcondens()}
#'   fits.
#' @param resp_threshold A number in \eqn{[0, 1]}; responsibilities below it
#'   are set to 0.
#'
#' @return A list with components:
#' \describe{
#'   \item{resp}{List of \code{TT} \eqn{M_t \times K} matrices of
#'     responsibilities, after thresholding.}
#'   \item{weight}{List of \code{TT} \eqn{M_t \times K} matrices of posterior
#'     weights, \code{resp} times \code{bin_mass}.}
#'   \item{idx}{List of \code{TT} \eqn{M_t \times K} logical matrices;
#'     \code{TRUE} where \code{resp} is positive.}
#' }
#'
#' @examples
#' \dontrun{
#' TT <- 4; K <- 3; p <- 2
#' # Residuals (5 bins, K components) and bin weights
#' residuals <- lapply(1:TT, function(t) matrix(rnorm(5 * K), ncol = K))
#' bin_mass  <- lapply(residuals, function(m) runif(nrow(m)))
#' densities <- replicate(K,
#'   modified_logcondens(rnorm(100), w = rep(1/100, 100)),
#'   simplify = FALSE
#' )
#' alpha <- matrix(rnorm(K * (p + 1)), nrow = K)
#' e_res <- estep_lcd(
#'   X              = matrix(rnorm(TT * p), nrow = TT),
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
      dens <- suppressWarnings(
        logcondens::evaluateLogConDens(
          residuals[[t]][, k],
          densities[[k]],
          which = 2 # fixT: which = 2 fills column 3 (density) only
        )[, 3]
      )
      lik_mat[, k] <- dens * pi_mat[t, k]
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
