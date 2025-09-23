# Generated from create-LCDmix.Rmd: do not edit by hand

#' Estimate component densities via log-concave density estimation (M-step)
#'
#' @description
#' For each mixture component \(k\), gathers the residuals across all time points
#' for bins assigned to component \(k\), weighs them by their posterior weights,
#' and fits a log-concave density using \code{modified_logcondens()}.
#'
#' @param residuals A list of length \code{TT}, where each element is an \eqn{M_t \times K}
#'   matrix of residuals at time point \code{t}, with \eqn{M_t \le n\_bins} bins.
#' @param weights A list of length \code{TT}, where each element is an \eqn{M_t \times K}
#'   matrix of posterior weights (e.g., responsibilities) corresponding to \code{residuals}.
#' @param idx A list of length \code{TT}, where each element is an \eqn{M_t \times K}
#'   logical or integer matrix.  \code{idx[[t]][i,k]} indicates whether the \(i\)th bin
#'   at time \(t\) contributes to component \(k\).
#'
#' @return A list of length \code{K}, where element \code{k} is the output of
#'   \code{modified_logcondens()}—the estimated log-concave density for component \(k\).
#'
#' @examples
#' \dontrun{
#' TT <- 3; K <- 2
#' # Simulate residuals and weights (5 bins × 2 components)
#' residuals <- lapply(1:TT, function(t) matrix(rnorm(5 * K), ncol = K))
#' weights   <- lapply(residuals, function(m) abs(m))  # just for demo
#' # Include all bins for both components
#' idx <- lapply(residuals, function(m) matrix(TRUE, nrow = nrow(m), ncol = ncol(m)))
#' densities <- mstep_g(residuals, weights, idx)
#' }
#' @export
mstep_g <- function(
  residuals,
  weights,
  idx
) {
  TT <- length(weights)
  K  <- ncol(idx[[1]])
  densities <- vector("list", K)
  
  for (k in seq_len(K)) {
    # Collect all residuals and weights for component k
    res_k <- unlist(lapply(seq_len(TT),
                           function(t) residuals[[t]][idx[[t]][, k], k]))
    w_k   <- unlist(lapply(seq_len(TT),
                           function(t) weights[[t]][idx[[t]][, k], k]))
    
    # Unique residual values & aggregated weights
    uniq_res <- unique(res_k)
    if (length(uniq_res) < 5) {
      message("Only ", length(uniq_res), " unique points for component ", k)
    }
    uniq_w   <- sapply(uniq_res, function(val) sum(w_k[res_k == val]))
    
    # Fit log-concave density
    densities[[k]] <- modified_logcondens(
      x     = uniq_res,
      w     = uniq_w / sum(uniq_w),
      print = FALSE
    )
  }
  
  return(densities)
}
