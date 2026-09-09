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
#' @param dedup_tol Numeric; residuals closer together than this are merged
#'   before fitting. Binning produces residuals differing only at machine
#'   precision, which prevents the log-concave active-set search from reaching
#'   its maximizer. Set to \code{0} to disable merging. Default \code{1e-10}.
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
  idx,
  dedup_tol = 1e-10   # NEW - appended. Set to 0 to restore the old behaviour.
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

    if (dedup_tol > 0) {
      # Merge residuals separated by less than dedup_tol. unique() treats values
      # ~1e-16 apart as distinct knot candidates, which prevents the active-set
      # search in modified_logcondens() from reaching its maximiser.
      ord <- order(res_k)
      rs  <- res_k[ord]
      ws  <- w_k[ord]
      grp <- cumsum(c(TRUE, diff(rs) > dedup_tol))
      uniq_res <- as.numeric(tapply(rs, grp, mean))
      uniq_w   <- as.numeric(tapply(ws, grp, sum))

      # ENDPOINT PRESERVATION -- do not drop these two lines.
      # The log-concave MLE's support is exactly [min(x), max(x)]. Taking the
      # cluster MEAN at the extremes moves the support endpoint inward, leaving
      # the true extreme residual outside its own component's support. estep_lcd()
      # then finds a bin with zero density under every component, divides 0/0,
      # and glmnet dies an iteration or two later. Traced on sim-61 cell 1-1-2-1:
      # the offending bin sat 8.882e-16 outside the boundary. Without these two
      # lines 5/40 grid cells fail; with them, 0/40.
      uniq_res[1]                <- min(rs)
      uniq_res[length(uniq_res)] <- max(rs)
    } else {
      # Original behaviour, kept for dedup_tol = 0
      uniq_res <- unique(res_k)
      uniq_w   <- sapply(uniq_res, function(val) sum(w_k[res_k == val]))
    }

    if (length(uniq_res) < 5) {
      message("Only ", length(uniq_res), " unique points for component ", k)
    }

    # Fit log-concave density
    densities[[k]] <- modified_logcondens(
      x     = uniq_res,
      w     = uniq_w / sum(uniq_w),
      print = FALSE
    )
  }

  return(densities)
}
