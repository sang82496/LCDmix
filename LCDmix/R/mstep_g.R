# Generated from create-LCDmix.Rmd: do not edit by hand

#' Update the component densities (log-concave M-step)
#'
#' @description
#' For each component \eqn{k}:
#'
#' 1. Pools the residuals of the bins that are active for \eqn{k} over all
#'    time points, with their posterior weights.
#' 2. Merges residuals closer than \code{dedup_tol} and sums their weights.
#'    The smallest and the largest residual are kept exactly, because the
#'    support of the estimate is their range.
#' 3. Normalizes the weights to sum to 1 and fits the weighted log-concave
#'    maximum likelihood estimate with \code{modified_logcondens()}.
#'
#' Stops with an error whose message starts with "degenerate component" when
#' a component has fewer than two distinct residuals; \code{count_degenerate()}
#' counts these failures in saved cross-validation results. A message is
#' printed when a component has fewer than five distinct residuals.
#'
#' @param residuals A list of length \code{TT}; each element is an
#'   \eqn{M_t \times K} matrix of residuals.
#' @param weights A list of length \code{TT}; each element is an
#'   \eqn{M_t \times K} matrix of posterior weights.
#' @param idx A list of length \code{TT}; each element is an
#'   \eqn{M_t \times K} logical matrix of active bins.
#' @param dedup_tol Residuals closer together than this are merged before
#'   fitting. Residuals that differ only by rounding error stop the active-set
#'   search of \code{modified_logcondens()} before it reaches the maximizer.
#'   Use \code{0} to switch merging off.
#'
#' @return A list of length \code{K}; element \eqn{k} is the
#'   \code{modified_logcondens()} fit for component \eqn{k}. Its support is the
#'   range of the pooled residuals of component \eqn{k}.
#'
#' @examples
#' \dontrun{
#' TT <- 3; K <- 2
#' # Residuals and weights for 5 bins and 2 components
#' residuals <- lapply(1:TT, function(t) matrix(rnorm(5 * K), ncol = K))
#' weights   <- lapply(residuals, function(m) abs(m))  # for the example only
#' # All bins active for both components
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
   # NEW: a component that owns no bins cannot be fitted
   if (length(res_k) == 0L) {
     stop(sprintf("degenerate component: component %d has %d distinct residual(s) (%d before merging); at least 2 are needed",
                  k, 0L, 0L), call. = FALSE)
   }
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
    
   # NEW: a log-concave MLE needs at least two distinct points; with one point
   # the support has zero width and modified_logcondens() fails at tmp[, 1]
   if (length(uniq_res) <= 1L) {
     stop(sprintf("degenerate component: component %d has %d distinct residual(s) (%d before merging); at least 2 are needed",
                  k, length(uniq_res), length(res_k)), call. = FALSE)
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
