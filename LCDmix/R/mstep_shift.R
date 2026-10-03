# Generated from create-LCDmix.Rmd: do not edit by hand

#' Update the component intercepts (centering step)
#'
#' @description
#' For each component \eqn{k}, sets the intercept to the weighted mean of
#' \eqn{y - X_t^\top \theta_k} over the active bins:
#' \deqn{\theta_{0k} = \frac{\sum_t \sum_{i \in I_{tk}} w_{tik} (y_{ti} - X_t^\top \theta_k)}{\sum_t \sum_{i \in I_{tk}} w_{tik}}.}
#' Afterwards the weighted mean residual of each component is 0.
#'
#' This is the exact maximizer of \eqn{Q} over the intercept only when the
#' component density is Gaussian. For a general log-concave density it centers
#' the residuals but does not maximize \eqn{Q}, which is why the algorithm is a
#' generalized EM algorithm and not an ECM algorithm. \code{iteration()}
#' checks the ascent of the whole iteration.
#'
#' @param Y_bin A list of length \code{TT}; each element is an
#'   \eqn{M_t \times 1} matrix of bin centers at time \eqn{t}.
#' @param X A numeric \eqn{TT \times p} covariate matrix.
#' @param weights A list of length \code{TT}; each element is an
#'   \eqn{M_t \times K} matrix of posterior weights (responsibility times bin
#'   weight), as returned in \code{weight} by \code{estep_lcd()}.
#' @param idx A list of length \code{TT}; each element is an
#'   \eqn{M_t \times K} logical matrix; \code{idx[[t]][i, k]} is \code{TRUE}
#'   if bin \eqn{i} at time \eqn{t} is active for component \eqn{k}.
#' @param slopes A list of length \code{K}; each element is the slope vector
#'   \eqn{\theta_k} of length \eqn{p}.
#'
#' @return A list of length \code{K}; element \eqn{k} is the new intercept
#'   \eqn{\theta_{0k}}. It is \code{NaN} if component \eqn{k} has no active bin.
#'
#' @examples
#' \dontrun{
#' TT <- 3; K <- 2; p <- 2; n_bins <- 5
#' Y_bin   <- lapply(1:TT, function(t) matrix(rnorm(n_bins), ncol = 1))
#' X       <- matrix(rnorm(TT * p), nrow = TT, ncol = p)
#' weights <- lapply(Y_bin, function(m) matrix(runif(nrow(m) * K), ncol = K))
#' idx     <- lapply(weights, function(w) w > 0.1)
#' slopes  <- replicate(K, runif(p), simplify = FALSE)
#' intercepts <- mstep_shift(Y_bin, X, weights, idx, slopes)
#' }
#' @export
mstep_shift <- function(
  Y_bin,
  X,
  weights,
  idx,
  slopes
) {
  K_comp <- length(slopes)
  TT     <- length(Y_bin)
  intercepts <- vector("list", K_comp)

  for (k in seq_len(K_comp)) {
    num <- 0
    den <- 0
    for (t in seq_len(TT)) {
      # select which bins at time t belong to component k
      idx_tk <- idx[[t]][, k]
      # corresponding biomass weights for those bins
      w_tk   <- weights[[t]][idx_tk, k]
      # sum of binned responses in those bins
      resp_sum <- as.numeric(w_tk %*% Y_bin[[t]][idx_tk, , drop = FALSE])
      # adjustment from current slopes
      slope_adj <- sum(w_tk) * sum(X[t, ] * slopes[[k]])
      num <- num + (resp_sum - slope_adj)
      den <- den + sum(w_tk)
    }
    intercepts[[k]] <- num / den
  }

  return(intercepts)
}
