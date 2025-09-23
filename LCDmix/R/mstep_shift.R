# Generated from create-LCDmix.Rmd: do not edit by hand

#' Update component intercepts in the M-step using binned data
#'
#' @description
#' Computes the updated intercepts \eqn{\theta_{0k}} for each mixture component \eqn{k}
#' by solving the weighted shift equations on binned responses. Here \code{Y_bin} and
#' \code{weights} both have at most \code{n_bins} rows per time point (empty bins removed).
#'
#' @param Y_bin A list of length \code{TT}; each element is an \eqn{M_t \times 1} matrix
#'   of binned response values for time point \eqn{t}, where \eqn{M_t \le n\_bins}.
#' @param X A numeric matrix of size \eqn{TT \times p}, where each row \code{X[t, ]} is
#'   the covariate vector at time \eqn{t}.
#' @param weights A list of length \code{TT}; each element is a numeric vector of length
#'   \eqn{M_t} giving the total biomass per bin (i.e.\ the \code{bin_mass} output
#'   from \code{binning()}, with empty bins removed).
#' @param idx A list of length \code{TT}; each element is an \eqn{M_t \times K}
#'   logical or integer matrix where \code{idx[[t]][i, k]} is \code{TRUE} if the
#'   \(i\)th bin at time \emph{t} contributes to component \emph{k}.
#' @param slopes A list of length \code{K}; each element is a numeric vector of length
#'   \eqn{p} giving the current slope parameters \eqn{\boldsymbol{\theta}_k}.
#'
#' @return A list of length \code{K}, where element \code{k} is the updated intercept
#'   \eqn{\theta_{0k}} for component \eqn{k}.
#'
#' @examples
#' \dontrun{
#' TT      <- 3; K <- 2; p <- 2; n_bins <- 5
#' # Simulate binned responses and bin_mass
#' Y_bin   <- lapply(1:TT, function(t) matrix(rnorm(n_bins), ncol = 1))
#' bin_mass<- lapply(Y_bin, function(m) runif(nrow(m)))
#' # Dummy X and idx
#' X       <- matrix(rnorm(TT * p), nrow = TT, ncol = p)
#' idx <- lapply(Y_bin, function(m) matrix(sample(c(TRUE,FALSE),
#'                       length(m) * K, replace = TRUE),
#'                       ncol = K))
#' slopes  <- replicate(K, runif(p), simplify = FALSE)
#' # Compute updated intercepts
#' intercepts <- mstep_shift(Y_bin, X, bin_mass, idx, slopes)
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
