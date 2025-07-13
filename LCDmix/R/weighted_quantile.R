# Generated from create-LCDmix.Rmd: do not edit by hand

#' Compute the weighted quantile of a numeric sample
#'
#' @description
#' Returns the \code{prob}-th weighted quantile of the values \code{x} with
#' corresponding nonnegative weights \code{w}. The observations are sorted by
#' \code{x}, and the smallest \code{x} such that the cumulative weight
#' reaches \code{prob * sum(w)} is returned.
#'
#' @param x Numeric vector of data values.
#' @param w Numeric vector of nonnegative weights of the same length as \code{x}.
#' @param prob Numeric scalar in \[0,1\]; the desired quantile level (e.g.\ 0.05 for
#'   the 5th percentile). Default: \code{0.05}.
#'
#' @return Numeric scalar: the weighted \code{prob}-quantile of \code{x}.
#'
#' @examples
#' \dontrun{
#' x   <- c(10, 20, 30, 40, 50)
#' w   <- c(1, 1, 1, 1, 6)     # more weight on the largest value
#' q05 <- weighted_quantile(x, w, prob = 0.05)  # should be 10
#' q95 <- weighted_quantile(x, w, prob = 0.95)  # should be 50
#' }
#' @export
weighted_quantile <- function(
  x,
  w,
  prob = 0.05
) {
  #— Validate inputs —#
  if (!is.numeric(x) || !is.numeric(w)) {
    stop("`x` and `w` must both be numeric vectors")
  }
  if (length(x) != length(w)) {
    stop("`x` and `w` must have the same length")
  }
  if (any(w < 0)) {
    stop("`w` must be nonnegative")
  }
  if (!is.numeric(prob) || length(prob) != 1 || prob < 0 || prob > 1) {
    stop("`prob` must be a single numeric value between 0 and 1")
  }

  #— Sort by x —#
  o     <- order(x)
  x_ord <- x[o]
  w_ord <- w[o]

  #— Compute cumulative weights and threshold —#
  cum_w     <- cumsum(w_ord)
  total_w   <- sum(w_ord)
  threshold <- prob * total_w

  #— Find first index where cumulative weight ≥ threshold —#
  idx <- which(cum_w >= threshold)[1]

  #— Return corresponding quantile value —#
  return(x_ord[idx])
}
