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
  # Sort by x
  o     <- order(x)
  sorted_x <- x[o]
  sorted_w <- w[o]

  #  Compute cumulative weights and threshold
  cum_w     <- cumsum(sorted_w)
  threshold <- prob * sum(sorted_w)

  # Find first index where cumulative weight ≥ threshold
  idx <- which(cum_w >= threshold)[1]
  if (idx == 1) {return(idx)}
  diff_x = sorted_x[idx] - sorted_x[idx-1]
  diff_w = sorted_w[idx] - sorted_w[idx-1]
  frac = (cum_w[idx]-threshold)/diff_w
  return(sorted_x[idx] - frac*diff_x)
}
