# Generated from create-LCDmix.Rmd: do not edit by hand

#' Weighted quantile with linear interpolation
#'
#' @description
#' Returns the \code{prob}-quantile of \code{x} with weights \code{w}. Points
#' with \code{NA} in \code{x}, or with a weight that is 0, negative or not
#' finite, are dropped. After sorting, let \eqn{C_j} be the cumulative weight
#' and \eqn{c} = \code{prob} times the total weight, and let \eqn{j} be the
#' first index with \eqn{C_j \ge c}. The result is
#' \deqn{x_{(j-1)} + \frac{c - C_{j-1}}{w_{(j)}} (x_{(j)} - x_{(j-1)}),}
#' kept within \eqn{[x_{(j-1)}, x_{(j)}]}. If \eqn{j = 1}, the smallest
#' \code{x} is returned; if a neighbor is not finite (for example
#' \code{-Inf}), \eqn{x_{(j)}} is returned without interpolation.
#' \code{eval_lcd()} uses this function for the median and for the trimming
#' threshold.
#'
#' @param x Numeric vector of values; \code{-Inf} is allowed.
#' @param w Numeric vector of weights of the same length as \code{x}.
#' @param prob The quantile level, a number in \eqn{[0, 1]}.
#'
#' @return A single number, or \code{NA} if no point is left after dropping.
#'
#' @examples
#' \dontrun{
#' x <- c(10, 20, 30, 40, 50)
#' w <- c(1, 1, 1, 1, 6)                 # more weight on the largest value
#' weighted_quantile(x, w, prob = 0.05)  # 10
#' weighted_quantile(x, w, prob = 0.5)   # 41.67
#' weighted_quantile(x, w, prob = 0.95)  # 49.17
#' }
#' @export
weighted_quantile <- function(
  x,
  w,
  prob = 0.05
) {
  ok <- !is.na(x) & is.finite(w) & w > 0          # zero-weight points cannot move a quantile
  x  <- x[ok]; w <- w[ok]
  if (!length(x)) return(NA_real_)

  o        <- order(x)
  sorted_x <- x[o]
  sorted_w <- w[o]
  cum_w     <- cumsum(sorted_w)
  threshold <- prob * cum_w[length(cum_w)]

  idx <- which(cum_w >= threshold)[1]
  if (idx == 1) return(sorted_x[1])                         # was: return(idx)
  lo <- sorted_x[idx - 1]; hi <- sorted_x[idx]
  if (!is.finite(lo) || !is.finite(hi)) return(hi)          # -Inf neighbor: no interpolation
  frac <- (threshold - cum_w[idx - 1]) / sorted_w[idx]      # was: / (sorted_w[idx] - sorted_w[idx-1])
  min(hi, max(lo, lo + frac * (hi - lo)))                   # clamp: floating overshoot past hi
}
