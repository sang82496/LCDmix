# Generated from create-LCDmix.Rmd: do not edit by hand

#' Compute a weighted histogram from a log‐concave density object
#'
#' @description
#' Given an object returned by \code{modified_logcondens()}, which contains the
#' sample points \code{xn} and their associated weights \code{w}, this function
#' constructs a normalized histogram (as a density estimate) over \code{n_bins}
#' equal‐width bins spanning the range of \code{xn}.
#'
#' @param g A list (output of \code{modified_logcondens()}) containing:
#'   \describe{
#'     \item{\code{xn}}{Numeric vector of sample points.}
#'     \item{\code{w}}{Numeric vector of nonnegative weights corresponding to each element of \code{xn}.}
#'   }
#' @param n_bins Integer; number of equal‐width bins to use. Default: \code{30}.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{midpoints}}{Numeric vector of length \code{n_bins}; the center of each bin.}
#'   \item{\code{densities}}{Numeric vector of length \code{n_bins}; the estimated density
#'     (normalized histogram) at each midpoint, such that the total area integrates to 1.}
#' }
#'
#' @examples
#' \dontrun{
#' # Suppose g_k is one component from mstep_g()
#' w_hist <- weighted_hist(
#'   g = g_k,
#'   n_bins         = 50
#' )
#' plot(w_hist$midpoints, w_hist$densities, type = "h",
#'      xlab = "Residual value", ylab = "Density")
#' }
#' @export
weighted_hist <- function(
  g,
  n_bins = 30
) {
  # Extract sample points and weights
  x_vals <- g$xn
  w_vals <- g$w

  # Determine range and bin cutpoints
  min_x   <- min(x_vals)
  max_x   <- max(x_vals)
  breaks  <- seq(from = min_x, to = max_x, length.out = n_bins + 1)

  # Assign each sample to a bin (1 through n_bins)
  bin_idx <- findInterval(x_vals, breaks, rightmost.closed = TRUE)
  bin_idx <- factor(bin_idx, levels = seq_len(n_bins))

  # Sum weights within each bin
  w_sum <- tapply(w_vals, bin_idx, sum)
  w_sum[is.na(w_sum)] <- 0

  # Normalize to form a density estimate: ensure area under histogram = 1
  densities <- w_sum * n_bins / (sum(w_sum) * (max_x - min_x))

  # Compute bin midpoints
  midpoints <- (breaks[-1] + breaks[-length(breaks)]) / 2

  return(list(
    midpoints = midpoints,
    densities = as.numeric(densities)
  ))
}
