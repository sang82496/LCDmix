# Generated from create-LCDmix.Rmd: do not edit by hand

#' Weighted histogram of the data behind a log-concave fit
#'
#' @description
#' Builds a histogram density of the weighted points that a
#' \code{modified_logcondens()} fit was computed from (for a component of
#' \code{mstep_g()}: the pooled residuals and their weights), on
#' \code{n_bins} equal-width bins over their range. It uses \code{g$xn} and
#' \code{g$w}, which are in the same order only when the fit was made with
#' weights given, as \code{mstep_g()} does.
#'
#' @param g A \code{modified_logcondens()} fit, with fields \code{xn} (points)
#'   and \code{w} (their weights).
#' @param n_bins Number of equal-width bins.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{midpoints}}{Numeric vector of length \code{n_bins}; the bin
#'     centers.}
#'   \item{\code{densities}}{Numeric vector of length \code{n_bins}; the
#'     histogram heights, scaled so that the total area is 1.}
#' }
#'
#' @examples
#' \dontrun{
#' # g_k is one component of mstep_g()
#' w_hist <- weighted_hist(g = g_k, n_bins = 50)
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
