# Generated from create-LCDmix.Rmd: do not edit by hand

#' Integration grid: the central range of the true density, widened to cover
#' the estimated support.
#' @export
tv_grid <- function(qfun, est_support = NULL, p = 5e-4) {
  lo <- qfun(p); hi <- qfun(1 - p)
  if (!is.null(est_support)) {
    lo <- min(lo, est_support[1]); hi <- max(hi, est_support[2])
  }
  pad <- 0.02 * (hi - lo)
  return(c(lo - pad, hi + pad))
}
