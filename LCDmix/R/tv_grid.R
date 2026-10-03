# Generated from create-LCDmix.Rmd: do not edit by hand

#' Integration range for tv_density()
#'
#' @description
#' Takes the central range of the true distribution, from its
#' \code{p}-quantile to its \code{(1 - p)}-quantile, widens it to cover
#' \code{est_support} if that is given, and adds 2 percent of the width on
#' each side.
#'
#' @param qfun Quantile function of the true distribution.
#' @param est_support Optional \code{c(lower, upper)}, the support of the
#'   estimated density.
#' @param p Probability left out in each tail.
#'
#' @return \code{c(lo, hi)}.
#'
#' @export
tv_grid <- function(qfun, est_support = NULL, p = 5e-4) {
  lo <- qfun(p); hi <- qfun(1 - p)
  if (!is.null(est_support)) {
    lo <- min(lo, est_support[1]); hi <- max(hi, est_support[2])
  }
  pad <- 0.02 * (hi - lo)
  return(c(lo - pad, hi + pad))
}
