# Generated from create-LCDmix.Rmd: do not edit by hand

#' Total variation distance between two densities on an interval
#'
#' @description
#' Computes \eqn{\frac{1}{2} \int_{lo}^{hi} |\hat f(x) - f(x)|\, dx} with the
#' trapezoid rule on \code{n} equally spaced points. The value lies in
#' \eqn{[0, 1]} when both functions are densities and the interval holds
#' almost all of their mass.
#'
#' @param fhat,ftrue Density functions that accept a vector.
#' @param lo,hi Integration limits.
#' @param n Number of grid points.
#'
#' @return A single number.
#'
#' @export
tv_density <- function(fhat, ftrue, lo, hi, n = 2001) {
  xs <- seq(lo, hi, length.out = n)
  d  <- abs(fhat(xs) - ftrue(xs))
  return(0.5 * sum(diff(xs) * (utils::head(d, -1) + utils::tail(d, -1)) / 2))
}
