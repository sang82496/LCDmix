# Generated from create-LCDmix.Rmd: do not edit by hand

#' Overlap of a density and its shifted copy
#'
#' @description
#' Computes \eqn{\int \min\{f(v), f(v - \delta)\}\, dv}, the overlap coefficient
#' of \eqn{f} and \eqn{f} shifted by \eqn{\delta}. In the simulation both
#' components share one error density and their means differ by a shift, so
#' this measures how much the components overlap: 1 when they are identical,
#' 0 when they do not overlap. If the integration fails, it is repeated with
#' the infinite limits replaced by -40 and \code{40 + delta}.
#'
#' @param f A density function that accepts a vector.
#' @param delta The shift.
#' @param lo,hi Integration limits.
#'
#' @return The overlap, a number in \eqn{[0, 1]}.
#'
#' @export
ovl <- function(f, delta, lo = -Inf, hi = Inf) {
  integrand <- function(v) pmin(f(v), f(v - delta))
  out <- tryCatch(
    stats::integrate(integrand, lo, hi, subdivisions = 2000L)$value,
    error = function(e) NA_real_)
  if (is.na(out)) {
    ## Fall back to wide finite limits when the infinite range fails.
    lo2 <- if (is.finite(lo)) lo else -40
    hi2 <- if (is.finite(hi)) hi else  40 + delta
    out <- stats::integrate(integrand, lo2, hi2, subdivisions = 2000L)$value
  }
  return(out)
}
