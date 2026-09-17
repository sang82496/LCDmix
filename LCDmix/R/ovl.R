# Generated from create-LCDmix.Rmd: do not edit by hand

#' Overlap coefficient between a density and its shift by delta.
#' Both components share one error density, so this is the whole story.
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
