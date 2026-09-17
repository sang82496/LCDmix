# Generated from create-LCDmix.Rmd: do not edit by hand

#' Half the L1 integral between two densities, so the value lies in [0, 1].
#' @export
tv_density <- function(fhat, ftrue, lo, hi, n = 2001) {
  xs <- seq(lo, hi, length.out = n)
  d  <- abs(fhat(xs) - ftrue(xs))
  return(0.5 * sum(diff(xs) * (utils::head(d, -1) + utils::tail(d, -1)) / 2))
}
