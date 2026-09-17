# Generated from create-LCDmix.Rmd: do not edit by hand

#' Estimated error density of one LCDmix component, zero outside the support.
#'
#' NOTE: evaluateLogConDens() returns a matrix, not a list, so the work order's
#' `$log.density` does not work. Column 3 is the density, which is what the
#' package's own dens_est_fun() uses.
lcd_err_density <- function(g) {
  sup <- range(g$x)
  return(function(v) {
    out <- numeric(length(v))
    inside <- v >= sup[1] & v <= sup[2]
    if (any(inside)) {
      val <- tryCatch(logcondens::evaluateLogConDens(v[inside], g)[, 3],
                      error = function(e) rep(NA_real_, sum(inside)))
      out[inside] <- val
    }
    out[!is.finite(out)] <- 0
    return(out)
  })
}
