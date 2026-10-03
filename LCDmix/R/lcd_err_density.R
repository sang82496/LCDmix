# Generated from create-LCDmix.Rmd: do not edit by hand

#' Estimated error density of one LCDmix component
#'
#' @description
#' Returns a function that evaluates a fitted log-concave density at residual
#' values: \code{logcondens::evaluateLogConDens()} (column 3, the density)
#' inside the support, and 0 outside it or wherever the evaluation fails.
#'
#' @param g One element of \code{g_new}, a \code{modified_logcondens()} fit.
#'
#' @return A function of a numeric vector \code{v} that returns the densities
#'   at \code{v}.
#' @keywords internal
lcd_err_density <- function(g) {
  sup <- range(g$x)
  return(function(v) {
    out <- numeric(length(v))
    inside <- v >= sup[1] & v <= sup[2]
    if (any(inside)) {
      val <- tryCatch(logcondens::evaluateLogConDens(v[inside], g, which = 2)[, 3],   # fixT: column 3 = density
                      error = function(e) rep(NA_real_, sum(inside)))
      out[inside] <- val
    }
    out[!is.finite(out)] <- 0
    return(out)
  })
}
