# Generated from create-LCDmix.Rmd: do not edit by hand

#' Piecewise-linear log-density with linear extension beyond the support
#'
#' @description
#' Builds an evaluator for the fitted log-concave log-density
#' \eqn{\hat g_k}: linear between knots, and continued beyond
#' \eqn{[L_k, U_k]} with the slopes of the two boundary pieces. The function
#' is then finite, concave and piecewise linear everywhere. This is
#' deliberate: the quasi-Newton arm of the LP ablation needs a finite
#' objective outside the support. With the \code{-Inf} convention,
#' \code{optim} stops at its first step outside the support and the
#' comparison shows nothing.
#'
#' @param density_k A \code{modified_logcondens()} fit (fields \code{x},
#'   \code{phi} and \code{IsKnot} are used). It needs at least two knots.
#'
#' @return A list with \code{value(u)} (the log-density at \code{u}),
#'   \code{slope(u)} (its slope at \code{u}; the right-hand slope at a knot)
#'   and the support endpoints \code{L} and \code{U}.
#' @keywords internal
make_logdens_ext <- function(density_k) {
  knot  <- as.logical(density_k$IsKnot)
  x_m   <- density_k$x[knot]
  phi_m <- density_k$phi[knot]

  M <- length(x_m)
  if (M < 2L) stop("make_logdens_ext(): fitted density has fewer than two knots.")

  L <- x_m[1L]
  U <- x_m[M]

  b_left  <- (phi_m[2L] - phi_m[1L]) / (x_m[2L] - x_m[1L])
  b_right <- (phi_m[M]  - phi_m[M - 1L]) / (x_m[M] - x_m[M - 1L])

  value <- function(u) {
    out <- numeric(length(u))
    lo  <- u < L
    hi  <- u > U
    mid <- !lo & !hi
    if (any(mid)) out[mid] <- stats::approx(x_m, phi_m, xout = u[mid])$y
    if (any(lo))  out[lo]  <- phi_m[1L] + b_left  * (u[lo] - L)
    if (any(hi))  out[hi]  <- phi_m[M]  + b_right * (u[hi] - U)
    out
  }

  # Local slope; at a kink the right-hand slope is returned (a valid
  # subgradient element, and the point of the experiment is that these exist).
  slope <- function(u) {
    out <- numeric(length(u))
    lo  <- u < L
    hi  <- u > U
    mid <- !lo & !hi
    if (any(lo)) out[lo] <- b_left
    if (any(hi)) out[hi] <- b_right
    if (any(mid)) {
      j <- findInterval(u[mid], x_m, rightmost.closed = TRUE, all.inside = TRUE)
      out[mid] <- (phi_m[j + 1L] - phi_m[j]) / (x_m[j + 1L] - x_m[j])
    }
    out
  }

  list(value = value, slope = slope, L = L, U = U)
}
