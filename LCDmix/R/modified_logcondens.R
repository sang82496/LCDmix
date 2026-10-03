# Generated from create-LCDmix.Rmd: do not edit by hand

#' Weighted log-concave maximum likelihood estimate with safeguards
#'
#' @description
#' A modified copy of \code{logcondens::activeSetLogCon()}, the active-set
#' algorithm for the log-concave maximum likelihood estimate. Knots are added
#' one at a time until the estimate converges. Two changes stop the knot
#' search early and keep the last good fit instead of failing:
#'
#' - When the line search stalls. Continuing has been seen to drive the
#'   log-density to numerical divergence a few iterations later.
#' - When \code{logcondens::LocalMLE()} fails or returns non-finite values
#'   (fixW). This happens when a support point has a very small weight: its
#'   log-density falls below about -745, \code{exp()} returns 0, and the Newton
#'   step is \code{NaN}.
#'
#' The field \code{refine_stop} records which of the two, if any, applied.
#'
#' @param x Numeric vector of data points. It need not be sorted.
#' @param xgrid Optional grid for \code{logcondens::preProcess()}; allowed only
#'   when \code{w = NA}.
#' @param print Logical; if \code{TRUE}, prints the log-likelihood and the
#'   number of knots at each step.
#' @param w \code{NA}, to compute weights with \code{logcondens::preProcess()},
#'   or a vector of nonnegative weights, one per element of \code{x}, that
#'   sums to 1 (\code{mstep_g()} passes normalized weights).
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{xn}}{The input \code{x}, sorted.}
#'   \item{\code{x}}{The sorted support points of the fit (equal to \code{xn}
#'     when \code{w} is given).}
#'   \item{\code{w}}{The weights, in the order of \code{x}.}
#'   \item{\code{phi}}{The log-density at \code{x}. It is linear between
#'     consecutive points and \code{-Inf} outside \code{range(x)}.}
#'   \item{\code{IsKnot}}{0/1 vector; 1 where \code{x} is a knot of \code{phi}.}
#'   \item{\code{L}}{Final value of the log-likelihood.}
#'   \item{\code{Fhat}}{Distribution function at \code{x}.}
#'   \item{\code{H}}{Directional derivatives used to choose the next knot.}
#'   \item{\code{n}}{\code{length(xn)}.}
#'   \item{\code{m}}{\code{length(x)}, the number of support points.}
#'   \item{\code{knots}}{\code{x[IsKnot == 1]}.}
#'   \item{\code{mode}}{The value(s) of \code{x} where \code{phi} is largest.}
#'   \item{\code{sig}}{Standard deviation estimate of \code{x} (weighted when
#'     \code{w} is given).}
#'   \item{\code{refine_stop}}{\code{"none"}, \code{"line search stalled"} or
#'     \code{"LocalMLE failed"}.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(42)
#' x <- c(rnorm(100, mean = -2), rnorm(150, mean = 3))
#'
#' # Weights computed by logcondens::preProcess()
#' res1 <- modified_logcondens(x)
#' plot(res1$x, res1$phi, type = "l", xlab = "x", ylab = "log-density")
#'
#' # Equal weights that sum to 1
#' res2 <- modified_logcondens(x, w = rep(1 / length(x), length(x)))
#' lines(res2$x, res2$phi, col = "blue")
#' res2$refine_stop
#' }
#' @export
modified_logcondens <- function(x, xgrid = NULL, print = FALSE, w = NA){
    prec <- 1e-10
    xn <- sort(x)
    if ((!identical(xgrid, NULL) & (!identical(w, NA)))) {
        stop("If w != NA then xgrid must be NULL!\n")
    }
    if (identical(w, NA)) {
        tmp <- logcondens::preProcess(x, xgrid = xgrid)
        x <- tmp$x
        w <- tmp$w
        sig <- tmp$sig
    }
    if (!identical(w, NA)) {
        tmp <- cbind(x, w)
        tmp <- tmp[order(x), ]
        x <- tmp[, 1]
        w <- tmp[, 2]
        est.m <- sum(w * x)
        est.sd <- sum(w * (x - est.m)^2)
        est.sd <- sqrt(est.sd * length(x)/(length(x) - 1))
        sig <- est.sd
    }
    n <- length(x)
    phi <- logcondens::LocalNormalize(x, 1:n * 0)
    IsKnot <- 1:n * 0
    IsKnot[c(1, n)] <- 1
    res1 <- logcondens::LocalMLE(x = x, w = w, IsKnot = IsKnot, phi_o = phi, 
        prec = prec)
    phi <- res1$phi
    L <- res1$L
    conv <- res1$conv
    H <- res1$H
    iter1 <- 1
    ## fixW: logcondens::LocalMLE() can return NaN, or stop with "missing value
    ## where TRUE/FALSE needed", when a support point has a very small weight:
    ## its log-density falls below about -745, exp() gives 0, and the Newton
    ## step is NaN. lmle_safe() returns NULL in that case, and the search below
    ## then keeps the last good fit instead of failing.
    lmle_safe <- function(IsKnot, phi) {
      r <- tryCatch(logcondens::LocalMLE(x, w, IsKnot, phi, prec), error = function(e) NULL)
      if (is.null(r) || !all(is.finite(r$phi)) || !all(is.finite(r$H)) || !all(is.finite(r$conv))) {
        return(NULL)
      }
      return(r)
    }
    refine_stop <- "none"                                  # fixW: why the search stopped early, if it did
    while ((iter1 < 500) & (max(H) > prec * mean(abs(H)))) {
        IsKnot_old <- IsKnot
        phi_old <- phi; L_old <- L; H_old <- H; conv_old <- conv   # fixW: the state to return to
        iter1 <- iter1 + 1
        tmp <- max(H)
        k <- (1:n) * (H == tmp)
        k <- min(k[k > 0])
        IsKnot[k] <- 1
        res2 <- lmle_safe(IsKnot, phi)                     # fixW: was logcondens::LocalMLE(x, w, IsKnot, phi, prec)
        if (is.null(res2)) {                               # fixW: discard the attempted knot, keep the last good fit
          IsKnot <- IsKnot_old; phi <- phi_old; L <- L_old; H <- H_old; conv <- conv_old
          refine_stop <- "LocalMLE failed"
          break
        }
        phi_new <- res2$phi
        L <- res2$L
        conv_new <- res2$conv
        H <- res2$H
        stalled <- FALSE                                  # NEW
        while ((max(conv_new) > prec * max(abs(conv_new)))) {
            JJ <- (1:n) * (conv_new > 0)
            JJ <- JJ[JJ > 0]
            if (length(JJ) == 1 && conv[JJ] == conv_new[JJ]){
              # Line search stalled: adopting phi_new and continuing has been
              # shown to drive phi into numerical divergence a few outer
              # iterations later (e.g. from a sane range like [-10,-1] out to
              # [-230121,-230112]), at which point logcondens::LocalMLE()
              # overflows internally and crashes with "NAs are not allowed
              # in subscripted assignments". Stop refining here instead of
              # letting that happen.
              stalled <- TRUE                              # NEW
              refine_stop <- "line search stalled"         # fixW
              break
              } 
            tmp <- conv[JJ]/(conv[JJ] - conv_new[JJ])
            lambda <- min(tmp)
            KK <- (1:length(JJ)) * (tmp == lambda)
            KK <- KK[KK > 0]
            IsKnot[JJ[KK]] <- 0
            phi <- (1 - lambda) * phi + lambda * phi_new
            conv <- pmin(c(logcondens::LocalConvexity(x, phi), 0))
            res3 <- lmle_safe(IsKnot, phi)                 # fixW: was logcondens::LocalMLE(x, w, IsKnot, phi, prec)
            if (is.null(res3)) {                           # fixW
              stalled <- TRUE
              refine_stop <- "LocalMLE failed"
              break
            }
            phi_new <- res3$phi
            L <- res3$L
            conv_new <- res3$conv
            H <- res3$H
        }
        if (stalled) {                                     # NEW
          IsKnot <- IsKnot_old                              # NEW: discard this
          phi <- phi_old; L <- L_old; H <- H_old; conv <- conv_old   # fixW: restore the whole state, not only the knots
          break                                             # NEW: outer iteration's
        }                                                   # NEW: attempted knot
        phi <- phi_new
        conv <- conv_new
        if (sum(IsKnot != IsKnot_old) == 0) {
            break
        }
        if (print == TRUE) {
            print(paste("iter1 = ", iter1 - 1, " / L = ", round(L, 
                4), " / max(H) = ", round(max(H), 4), " / #knots = ", 
                sum(IsKnot), sep = ""))
        }
    }
    Fhat <- logcondens::LocalF(x, phi)
    res <- list(xn = xn, x = x, w = w, phi = as.vector(phi), 
        IsKnot = IsKnot, L = L, Fhat = as.vector(Fhat), H = as.vector(H), 
        n = length(xn), m = n, knots = x[IsKnot == 1], mode = x[phi == 
            max(phi)], sig = sig,
        refine_stop = refine_stop)                         # fixW: "none", "line search stalled" or "LocalMLE failed"
    return(res)
}
