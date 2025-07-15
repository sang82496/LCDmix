# Generated from create-LCDmix.Rmd: do not edit by hand

#' Modified log‐concave density estimation with optional preprocessing
#'
#' @description
#' Fits a univariate log‐concave density to data \code{x} using the
#' \pkg{logcondens} routines.  If \code{w} is \code{NA}, a preprocessing step
#' chooses an appropriate weight vector and grid; otherwise \code{x} and \code{w}
#' are assumed sorted and used directly.  The estimator is refined by adding
#' knots until convergence.
#'
#' @param x Numeric vector of data points at which to estimate a log‐concave density.
#' @param xgrid Optional numeric vector of grid points for preprocessing.  If
#'   provided, \code{w} must be \code{NA}.  Default: \code{NULL}.
#' @param print Logical; if \code{TRUE}, iteration progress and likelihood values
#'   are printed.  Default: \code{FALSE}.
#' @param w Numeric vector of nonnegative weights for each \code{x}, or \code{NA}
#'   to compute weights automatically.  Default: \code{NA}.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{xn}}{Sorted original \code{x} values.}
#'   \item{\code{x}}{Processed \code{x} after any preprocessing.}
#'   \item{\code{w}}{Weights corresponding to \code{x}.}
#'   \item{\code{phi}}{Numeric vector of estimated log‐density values at \code{x}.}
#'   \item{\code{IsKnot}}{Integer or logical vector indicating which points are knots.}
#'   \item{\code{L}}{Final log‐likelihood value.}
#'   \item{\code{Fhat}}{Estimated CDF values at \code{x}.}
#'   \item{\code{H}}{Numeric vector of directional derivatives (used in knot selection).}
#'   \item{\code{n}}{Number of original data points (\code{length(xn)}).}
#'   \item{\code{m}}{Same as \code{n}.}
#'   \item{\code{knots}}{Vector of \code{x} values selected as knots.}
#'   \item{\code{mode}}{Value of \code{x} at which \code{phi} is maximized.}
#'   \item{\code{sig}}{Estimated standard deviation used for preprocessing.}
#' }
#'
#' @examples
#' \dontrun{
#' # Example: estimate log‐concave density of a mixture sample
#' set.seed(42)
#' x1 <- rnorm(100, mean = -2)
#' x2 <- rnorm(150, mean =  3)
#' x  <- c(x1, x2)
#'
#' # Default call (weights computed automatically)
#' res1 <- modified_logcondens(x)
#' plot(res1$xn, res1$phi, type = "l", xlab = "x", ylab = "log‐density")
#'
#' # Provide custom weights (e.g., uniform)
#' w <- rep(1, length(x))
#' res2 <- modified_logcondens(x, w = w)
#' lines(res2$xn, res2$phi, col = "blue")
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
    while ((iter1 < 500) & (max(H) > prec * mean(abs(H)))) {
        IsKnot_old <- IsKnot
        iter1 <- iter1 + 1
        tmp <- max(H)
        k <- (1:n) * (H == tmp)
        k <- min(k[k > 0])
        IsKnot[k] <- 1
        res2 <- logcondens::LocalMLE(x, w, IsKnot, phi, prec)
        phi_new <- res2$phi
        L <- res2$L
        conv_new <- res2$conv
        H <- res2$H
        while ((max(conv_new) > prec * max(abs(conv_new)))) {
            JJ <- (1:n) * (conv_new > 0)
            JJ <- JJ[JJ > 0]
            if (length(JJ) == 1 && conv[JJ] == conv_new[JJ]){ # inserted
              print('break')
              break
              } 
            tmp <- conv[JJ]/(conv[JJ] - conv_new[JJ])
            lambda <- min(tmp)
            KK <- (1:length(JJ)) * (tmp == lambda)
            KK <- KK[KK > 0]
            IsKnot[JJ[KK]] <- 0
            phi <- (1 - lambda) * phi + lambda * phi_new
            conv <- pmin(c(logcondens::LocalConvexity(x, phi), 0))
            res3 <- logcondens::LocalMLE(x, w, IsKnot, phi, prec)
            phi_new <- res3$phi
            L <- res3$L
            conv_new <- res3$conv
            H <- res3$H
        }
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
            max(phi)], sig = sig)
    return(res)
}
