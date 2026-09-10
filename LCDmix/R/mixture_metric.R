# Generated from create-LCDmix.Rmd: do not edit by hand

#' Compute biomass‐weighted L1 (total‐variation) distance between estimated and true mixtures
#'
#' @description
#' For each time point \code{t}, evaluates the estimated mixture density 
#' (either from an LCDmix fit or a flowmix fit) and the true mixture density 
#' at the observed bin midpoints, computes the \(L^1\) distance, and then 
#' returns both the per‐time distances and the overall weighted average 
#' distance (weighted by total biomass at each \code{t}).
#'
#' @param Y_bin List of length \code{TT}; each element is a numeric vector of bin midpoints for time \code{t}.
#' @param X Numeric \code{TT × p} matrix of covariates (rows = time points).
#' @param bin_mass List of length \code{TT}; each element is a numeric vector of biomass weights aligned with \code{Y_bin[[t]]}.
#' @param est_res A fitted model object:
#'   \describe{
#'     \item{LCDmix fits}{A list containing \code{alpha_new}, \code{theta0_new}, \code{theta_new}, \code{g_new}.}
#'     \item{flowmix fits}{A list containing \code{alpha}, \code{mn}, and \code{sigma}.}
#'   }
#' @param true_res A list with components:
#'   \describe{
#'     \item{\code{prob}}{\code{TT × K} numeric matrix of true mixing proportions \eqn{\pi_{t,k}}.}
#'     \item{\code{dens_true}}{List of length \code{TT}, each an \code{n_t × K} numeric matrix of the true component densities evaluated at \code{Y_bin[[t]]}.}
#'   }
#' @param rescale Logical; if \code{TRUE} (default) multiply the midpoint
#' Riemann sum by the bin width \eqn{h}, returning the integrated absolute
#'   error \eqn{\int |\hat f - f|\,dy \in [0, 2]}. If \code{FALSE}, return the
#'   unscaled sum \eqn{\sum_b |\hat f_b - f_b|}, which is \eqn{1/h} times that
#'   and is what versions of this function before the rescale returned. Ignored
#'   with a warning when the bin midpoints do not lie on an equally spaced
#'   lattice (e.g. data binned with \code{n_bins = 0}).
#'
#' @return A list with elements:
#' \describe{
#'   \item{\code{per_time}}{Numeric vector of length \code{TT}, the \(L^1\) distance at each time point.}
#'   \item{\code{weighted}}{Scalar, the biomass‐weighted average of those distances.}
#' }
#' @return Scalar: the biomass-weighted average over time of the \eqn{L^1}
#'   distance between the estimated and true mixture densities. Carries
#'   attributes \code{"h"} (the bin width used, \code{NA} if not rescaled) and
#'   \code{"per_time"} (the length-\code{TT} vector of per-time distances).
#'   
#' @examples
#' \dontrun{
#' sim <- generate_skewed_data(seed = 42)
#' Y_bin    <- sim$ylist
#' X        <- sim$X
#' bin_mass <- sim$countslist
#'
#' # Fit an LCDmix model
#' fit <- main(Y = Y_bin, X = X, biomass = bin_mass, K = sim$numclust)
#' true_res <- list(prob = sim$prob, dens_true = sim$dens_true)
#'
#' # Compute metric
#' dist <- mixture_metric(Y_bin, X, bin_mass, fit$iter, true_res)
#' plot(dist$per_time, type = "b", xlab = "Time", ylab = "L1 distance")
#' print(dist$weighted)
#' }
#' @export
mixture_metric <- function(
  sim,
  est_res,
  rescale = TRUE                      # NEW - appended
) {
  Y_bin     = sim$Y_bin
  X         = sim$X
  bin_mass  = sim$bin_mass
  TT        = length(Y_bin)
  pi_true = pi_k(X, t(sim$alpha))
  per_time <- numeric(TT)
  
  if (!is.null(est_res$alpha_new)) { # if LCDmix 
    pi_est <- pi_k(X, est_res$alpha_new)
  } else { # if flowmix
    pi_est <- pi_k(X, est_res$alpha)
  }
  
  for (t in seq_len(TT)) {
    dens_est  <- dens_est_fun(est_res, t, Y_bin[[t]], X)
    mix_est   <- dens_est %*% pi_est[t, ]
    dens_true <- dens_true_fun(sim, t, Y_bin[[t]])
    mix_true  <- dens_true %*% pi_true[t, ]
    per_time[t] <- sum(abs(mix_est - mix_true))
  }
  
  ## --- NEW: recover the bin width and convert the sum into an integral -------
  ## All time points share one global equal-width grid (see binning()), so the
  ## smallest gap between distinct midpoints IS the bin width. Empty bins make
  ## some gaps 2h, 3h, ..., hence the integer-multiple test rather than an
  ## all-gaps-equal test.
  h <- NA_real_
  if (isTRUE(rescale)) {
    grid <- sort(unique(as.numeric(unlist(Y_bin))))
    if (length(grid) < 2L) {
      warning("mixture_metric(): fewer than two distinct bin midpoints; ",
              "returning the unscaled sum")
    } else {
      d  <- diff(grid)
      h0 <- min(d)
      if (max(abs(d / h0 - round(d / h0))) > 1e-6) {
        warning("mixture_metric(): bin midpoints are not on an equally spaced ",
                "lattice (n_bins = 0?); returning the unscaled sum")
      } else {
        h        <- h0
        per_time <- h * per_time
      }
    }
  }
  ## --------------------------------------------------------------------------

  w_t      <- vapply(bin_mass, sum, numeric(1))
  metric   <- sum(w_t * per_time) / sum(w_t)

  attr(metric, "h")        <- h
  attr(metric, "per_time") <- per_time
  return(metric)
}
