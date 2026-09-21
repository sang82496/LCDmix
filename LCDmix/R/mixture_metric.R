# Generated from create-LCDmix.Rmd: do not edit by hand

#' Biomass-weighted L1 distance between the estimated and true mixture densities
#'
#' @description
#' For each time point t, computes the L1 distance between the estimated
#' conditional density and the true one,
#' \deqn{D_t = \int |\hat p_t(y) - p_t(y)| \, dy,}
#' and returns the average of \eqn{D_t} weighted by the total biomass at each
#' time point. The value lies in \eqn{[0, 2]}.
#'
#' The integral is computed on a fine grid, not at the bin midpoints. Each
#' log-concave component density is evaluated exactly: its log-density is
#' linear between its sorted support points and it is zero outside the
#' support. The grid holds every knot and both ends of every support, so the
#' kinks and the jumps to zero are resolved, and the trapezoid rule integrates
#' the absolute difference. flowmix components are Gaussian and need no knots.
#'
#' @param sim     A simulated dataset from \code{gen_simul_data()}.
#' @param est_res A fitted model: the \code{iter} element of an LCDmix fit, or a
#'   flowmix best result (\code{alpha}, \code{mn}, \code{sigma}).
#' @param rescale Kept for backward compatibility; it has no effect. The value
#'   is always the integral.
#' @param n_grid  Number of equally spaced grid points per time point, before
#'   the knots and support ends are added. Default 4000.
#'
#' @return A single number, the biomass-weighted mean of \eqn{D_t}, with
#'   attributes \code{per_time} (the length-TT vector of \eqn{D_t}),
#'   \code{mass_est} (the range over time points of the integral of the
#'   estimated density, which should be close to 1), and \code{h} (the bin
#'   width of the data, reported for scripts that read it).
#' @export
mixture_metric <- function(
  sim,
  est_res,
  rescale = TRUE,                     # kept for backward compatibility; no effect
  n_grid  = 4000                      # NEW - appended
) {
  if (!isTRUE(rescale))
    warning("mixture_metric(): 'rescale' has no effect; the value is always the integral")
  Y_bin    <- sim$Y_bin
  X        <- sim$X
  bin_mass <- sim$bin_mass
  TT       <- length(Y_bin)
  pi_true  <- pi_k(X, t(sim$alpha))
  is_lcd   <- !is.null(est_res$alpha_new)

  if (is_lcd) {                       # LCDmix
    pi_est <- pi_k(X, est_res$alpha_new)
    K      <- length(est_res$g_new)
    ## A log-concave MLE has a log-density that is linear between its sorted
    ## support points, so linear interpolation of phi is exact. It matches
    ## logcondens::evaluateLogConDens() to 4.4e-16 and is far faster.
    dens_k <- lapply(est_res$g_new, function(g) {
      x <- g$x; ph <- g$phi
      return(function(v) {
        out <- numeric(length(v))
        inside <- v >= x[1] & v <= x[length(x)]
        if (any(inside))
          out[inside] <- exp(stats::approx(x, ph, v[inside], ties = "ordered")$y)
        return(out)
      })
    })
    knots_k <- lapply(est_res$g_new, function(g) {
      x <- if (!is.null(g$IsKnot)) g$x[g$IsKnot == 1] else g$x
      return(sort(unique(c(x, range(g$x)))))
    })
    ends_k <- lapply(est_res$g_new, function(g) range(g$x))
  } else {                            # flowmix
    pi_est <- pi_k(X, est_res$alpha)
    K      <- dim(est_res$mn)[3]
    sd_k   <- sqrt(as.numeric(est_res$sigma))
  }

  per_time <- numeric(TT)
  mass_est <- numeric(TT)
  for (t in seq_len(TT)) {
    if (is_lcd) {
      mu <- vapply(seq_len(K), function(k)
        est_res$theta0_new[[k]] + sum(X[t, ] * est_res$theta_new[[k]]), numeric(1))
      ends  <- unlist(lapply(seq_len(K), function(k) mu[k] + ends_k[[k]]))
      extra <- c(unlist(lapply(seq_len(K), function(k) mu[k] + knots_k[[k]])),
                 ends - 1e-9, ends + 1e-9)       # resolve the kinks and the jumps to zero
    } else {
      mu    <- est_res$mn[t, 1, ]
      ends  <- c(mu - 8 * sd_k, mu + 8 * sd_k)
      extra <- numeric(0)
    }
    lo <- min(min(Y_bin[[t]]) - 4, ends)
    hi <- max(max(Y_bin[[t]]) + 4, ends)
    z  <- sort(unique(c(seq(lo, hi, length.out = n_grid), extra)))
    z  <- z[z >= lo & z <= hi]

    if (is_lcd) {
      p_est <- rowSums(vapply(seq_len(K), function(k)
        pi_est[t, k] * dens_k[[k]](z - mu[k]), numeric(length(z))))
    } else {
      p_est <- rowSums(vapply(seq_len(K), function(k)
        pi_est[t, k] * stats::dnorm(z, mu[k], sd_k[k]), numeric(length(z))))
    }
    p_true <- as.numeric(dens_true_fun(sim, t, z) %*% pi_true[t, ])

    dz <- diff(z)
    per_time[t] <- sum(dz * (utils::head(abs(p_est - p_true), -1) +
                             utils::tail(abs(p_est - p_true), -1)) / 2)
    mass_est[t] <- sum(dz * (utils::head(p_est, -1) + utils::tail(p_est, -1)) / 2)
  }

  w_t    <- vapply(bin_mass, sum, numeric(1))
  metric <- sum(w_t * per_time) / sum(w_t)
  grid   <- sort(unique(as.numeric(unlist(Y_bin))))
  attr(metric, "per_time") <- per_time
  attr(metric, "mass_est") <- range(mass_est)
  attr(metric, "h")        <- if (length(grid) > 1) min(diff(grid)) else NA_real_
  return(metric)
}
