# Generated from create-LCDmix.Rmd: do not edit by hand

#' Update the expert coefficients of all components
#'
#' @description
#' For each component \eqn{k}, updates the intercept and the slopes with
#' \code{mstep_theta_lp()} (\code{update = "lp"}, the method of the paper) or
#' \code{mstep_theta_optim()} (\code{update = "optim"}, the comparison arm of
#' the LP ablation), and records diagnostics that are computed the same way
#' for both. The returned intercept is the update's own intercept;
#' \code{iteration()} then replaces it with \code{mstep_shift()}.
#'
#' @inheritParams mstep_theta_lp
#' @param densities A list of length \code{K} of \code{modified_logcondens()}
#'   fits, one per component.
#' @param intercepts A list of length \code{K} of current intercepts.
#' @param slopes A list of length \code{K} of current slope vectors.
#' @param update \code{"lp"} or \code{"optim"}; see \code{main()}.
#' @param maxdev \code{NULL} or a positive number; see
#'   \code{mstep_theta_lp()}. Only \code{update = "lp"} supports it; with
#'   \code{"optim"} a non-\code{NULL} value is an error.
#'
#' @return A list with components:
#' \describe{
#'   \item{theta0}{List of \code{K} updated intercepts.}
#'   \item{theta}{List of \code{K} updated slope vectors.}
#'   \item{diag}{List of \code{K} diagnostic lists, each with \code{arm}
#'     (the value of \code{update}), \code{seconds} (time of the update),
#'     \code{n_active} (number of active bins), \code{n_outside} (active bins
#'     whose new residual lies outside the support of the density),
#'     \code{convergence}, \code{obj_start}, \code{obj_end} (from
#'     \code{mstep_theta_optim()}; \code{NA} for the LP), and
#'     \code{maxdev_n_relaxed}, \code{maxdev_max_excess} (from
#'     \code{mstep_theta_lp()} with \code{maxdev}; \code{NA} otherwise).}
#' }
#'
#' @examples
#' \dontrun{
#' TT <- 6; K <- 2; p <- 2
#' X          <- matrix(rnorm(TT * p), nrow = TT, ncol = p)
#' Y_bin      <- lapply(1:TT, function(t) matrix(sort(rnorm(20)), ncol = 1))
#' intercepts <- list(-0.5, 0.5)
#' slopes     <- replicate(K, rep(0, p), simplify = FALSE)
#' resi       <- comp_resi(Y_bin, X, intercepts, slopes)
#' weights    <- lapply(Y_bin, function(m) matrix(runif(nrow(m) * K), ncol = K))
#' idx        <- lapply(weights, function(w) w > 0)
#' # The densities must be fitted to the current residuals, or the
#' # linear program has no feasible point
#' densities  <- mstep_g(resi, weights, idx)
#' result <- mstep_theta(
#'   Y_bin, X, weights, resi,
#'   densities, idx,
#'   intercepts, slopes,
#'   lambda_theta = 1e-3
#' )
#' str(result$theta)
#' }
#' @export
mstep_theta <- function(
  Y_bin,
  X,
  weights,
  residuals,
  densities,
  idx,
  intercepts,
  slopes,
  lambda_theta,
  lp_time_limit = 3600,
  update        = c("lp", "optim"),  # NEW; default preserves current behavior
  maxdev        = NULL               # NEW (fixP) - appended last
) {
  update <- match.arg(update)
  if (!is.null(maxdev) && update != "lp")                       # NEW (fixP)
    stop("mstep_theta(): maxdev is implemented only for update = \"lp\".")

  K <- length(densities)
  theta0_new <- vector("list", K)
  theta_new  <- vector("list", K)
  diag_list  <- vector("list", K)

  for (k in seq_len(K)) {

    t0 <- proc.time()[["elapsed"]]

    tmp <- switch(
      update,
      lp = mstep_theta_lp(
        Y_bin = Y_bin, X = X, weights = weights, residuals = residuals,
        density_k = densities[[k]], idx = idx,
        intercept_k = intercepts[[k]], slopes_k = slopes[[k]],
        lambda_theta = lambda_theta, component = k,
        lp_time_limit = lp_time_limit, maxdev = maxdev   # fixP
      ),
      optim = mstep_theta_optim(
        Y_bin = Y_bin, X = X, weights = weights, residuals = residuals,
        density_k = densities[[k]], idx = idx,
        intercept_k = intercepts[[k]], slopes_k = slopes[[k]],
        lambda_theta = lambda_theta, component = k,
        lp_time_limit = lp_time_limit
      )
    )

    elapsed <- proc.time()[["elapsed"]] - t0

    theta0_new[[k]] <- tmp$theta0_k
    theta_new[[k]]  <- tmp$theta_k

    ## ---- diagnostics, computed identically for BOTH arms --------------------
    ## n_outside must be measured the same way in each arm or the headline
    ## metric is not comparable. mstep_theta_optim returns its own, but we
    ## recompute here so the LP arm gets the same number by the same route.
    g_ext <- make_logdens_ext(densities[[k]])
    u_new <- numeric(0)
    for (t in seq_len(length(Y_bin))) {
      idx_tk <- idx[[t]][, k]
      if (any(idx_tk)) {
        u_new <- c(u_new,
                   Y_bin[[t]][idx_tk, 1] - tmp$theta0_k -
                     sum(X[t, ] * tmp$theta_k))
      }
    }

    diag_list[[k]] <- list(
      arm         = update,
      seconds     = elapsed,
      n_active    = length(u_new),
      n_outside   = sum(u_new < g_ext$L | u_new > g_ext$U),
      convergence = if (is.null(tmp$convergence)) NA_integer_ else tmp$convergence,
      obj_start   = if (is.null(tmp$obj_start))   NA_real_    else tmp$obj_start,
      obj_end     = if (is.null(tmp$obj_end))     NA_real_    else tmp$obj_end,
      # NEW (fixP): how often the maxdev guard was needed in this LP call
      maxdev_n_relaxed  = if (is.null(tmp$maxdev_n_relaxed))  NA_integer_ else tmp$maxdev_n_relaxed,
      maxdev_max_excess = if (is.null(tmp$maxdev_max_excess)) NA_real_    else tmp$maxdev_max_excess
    )
  }

  list(
    theta0 = theta0_new,
    theta  = theta_new,
    diag   = diag_list          # NEW; ignored by existing callers
  )
}
