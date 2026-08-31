# Generated from create-LCDmix.Rmd: do not edit by hand

#' M-step update of intercept and slopes by quasi-Newton (no LP, no constraints)
#'
#' @description
#' Comparison arm for the LP ablation. Maximizes the same weighted
#' log-likelihood with the same L1 penalty as \code{mstep_theta_lp()}, but with
#' \code{optim(method = "L-BFGS-B")} and without the log-concavity epigraph
#' constraints or the support-feasibility bounds.
#'
#' The slope vector is split as \eqn{\theta = \theta^+ - \theta^-} with both
#' parts nonnegative, mirroring the LP's variable structure so that the two arms
#' differ in solver rather than in problem statement, and so the non-smooth L1
#' penalty is handled the same way in both.
#'
#' @inheritParams mstep_theta_lp
#' @param maxit Integer; \code{optim} iteration cap. Default 500.
#' @param use_gradient Logical; supply the analytic (sub)gradient. Default TRUE.
#'
#' @return A list with \code{theta0_k}, \code{theta_k}, and diagnostics
#'   \code{convergence}, \code{obj_start}, \code{obj_end}, \code{n_outside},
#'   \code{counts}. \code{mstep_theta()} reads only the first two, so the
#'   diagnostics are free to carry.
#' @export
mstep_theta_optim <- function(
  Y_bin,
  X,
  weights,
  residuals,
  density_k,
  idx,
  intercept_k,
  slopes_k,
  lambda_theta,
  component,
  lp_time_limit = NULL,   # accepted and ignored; keeps the signature identical
  maxit         = 500L,
  use_gradient  = TRUE
) {
  TT <- length(Y_bin)
  p  <- ncol(X)

  ## ---- collect the same bins the LP would use -------------------------------
  w_k <- numeric(0)
  Y_k <- numeric(0)
  X_k <- matrix(nrow = 0, ncol = p)

  for (t in seq_len(TT)) {
    idx_tk <- idx[[t]][, component]
    if (any(idx_tk)) {
      w_k <- c(w_k, weights[[t]][idx_tk, component])
      Y_k <- c(Y_k, Y_bin[[t]][idx_tk, 1])
      X_k <- rbind(X_k,
                   matrix(rep(X[t, ], sum(idx_tk)), nrow = sum(idx_tk), byrow = TRUE))
    }
  }

  if (length(Y_k) == 0L) {
    return(list(theta0_k = intercept_k, theta_k = slopes_k,
                convergence = NA_integer_, obj_start = NA_real_,
                obj_end = NA_real_, n_outside = 0L, counts = c(NA, NA)))
  }

  ## ---- objective ------------------------------------------------------------
  g       <- make_logdens_ext(density_k)
  N_total <- sum(unlist(weights))          # matches mstep_theta_lp L1658
  pen     <- N_total * lambda_theta

  unpack <- function(par) {
    list(theta0 = par[1L],
         theta  = par[2L:(p + 1L)] - par[(p + 2L):(2L * p + 1L)])
  }

  # optim minimizes, so return the negated objective
  negobj <- function(par) {
    q <- unpack(par)
    u <- Y_k - q$theta0 - as.vector(X_k %*% q$theta)
    -(sum(w_k * g$value(u)) - pen * sum(par[-1L]))
  }

  # d/dtheta0 g(u) = -slope(u);  d/dtheta_j g(u) = -slope(u) * X_j
  # d/dtheta^+_j ||theta||_1 = d/dtheta^-_j ||theta||_1 = 1
  neggrad <- function(par) {
    q  <- unpack(par)
    u  <- Y_k - q$theta0 - as.vector(X_k %*% q$theta)
    s  <- g$slope(u)
    ws <- w_k * s
    d_theta0 <- -sum(ws)
    d_theta  <- -as.vector(crossprod(X_k, ws))
    -c(d_theta0, d_theta - pen, -d_theta - pen)
  }

  ## ---- start at the current iterate, as the LP effectively does -------------
  par0 <- c(intercept_k, pmax(slopes_k, 0), pmax(-slopes_k, 0))
  lower <- c(-Inf, rep(0, 2L * p))

  obj_start <- -negobj(par0)

  fit <- tryCatch(
    stats::optim(
      par     = par0,
      fn      = negobj,
      gr      = if (use_gradient) neggrad else NULL,
      method  = "L-BFGS-B",
      lower   = lower,
      control = list(maxit = maxit)
    ),
    error = function(e) list(par = par0, value = -obj_start,
                             convergence = 99L, counts = c(NA, NA),
                             message = conditionMessage(e))
  )

  q <- unpack(fit$par)

  ## ---- diagnostic: how far off the support did this update push us? ---------
  u_new     <- Y_k - q$theta0 - as.vector(X_k %*% q$theta)
  n_outside <- sum(u_new < g$L | u_new > g$U)

  list(
    theta0_k    = q$theta0,
    theta_k     = q$theta,
    convergence = fit$convergence,
    obj_start   = obj_start,
    obj_end     = -fit$value,
    n_outside   = n_outside,
    counts      = fit$counts
  )
}
