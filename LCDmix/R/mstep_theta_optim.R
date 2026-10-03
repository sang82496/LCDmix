# Generated from create-LCDmix.Rmd: do not edit by hand

#' Update the coefficients of one expert by quasi-Newton (LP ablation arm)
#'
#' @description
#' Comparison arm of the LP ablation. Maximizes the same weighted
#' log-likelihood with the same L1 penalty as \code{mstep_theta_lp()}, but
#' with \code{optim(method = "L-BFGS-B")} and without the support constraints.
#' The log-density is evaluated with \code{make_logdens_ext()}, which
#' continues it linearly beyond the support.
#'
#' The slope vector is split as \eqn{\theta = \theta^+ - \theta^-} with both
#' parts nonnegative, as in the linear program, so that the two arms differ in
#' the solver and not in the problem, and the L1 penalty is handled the same
#' way in both.
#'
#' @inheritParams mstep_theta_lp
#' @param residuals Not used; an argument so that both updates take the same
#'   arguments.
#' @param intercept_k Current intercept; the starting value of \code{optim}.
#' @param slopes_k Current slope vector; the starting value of \code{optim}.
#' @param lp_time_limit Not used; an argument so that both updates take the
#'   same arguments.
#' @param maxit Maximum number of \code{optim} iterations.
#' @param use_gradient Logical; supply the analytic (sub)gradient to
#'   \code{optim}.
#'
#' @return A list with \code{theta0_k}, \code{theta_k} and the diagnostics
#'   \code{convergence} (code from \code{optim}; 99 if \code{optim} stopped
#'   with an error), \code{obj_start}, \code{obj_end} (objective before and
#'   after), \code{n_outside} (active bins whose new residual lies outside the
#'   support) and \code{counts} (from \code{optim}). If the component has no
#'   active bin, the current values are returned with \code{NA} diagnostics.
#'   \code{mstep_theta()} uses only \code{theta0_k} and \code{theta_k} for the
#'   fit.
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
