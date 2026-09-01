# Generated from create-LCDmix.Rmd: do not edit by hand

#' Run EM-style iterations for log-concave mixture-of-experts model
#'
#' @description
#' Performs up to \code{max_iter} iterations of the EM algorithm:
#' 1. E-step via \code{estep_lcd()}  
#' 2. M-steps:
#'    - Update mixture weights with \code{mstep_alpha()}  
#'    - Update regression slopes/intercepts via \code{mstep_theta()} and \code{mstep_shift()}  
#'    - Update log‐concave component densities with \code{mstep_g()}  
#' Tracks the surrogate log‐likelihood \code{Q} and stops early if the relative increase falls below \code{iter_eta}, or if it decreases.
#'
#' @param Y_bin List of length \code{TT}; each element is an \eqn{M_t \times 1} matrix of binned responses.
#' @param X Numeric \eqn{TT \times p} covariate matrix (rows = time points).
#' @param bin_mass List of length \code{TT}; each element is a numeric vector of length \eqn{M_t} of biomass weights.
#' @param init_res List returned by \code{initialization()}, containing:\cr
#'   \code{idx_init}, \code{resp_init}, \code{weight_init}, \code{resi_init},\cr
#'   \code{alpha_init}, \code{theta0_init}, \code{theta_init}, \code{g_init}, \code{Q_every}.
#' @param lambda_alpha Nonnegative numeric L1 penalty on mixture‐weight coefficients.
#' @param lambda_theta Nonnegative numeric L1 penalty on regression slopes.
#' @param iter_eta Numeric; relative change threshold for stopping. Default: \code{1e-6}.
#' @param max_iter Integer; maximum number of EM iterations. Default: \code{30}.
#' @param resp_threshold Numeric in [0,1]; responsibilities below this are set to zero. Default: \code{1e-3}.
#'
#' @return A list with components:
#' \describe{
#'   \item{idx_new}{Final list of \eqn{TT} logical matrices of “active” bins per component.}
#'   \item{resp_new}{Final list of \eqn{TT} responsibility matrices.}
#'   \item{weight_new}{Final list of \eqn{TT} posterior‐weight matrices.}
#'   \item{resi_new}{Final list of \eqn{TT} residual matrices.}
#'   \item{alpha_new}{Final \eqn{K \times (p+1)} mixture‐weight matrix.}
#'   \item{theta0_new}{List of length \eqn{K} of final intercepts.}
#'   \item{theta_new}{List of length \eqn{K} of final slope vectors.}
#'   \item{g_new}{List of length \eqn{K} of final log‐concave densities.}
#'   \item{lambda_alpha}{Echo of input penalty on α.}
#'   \item{lambda_theta}{Echo of input penalty on θ.}
#'   \item{Q}{Vector of surrogate log‐likelihood values at each iteration.}
#'   \item{Q_every}{Same as \code{Q}, for compatibility with initialization.}
#' }
#'
#' @export
iteration <- function(
  Y_bin,
  X,
  bin_mass,
  init_res,
  lambda_alpha,
  lambda_theta,
  iter_eta       = 1e-6,
  max_iter       = 30,
  resp_threshold = 1e-3,
  calc_Q_every   = FALSE,
  debug          = FALSE,
  lp_time_limit  = 3600,
  update         = c("lp", "optim")     # NEW - must be LAST
) {
  update <- match.arg(update)           # NEW
  TT <- nrow(X)
  p  <- ncol(X)
  K  <- length(init_res$g_init)
  
  # Unpack initial values
  idx_old     <- init_res$idx_init
  resp_old    <- init_res$resp_init
  weight_old  <- init_res$weight_init
  resi_old    <- init_res$resi_init
  alpha_old   <- init_res$alpha_init
  theta0_old  <- init_res$theta0_init
  theta_old   <- init_res$theta_init
  g_old       <- init_res$g_init
  Q           <- init_res$Q
  Q_every     <- init_res$Q_every
  n_outside_every <- integer(0)      # NEW: parallel to Q_every
  theta_diag <- list()
  lp_check_every  <- list()                       # NEW
  
  last_state <- list(
    idx_old     = idx_old,
    resp_old    = resp_old,
    weight_old  = weight_old,
    resi_old    = resi_old,
    alpha_old   = alpha_old,
    theta0_old  = theta0_old,
    theta_old   = theta_old, 
    g_old       = g_old,
    Q           = Q,
    Q_every     = Q_every,
    theta_diag  = theta_diag,
    n_outside_every = n_outside_every,
    iter_num    = 0)
  
  idx_new    <- idx_old
  resp_new   <- resp_old
  weight_new <- weight_old
  resi_new   <- resi_old
  alpha_new  <- alpha_old
  theta0_new <- theta0_old
  theta_new  <- theta_old
  g_new      <- g_old
  i          <- 0L
  
  res <- tryCatch({
  for (i in seq_len(max_iter)){
    
    #— E‐step —#
    Estep   <- estep_lcd(
      X               = X,
      bin_mass        = bin_mass,
      residuals       = resi_old,
      alpha           = alpha_old,
      densities       = g_old,
      resp_threshold  = resp_threshold
    )
    idx_new    <- Estep$idx
    resp_new   <- Estep$resp
    weight_new <- Estep$weight
#     Q(Theta^(m) | Theta^(m)) -- the reference point for the ascent test.
    if (calc_Q_every) {
     Q_new <- comp_Q(X, g_old, resi_old, theta_old, alpha_old, idx_new,
                     weight_new, lambda_alpha, lambda_theta,
                     Y_bin = Y_bin, intercepts = theta0_old)          # NEW
     Q_every         <- c(Q_every, Q_new)
     n_outside_every <- c(n_outside_every, attr(Q_new, "n_outside"))  # NEW
    }
    message("✔ E‐step complete")
    
    #— M‐step α —#
    alpha_new <- mstep_alpha(
      X                    = X,
      weights              = weight_new,
      idx                  = idx_new,
      lambda_alpha         = lambda_alpha
    )
    if (calc_Q_every) {
     Q_new <- comp_Q(X, g_old, resi_old, theta_old, alpha_new, idx_new,
                     weight_new, lambda_alpha, lambda_theta,
                     Y_bin = Y_bin, intercepts = theta0_old)          # NEW
     Q_every         <- c(Q_every, Q_new)
     n_outside_every <- c(n_outside_every, attr(Q_new, "n_outside"))  # NEW
    }
    message("✔ Updated α")
    
    #— M‐step θ via LP + shift —#
    #     THIS IS FIX 1. theta_lp$theta0 is the coefficient update's own intercept,
    #     before mstep_shift() overwrites it at L1147. Naming it explicitly keeps
    #     the distinction from being lost to a later edit.
    theta_lp <- mstep_theta(
      Y_bin = Y_bin, X = X, weights = weight_new, residuals = resi_old,
      densities = g_old, idx = idx_old, intercepts = theta0_old,
      slopes = theta_old, lambda_theta = lambda_theta,
      lp_time_limit = lp_time_limit,
      update = update                                   # NEW, for the ablation
    )
    theta0_lp  <- theta_lp$theta0    # NEW: the update's own intercept, pre-shift
    theta0_new <- theta0_lp
    theta_new  <- theta_lp$theta
    theta_diag[[i]] <- theta_lp$diag                    # NEW, for the ablation
    
    if (calc_Q_every) {
      lp_check <- lapply(seq_len(K), function(k) {
        g_ext <- make_logdens_ext(g_old[[k]])
        u <- unlist(lapply(seq_len(TT), function(t) {
          ii <- idx_old[[t]][, k]                     # the bin set the LP was given
          if (!any(ii)) return(numeric(0))
          Y_bin[[t]][ii, 1] - theta0_lp[[k]] - sum(X[t, ] * theta_new[[k]])
        }))
        tol <- 1e-8
        c(n_out = sum(u < g_ext$L - tol | u > g_ext$U + tol),
          max_over = max(0, g_ext$L - min(u), max(u) - g_ext$U))
      })
      lp_check_every[[i]] <- do.call(rbind, lp_check)
    }
    message("✔ Updated θ via LP")
    
    # Shift intercepts analytically
    theta0_new <- mstep_shift(
      Y_bin              = Y_bin,
      X                  = X,
      weights            = weight_new,
      idx                = idx_new,
      slopes             = theta_new
    )
    resi_new <- comp_resi(
      Y_bin   = Y_bin,
      X       = X,
      intercepts = theta0_new,
      slopes     = theta_new
    )
    if (calc_Q_every) {
     Q_new <- comp_Q(X, g_old, resi_new, theta_new, alpha_new, idx_new,
                     weight_new, lambda_alpha, lambda_theta,
                     Y_bin = Y_bin, intercepts = theta0_new)          # NEW
     Q_every         <- c(Q_every, Q_new)
     n_outside_every <- c(n_outside_every, attr(Q_new, "n_outside"))  # NEW
    }
    message("✔ Centered the intercepts")
    
    #— M‐step g —#
    g_new <- mstep_g(
      residuals       = resi_new,
      weights         = weight_new,
      idx             = idx_new
    )
    message("✔ Updated g")
    
    # Surrogate log-likelihood
    Q_new <- comp_Q(X, g_new, resi_new, theta_new, alpha_new, idx_new,
                   weight_new, lambda_alpha, lambda_theta)
    Q_every         <- c(Q_every, Q_new)
    Q               <- c(Q, Q_new)
    n_outside_every <- c(n_outside_every, attr(Q_new, "n_outside"))   # NEW, safe

    message("✔ Q[i] = ", round(Q_new, 6))
    
    
    # Check convergence or decrease
    inc <- (Q[i + 1] - Q[i]) / abs(Q[i])
    if (inc < 0) {
      message("⚠ Q decreased at iteration ", i,  "; reverting to previous iteration")
      idx_new    <- idx_old
      resp_new   <- resp_old
      weight_new <- weight_old
      resi_new   <- resi_old
      alpha_new  <- alpha_old
      theta0_new <- theta0_old
      theta_new  <- theta_old
      g_new      <- g_old
      break
    }
    if (inc <= iter_eta || i == max_iter) {
      message("Converged at iteration ", i, " (inc = ", round(inc, 6), ")")
      break
    }
    
    # Prepare for next iteration
    idx_old    <- idx_new
    resp_old   <- resp_new
    weight_old <- weight_new
    resi_old   <- resi_new
    alpha_old  <- alpha_new
    theta0_old <- theta0_new
    theta_old  <- theta_new
    g_old      <- g_new
    
    # Store the current parameters
    last_state <- list(
      idx_old     = idx_old,
      resp_old    = resp_old,
      weight_old  = weight_old,
      resi_old    = resi_old,
      alpha_old   = alpha_old,
      theta0_old  = theta0_old,
      theta_old   = theta_old, 
      g_old       = g_old,
      Q           = Q,
      Q_every     = Q_every,
      n_outside_every = n_outside_every,
      theta_diag      = theta_diag,
      lp_check_every  = lp_check_every,               # NEW
      iter_num        = i
      )
  }
    list(final = last_state,
         error = NULL)
    }, error = function(e){
    # on *any* error inside the loop:
    if (debug) {
      # return the last successful state + the error message
      list(
        final       = last_state,
        error       = conditionMessage(e),
        failed_iter = i
      )
    } else {
      stop(e)
    }}
  )
  return(list(
    idx_new     = idx_new,
    resp_new    = resp_new,
    weight_new  = weight_new,
    resi_new    = resi_new,
    alpha_new   = alpha_new,
    theta0_new  = theta0_new,
    theta_new   = theta_new, 
    lambda_alpha = lambda_alpha,
    lambda_theta = lambda_theta,
    g_new       = g_new,
    Q           = Q,
    Q_every     = Q_every,
    n_outside_every = n_outside_every,
    theta_diag      = theta_diag,
    lp_check_every  = lp_check_every,               # NEW
    error           = res$error,           # NEW
    failed_iter     = res$failed_iter,      # NEW
    iter_num    = i))
}
