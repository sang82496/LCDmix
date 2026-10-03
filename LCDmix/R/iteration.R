# Generated from create-LCDmix.Rmd: do not edit by hand

#' Run the EM-type iterations of the log-concave mixture-of-experts model
#'
#' @description
#' Starts from \code{init_res} and repeats, for \eqn{m = 1, 2, \ldots}:
#'
#' 1. E-step with \code{estep_lcd()}. Then
#'    \eqn{Q_{ref} = Q(\Theta^{(m)} \mid \Theta^{(m)})}, the objective at the
#'    current parameters under the new responsibilities (checkpoint A).
#' 2. Gate update with \code{mstep_alpha()} (checkpoint B).
#' 3. Expert coefficient update with \code{mstep_theta()} (checkpoint C). It
#'    uses the active bins of the previous E-step (\code{idx_old}) on purpose:
#'    the support bounds of the linear program must match the density that was
#'    fitted to those residuals, or the current point can be infeasible.
#' 4. Intercept centering with \code{mstep_shift()} (checkpoint D).
#' 5. Density update with \code{mstep_g()}.
#' 6. \eqn{Q_E = Q(\Theta^{(m+1)} \mid \Theta^{(m)})} with \code{comp_Q()}
#'    (checkpoint E).
#'
#' Stopping rule: let \eqn{inc = (Q_E - Q_{ref}) / |Q_{ref}|}. If
#' \eqn{inc < 0}, the iteration is discarded, the parameters of the previous
#' iteration are kept, and the loop stops. Otherwise the loop stops when
#' \eqn{inc \le} \code{iter_eta} or \eqn{m =} \code{max_iter}, with the new
#' parameters.
#'
#' @param init_res The list returned by \code{initialization()}.
#' @param iter_eta Stopping tolerance on the relative ascent \eqn{inc}; see the
#'   stopping rule above.
#' @inheritParams main
#' @inheritParams initialization
#'
#' @details
#' Note for developers: \code{main()} calls \code{iteration()} with positional
#' arguments. Add new arguments only at the end of the signature.
#'
#' @return A list with components:
#' \describe{
#'   \item{idx_new, resp_new, weight_new}{E-step output that belongs to the
#'     returned parameters.}
#'   \item{resi_new}{Residuals at the returned parameters.}
#'   \item{alpha_new}{\eqn{K \times (p+1)} gate coefficients; row 1 is 0.}
#'   \item{theta0_new, theta_new}{Lists of \code{K} intercepts and \code{K}
#'     slope vectors.}
#'   \item{g_new}{List of \code{K} log-concave fits.}
#'   \item{lambda_alpha, lambda_theta}{The penalties, returned unchanged.}
#'   \item{Q}{The starting \eqn{Q}, then \eqn{Q_E} of each iteration,
#'     including a last iteration that was discarded because \eqn{Q}
#'     decreased.}
#'   \item{Q_every}{With \code{calc_Q_every = FALSE}, the same as \code{Q}.
#'     With \code{TRUE}, \eqn{1 + 5 \times} \code{iter_num} values: the
#'     starting \eqn{Q}, then checkpoints A to E of each iteration.
#'     \code{Q_every[5 * m + 1] - Q_every[5 * m - 3]} is the ascent of
#'     iteration \eqn{m}.}
#'   \item{n_outside_every}{Attribute \code{n_outside} of \code{comp_Q()} at
#'     each recorded checkpoint: 5 values per iteration with
#'     \code{calc_Q_every = TRUE} (no starting value), otherwise 1.}
#'   \item{theta_diag}{List with one element per iteration: the \code{diag}
#'     list returned by \code{mstep_theta()}.}
#'   \item{lp_check_every}{With \code{calc_Q_every = TRUE}, one
#'     \eqn{K \times 2} matrix per iteration with columns \code{n_out} (bins
#'     given to the \eqn{\theta} update whose new residual lies outside the
#'     support \eqn{[L_k, U_k]}, tolerance 1e-8) and \code{max_over} (largest
#'     distance outside). Empty otherwise.}
#'   \item{error, failed_iter}{\code{NULL}, or with \code{debug = TRUE} the
#'     error message and the iteration at which the loop failed.}
#'   \item{iter_num}{Number of the last iteration that was started. If it was
#'     discarded, the returned parameters come from iteration
#'     \code{iter_num - 1}.}
#' }
#' After an error with \code{debug = TRUE}, the parameter fields hold the
#' values at the moment of the error, which can be partly updated by the
#' failed iteration.
#'
#' @export
iteration <- function(
  Y_bin,
  X,
  bin_mass,
  init_res,
  lambda_alpha,
  lambda_theta,
  iter_eta       = 1e-4,
  max_iter       = 30,
  resp_threshold = 1e-3,
  calc_Q_every   = FALSE,
  debug          = FALSE,
  lp_time_limit  = 3600,
  update         = c("lp", "optim"),    # NEW
  maxdev         = NULL                 # NEW (fixP) - appended last
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
    ## Q(Theta^(m) | Theta^(m)) -- the reference point for the ascent test.
    ## Computed unconditionally now: the stopping rule below needs it.
    Q_ref <- comp_Q(X, g_old, resi_old, theta_old, alpha_old, idx_new,
                    weight_new, lambda_alpha, lambda_theta,
                    Y_bin = Y_bin, intercepts = theta0_old)
    if (calc_Q_every) {
     Q_every         <- c(Q_every, Q_ref)
     n_outside_every <- c(n_outside_every, attr(Q_ref, "n_outside"))
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
    
    # M‐step θ via LP + shift #
    #     THIS IS FIX 1. theta_lp$theta0 is the coefficient update's own intercept,
    #     before mstep_shift() overwrites it at L1147. Naming it explicitly keeps
    #     the distinction from being lost to a later edit.
    theta_lp <- mstep_theta(
      Y_bin = Y_bin, X = X, weights = weight_new, residuals = resi_old,
      densities = g_old, idx = idx_old, intercepts = theta0_old,
      slopes = theta_old, lambda_theta = lambda_theta,
      lp_time_limit = lp_time_limit,
      update = update,                                  # NEW, for the ablation
      maxdev = maxdev                                   # NEW (fixP)
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
     Q_new <- comp_Q(X, g_old, resi_old, theta_new, alpha_new, idx_new,
                     weight_new, lambda_alpha, lambda_theta,
                     Y_bin = Y_bin, intercepts = theta0_lp)           # NEW
     Q_every         <- c(Q_every, Q_new)
     n_outside_every <- c(n_outside_every, attr(Q_new, "n_outside"))  # NEW
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
    Q_E <- Q_new                      # NEW: explicit alias, checkpoint E
    Q_every         <- c(Q_every, Q_new)
    Q               <- c(Q, Q_new)
    n_outside_every <- c(n_outside_every, attr(Q_new, "n_outside"))

    message("✔ Q[i] = ", round(Q_new, 6))
    
    
    # Check convergence or decrease
    inc <- (as.numeric(Q_E) - as.numeric(Q_ref)) / abs(as.numeric(Q_ref))
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
  }
    list(error = NULL)
    }, error = function(e){
    # on *any* error inside the loop:
    if (debug) {
      # return the error message and the iteration at which it happened
      list(
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
