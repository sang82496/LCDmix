# Generated from create-LCDmix.Rmd: do not edit by hand

#' Run EM‐style iterations for log‐concave mixture‐of‐experts model
#'
#' @description
#' Performs up to \code{max_iter} iterations of the EM algorithm:
#' 1. E‐step via \code{e_step_log_concave()}  
#' 2. M‐steps:
#'    - Update mixture weights with \code{mstep_update_alpha()}  
#'    - Update regression slopes/intercepts via \code{mstep_update_theta_log_concave()} and \code{mstep_update_intercepts()}  
#'    - Update log‐concave component densities with \code{mstep_estimate_log_concave_densities()}  
#' Tracks the surrogate log‐likelihood \code{Q} and stops early if the relative increase falls below \code{iter_eta}, or if it decreases.
#'
#' @param Y_bin List of length \code{TT}; each element is an \eqn{M_t \times 1} matrix of binned responses.
#' @param X Numeric \eqn{TT \times p} covariate matrix (rows = time points).
#' @param bin_mass List of length \code{TT}; each element is a numeric vector of length \eqn{M_t} of biomass weights.
#' @param init_res List returned by \code{initialize_model()}, containing:\cr
#'   \code{idx_init}, \code{resp_init}, \code{weight_init}, \code{resi_init},\cr
#'   \code{alpha_init}, \code{theta0_init}, \code{theta_init}, \code{g_init}, \code{Q_every}.
#' @param lambda_alpha Nonnegative numeric L1 penalty on mixture‐weight coefficients.
#' @param lambda_theta Nonnegative numeric L1 penalty on regression slopes.
#' @param iter_eta Numeric; relative change threshold for stopping. Default: \code{1e-6}.
#' @param max_iter Integer; maximum number of EM iterations. Default: \code{30}.
#' @param resp_threshold Numeric in [0,1]; responsibilities below this are set to zero. Default: \code{1e-3}.
#' @param maxdev Numeric or \code{NULL}; optional max‐deviation constraint. Default: \code{NULL}.
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
  debug          = FALSE,
  maxdev         = NULL
) {
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
    Q_every     = Q_every
  )
  
  res <- tryCatch({
  for (i in seq_len(max_iter)){
    
    #— E‐step —#
    Estep   <- e_step_log_concave(
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
    Q_new    <- compute_surrogate_loglikelihood(
      X                    = X,
      densities            = g_old,
      residuals            = resi_old,
      slopes               = theta_old,
      alpha                = alpha_old,
      indices              = idx_new,
      responsibilities     = weight_new,
      lambda_alpha         = lambda_alpha,
      lambda_theta         = lambda_theta
    )
    Q_every <- c(Q_every, Q_new)
    message("✔ E‐step complete")
    
    #— M‐step α —#
    alpha_new <- mstep_update_alpha(
      X                    = X,
      posterior_weights    = weight_new,
      responsibility_mask  = idx_new,
      lambda_alpha         = lambda_alpha
    )
    Q_new    <- compute_surrogate_loglikelihood(
      X                    = X,
      densities            = g_old,
      residuals            = resi_old,
      slopes               = theta_old,
      alpha                = alpha_new,
      indices              = idx_new,
      responsibilities     = weight_new,
      lambda_alpha         = lambda_alpha,
      lambda_theta         = lambda_theta
    )
    Q_every <- c(Q_every, Q_new)
    message("✔ Updated α")
    
    #— M‐step θ via LP + shift —#
    theta_lp <- mstep_update_theta_log_concave(
      Y_bin               = Y_bin,
      X                   = X,
      posterior_weights   = weight_new,
      residuals           = resi_old,
      densities           = g_old,
      responsibility_mask = idx_old,
      intercepts          = theta0_old,
      slopes              = theta_old,
      lambda_theta        = lambda_theta,
      maxdev              = maxdev
    )
    theta0_new <- theta_lp$theta0
    theta_new  <- theta_lp$theta
    Q_new    <- compute_surrogate_loglikelihood(
      X                    = X,
      densities            = g_old,
      residuals            = resi_old,
      slopes               = theta_new,
      alpha                = alpha_new,
      indices              = idx_new,
      responsibilities     = weight_new,
      lambda_alpha         = lambda_alpha,
      lambda_theta         = lambda_theta
    )
    Q_every <- c(Q_every, Q_new)
    message("✔ Updated θ via LP")
    
    # Shift intercepts analytically
    theta0_new <- mstep_update_intercepts(
      Y_bin              = Y_bin,
      X                  = X,
      weights            = weight_new,
      indices            = idx_new,
      slopes             = theta_new
    )
    resi_new <- compute_residuals(
      Y_bin   = Y_bin,
      X       = X,
      intercepts = theta0_new,
      slopes     = theta_new
    )
    Q_new    <- compute_surrogate_loglikelihood(
      X                    = X,
      densities            = g_old,
      residuals            = resi_new,
      slopes               = theta_new,
      alpha                = alpha_new,
      indices              = idx_new,
      responsibilities     = weight_new,
      lambda_alpha         = lambda_alpha,
      lambda_theta         = lambda_theta
    )
    Q_every <- c(Q_every, Q_new)
    
    #— M‐step g —#
    g_new <- mstep_estimate_log_concave_densities(
      residuals           = resi_new,
      weights             = weight_new,
      indices             = idx_new
    )
    
    # Surrogate log‐likelihood
    Q_new    <- compute_surrogate_loglikelihood(
      X                    = X,
      densities            = g_new,
      residuals            = resi_new,
      slopes               = theta_new,
      alpha                = alpha_new,
      indices              = idx_new,
      responsibilities     = weight_new,
      lambda_alpha         = lambda_alpha,
      lambda_theta         = lambda_theta
    )
    Q_every <- c(Q_every, Q_new)
    Q       <- c(Q, Q_new)
    message("✔ Q[i] = ", Q_new)
    
    
    # Check convergence or decrease
    inc <- (Q[i + 1] - Q[i]) / abs(Q[i])
    if (inc < 0) {
      message("⚠ Q decreased; reverting to previous iteration")
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
      message("Converged at iteration ", i, " (inc = ", inc, ")")
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
    i           = i)
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
  return(res)
}
