# Generated from create-LCDmix.Rmd: do not edit by hand

#' Evaluate fitted LCD mixture model on new data
#'
#' @description
#' Computes individual log‐likelihoods and a trimmed weighted average log‐likelihood
#' for a fitted log‐concave mixture‐of‐experts model on new observations.
#'
#' @param model A list returned by \code{main()} or \code{iteration()}, containing:
#'   \describe{
#'     \item{\code{g_new}}{List of length \eqn{K} of log‐concave density objects.}
#'     \item{\code{theta_new}}{List of length \eqn{K} of slope vectors \(\theta_k\).}
#'     \item{\code{theta0_new}}{List of length \eqn{K} of intercepts \(\theta_{0k}\).}
#'     \item{\code{alpha_new}}{Numeric \eqn{K \times (p+1)} matrix of mixture parameters.}
#'     \item{\code{lambda_alpha}}{Numeric; L1 penalty used on mixture‐weight coefficients.}
#'     \item{\code{lambda_theta}}{Numeric; L1 penalty used on slope coefficients.}
#'   }
#' @param Y A list of length \code{TT}, where each element is an \eqn{n_t}-row vector or single‐column matrix of observed responses at time \(t\).
#' @param X A numeric \eqn{TT \times p} matrix of covariates (rows = time points).
#' @param biomass A list of length \code{TT}, where each element is a numeric vector of biomass weights for each observation at time \(t\).
#' @param trim_prob Numeric in \[0,0.5\]; proportion of lowest‐likelihood observations to trim when computing the average. Default: \code{0.05}.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{loglike}}{Numeric vector of individual log‐likelihoods \(\log\bigl(\sum_{k} \pi_{t,k} f_k(r_{t,i,k})\bigr)\).}
#'   \item{\code{weights}}{Numeric vector of biomass weights corresponding to each loglike entry.}
#'   \item{\code{prop_infinite}}{Proportion of \code{loglike} values that are \code{Inf}.}
#'   \item{\code{trimmed_loglike}}{Numeric vector of log‐likelihoods above the trimming threshold.}
#'   \item{\code{trimmed_weights}}{Numeric vector of weights corresponding to \code{trimmed_loglike}.}
#'   \item{\code{trimmed_avg_loglike}}{Numeric scalar; weighted average of \code{trimmed_loglike} minus L1 penalties.}
#' }
#'
#' @examples
#' \dontrun{
#' # Assume `fit` is returned by main()
#' eval_res <- evaluate_lcd_model(
#'   model      = fit,
#'   Y          = Y_test,
#'   X          = X_test,
#'   biomass    = biomass_test,
#'   trim_prob  = 0.05
#' )
#' str(eval_res)
#' }
#' @export
evaluate_lcd_model <- function(
  model,
  Y,
  X,
  biomass,
  trim_prob = 0.05
) {
  # Unpack fitted parameters
  densities     <- model$g_new
  slopes_list   <- model$theta_new
  intercepts    <- model$theta0_new
  alpha_mat     <- model$alpha_new
  lambda_alpha  <- model$lambda_alpha
  lambda_theta  <- model$lambda_theta
  
  TT     <- nrow(X)
  K_comp <- length(densities)
  
  # Flatten biomass weights for trimming
  weights_all <- unlist(biomass)
  
  # Compute mixture probabilities π_{t,k}
  pi_mat <- pi_k(X, alpha_mat)
  
  # Initialize vector for individual log‐likelihoods
  loglike_vals <- numeric(0)
  
  # Loop over time points and observations
  for (t in seq_len(TT)) {
    y_t   <- as.numeric(Y[[t]])
    n_t   <- length(y_t)
    # Compute component‐wise densities f_k(r_{t,i,k})
    comp_dens <- matrix(0, nrow = n_t, ncol = K_comp)
    for (k in seq_len(K_comp)) {
      # Residuals under component k
      pred_tk      <- intercepts[[k]] + sum(X[t, ] * slopes_list[[k]])
      resi_tk_vec  <- y_t - pred_tk
      # Evaluate log‐concave density; column 3 = density
      dens_vals    <- suppressWarnings(
        logcondens::evaluateLogConDens(resi_tk_vec, densities[[k]])[, 3]
      )
      comp_dens[, k] <- dens_vals * pi_mat[t, k]
    }
    # Mixture density and log
    mix_dens  <- rowSums(comp_dens)
    loglike_t <- log(mix_dens)
    loglike_vals <- c(loglike_vals, loglike_t)
  }
  
  # Proportion of infinite log‐likelihoods
  prop_inf <- mean(is.infinite(loglike_vals))
  
  # Determine trimming threshold on loglike
  threshold <- weighted_quantile(loglike_vals, weights_all, prob = trim_prob)
  keep_idx  <- loglike_vals > threshold
  
  trimmed_ll     <- loglike_vals[keep_idx]
  trimmed_w      <- weights_all[keep_idx]
  avg_trim_ll    <- sum(trimmed_ll * trimmed_w) / sum(trimmed_w) -
                    lambda_alpha * sum(abs(alpha_mat[, -1])) -
                    lambda_theta * sum(abs(unlist(slopes_list)))
  
  return(list(
    loglike             = loglike_vals,
    weights             = weights_all,
    prop_infinite       = prop_inf,
    trimmed_loglike     = trimmed_ll,
    trimmed_weights     = trimmed_w,
    trimmed_avg_loglike = avg_trim_ll
  ))
}
