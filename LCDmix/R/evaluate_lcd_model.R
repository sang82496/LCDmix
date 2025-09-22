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
#' @param trim_prob Numeric in \[0,0.5\]; proportion of lowest-likelihood observations to trim when computing the average. Default: \code{0.01}.
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
#'   trim_prob  = 0.01
#' )
#' str(eval_res)
#' }
#' @export
evaluate_lcd_model <- function(
  model,
  Y_test,
  X_test,
  biomass,
  trim_prob = 0.01
) {
  # Unpack fitted parameters
  densities     <- model$g_new
  slopes        <- model$theta_new
  intercepts    <- model$theta0_new
  alpha         <- model$alpha_new
  lambda_alpha  <- model$lambda_alpha
  lambda_theta  <- model$lambda_theta
  
  TT <- nrow(X_test)
  K  <- nrow(alpha)
  
  # Flatten biomass weights for trimming
  weights <- unlist(biomass)
  
  # Compute mixture probabilities π_{t,k}
  pi_mat <- pi_k(X_test, alpha)
  
  # Initialize vector for individual log‐likelihoods
  loglikes <- vector("list", TT)
  
  # Loop over time points and observations
  for (t in seq_len(TT)) {
    y_t   <- as.numeric(Y_test[[t]])
    n_t   <- length(y_t)
    lt    <- matrix(NA_real_, n_t, K)
    # Compute component‐wise densities f_k(r_{t,i,k})
    for (k in seq_len(K)) {
      pred_tk      <- intercepts[[k]] + sum(X_test[t, ] * slopes[[k]])
      resi_tk      <- y_t - pred_tk
      # Evaluate log‐concave density; column 3 = density
      dens_vals    <- suppressWarnings(
        logcondens::evaluateLogConDens(resi_tk, densities[[k]])[, 3]
      )
      lt[, k] <- dens_vals * pi_mat[t, k]
    }
    # Mixture density and log
    mix_dens      <- rowSums(lt)
    loglikes[[t]] <- log(mix_dens)
  }
  loglikes <- unlist(loglikes)
  
  # Proportion of infinite log‐likelihoods
  finite_mask <- is.finite(loglikes)
  prop_inf    <- mean(!finite_mask)
  
  # if all -Inf
  if (!any(finite_mask)) {
    return(list(
      loglikes            = loglikes,
      weights             = weights,
      prop_inf            = 1,
      finite_loglik       = -Inf,
      trimmed_ll          = -Inf,
      trimmed_w           = 0,
      trimmed_loglik      = -Inf
    ))
  }
    
  # Finite loglikelihood
  base_ll <- sum(weights[finite_mask] * loglikes[finite_mask]) / sum(weights[finite_mask])
  
  # Trimmed loglikelihood
  threshold <- weighted_quantile(loglikes, weights, prob = trim_prob)
  keep_idx  <- loglikes > threshold
  trimmed_ll     <- loglikes[keep_idx]
  trimmed_w      <- weights[keep_idx]
  base_trimmed   <- sum(trimmed_ll * trimmed_w) / sum(trimmed_w)
  
  # Penalty 
  l1_alpha  <- sum(abs(alpha[, -1]))
  l1_theta  <- sum(abs(unlist(slopes)))
  pen_total <- lambda_alpha * l1_alpha + lambda_theta * l1_theta
  
  finite_loglik     <- base_ll - pen_total
  trimmed_loglik    <- base_trimmed - pen_total
  
  return(list(
    loglikes            = loglikes,
    weights             = weights,
    prop_inf            = prop_inf,
    finite_loglik       = finite_loglik,
    trimmed_ll          = trimmed_ll,
    trimmed_w           = trimmed_w,
    trimmed_loglik      = trimmed_loglik
  ))
}
