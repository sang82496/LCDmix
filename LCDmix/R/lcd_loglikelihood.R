# Generated from create-LCDmix.Rmd: do not edit by hand

#' Penalized and trimmed log-likelihood for LCDmix
#'
#' @description
#' Computes the penalized average log-likelihood for an LCDmix model and a
#' quantile-trimmed variant. Quantile trimming removes a user-specified
#' fraction of the lowest per-observation log-likelihoods (after the
#' log-sum-exp over components) before averaging. L1 penalties are applied to
#' \eqn{\alpha} (excluding its intercept column) and to all slope parameters
#' \eqn{\theta}.
#'
#' @param model A fitted LCDmix iteration/state object containing:
#' \itemize{
#' \item \code{g_new}: list of length \eqn{K} with log-concave fits
#' for each component; compatible with
#' \code{logcondens::evaluateLogConDens()}.
#' \item \code{theta_new}: slopes used for the \eqn{L_1} penalty.
#' \item \code{alpha_new}: \eqn{K \times (p+1)} mixture-weight coefficients
#' (first column is the intercept; not penalized).
#' \item \code{resi_new}: list of length \eqn{TT}; element \eqn{t} is an
#' \eqn{n_t \times K} matrix of residuals at time \eqn{t}.
#' \item \code{lambda_alpha}, \code{lambda_theta}: nonnegative penalty
#' weights for \eqn{\alpha} and \eqn{\theta}.
#' }
#' @param X A \eqn{TT \times p} covariate matrix (rows are time points).
#' @param biomass List of length \eqn{TT}; element \eqn{t} is a numeric vector
#' of bin masses (weights) of length \eqn{n_t}. (Required for weighted
#' averaging.)
#' @param trim_prob Numeric in [0,1). Fraction of observations (by weight)
#' to trim from the lower tail of the per-observation log-likelihoods
#' (after log-sum-exp across components). Default: \code{0.01}.
#'
#' @details
#' For each observation, the function forms
#' \deqn{\log \sum_{k=1}^K \pi_{t,k}(X), f_k(\text{residual}{t,k}),}
#' computed via a stable rowwise log-sum-exp. The untrimmed average is
#' taken over all rows; if any row is \code{-Inf}, the untrimmed average is set
#' to \code{-Inf}. The trimmed average discards the lowest
#' \code{trim_prob} fraction (by weight) of the per-row log-likelihoods and then
#' averages the remainder. Both figures are then reduced by the L1 penalties:
#' \deqn{\lambda\alpha |\alpha_{\text{no-int}}|1 +
#' \lambda\theta |\theta|_1.}
#'
#' Warnings from \code{logcondens::evaluateLogConDens()} when evaluating outside
#' support are suppressed (these yield \code{-Inf} log-densities).
#'
#' @return A list with:
#' \describe{
#' \item{\code{loglik}}{Penalized average log-likelihood using all rows
#' (\code{-Inf} if any row is non-finite before penalization).}
#' \item{\code{trim_loglik}}{Penalized quantile-trimmed average
#' log-likelihood (lower \code{trim_prob} fraction removed).}
#' \item{\code{prop_inf}}{Proportion of rows with non-finite per-row
#' log-likelihoods.}
#' }
#'
#' @examples
#' \dontrun{
#' out <- lcd_loglikelihood(model = fit$iter, X = X, biomass = bin_mass, trim_prob = 0.05)
#' out$loglik
#' out$trim_loglik
#' out$prop_inf
#' }
#'
#' @seealso \code{\link{pi_k}}, \code{logcondens::evaluateLogConDens}
#' @export
lcd_loglikelihood <- function(
  model,
  X,
  biomass = NULL,
  trim_prob = 0.01
) {
  densities    <- model$g_new             # list length K (log-concave fits)
  slopes       <- model$theta_new         # penalized slopes (any shape)
  alpha        <- model$alpha_new         # K × (p+1)
  resi         <- model$resi_new          # list length TT; each n_t × K residuals
  lambda_alpha <- model$lambda_alpha
  lambda_theta <- model$lambda_theta

  TT <- nrow(X)
  K  <- nrow(alpha)

  # 1) π_{t,k} and log π
  pi_mat <- pi_k(X, alpha)                # TT × K
  log_pi <- log(pi_mat)

  # 2) Build per-observation log-terms (N × K after rbind)
  log_terms_list <- vector("list", TT)
  for (t in seq_len(TT)) {
    res_t <- resi[[t]]                    # n_t × K
    n_t   <- nrow(res_t)
    lt    <- matrix(NA_real_, n_t, K)
    for (k in seq_len(K)) {
      # evaluateLogConDens(...)[,2] gives log-density; out-of-support → -Inf (with warnings)
      ldens <- suppressWarnings(
        logcondens::evaluateLogConDens(res_t[, k], densities[[k]])[, 2]
      )
      lt[, k] <- log_pi[t, k] + ldens
    }
    log_terms_list[[t]] <- lt
  }
  log_terms <- do.call(rbind, log_terms_list)   # N × K

  # 3) Weights
  weights <- unlist(biomass)

  # 4) Stable rowwise log-sum-exp
  any_finite_comp <- rowSums(is.finite(log_terms)) > 0
  loglikes        <- rep(-Inf, nrow(log_terms))
  if (any(any_finite_comp)) {
    M <- apply(log_terms[any_finite_comp, , drop = FALSE], 1, max)
    loglikes[any_finite_comp] <-
      M + log(rowSums(exp(log_terms[any_finite_comp, , drop = FALSE] - M)))
  }

  # 5) Base and trimmed summaries
  finite_mask <- is.finite(loglikes)
  prop_inf    <- mean(!finite_mask)

  ## actual (untrimmed): if any non-finite exists → -Inf
  if (any(!finite_mask)) {
    base_ll <- -Inf
  } else {
    base_ll <- sum(weights * loglikes) / sum(weights)
  }

  ## trimmed:
  if (any(finite_mask)) {
    threshold  <- weighted_quantile(loglikes, weights, prob = trim_prob)
    keep_idx   <- loglikes > threshold
    trimmed_ll <- loglikes[keep_idx]
    trimmed_w  <- weights[keep_idx]
    
    trim_ll_base <- sum(trimmed_ll * trimmed_w) / sum(trimmed_w)
  } else {
    trim_ll_base <- -Inf
  }
  
  # 6) L1 penalties (α intercept column is the first; K >= 2 assumed)
  l1_alpha  <- sum(abs(alpha[, -1]))
  l1_theta  <- sum(abs(unlist(slopes)))
  pen_total <- lambda_alpha * l1_alpha + lambda_theta * l1_theta

  loglik      <- base_ll      - pen_total
  trim_loglik <- trim_ll_base - pen_total

  return(list(
    loglik       = loglik,     # penalized, using all rows; -Inf if any -Inf present
    trim_loglik  = trim_loglik,    # penalized, trimmed
    prop_inf     = prop_inf
  ))
}
