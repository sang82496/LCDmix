# Generated from create-LCDmix.Rmd: do not edit by hand

#' Evaluate held-out (penalized) log-likelihood for LCDmix
#'
#' Computes per-observation held-out log-likelihoods for a fitted LCDmix
#' \code{model} on test data \code{(Y_test, X_test)} and returns the untrimmed
#' and trimmed **penalized** weighted means. Trimming is performed over
#' \emph{all rows, including \code{-Inf} rows}, so that models with different
#' proportions of finite rows remain comparable. The trimmed threshold is a
#' weighted quantile of the log-likelihoods using \code{weighted_quantile()}.
#'
#' @param model List; fitted LCDmix object containing at least:
#'   \code{g_new} (list of K log-concave density fits),
#'   \code{theta_new} (list/array of slopes),
#'   \code{theta0_new} (list of intercepts),
#'   \code{alpha_new} (K × (p+1) gating coefficients),
#'   \code{lambda_alpha}, \code{lambda_theta} (penalties).
#' @param Y_test List of length \eqn{TT}; test responses for each time point
#'   (\eqn{t = 1,\dots,TT}). Element \eqn{t} is a numeric vector of length
#'   \eqn{n_t}.
#' @param X_test Numeric matrix \eqn{TT \times p}; covariates aligned by time.
#'   Row \eqn{t} is used with \code{Y_test[[t]]}.
#' @param biomass_test List of length \eqn{TT}; per-time weights (flattened to
#'   a vector of length \eqn{N = \sum_t n_t} via \code{unlist(biomass_test)}).
#' @param trim_prob Numeric in \eqn{[0,1)}; fraction of weight to trim based on
#'   the weighted quantile of \emph{all} log-likelihoods (including \code{-Inf}).
#'
#' @details
#' For each time \eqn{t} and component \eqn{k}, residuals are
#' \eqn{r_{t,i,k} = y_{t,i} - (\theta_{0k} + x_t^\top \theta_k)}.
#' Densities are evaluated via
#' \code{logcondens::evaluateLogConDens(r_{t,i,k}, g_k)[,3]} (density column).
#' Mixture densities are \eqn{\sum_k \pi_{t,k} f_k(r_{t,i,k})}, with
#' \eqn{\pi_{t,k}} from \code{pi_k(X_test, alpha)}; per-row log-likelihoods are
#' the \eqn{\log} of those mixture densities. Penalization subtracts
#' \eqn{\lambda_\alpha \lVert \alpha_{\cdot,-1}\rVert_1 + \lambda_\theta \sum_k \lVert \theta_k\rVert_1}.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{prop_inf}}{Proportion of rows where the per-row log-likelihood is not finite.}
#'   \item{\code{finite_loglik}}{Untrimmed \emph{penalized} weighted mean over finite rows; \code{-Inf} if any \code{-Inf} present.}
#'   \item{\code{trimmed_loglik}}{Trimmed \emph{penalized} weighted mean (trim over all rows, ties kept with \code{>=}).}
#' }
#'
#' @seealso \code{\link{pi_k}}, \code{\link{weighted_quantile}},
#'   \code{\link[logcondens]{evaluateLogConDens}}
#'
#' @examples
#' \dontrun{
#' res <- eval_lcd(
#'   model        = fit$iter,
#'   Y_test       = Y_bin[test_idx],
#'   X_test       = X[test_idx, , drop = FALSE],
#'   biomass_test = bin_mass[test_idx],
#'   trim_prob    = 0.03
#' )
#' res$trimmed_loglik
#' }
#'
#' @export
eval_lcd <- function(
  model,
  Y_test,
  X_test,
  biomass_test,
  trim_prob = 0.03
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

  # Flatten biomass weights
  weights <- unlist(biomass_test)

  # Mixture probabilities π_{t,k}
  pi_mat <- pi_k(X_test, alpha)

  # Per-observation log-likelihoods
  loglikes <- vector("list", TT)
  for (t in seq_len(TT)) {
    y_t <- as.numeric(Y_test[[t]])
    n_t <- length(y_t)
    lt  <- matrix(NA_real_, n_t, K)

    # Component-wise densities f_k(r_{t,i,k})
    for (k in seq_len(K)) {
      pred_tk   <- intercepts[[k]] + sum(X_test[t, ] * slopes[[k]])
      resi_tk   <- y_t - pred_tk
      # Column 3 = density; out-of-support → 0 (log -> -Inf)
      dens_vals <- suppressWarnings(
        logcondens::evaluateLogConDens(resi_tk, densities[[k]])[, 3]
      )
      lt[, k] <- dens_vals * pi_mat[t, k]
    }

    # Mixture density and log
    mix_dens      <- rowSums(lt)
    loglikes[[t]] <- log(mix_dens)
  }
  loglikes <- unlist(loglikes)

  # Proportion of non-finite rows
  finite_mask <- is.finite(loglikes)
  #prop_inf    <- mean(!finite_mask)
  prop_inf    <- sum(weights[!finite_mask])/sum(weights)

  # All -Inf → return early
  if (!any(finite_mask)) {
    return(list(
      prop_inf       = 1,
      finite_loglik  = -Inf,
      trimmed_loglik = -Inf
    ))
  }

  # Untrimmed weighted average over finite rows
  base_ll <- sum(weights[finite_mask] * loglikes[finite_mask]) /
             sum(weights[finite_mask])

  # Trim over all rows (including -Inf), with ties kept (>=)
  threshold    <- weighted_quantile(loglikes, weights, prob = trim_prob)
  keep_idx     <- loglikes >= threshold
  trimmed_ll   <- loglikes[keep_idx]
  trimmed_w    <- weights[keep_idx]
  base_trimmed <- sum(trimmed_ll * trimmed_w) / sum(trimmed_w)

  # L1 penalties (exclude α intercept column)
  l1_alpha  <- sum(abs(alpha[, -1]))
  l1_theta  <- sum(abs(unlist(slopes)))
  pen_total <- lambda_alpha * l1_alpha + lambda_theta * l1_theta

  return(list(
    prop_inf       = prop_inf,
    finite_loglik  = base_ll       - pen_total,
    trimmed_loglik = base_trimmed  - pen_total
  ))
}
