# Generated from create-LCDmix.Rmd: do not edit by hand

#' Penalized (and trimmed) log-likelihood for LCDmix
#'
#' @description
#' Computes the penalized average log-likelihood of the LCDmix model on a
#' binned dataset, and also returns a trimmed version that excludes only
#' non-finite rows (i.e., \code{-Inf}) from the average. L1 penalties are applied
#' to \eqn{\alpha} (excluding the intercept column) and to all slope parameters
#' \eqn{\theta}.
#'
#' @param model A fitted LCDmix model object (typically an iteration/state list)
#' containing at least:
#' \itemize{
#' \item \code{g_new}: list of length \eqn{K} with log-concave density fits
#' usable by \code{logcondens::evaluateLogConDens()}.
#' \item \code{theta_new}: slopes for the mean functions (used for the
#' \eqn{L_1} penalty).
#' \item \code{alpha_new}: \eqn{K \times (p+1)} matrix of mixture-weight
#' coefficients (first column is the intercept; not penalized).
#' \item \code{resi_new}: list of length \eqn{TT}; element \eqn{t} is an
#' \eqn{n_t \times K} matrix of residuals for time \eqn{t}.
#' \item \code{lambda_alpha}, \code{lambda_theta}: nonnegative penalty
#' weights for \eqn{\alpha} and \eqn{\theta}.
#' }
#' @param X A \eqn{TT \times p} covariate matrix (rows are time points).
#' @param bin_mass Optional list of length \eqn{TT}; element \eqn{t} is a numeric
#' vector of bin masses (weights) of length \eqn{n_t}. If \code{NULL}, unit
#' weights are used.
#'
#' @details
#' Let \eqn{\pi_{t,k}} be mixture weights given by \code{pi_k(X, alpha_new)}
#' and \eqn{f_k} the component densities given by \code{g_new}. For each
#' observation (bin midpoint) the function forms
#' \eqn{\log \sum_k \pi_{t,k}, f_k(\text{residual}{t,k})} via a stable
#' log-sum-exp.
#'
#' The actual penalized log-likelihood (\code{loglik}) is the
#' weighted average over all rows; if any row is non-finite
#' (\code{-Inf}), \code{loglik} is \code{-Inf} before penalization.
#'
#' The trimmed penalized log-likelihood (\code{trim_ll}) is the same weighted
#' average but excluding only the non-finite rows; no additional quantile
#' trimming is performed here.
#'
#' L1 penalty:
#' \deqn{\lambda\alpha |\alpha_{\text{no-int}}|1 ;+; \lambda\theta |\theta|_1,}
#' where the intercept column of \eqn{\alpha} is excluded from the \eqn{L_1}
#' norm.
#'
#' The function suppresses warnings from \code{logcondens::evaluateLogConDens()}
#' when evaluating outside the support (these yield \code{-Inf} log-densities).
#'
#' @return A list with:
#' \describe{
#' \item{\code{loglik}}{Penalized average log-likelihood using all rows; \code{-Inf}
#' if any row is non-finite (before penalties).}
#' \item{\code{trim_ll}}{Penalized average log-likelihood after dropping only
#' non-finite rows.}
#' \item{\code{prop_inf}}{Proportion of rows with non-finite log-likelihoods.}
#' }
#'
#' @examples
#' \dontrun{
#' ll <- lcd_loglikelihood(model = fit$iter, X = X, bin_mass = bin_mass)
#' ll$loglik
#' ll$trim_ll
#' ll$prop_inf
#' }
#'
#' @seealso \code{\link{pi_k}}, \code{logcondens::evaluateLogConDens}
#' @export
lcd_loglikelihood <- function(
  model,
  X,
  bin_mass = NULL
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
  weights <- unlist(bin_mass)

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

  ## trimmed: drop non-finite rows only (no quantile trimming)
  if (any(finite_mask)) {
    trim_ll_base <- sum(weights[finite_mask] * loglikes[finite_mask]) /
                    sum(weights[finite_mask])
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
    trim_loglik  = trim_loglik,    # penalized, dropping only non-finite rows
    prop_inf     = prop_inf
  ))
}
