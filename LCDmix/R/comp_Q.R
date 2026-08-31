# Generated from create-LCDmix.Rmd: do not edit by hand

#' Compute the surrogate log-likelihood Q
#'
#' @description
#' Evaluates the (normalized) surrogate log-likelihood
#' \deqn{Q = \frac{1}{N}\sum_{t}\sum_{k}\sum_{i \in \mathcal{I}_{tk}}
#'   w_{tik}\left[\log \hat f_k(u_{tik}) + \log \pi_k(X_t)\right]
#'   - \lambda_\alpha \|\alpha\|_1 - \lambda_\theta \|\theta\|_1 .}
#'
#' Residuals may either be supplied directly through \code{residuals}, or
#' recomputed internally from \code{(intercepts, slopes)} when both
#' \code{Y_bin} and \code{intercepts} are given. The second form exists because
#' \code{slopes} otherwise enters only through the L1 penalty: passing stale
#' residuals alongside updated slopes measures the penalty change and nothing
#' else, which makes the M-step diagnostic in \code{iteration()} vacuous.
#'
#' @param X A numeric \eqn{TT \times p} covariate matrix (rows = time points).
#' @param densities A list of length \eqn{K} of \code{modified_logcondens()} objects.
#' @param residuals A list of length \eqn{TT}, each an \eqn{M_t \times K} matrix of
#'   residuals. Ignored when both \code{Y_bin} and \code{intercepts} are supplied.
#' @param slopes A list of length \eqn{K} of slope vectors.
#' @param alpha A numeric \eqn{K \times (p+1)} matrix of gating parameters.
#' @param idx A list of length \eqn{TT}, each an \eqn{M_t \times K} logical matrix.
#' @param resp A list of length \eqn{TT}, each an \eqn{M_t \times K} matrix of
#'   posterior weights \eqn{w_{tik}}.
#' @param lambda_alpha Nonnegative numeric L1 penalty on non-intercept columns of \code{alpha}.
#' @param lambda_theta Nonnegative numeric L1 penalty on \code{slopes}.
#' @param Y_bin Optional list of length \eqn{TT} of binned responses. Supply
#'   together with \code{intercepts} to recompute residuals internally.
#' @param intercepts Optional list of length \eqn{K} of intercepts. Supply
#'   together with \code{Y_bin} to recompute residuals internally.
#'
#' @return A single numeric: the normalized surrogate log-likelihood minus the
#'   L1 penalties. Carries three attributes, used by the LP ablation and ignored
#'   by ordinary arithmetic:
#'   \describe{
#'     \item{\code{n_eval}}{Number of active bins evaluated.}
#'     \item{\code{n_outside}}{Number of those whose residual fell outside the
#'       fitted support, so the log-density was \code{-Inf} and the bin was
#'       dropped from the sum.}
#'     \item{\code{mass_outside}}{Total posterior weight of the dropped bins.}
#'   }
#'   Bins outside the support are silently excluded from \eqn{Q}; a larger
#'   \code{n_outside} therefore means \eqn{Q} is computed over less of the data,
#'   and \eqn{Q} values with different \code{n_outside} are not comparable.
#'
#' @examples
#' \dontrun{
#' TT <- 3; p <- 2; K <- 2
#' X  <- matrix(rnorm(TT * p), nrow = TT)
#' densities <- replicate(K,
#'   modified_logcondens(rnorm(50), w = rep(1/50, 50)),
#'   simplify = FALSE
#' )
#' residuals <- lapply(1:TT, function(t) matrix(rnorm(5 * K), ncol = K))
#' idx    <- lapply(residuals, function(m) matrix(TRUE, nrow(m), ncol(m)))
#' resp   <- lapply(idx, function(ii) matrix(runif(length(ii)), nrow = nrow(ii)))
#' slopes <- replicate(K, rnorm(p), simplify = FALSE)
#' alpha  <- matrix(rnorm(K * (p + 1)), nrow = K)
#'
#' ## (a) residuals supplied directly -- unchanged behavior
#' Q_val <- comp_Q(
#'   X, densities, residuals, slopes, alpha,
#'   idx, resp,
#'   lambda_alpha = 1e-3, lambda_theta = 1e-3
#' )
#'
#' ## (b) residuals recomputed from (intercepts, slopes)
#' Y_bin      <- lapply(1:TT, function(t) matrix(rnorm(5), ncol = 1))
#' intercepts <- replicate(K, rnorm(1), simplify = FALSE)
#' Q_val2 <- comp_Q(
#'   X, densities, residuals, slopes, alpha,
#'   idx, resp,
#'   lambda_alpha = 1e-3, lambda_theta = 1e-3,
#'   Y_bin = Y_bin, intercepts = intercepts
#' )
#' attr(Q_val2, "n_outside")
#' }
#' @export
comp_Q <- function(
  X,
  densities,
  residuals,
  slopes,
  alpha,
  idx,
  resp,
  lambda_alpha,
  lambda_theta,
  Y_bin      = NULL,   # NEW, appended
  intercepts = NULL    # NEW, appended
) {
  # Number of time points
  TT     <- nrow(X)
  K_comp <- length(densities)

  ## ---- NEW: recompute residuals when both new arguments are supplied -------
  ## comp_resi() is defined earlier in the file, so it is always available here.
  if (!is.null(Y_bin) && !is.null(intercepts)) {
    if (length(intercepts) != K_comp || length(slopes) != K_comp) {
      stop("comp_Q(): `intercepts` and `slopes` must both have length ",
           K_comp, " to recompute residuals.")
    }
    residuals <- comp_resi(Y_bin, X, intercepts, slopes)
  } else if (xor(is.null(Y_bin), is.null(intercepts))) {
    warning("comp_Q(): `Y_bin` and `intercepts` must be supplied together; ",
            "falling back to the `residuals` argument.")
  }

  # Total "effective sample size"
  N_total <- sum(unlist(resp))

  # Mixture probabilities: TT x K matrix
  pi_mat <- pi_k(X, alpha)

  total_ll     <- 0
  n_eval       <- 0L   # NEW: active bins evaluated
  n_outside    <- 0L   # NEW: of those, how many fell outside the support
  mass_outside <- 0     # NEW: their total posterior weight

  for (k in seq_len(K_comp)) {
    for (t in seq_len(TT)) {
      # Select bins for component k at time t
      idx_tk  <- idx[[t]][, k]
      w_tk    <- resp[[t]][idx_tk, k]
      resi_tk <- residuals[[t]][idx_tk, k]

      if (length(w_tk) == 0) next

      # Evaluate log-density for each residual under component k
      log_dens <- suppressWarnings(
        logcondens::evaluateLogConDens(resi_tk, densities[[k]])[, 2]
      )
      finite_mask <- is.finite(log_dens)

      # NEW: record how much of the data this Q is actually computed over
      n_eval       <- n_eval + length(w_tk)
      n_outside    <- n_outside + sum(!finite_mask)
      mass_outside <- mass_outside + sum(w_tk[!finite_mask])

      # Accumulate weighted log-density and mixture log-prob terms
      total_ll <- total_ll +
        sum(w_tk[finite_mask] * log_dens[finite_mask]) +
        sum(w_tk[finite_mask]) * log(pi_mat[t, k])
    }
  }

  # Normalize and subtract L1 penalties
  l1_alpha  <- sum(abs(alpha[, -1]))
  l1_theta  <- sum(abs(unlist(slopes)))
  pen_total <- lambda_alpha * l1_alpha + lambda_theta * l1_theta

  Q_val <- total_ll / N_total - pen_total

  ## ---- NEW: diagnostics ride along as attributes; return stays scalar ------
  attr(Q_val, "n_eval")       <- n_eval
  attr(Q_val, "n_outside")    <- n_outside
  attr(Q_val, "mass_outside") <- mass_outside

  return(Q_val)
}
