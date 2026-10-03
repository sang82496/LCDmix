# Generated from create-LCDmix.Rmd: do not edit by hand

#' Compute the surrogate objective Q
#'
#' @description
#' Evaluates the normalized surrogate objective
#' \deqn{Q = \frac{1}{N}\sum_{t}\sum_{k}\sum_{i \in I_{tk}} w_{tik}\left[\log \hat f_k(u_{tik}) + \log \pi_{tk}\right] - \lambda_\alpha \sum_k \sum_{j \ge 1} |\alpha_{kj}| - \lambda_\theta \sum_k \|\theta_k\|_1,}
#' where \eqn{I_{tk}} is the set of active bins (\code{idx}), \eqn{w_{tik}} are
#' the posterior weights (\code{resp}), \eqn{N} is the sum of all entries of
#' \code{resp}, \eqn{\hat f_k} is the fitted density of component \eqn{k},
#' \eqn{u_{tik}} is the residual and \eqn{\pi_{tk}} is the gate probability.
#' The intercept column of \code{alpha} is not penalized. The log-density is
#' evaluated at the residual of the bin center; it is not integrated over the
#' bin (\code{eval_lcd()} scores bins by their probability).
#'
#' Residuals may be supplied directly through \code{residuals}, or recomputed
#' from \code{(intercepts, slopes)} when both \code{Y_bin} and
#' \code{intercepts} are given. The second form exists because \code{slopes}
#' otherwise enters only through the L1 penalty: passing old residuals with
#' new slopes measures the change in the penalty and nothing else, which makes
#' the checkpoint diagnostics in \code{iteration()} meaningless.
#'
#' @param X A numeric \eqn{TT \times p} covariate matrix.
#' @param densities A list of length \eqn{K} of \code{modified_logcondens()}
#'   fits.
#' @param residuals A list of length \eqn{TT}, each an \eqn{M_t \times K}
#'   matrix of residuals. Ignored when both \code{Y_bin} and
#'   \code{intercepts} are supplied.
#' @param slopes A list of length \eqn{K} of slope vectors.
#' @param alpha A numeric \eqn{K \times (p+1)} matrix of gate coefficients.
#' @param idx A list of length \eqn{TT}, each an \eqn{M_t \times K} logical
#'   matrix of active bins.
#' @param resp A list of length \eqn{TT}, each an \eqn{M_t \times K} matrix of
#'   posterior weights \eqn{w_{tik}}.
#' @param lambda_alpha Nonnegative L1 penalty on the non-intercept columns of
#'   \code{alpha}.
#' @param lambda_theta Nonnegative L1 penalty on \code{slopes}.
#' @param Y_bin Optional list of length \eqn{TT} of bin centers. Supply it
#'   together with \code{intercepts} to recompute the residuals.
#' @param intercepts Optional list of length \eqn{K} of intercepts. Supply it
#'   together with \code{Y_bin} to recompute the residuals.
#'
#' @details
#' Note for developers: \code{iteration()} and \code{initialization()} call
#' \code{comp_Q()} with positional arguments. Add new arguments only at the
#' end of the signature, and keep the return value a single number.
#'
#' @return A single number: the normalized surrogate log-likelihood minus the
#'   L1 penalties. It carries four attributes, which ordinary arithmetic
#'   ignores:
#'   \describe{
#'     \item{\code{n_eval}}{Number of active bins evaluated.}
#'     \item{\code{n_outside}}{Number of those whose residual fell outside the
#'       support of the fitted density, so the log-density was \code{-Inf}
#'       and the bin was dropped from the sum.}
#'     \item{\code{mass_outside}}{Total posterior weight of the dropped bins.}
#'     \item{\code{max_over}}{Largest distance by which a dropped residual lies
#'       outside the support (0 if no bin was dropped).}
#'   }
#'   Bins outside the support are excluded from \eqn{Q} without an error. A
#'   larger \code{n_outside} therefore means that \eqn{Q} is computed over
#'   less of the data, and \eqn{Q} values with different \code{n_outside} are
#'   not comparable.
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
#' ## (a) residuals supplied directly
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
  max_over     <- 0     # ADD

  for (k in seq_len(K_comp)) {
    for (t in seq_len(TT)) {
      # Select bins for component k at time t
      idx_tk  <- idx[[t]][, k]
      w_tk    <- resp[[t]][idx_tk, k]
      resi_tk <- residuals[[t]][idx_tk, k]

      if (length(w_tk) == 0) next

      # Evaluate log-density for each residual under component k
      log_dens <- suppressWarnings(
      logcondens::evaluateLogConDens(resi_tk, densities[[k]], which = 1)[, 2]   # fixT: which = 1 fills column 2 (log-density) only
      )
      finite_mask <- is.finite(log_dens)
      
      if (any(!finite_mask)) {
        rng <- range(densities[[k]]$x)
        over <- pmax(rng[1] - resi_tk[!finite_mask], resi_tk[!finite_mask] - rng[2])
        max_over <- max(max_over, over)
      }

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

  Q_val <- unname(total_ll / N_total - pen_total)

  ## ---- NEW: diagnostics ride along as attributes; return stays scalar ------
  attr(Q_val, "n_eval")       <- n_eval
  attr(Q_val, "n_outside")    <- n_outside
  attr(Q_val, "mass_outside") <- mass_outside
  attr(Q_val, "max_over")     <- max_over     # ADD

  return(Q_val)
}
