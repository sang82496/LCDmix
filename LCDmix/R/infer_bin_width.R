# Generated from create-LCDmix.Rmd: do not edit by hand

#' Width of one bin on the binning grid.
#'
#' NOTE: binning() builds equally spaced cutpoints, so the smallest gap between
#' distinct bin centers is the bin width. Empty grid positions leave gaps that
#' are multiples of it, which is why the minimum is used and not the mean.
infer_bin_width <- function(Y_test) {
  centers <- sort(unique(as.numeric(unlist(Y_test))))
  if (length(centers) < 2L) {
    return(NA_real_)
  }
  gaps  <- diff(centers)
  width <- min(gaps)
  if (any(abs(gaps / width - round(gaps / width)) > 1e-6)) {
    warning("infer_bin_width(): the bin centers are not on an equally spaced grid")
  }
  return(width)
}
#' Evaluate held-out (penalized) log-likelihood for LCDmix
#'
#' Computes per-observation held-out log-likelihoods for a fitted LCDmix
#' \code{model} on test data \code{(Y_test, X_test)} and returns the untrimmed
#' and trimmed **penalized** weighted means. Trimming is performed over
#' \emph{all rows, including \code{-Inf} rows}, so that models with different
#' proportions of finite rows remain comparable. The trimmed threshold is a
#' weighted quantile of the log-likelihoods using \code{weighted_quantile()}.
#' The data are binned, so a bin is an interval. Scoring it by the density at
#' its center rewards a fitted density that is peaked inside the bin, which is
#' why the bin probability is used instead.
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
#' @param bin_width Numeric; the width of a bin. \code{NULL} (the default)
#'   infers it from \code{Y_test}, which is correct for the equally spaced grid
#'   that \code{binning()} builds.
#' @param score Either \code{"bin_prob"} (the default) to score each bin by its
#'   probability divided by the bin width, that is the average density over the
#'   bin, or \code{"center"} for the density at the bin center, which was the
#'   behavior before this fix.
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
  trim_prob = 0.03,
  bin_width = NULL,                       # NEW, appended
  score     = c("bin_prob", "center")     # NEW, appended
) {
  score <- match.arg(score)

  # Unpack fitted parameters
  densities     <- model$g_new
  slopes        <- model$theta_new
  intercepts    <- model$theta0_new
  alpha         <- model$alpha_new
  lambda_alpha  <- model$lambda_alpha
  lambda_theta  <- model$lambda_theta

  TT <- nrow(X_test)
  K  <- nrow(alpha)

  # NEW: the width of a bin. The data are binned, so a bin is an interval and
  # its probability, not the density at its center, is what the model predicts.
  if (score == "bin_prob" && is.null(bin_width)) {
    bin_width <- infer_bin_width(Y_test)
  }
  if (score == "bin_prob" && (is.na(bin_width) || bin_width <= 0)) {
    warning("eval_lcd(): could not infer the bin width; scoring at bin centers instead")
    score <- "center"
  }

  # Flatten biomass weights
  weights <- unlist(biomass_test)

  # Mixture probabilities pi_{t,k}
  pi_mat <- pi_k(X_test, alpha)

  # Per-observation log-likelihoods
  loglikes <- vector("list", TT)
  for (t in seq_len(TT)) {
    nt <- length(Y_test[[t]])
    lt <- matrix(NA_real_, nt, K)
    # Component-wise values f_k(r_{t,i,k})
    for (k in seq_len(K)) {
      pred_tk <- intercepts[[k]] + sum(X_test[t, ] * slopes[[k]])
      resid   <- as.numeric(Y_test[[t]]) - pred_tk
      if (score == "bin_prob") {
        # NEW: average density over the bin, from the exact distribution
        # function (column 4). It is 0 below the support and 1 above it, so a
        # bin outside the support still gives 0 and then -Inf, as before.
        upper <- suppressWarnings(
          logcondens::evaluateLogConDens(resid + bin_width / 2, densities[[k]], which = 3)[, 4]   # fixT: which = 3 fills column 4 (CDF) only
        )
        lower <- suppressWarnings(
          logcondens::evaluateLogConDens(resid - bin_width / 2, densities[[k]], which = 3)[, 4]   # fixT
        )
        dens_vals <- pmax(upper - lower, 0) / bin_width
      } else {
        # Column 3 = density; out-of-support -> 0 (log -> -Inf)
        dens_vals <- suppressWarnings(
          logcondens::evaluateLogConDens(resid, densities[[k]], which = 2)[, 3]   # fixT: which = 2 fills column 3 (density) only
        )
      }
      lt[, k] <- dens_vals * pi_mat[t, k]
    }

    # Mixture value and log
    mix_dens      <- rowSums(lt)
    loglikes[[t]] <- log(mix_dens)
  }
  loglikes <- unlist(loglikes)

  # Proportion of non-finite rows
  finite_mask <- is.finite(loglikes)
  prop_inf    <- sum(weights[!finite_mask]) / sum(weights)

  # All -Inf -> return early
  if (!any(finite_mask)) {
    return(list(
      prop_inf        = 1,
      loglik          = -Inf,
      finite_loglik   = -Inf,
      med_loglik      = -Inf,
      trimmed_loglik  = -Inf,
      penalty         = NA_real_,
      sum_w           = NA_real_,
      sum_trimmed_w   = NA_real_
    ))
  }

  # Untrimmed weighted average
  base_ll <- sum(weights * loglikes) / sum(weights)

  # Untrimmed weighted average over finite rows
  base_finite <- sum(weights[finite_mask] * loglikes[finite_mask]) /
    sum(weights[finite_mask])

  # Median log-likelihood
  base_med <- weighted_quantile(loglikes, weights, prob = 0.5)

  # Trim over all rows, with ties kept (>=)
  threshold    <- weighted_quantile(loglikes, weights, prob = trim_prob)
  keep_idx     <- loglikes >= threshold
  base_trimmed <- sum(loglikes[keep_idx] * weights[keep_idx]) / sum(weights[keep_idx])

  # L1 penalties (exclude the alpha intercept column)
  l1_alpha <- sum(abs(alpha[, -1]))
  l1_theta <- sum(abs(unlist(slopes)))
  penalty  <- lambda_alpha * l1_alpha + lambda_theta * l1_theta

  return(list(
    prop_inf            = prop_inf,
    loglik              = base_ll - penalty,
    finite_loglik       = base_finite - penalty,
    med_loglik          = base_med - penalty,
    trimmed_loglik      = base_trimmed - penalty,
    penalty             = penalty,
    sum_w               = sum(weights),
    sum_trimmed_w       = sum(weights[keep_idx])
  ))
}
