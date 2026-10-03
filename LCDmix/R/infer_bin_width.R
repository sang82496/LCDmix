# Generated from create-LCDmix.Rmd: do not edit by hand

#' Width of one bin on the binning grid
#'
#' @description
#' \code{binning()} builds equally spaced cut points, so the smallest gap
#' between distinct bin centers is the bin width. Empty bins leave gaps that
#' are multiples of it, which is why the minimum is used and not the mean.
#' Gives a warning if a gap is not a whole multiple of the minimum (relative
#' tolerance 1e-6).
#'
#' @param Y_test A list of bin centers (vectors or one-column matrices), as in
#'   \code{Y_bin}.
#'
#' @return The bin width, or \code{NA} if there are fewer than two distinct
#'   bin centers.
#' @keywords internal
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
#' Log-likelihood of a fitted LCDmix model on binned data
#'
#' @description
#' Computes, for every bin \eqn{i} at every time point \eqn{t}, the log of the
#' fitted mixture density,
#' \deqn{\ell_{ti} = \log \sum_k \pi_{tk} \hat f_{k,ti},}
#' and returns weighted means of \eqn{\ell_{ti}} (weights
#' \code{biomass_test}) minus the L1 penalty of the model. With
#' \code{score = "bin_prob"} (the default), \eqn{\hat f_{k,ti}} is the
#' probability of the bin under component \eqn{k} divided by the bin width
#' \eqn{h}, \eqn{(\hat F_k(u + h/2) - \hat F_k(u - h/2)) / h}, where \eqn{u} is
#' the residual of the bin center. With \code{score = "center"} it is the
#' density at the bin center. The data are binned, so a bin is an interval;
#' scoring by the density at the center rewards a component that is peaked
#' inside one bin, which is why the bin probability is the default. A bin
#' outside the support of every component gets \eqn{\ell_{ti} = -\infty}.
#'
#' \code{main()} calls this function on the training data (the result is
#' \code{fit$L}); \code{cv_lcd_onejob()} calls it on the held-out time points.
#'
#' @param model A fitted model: the \code{iter} element of a \code{main()} fit.
#'   The fields \code{g_new}, \code{theta0_new}, \code{theta_new},
#'   \code{alpha_new}, \code{lambda_alpha} and \code{lambda_theta} are used.
#' @param Y_test A list of length \eqn{TT}; element \eqn{t} holds the bin
#'   centers at time \eqn{t} (a vector or a one-column matrix).
#' @param X_test A numeric \eqn{TT \times p} covariate matrix; row \eqn{t}
#'   goes with \code{Y_test[[t]]}.
#' @param biomass_test A list of length \eqn{TT}; the weight of each bin.
#' @param trim_prob A number in \eqn{[0, 1)}: the fraction of the total weight
#'   trimmed from the bottom for \code{trimmed_loglik}. The threshold is the
#'   weighted \code{trim_prob}-quantile (\code{weighted_quantile()}) of all
#'   \eqn{\ell_{ti}}, \code{-Inf} included, so models with different numbers of
#'   \code{-Inf} bins are trimmed in the same way.
#' @param bin_width The width of a bin. \code{NULL} (the default) infers it
#'   from \code{Y_test} with \code{infer_bin_width()}, which is correct for the
#'   equally spaced grid that \code{binning()} builds. If it cannot be
#'   inferred, the function warns and scores at the bin centers.
#' @param score \code{"bin_prob"} (the default) to score each bin by its
#'   probability divided by the bin width, that is, the average density over
#'   the bin; \code{"center"} for the density at the bin center.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{prop_inf}}{Share of the total weight in bins with
#'     \eqn{\ell_{ti} = -\infty}.}
#'   \item{\code{loglik}}{Weighted mean of \eqn{\ell_{ti}} over all bins, minus
#'     the penalty; \code{-Inf} if any bin has \eqn{\ell_{ti} = -\infty}.}
#'   \item{\code{finite_loglik}}{The same over the bins with finite
#'     \eqn{\ell_{ti}} only.}
#'   \item{\code{med_loglik}}{Weighted median of \eqn{\ell_{ti}}, minus the
#'     penalty.}
#'   \item{\code{trimmed_loglik}}{Weighted mean of \eqn{\ell_{ti}} over the
#'     bins at or above the trimming threshold, minus the penalty.}
#'   \item{\code{penalty}}{\eqn{\lambda_\alpha \sum_k \sum_{j \ge 1} |\alpha_{kj}| + \lambda_\theta \sum_k \|\theta_k\|_1},
#'     with the penalties stored in \code{model}. Add it back to get an
#'     unpenalized value.}
#'   \item{\code{sum_w}}{Total weight.}
#'   \item{\code{sum_trimmed_w}}{Weight of the bins kept after trimming.}
#' }
#' If every bin has \eqn{\ell_{ti} = -\infty}, \code{prop_inf} is 1, the four
#' log-likelihoods are \code{-Inf}, and \code{penalty}, \code{sum_w} and
#' \code{sum_trimmed_w} are \code{NA}.
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
#' res$trimmed_loglik + res$penalty   # unpenalized trimmed log-likelihood
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
