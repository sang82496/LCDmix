# Generated from create-LCDmix.Rmd: do not edit by hand

#' Simulate one dataset from the two-component simulation design
#'
#' @description
#' Draws binned data from a mixture of two regressions over \code{TT} time
#' points:
#'
#' - Covariates: \code{X = (par, ramp, noise1, ..., noise(p-2))}. \code{par}
#'   is the first \code{TT} values of the series in \code{simul_helper.rds},
#'   scaled to mean 0 and standard deviation 1. \code{ramp} rises linearly
#'   from 0 to 1 over the first half of the time points and stays at 1 over
#'   the second half; it is not standardized. The noise covariates are
#'   independent \eqn{N(0, 1)}.
#' - Expert means: component 1 has intercept 0 and slope \code{theta_par} on
#'   \code{par}; component 2 has intercept \code{gap} and slope
#'   \code{-theta_par} on \code{par}. All other slopes are 0.
#' - Gate: \eqn{P(Z_t = 2) = \mathrm{logit}^{-1}(\mathrm{logit}(0.10) + (\mathrm{logit}(0.30) - \mathrm{logit}(0.10))\, ramp_t)},
#'   which rises from 0.10 at \code{ramp = 0} to 0.30 at \code{ramp = 1}. All
#'   other gate coefficients are 0.
#' - Errors: the same mean-zero distribution in both components, chosen by
#'   \code{noisetype}: \code{"gaussian"}, \eqn{N(0, 1)}; \code{"skewed"},
#'   skew-normal with shape \code{skew_alpha}, scaled to variance 1 and
#'   shifted to mean 0; \code{"heavytail"}, \eqn{t} with \code{df} degrees of
#'   freedom, scaled to variance 1; \code{"laplace"}, Laplace with scale 1
#'   (variance 2); \code{"exponential"}, \eqn{Exp(1) - 1}.
#' - Each time point has \code{nt} observations with weight 1, binned with
#'   \code{binning()} into \code{B} equal-width bins over the pooled range.
#'
#' @param sim_seed Seed for \code{set.seed()}; \code{NULL} leaves the random
#'   number generator as it is.
#' @param nt Number of observations per time point (a multiple of 5).
#' @param TT Number of time points. Use an even number, so that the two halves
#'   of \code{ramp} have equal length.
#' @param theta_par Slope on \code{par}: \code{+theta_par} for component 1 and
#'   \code{-theta_par} for component 2.
#' @param p Number of covariates, at least 3.
#' @param B Number of bins.
#' @param noisetype One of \code{"gaussian"}, \code{"skewed"},
#'   \code{"heavytail"}, \code{"laplace"} and \code{"exponential"}.
#' @param df Degrees of freedom for \code{"heavytail"} (at least 3).
#' @param skew_alpha Shape parameter for \code{"skewed"}.
#' @param gap Intercept of component 2 (component 1 has intercept 0).
#' @param sim_helper_dir Directory that holds \code{simul_helper.rds}.
#'
#' @return A list with components:
#' \describe{
#'   \item{Y_bin, X, bin_mass}{The data, in the form \code{main()} takes with
#'     \code{binned = TRUE}.}
#'   \item{noisetype, df, skew_alpha}{The arguments, returned unchanged.}
#'   \item{mnmat}{\eqn{TT \times 2} matrix of the true component means.}
#'   \item{prob}{\eqn{TT \times 2} matrix of the true gate probabilities.}
#'   \item{alpha, theta}{True gate and expert coefficients, each a
#'     \eqn{(p+1) \times 2} matrix with rows intercept, \code{par},
#'     \code{ramp}, \code{noise1}, ... . \code{alpha} is the transpose of the
#'     \eqn{K \times (p+1)} layout used by the rest of the package, so the
#'     gate probabilities are \code{pi_k(X, t(alpha))}.}
#'   \item{omega, mn_shift}{Scale and mean shift of the skew-normal
#'     (\code{"skewed"} only; \code{NULL} otherwise).}
#'   \item{variance}{Variance of the \eqn{t} distribution before scaling
#'     (\code{"heavytail"} only; \code{NULL} otherwise).}
#' }
#'
#' @export
gen_simul_data <- function(
  sim_seed       = NULL,
  nt             = 1000,
  TT             = 100,
  theta_par      = 0.5,
  p              = 10,
  B              = 30,
  noisetype      = 'gaussian',
  df             = NULL,
  skew_alpha     = NULL,
  gap            = 4,
  sim_helper_dir = "."
) {
  ## Setup and basic checks
  assertthat::assert_that(nt %% 5 == 0)
  assertthat::assert_that(noisetype %in% c('heavytail', 'skewed', 'laplace', 'exponential', 'gaussian'))
  K = 2
  stopifnot(p >= 3)
  if (noisetype == 'heavytail') {
    assertthat::assert_that(!is.null(df))
    assertthat::assert_that(df >= 3)
  } else if (noisetype == 'skewed') {
    assertthat::assert_that(!is.null(skew_alpha))
  }
  
  if(!is.null(sim_seed)) set.seed(sim_seed)
  ntlist = rep(nt, TT)

  ## Generate covariate
  par = readRDS(file.path(sim_helper_dir, "simul_helper.rds"))

  Xrest = do.call(cbind, lapply(1:(p-2), function(ii) rnorm(TT)) )
  X2 = c(seq(0, 1, length = TT/2), rep(1, TT/2))   # fixV: a ramp from 0 to 1, then flat (was 0.8 to 1)
  X = cbind(scale(par[1:TT]), X2, Xrest)
  colnames(X) = c("par", "ramp", paste0("noise", 1:(p-2)))   # fixV: "cp" renamed to "ramp"

  ## theta coefficients
  theta = matrix(0, ncol = K, nrow = p+1)
  theta[0+1,1] = 0
  theta[1+1,1] = theta_par
  theta[0+1,2] = gap
  theta[1+1,2] = -theta_par
  colnames(theta) = paste0("clust", 1:K)
  rownames(theta) = c("intercept", "par", "ramp", paste0("noise", 1:(p-2)))

  ## alpha coefficients
  alpha = matrix(0, ncol = K, nrow = p+1)
  ## fixV: P(component 2) is 10% at ramp = 0 and 30% at ramp = 1
  ## (was -10 and 10 + log(1/4): 4.3% at X2 = 0.8 and 20% at X2 = 1)
  alpha[0+1, 2] = stats::qlogis(0.10)
  alpha[2+1, 2] = stats::qlogis(0.30) - stats::qlogis(0.10)

  colnames(alpha) = paste0("clust", 1:K)
  rownames(alpha) = c("intercept", "par", "ramp", paste0("noise", 1:(p-2)))

  ## Generate means and probabilities
  mnmat = cbind(1, X) %*% theta
  prob = pi_k(X, t(alpha))
  
  ## Samples |nt| memberships out of (1:K) according to the probs in prob.
  ## Data is a probabilistic mixture from these two means, over time.
  omega = NULL
  mn_shift = NULL
  variance = NULL
  if (noisetype == 'heavytail') {
     variance = df / (df - 2)
  } else if (noisetype == 'skewed') {
      omega = sqrt(1/(1 - 2 * (1/pi) * skew_alpha^2 / (1 + skew_alpha^2)))
      mn_shift = omega * skew_alpha * (1 / sqrt(1+skew_alpha^2)) * sqrt(2/pi)
  }
  ylist = lapply(1:TT, function(tt){
     draws = sample(1:K, size = ntlist[tt], replace = TRUE,
                    prob = c(prob[tt,1], prob[tt,2]))
     mns = mnmat[tt,]
     means = mns[draws]
     ## Add noise to obtain data points.
     if (noisetype == 'heavytail') {
       noise = rt(ntlist[tt], df = df) / sqrt(variance)
     } else if (noisetype == 'skewed'){
       noise = sn::rsn(ntlist[tt], xi = 0, omega = omega, alpha = skew_alpha) - mn_shift
     } else if (noisetype == 'laplace'){
       noise = VGAM::rlaplace(ntlist[tt], 0, 1)
     } else if (noisetype == 'exponential'){
       noise = rexp(ntlist[tt], 1) - 1
     } else { # gaussian
       noise = rnorm(ntlist[tt], 0, 1)
     } 
     datapoints = means + noise
     cbind(datapoints)
   })

  # Binning
  biomass = vector("list", TT)
  biomass = lapply(1:TT, function(t){biomass[[t]] = rep(1, length(ylist[[t]]))})
  binned = LCDmix::binning(ylist, biomass, n_bins = B)

  return(list(Y_bin = binned$Y_bin, 
              X = X,
              bin_mass = binned$bin_mass,
              ## The true generating model:
              noisetype = noisetype,
              mnmat = mnmat,
              prob = prob,
              alpha = alpha,
              theta = theta,
              skew_alpha = skew_alpha,
              omega = omega,
              mn_shift = mn_shift,
              df = df,
              variance = variance
              ))
}
