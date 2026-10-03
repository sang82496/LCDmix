# Generated from create-LCDmix.Rmd: do not edit by hand

#' Fit a log-concave mixture-of-experts model
#'
#' @description
#' Runs the whole fitting pipeline for one dataset:
#'
#' 1. If \code{binned = FALSE}, bins the responses with \code{binning()}.
#' 2. Computes starting values with \code{initialization()}, which fits a
#'    Gaussian mixture of regressions with \code{flowmix::flowmix()}.
#' 3. Runs the EM-type iterations with \code{iteration()}.
#' 4. Evaluates the fitted model on the training data with \code{eval_lcd()}:
#'    penalized log-likelihood, each bin scored by its probability.
#'
#' The flowmix start in step 2 is random. Call \code{set.seed()} first for a
#' reproducible fit.
#'
#' @param Y A list of length \code{TT}; \code{Y[[t]]} is a numeric vector or a
#'   one-column matrix of responses at time \eqn{t}. If \code{binned = TRUE},
#'   the bin centers, as returned in \code{Y_bin} by \code{binning()}.
#' @param X A numeric \eqn{TT \times p} matrix; row \code{X[t, ]} is the
#'   covariate vector at time \eqn{t}. It is used as given; see the Details
#'   of \code{main()} on covariate scale.
#' @param biomass A list of length \code{TT}; \code{biomass[[t]]} holds one
#'   nonnegative weight per element of \code{Y[[t]]} (per observation, or per
#'   bin if \code{binned = TRUE}).
#' @param binned Logical; \code{TRUE} if \code{Y} and \code{biomass} are
#'   already binned. Then they are used as given and \code{n_bins} is ignored.
#' @param n_bins Number of equal-width bins over the pooled range of \code{Y},
#'   used when \code{binned = FALSE}. \code{0} makes each distinct response
#'   value its own bin.
#' @param K Number of mixture components (experts).
#' @param lambda_alpha Nonnegative L1 penalty on the gate coefficients
#'   \eqn{\alpha}; the intercepts are not penalized.
#' @param lambda_theta Nonnegative L1 penalty on the expert slopes
#'   \eqn{\theta}; the intercepts are not penalized.
#' @param max_iter Maximum number of EM iterations.
#' @param iter_eta Stopping tolerance on the relative ascent of the surrogate
#'   objective in one iteration,
#'   \eqn{(Q(\Theta^{(m+1)} \mid \Theta^{(m)}) - Q(\Theta^{(m)} \mid \Theta^{(m)})) / |Q(\Theta^{(m)} \mid \Theta^{(m)})|}.
#'   The stopping rule is described in \code{iteration()}.
#' @param resp_threshold Responsibilities below this value are set to 0 in the
#'   E-step, and those bins are left out of the M-step for that component.
#' @param trim_prob Fraction of the total weight that \code{eval_lcd()} trims
#'   for the trimmed log-likelihood in \code{L}. It does not change the fit.
#' @param calc_Q_every Logical; if \code{TRUE}, \code{iteration()} records
#'   \eqn{Q} at five checkpoints per iteration (\code{Q_every}), the number of
#'   bins outside the support (\code{n_outside_every}) and the feasibility
#'   check of the \eqn{\theta} update (\code{lp_check_every}). This adds four
#'   \code{comp_Q()} calls per iteration; leave it \code{FALSE} for production
#'   runs.
#' @param debug Logical; if \code{TRUE}, an error inside the EM loop does not
#'   stop \code{main()}: the partial result is returned with the error message
#'   (see Value).
#' @param lp_time_limit Time limit in seconds for each linear program
#'   (Rsymphony) in the \eqn{\theta} update.
#' @param update \code{"lp"} (default) solves the \eqn{\theta} update as a
#'   linear program with \code{mstep_theta_lp()}. \code{"optim"} uses the
#'   unconstrained quasi-Newton update \code{mstep_theta_optim()}, which exists
#'   as the comparison arm of the LP ablation.
#' @param maxdev \code{NULL} (default, no constraint) or a positive number.
#'   Bounds the deviation of each component mean from its intercept,
#'   \eqn{|X_t^\top \theta_k| \le} \code{maxdev} at every time point \eqn{t},
#'   as in flowmix. It is passed to the flowmix initialization and to the LP
#'   update of \eqn{\theta}. Requires \code{update = "lp"}.
#'
#' @details
#' Covariate scale. The package does not standardize \code{X}. The gate update
#' calls \code{glmnet} with its default \code{standardize = TRUE}, so in the
#' original scale the gate penalty is
#' \eqn{\lambda_\alpha \sum_k \sum_j \mathrm{sd}(x_j) |\alpha_{kj}|}, while
#' \code{comp_Q()} and \code{eval_lcd()} use
#' \eqn{\lambda_\alpha \sum_k \sum_j |\alpha_{kj}|}. The two agree when every
#' column of \code{X} has standard deviation 1. flowmix behaves the same way.
#' Center and scale continuous covariates before fitting (0/1 indicators can
#' stay as they are), and give the same \code{X} to both methods.
#'
#' @return A list with components:
#' \describe{
#'   \item{Y_bin}{List of binned responses (\code{Y} itself if \code{binned = TRUE}).}
#'   \item{X}{Covariate matrix (unchanged).}
#'   \item{bin_mass}{List of the total weight in each bin.}
#'   \item{K}{Number of components.}
#'   \item{initial}{The list returned by \code{initialization()}.}
#'   \item{iter}{The list returned by \code{iteration()}: the fitted parameters
#'     (\code{alpha_new}, \code{theta0_new}, \code{theta_new}, \code{g_new}),
#'     the trace of \eqn{Q} and diagnostics.}
#'   \item{L}{The list returned by \code{eval_lcd()} on the training data.
#'     \code{L$loglik} is the penalized mean log-likelihood that
#'     \code{refit_lcd()} uses to choose among restarts.}
#' }
#' If \code{debug = TRUE} and the EM loop fails, a warning is given and the
#' list holds \code{Y_bin}, \code{X}, \code{bin_mass}, \code{initial},
#' \code{iter_partial} (the list returned by \code{iteration()}, with
#' \code{error} and \code{failed_iter} set) and \code{iter = NULL}.
#'
#' @examples
#' \dontrun{
#' # TT = 50 time points, p = 3 covariates
#' set.seed(123)
#' Y_list  <- lapply(1:50, function(t) matrix(rnorm(sample(20:50, 1)), ncol = 1))
#' biomass <- lapply(Y_list, function(y) runif(nrow(y), 0.5, 2))
#' X_mat   <- matrix(rnorm(50 * 3), nrow = 50, ncol = 3)
#' # Fit a 2-component mixture
#' result  <- main(
#'   Y       = Y_list,
#'   X       = X_mat,
#'   biomass = biomass,
#'   K       = 2
#' )
#' plot(result$iter$Q, type = "b")   # surrogate objective per iteration
#' result$L$loglik                   # penalized training log-likelihood
#' }
#' @export
main <- function(
  Y,
  X,
  biomass,
  binned         = FALSE,
  n_bins         = 40,
  K              = 2,
  lambda_alpha   = 1e-3,
  lambda_theta   = 1e-3,
  max_iter       = 30,
  iter_eta       = 1e-4,
  resp_threshold = 1e-3,
  trim_prob      = 0.03,
  calc_Q_every   = FALSE,
  debug          = FALSE,
  lp_time_limit  = 3600,
  update         = c("lp", "optim"),    # NEW
  maxdev         = NULL                 # NEW (fixP) - appended last
) {
  update <- match.arg(update)           # NEW
  #— Step 1: Binning (if needed) —#
  if (binned) {
    Y_bin    <- Y
    bin_mass <- biomass
  } else {
    bin_res  <- binning(Y, biomass, n_bins)
    Y_bin    <- bin_res$Y_bin
    bin_mass <- bin_res$bin_mass
  }
  message("✔ Binning complete")
  
  #— Step 2: Initialization via GMR —#
  init_res <- initialization(
    Y_bin,
    X,
    bin_mass,
    K,
    lambda_alpha,
    lambda_theta,
    resp_threshold,
    maxdev = maxdev                     # NEW (fixP)
  )
  message("✔ Initialization complete")
  
  #— Step 3: EM‐style iterations —#
  iter_res <- iteration(
    Y_bin,
    X,
    bin_mass,
    init_res,
    lambda_alpha,
    lambda_theta,
    iter_eta,
    max_iter,
    resp_threshold,
    calc_Q_every,
    debug,
    lp_time_limit,
    update = update,                    # NEW
    maxdev = maxdev                     # NEW (fixP)
  )
  
  if (debug && !is.null(iter_res$error)) {
    # iteration() caught the error: iter_res holds the message (error) and the iteration (failed_iter)
    warning("EM failed at iteration ", iter_res$failed_iter, 
            ": ", iter_res$error)
    return(list(
            Y_bin        = Y_bin,
            X            = X,
            bin_mass     = bin_mass,
            initial      = init_res,
            iter_partial = iter_res,
            iter         = NULL
          ))
  }
  message("✔ Iterations complete")
  
  #— Step 4: Compute loglikelihood —#
  L = eval_lcd(
    model         = iter_res,
    Y_test        = Y_bin,
    X_test        = X,
    biomass_test  = bin_mass,
    trim_prob     = trim_prob
  )
  message(paste0("✔ Calculating loglikelihood complete: L = ", round(L$loglik, 6)) )
  
  #— Return all key results —#
  return(list(
    Y_bin    = Y_bin,
    X        = X,
    bin_mass = bin_mass,
    K        = K,
    initial  = init_res,
    iter     = iter_res,
    L        = L
  ))
}
