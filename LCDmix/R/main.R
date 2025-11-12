# Generated from create-LCDmix.Rmd: do not edit by hand

#' Fit a log-concave mixture-of-experts model with optional binning
#'
#' @description
#' Runs the full pipeline for fitting a log-concave mixture-of-experts model:
#' 1. (Optional) Bin responses by biomass  
#' 2. Initialize via Gaussian mixture regression (GMR) with \code{flowmix}  
#' 3. Perform EM‐style iterations on mixture parameters  
#'
#' @param Y A list of length \code{TT}, where \code{Y[[t]]} is an \eqn{n_t}-row vector or single‐column matrix of responses at time \eqn{t}.
#' @param X A numeric matrix of dimension \eqn{TT \times p}, where each row \code{X[t, ]} is the covariate vector at time \eqn{t}.
#' @param biomass A list of length \code{TT}, where \code{biomass[[t]]} is a numeric vector of biomass weights per observation or bin at time \eqn{t}.
#' @param binned Logical; if \code{TRUE}, \code{Y} and \code{biomass} are assumed already binned. Default: \code{FALSE}.
#' @param n_bins Integer number of equal‐width bins if \code{binned = FALSE}. Default: \code{40}.
#' @param K Integer number of mixture components.
#' @param lambda_alpha Positive numeric L1 penalty on mixture‐weight coefficients. Default: \code{1e-3}.
#' @param lambda_theta Positive numeric L1 penalty on regression‐slope coefficients. Default: \code{1e-3}.
#' @param max_iter Integer maximum number of EM iterations. Default: \code{30}.
#' @param iter_eta Numeric step‐size (learning rate) for parameter updates. Default: \code{1e-3}.
#' @param resp_threshold Numeric threshold on responsibilities for soft‐assignment: any posterior probability below this value is treated as zero to improve numerical stability and computational speed. Default: \code{1e-3}.
#'
#' @return A list with components:
#' \describe{
#'   \item{Y_bin}{List of binned responses (or original \code{Y} if \code{binned = TRUE}).}
#'   \item{X}{Covariate matrix (unchanged).}
#'   \item{bin_mass}{List of biomass‐per‐bin weights.}
#'   \item{initial}{List returned by \code{initialization()}, containing starting parameters.}
#'   \item{iter}{List returned by \code{iteration()}, containing fitted parameters over EM iterations.}
#' }
#'
#' @examples
#' \dontrun{
#' # Simulate TT = 50 time points, p = 3 covariates
#' set.seed(123)
#' Y_list   <- lapply(1:50, function(t) matrix(rnorm(sample(20:50,1)), ncol = 1))
#' biomass  <- lapply(Y_list, function(y) runif(nrow(y), 0.5, 2))
#' X_mat    <- matrix(rnorm(50 * 3), nrow = 50, ncol = 3)
#' # Fit a 2‐component mixture
#' result   <- main(
#'   Y       = Y_list,
#'   X       = X_mat,
#'   biomass = biomass,
#'   K       = 2
#' )
#' plot(result$iter$logLik)
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
  iter_eta       = 1e-3,
  resp_threshold = 1e-3,
  trim_prob      = 0.03,
  calc_Q_every   = FALSE,
  debug          = FALSE
) {
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
    resp_threshold
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
    debug
  )
  
  if (debug && !is.null(iter_res$error)) {
    # here iter_res$final contains your last parameters
    warning("EM failed at iteration ", iter_res$failed_iter, 
            ": ", iter_res$error)
    return(list(
            Y_bin        = Y_bin,
            X            = X,
            bin_mass     = bin_mass,
            initial      = init_res,
            iter_partial = iter_res,
            iter         = NULL,
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
  message(paste0("✔ Calculating loglikelihood complete: trimmed L = ", round(L$trimmed_loglik, 6)) )
  
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
