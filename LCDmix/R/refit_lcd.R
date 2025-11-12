# Generated from create-LCDmix.Rmd: do not edit by hand

#' Refit LCDmix multiple times and select the best training fit
#'
#' Runs multiple refits of an LCDmix model on the full training data using a
#' grid of random seeds (one outer parallel cluster; no nested parallel). Each
#' refit is cached to \code{save_dir/refit_<seed>.rds} so the procedure is
#' fully resumable. The best fit is chosen by the largest training objective
#' \eqn{L} among the completed refits.
#'
#' @param Y_bin List of responses (one element per time point/bin), as used by \code{main()}.
#' @param X Numeric matrix of covariates with \code{nrow(X) == length(Y_bin)}.
#' @param bin_mass List of nonnegative weights aligned with \code{Y_bin}.
#' @param K Integer; number of mixture components.
#' @param opt_lambdas Numeric vector of length 2 giving \code{c(lambda_alpha, lambda_theta)}.
#' @param seeds Integer vector of seeds to run. If \code{NULL}, supply \code{cv_reps}.
#' @param cv_reps Integer; number of repeats used only when \code{seeds} is \code{NULL}.
#'   The seeds will be \code{1:cv_reps}.
#' @param max_iter Integer; maximum EM iterations. Default \code{30}.
#' @param iter_eta Numeric; convergence tolerance for the surrogate objective. Default \code{1e-3}.
#' @param resp_threshold Numeric in \eqn{[0,1]}; responsibilities below this are set to zero. Default \code{1e-3}.
#' @param trim_prob Numeric in \eqn{[0,1)}; trimming fraction used by \code{eval_lcd()} during fitting. Default \code{0.03}.
#' @param save_dir Character; directory to write/read cached refits (\code{refit_<seed>.rds}). Default \code{"./refits"}.
#' @param n_cores Integer or \code{"max"}; number of workers for the single outer cluster. Default \code{"max"}.
#' @param debug Logical; forwarded to \code{main()} for verbose diagnostics. Default \code{FALSE}.
#'
#' @details
#' Each seed triggers a call to \code{refit_onejob()}, which caches its result to
#' disk and returns a log string and the final training objective \eqn{L}.
#' Existing cache files are reused and not recomputed.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{logs}}{Character vector of log messages (one per seed plus a summary line).}
#'   \item{\code{refit_scores}}{Numeric vector of per-seed training objectives \eqn{L} (may contain \code{NA}).}
#'   \item{\code{best_fit}}{The list returned by \code{refit_onejob()} for the best seed
#'     (largest \eqn{L}); \code{NULL} if all refits failed.}
#' }
#'
#' @seealso \code{\link{refit_onejob}}
#'
#' @examples
#' \dontrun{
#' # Given Y_bin, X, bin_mass, and chosen penalties:
#' opt <- c(lambda_alpha = 1e-3, lambda_theta = 1e-3)
#' out <- refit_lcd(
#'   Y_bin = Y_bin, X = X, bin_mass = bin_mass, K = 2,
#'   opt_lambdas = opt,
#'   seeds = 1:10,
#'   save_dir = "refits",
#'   n_cores = 8
#' )
#' out$best_fit$L       # best training objective
#' out$best_fit$file    # path to the cached best refit
#' }
#'
#' @export
refit_lcd <- function(
  Y_bin, 
  X, 
  bin_mass, 
  K, 
  opt_lambdas,
  seeds = NULL, 
  cv_reps = NULL,
  max_iter = 30, 
  iter_eta = 1e-3, 
  resp_threshold = 1e-3, 
  trim_prob = 0.03,
  save_dir = "./refits", 
  n_cores = "max",
  trimmed = T
) {
  if (is.null(seeds) && is.null(cv_reps)) stop("`seeds` or `cv_reps` required")
  if (is.null(seeds)) seeds <- seq_len(cv_reps)
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  seeds        = as.integer(seeds)
  lambda_alpha = as.numeric(opt_lambdas[1])
  lambda_theta = as.numeric(opt_lambdas[2])

  n_workers <- if (identical(n_cores, "max")) parallel::detectCores(logical = FALSE) else as.integer(n_cores)
  cl <- parallel::makeCluster(n_workers)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

  parallel::clusterEvalQ(cl, {library(LCDmix); NULL })
  parallel::clusterExport(
    cl,
    varlist = c("Y_bin","X","bin_mass","K","lambda_alpha","lambda_theta",
                "max_iter","iter_eta","resp_threshold","trim_prob","save_dir", "trimmed"),
    envir = environment()
  )

  # One refit per task; TRUE on success, FALSE on fail (cached or fresh)
  res <- parallel::parLapply(cl, seeds, function(ii) {
    out_path <- file.path(save_dir, sprintf("refit_%d.rds", ii))

    # Cached result?
    if (file.exists(out_path)) {
      obj <- readRDS(out_path)
      return(!is.null(obj$fit))
    }

    # Run one refit (writes cache)
    res_ii <- refit_onejob(
      Y_bin = Y_bin, X = X, bin_mass = bin_mass, K = K,
      lambda_alpha = lambda_alpha, lambda_theta = lambda_theta,
      seed = ii, max_iter = max_iter, iter_eta = iter_eta,
      resp_threshold = resp_threshold, trim_prob = trim_prob,
      save_dir = save_dir, trimmed = trimmed)
    return(res_ii)
  })

  success <- unlist(res, use.names = FALSE)
  summary <- sprintf("Refit failures: %d/%d (%.1f%%)", sum(!success), length(success), 100 * sum(!success) / length(success))
  
  L_vec <- rep(NA_real_, length(seeds))
  for (i in seq_along(seeds)) {
    sd <- as.integer(seeds[i])
    file_name <- file.path(save_dir, sprintf("refit_%d.rds", sd))
    if (!file.exists(file_name)) next
    obj <- readRDS(file_name)
    if (!is.null(obj$fit)) {
      if (trimmed) {
        L_vec[i]  <- obj$fit$L$trimmed_loglik
      } else {
        L_vec[i]  <- obj$fit$L$med_loglik
      }
      
    }
  }
  
  if (all(!is.finite(L_vec))) {
    warning("No finite loglikelihood among refits")
    best_fit <- NULL
  } else {
      file_name <- file.path(save_dir, sprintf("refit_%d.rds", seeds[which.max(L_vec)]))
      obj <- readRDS(file_name)
      best_fit = obj$fit
    }
  return(list(summary = summary,
              L_vec   = L_vec,
              best_fit = best_fit))
}
