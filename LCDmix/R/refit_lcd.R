# Generated from create-LCDmix.Rmd: do not edit by hand

#' Refit the model with several seeds and keep the best fit
#'
#' @description
#' Runs \code{refit_onejob()} for every seed in parallel (one cluster, no
#' nested parallelism) on the full data with the chosen penalty pair. Each
#' refit is saved to \code{save_dir/refit_<seed>.rds}; a seed whose file
#' already exists is not run again, so an interrupted run can be resumed.
#' After the run, the saved files are read and the fit with the largest
#' \code{fit_L} (the penalized training log-likelihood) is returned.
#'
#' The workers load the installed LCDmix package, so install the current
#' version before running.
#'
#' @inheritParams main
#' @inheritParams initialization
#' @param opt_lambdas \code{c(lambda_alpha, lambda_theta)}, for example the
#'   \code{opt_lambdas} of \code{cv_lcd_summary()}.
#' @param seeds Seeds to run. If \code{NULL}, \code{1:cv_reps}.
#' @param cv_reps Number of seeds, used only when \code{seeds} is \code{NULL}.
#' @param save_dir Directory for the refit files; it is created if needed.
#' @param n_cores Number of worker processes, or \code{"max"} for all physical
#'   cores.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{summary}}{A text with the number of failed refits.}
#'   \item{\code{loglik}}{Numeric vector, \code{fit_L} of each seed
#'     (\code{NA} for a failed seed).}
#'   \item{\code{best_fit}}{The \code{main()} result of the seed with the
#'     largest \code{fit_L}; \code{NULL}, with a warning, if every refit
#'     failed.}
#' }
#'
#' @seealso \code{\link{refit_onejob}}, \code{\link{cv_lcd}}
#'
#' @examples
#' \dontrun{
#' out <- refit_lcd(
#'   Y_bin = Y_bin, X = X, bin_mass = bin_mass, K = 2,
#'   opt_lambdas = c(1e-3, 1e-3),
#'   seeds = 1:10,
#'   save_dir = "refits",
#'   n_cores = 8
#' )
#' out$loglik                  # penalized training log-likelihood per seed
#' out$best_fit$iter$theta_new # slopes of the best fit
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
  iter_eta = 1e-4, 
  resp_threshold = 1e-3, 
  trim_prob = 0.03,
  save_dir = "./refits", 
  n_cores = "max",
  lp_time_limit = 3600,
  update = c("lp", "optim"),    # NEW
  maxdev = NULL                 # NEW (fixP) - appended last
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
                "max_iter","iter_eta","resp_threshold","trim_prob","save_dir",
                "lp_time_limit", "update","maxdev"), # NEW (fixP)
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
      save_dir = save_dir, lp_time_limit = lp_time_limit, update = update,
      maxdev = maxdev)                                                     # NEW (fixP)
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
    if (!is.null(obj$fit_L)) {
      L_vec[i]  <- obj$fit_L
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
              loglik   = L_vec,
              best_fit = best_fit))
}
