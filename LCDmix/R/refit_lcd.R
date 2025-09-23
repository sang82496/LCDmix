# Generated from create-LCDmix.Rmd: do not edit by hand

#' Refit LCDmix multiple times and select the best training fit
#'
#' Given optimal penalties from CV, repeatedly refits \code{main()} on the
#' full (binned) dataset using different random seeds, records the final
#' training objective \code{L$loglik} for each repeat, caches results to disk,
#' and returns logs, per-repeat scores, and the best fit.
#'
#' @param Y_bin List of length \eqn{TT}; binned responses for each time point.
#'   Element \eqn{t} is a numeric vector of length \eqn{n_t}.
#' @param X Numeric matrix \eqn{TT \times p}; covariates aligned by time.
#' @param bin_mass List of length \eqn{TT}; per-time weights for each bin/observation.
#' @param K Integer; number of mixture components.
#' @param opt_lambdas Numeric length-2 vector \code{c(lambda_alpha, lambda_theta)}
#'   chosen from cross-validation.
#' @param max_iter Integer; maximum EM iterations per fit. Default \code{30}.
#' @param iter_eta Numeric; convergence threshold on relative change in the
#'   surrogate objective. Default \code{1e-3}.
#' @param resp_threshold Numeric in \eqn{[0,1]}; responsibilities below this are
#'   zeroed for stability. Default \code{1e-3}.
#' @param seeds Integer vector of random seeds (one per repeat). If \code{NULL},
#'   provide \code{cv_reps} to derive \code{seeds <- 1:cv_reps}. Default \code{NULL}.
#' @param trim_prob Numeric in \eqn{[0,1)}; trimming fraction passed to
#'   \code{main()} for any internal diagnostics (not used in the training
#'   objective). Default \code{0.01}.
#' @param save_dir Character; directory to cache per-repeat results as
#'   \code{"refit_<seed>.rds"}. Default \code{"./result"}.
#' @param n_cores Integer or \code{"max"}; number of parallel workers for
#'   \pkg{parallel}. \code{"max"} uses all physical cores minus one. Default \code{"max"}.
#' @param debug Logical; forwarded to \code{main()} to print extra diagnostics.
#'   Default \code{FALSE}.
#' @param cv_reps Integer; number of repeats used only when \code{seeds} is
#'   \code{NULL}. Default \code{NULL}.
#'
#' @details
#' Each repeat runs a full fit with \code{main()} and extracts the final training
#' objective \code{L$loglik}. Results are cached to enable resuming; cached runs
#' are skipped on subsequent calls. Parallelization uses a socket cluster and is
#' cleaned up on exit.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{logs}}{Character vector of per-repeat logs plus a summary line.}
#'   \item{\code{refit_scores}}{Numeric vector of final training log-likelihoods
#'     (\code{L$loglik}) across repeats (one per seed).}
#'   \item{\code{best_fit}}{The list for the best repeat (by \code{refit_scores}),
#'     containing \code{log_msg}, \code{L}, and \code{fit_try} (the fitted object).}
#' }
#'
#' @seealso \code{\link{cv_lcd}}, \code{\link{evaluate_lcd_model}}
#'
#' @examples
#' \dontrun{
#' best <- refit_lcd(
#'   Y_bin        = Y_bin,
#'   X            = X,
#'   bin_mass     = bin_mass,
#'   K            = 2,
#'   opt_lambdas  = c(1e-3, 1e-3),
#'   seeds        = 1:5,
#'   save_dir     = "refits"
#' )
#' best$refit_scores
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
  trim_prob = 0.01,
  save_dir = "./refits", 
  n_cores = "max", 
  debug = FALSE
) {
  if (is.null(seeds) && is.null(cv_reps)) stop("`seeds` or `cv_reps` required")
  if (is.null(seeds)) seeds <- seq_len(cv_reps)
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  lambda_alpha <- opt_lambdas[1]
  lambda_theta <- opt_lambdas[2]

  n_workers <- if (identical(n_cores, "max")) parallel::detectCores(logical = FALSE) else as.integer(n_cores)
  cl <- parallel::makeCluster(n_workers)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

  parallel::clusterEvalQ(cl, {library(LCDmix); NULL })
  parallel::clusterExport(
    cl,
    varlist = c("Y_bin","X","bin_mass","K","lambda_alpha","lambda_theta",
                "max_iter","iter_eta","resp_threshold","trim_prob","save_dir","debug"),
    envir = environment()
  )

  res_list <- parallel::parLapply(cl, seeds, function(sd) {
    refit_onejob(
      Y_bin = Y_bin, X = X, bin_mass = bin_mass, K = K,
      lambda_alpha = lambda_alpha, lambda_theta = lambda_theta,
      seed = sd, max_iter = max_iter, iter_eta = iter_eta,
      resp_threshold = resp_threshold, trim_prob = trim_prob,
      save_dir = save_dir, debug = debug
    )
  })

  logs <- vapply(res_list, `[[`, character(1), "log_msg")
  scores <- vapply(res_list, function(x) x$L, numeric(1))

  if (all(is.na(scores))) {
    logs <- c(logs, "All refits failed.")
    return(list(logs = logs, refit_scores = scores, best_fit = NULL))
  }

  best_idx <- which.max(scores)
  logs <- c(logs, paste0(
    "Completed ", length(seeds), " repeats; ",
    sum(is.na(scores)), " failures; best L = ", round(scores[best_idx], 6),
    " (seed index ", best_idx, ").\n"
  ))

  return(list(logs = logs, refit_scores = scores, best_fit = res_list[[best_idx]]))
}
