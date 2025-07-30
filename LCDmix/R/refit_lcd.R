# Generated from create-LCDmix.Rmd: do not edit by hand

#' Refit the LCD mixture model at optimal penalties with multiple repeats
#'
#' @description
#' Runs the full \code{main()} pipeline on the entire binned dataset multiple times
#' at the optimal L1 penalties to assess variability and select the best repeat.
#' Computations are dispatched in parallel.
#'
#' @param Y_bin A list of length \eqn{TT}, where each element is an \eqn{M_t \times 1}
#'   matrix of binned responses per time point.
#' @param X Numeric \eqn{TT \times p} covariate matrix (rows = time points).
#' @param bin_mass A list of length \eqn{TT}, where each element is a numeric vector
#'   of biomass weights per bin at time \(t\).
#' @param K Integer number of mixture components.
#' @param opt_lambdas Numeric vector of length 2: optimal \(\lambda_{\alpha}\) and
#'   \(\lambda_{\theta}\).
#' @param max_iter Integer maximum number of EM iterations. Default: 30.
#' @param iter_eta Numeric convergence threshold (relative change in Q). Default: 1e-3.
#' @param resp_threshold Numeric in [0,1]; responsibilities below this are zeroed. Default: 1e-3.
#' @param seeds Integer vector of random seeds (one per repeat).  Must supply either
#'   \code{seeds} or \code{cv_reps}, not both. Default: \code{NULL}.
#' @param trim_prob Numeric in [0,0.5]; fraction of lowest‐likelihood points to trim in evaluation. Default: 0.05.
#' @param save_dir Character; directory in which to save per‐repeat results. Default: "./result".
#' @param n_cores Integer or "max"; number of parallel workers. "max" uses all physical cores minus one. Default: "max".
#' @param cv_reps Integer number of repeated fits. Default: \code{NULL}.
#' @param maxdev Numeric or \code{NULL}; optional max‐deviation constraint. Default: \code{NULL}.
#' @param n_restarts Integer number of random restarts for \code{flowmix::flowmix()}. Default: 1.
#'
#' @return A list with components:
#' \describe{
#'   \item{logs}{Character vector of log messages (one per repeat, plus summary).}
#'   \item{refit_scores}{Numeric vector of final surrogate log‐likelihoods (\code{Q}) for each repeat.}
#'   \item{best_fit}{The \code{iteration} result list from the repeat with the highest final \code{Q}.}
#' }
#'
#' @export
refit_lcd <- function(
  Y_bin,
  X,
  bin_mass,
  K,
  opt_lambdas,
  max_iter       = 30,
  iter_eta       = 1e-3,
  resp_threshold = 1e-3,
  seeds          = NULL,
  trim_prob      = 0.05,
  save_dir       = "./result",
  n_cores        = "max",
  debug          = FALSE,
  cv_reps        = NULL,
  maxdev         = NULL,
  n_restarts     = 1
) {
  # Unpack optimal penalties
  lambda_alpha <- opt_lambdas[1]
  lambda_theta <- opt_lambdas[2]
  
  # Determine number of workers
  workers <- if (identical(n_cores, "max")) {
    parallel::detectCores(logical = FALSE) - 1
  } else {
    n_cores
  }
  cl <- parallel::makeCluster(workers)
  
  # Export necessary objects and functions
  parallel::clusterExport(
    cl,
    varlist = c(
      "main", "Y_bin", "X", "bin_mass", "K",
      "lambda_alpha", "lambda_theta",
      "n_restarts", "max_iter", "iter_eta",
      "maxdev", "resp_threshold", "save_dir"
    ),
    envir = environment()
  )
  
  # Perform repeated fits in parallel
  results <- parallel::parLapply(
    cl,
    seeds,
    function(seed_idx) {
      # start log for this seed
      log_msg <- paste0("▶ Refit seed ", seed_idx, "\n")
      set.seed(seed_idx)
      
      # capture ALL output from main(), including messages
      out_log <- capture.output({
        fit_try <- tryCatch(
          main(
            Y               = Y_bin,
            X               = X,
            biomass         = bin_mass,
            binned          = TRUE,
            n_bins          = 0,
            K               = K,
            lambda_alpha    = lambda_alpha,
            lambda_theta    = lambda_theta,
            max_iter        = max_iter,
            iter_eta        = iter_eta,
            resp_threshold  = resp_threshold,
            maxdev          = maxdev,
            n_restarts      = n_restarts
          ),
          error = function(e) {
            message("✖ Error on seed ", seed_idx, ": ", e$message)
            return(NA)
          }
        )
      }, type = "message")
      
      # append captured output to our log
      log_msg <- paste0(log_msg, paste(out_log, collapse = "\n"), "\n")
      
      # check success / failure
      if ((!is.list(fit_try))) {
        # failed
        return(list(
          log = paste0(log_msg, "✖ Fit failed on seed ", seed_idx, "\n"),
          Q   = NA,
          iter= NULL
        ))
      } else if (!is.null(fit_try$iter_partial)) {
        # Error at EM iteration
        return(list(
          log = paste0("Error at EM iteration ", fit_try$iter_partial$failed_iter, 
                          ": ", fit_try$error, "\n"),
          Q   = NA,
          iter= NULL
        ))
      }
    
    # if successful, record final Q and save
      final_Q  <- fit_try$iter$Q[length(fit_try$iter$Q)]
      save_path <- file.path(save_dir, paste0("refit_", seed_idx, ".Rdata"))
      save(final_Q, fit_try, file = save_path)
      log_msg <- paste0(log_msg,
                        "✔ Completed seed ", seed_idx,
                        "; final Q = ", round(final_Q, 4), "\n")
      return(list(
        log  = log_msg,
        Q    = final_Q,
        iter = fit_try
      ))
    }
  )
  parallel::stopCluster(cl)
  
  # Aggregate
  log_msgs     <- sapply(results, `[[`, "log")
  refit_scores <- sapply(results, `[[`, "Q")
  best_idx     <- which.max(refit_scores)
  best_fit     <- results[[best_idx]]
  
  # Summary
  summary_msg <- paste0(
    "Completed ", cv_reps, " repeats; ",
    sum(is.na(refit_scores)), " failures; ",
    "best Q = ", round(refit_scores[best_idx],4),
    " on repeat ", best_idx, "\n"
  )
  log_msgs <- c(log_msgs, summary_msg)
  
  return(list(
    logs         = log_msgs,
    refit_scores = refit_scores,
    best_fit     = best_fit
  ))
}
