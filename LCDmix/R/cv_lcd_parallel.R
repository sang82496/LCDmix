# Generated from create-LCDmix.Rmd: do not edit by hand

#' Cross‐validate the LCD mixture‐of‐experts model in parallel
#'
#' @description
#' Performs grid‐search cross‐validation over the L1 penalties \code{lambda_alpha}
#' and \code{lambda_theta} for the log‐concave mixture‐of‐experts model.  Creates
#' folds (with optional blocking), repeats each fold either by explicit seeds
#' or by a count of repeats, and evaluates models on held‐out data.  Computations
#' are dispatched in parallel, and per‐run results are saved.
#'
#' @param Y_bin List of length \eqn{TT}, where each element is an \eqn{M_t \times 1}
#'   matrix of binned responses per time point.
#' @param X Numeric \eqn{TT \times p} matrix of covariates (rows = time points).
#' @param bin_mass List of length \eqn{TT}, where each element is a numeric vector of
#'   biomass weights per bin at time \eqn{t}.
#' @param K Integer number of mixture components.
#' @param lambda_alpha_range Numeric vector of length 2 giving the min and max values
#'   for the \code{lambda_alpha} grid (on a log scale). Default: \code{c(1e-8, 1e-1)}.
#' @param lambda_theta_range Numeric vector of length 2 giving the min and max values
#'   for the \code{lambda_theta} grid (on a log scale). Default: \code{c(1e-8, 1e-1)}.
#' @param cv_gridsize Integer number of candidate lambdas per penalty. Default: 5.
#' @param max_iter Integer maximum number of EM iterations per model fit. Default: 30.
#' @param iter_eta Numeric convergence threshold (relative change in Q). Default: 1e-3.
#' @param resp_threshold Numeric in \[0,1\]; responsibilities below this are zeroed. Default: 1e-3.
#' @param nfold Integer number of cross‐validation folds. Default: 5.
#' @param seeds Integer vector of random seeds (one per repeat).  Must supply either
#'   \code{seeds} or \code{cv_reps}, not both. Default: \code{NULL}.
#' @param trim_prob Numeric fraction of lowest‐likelihood observations to trim when averaging. Default: 0.05.
#' @param save_dir Character; directory in which to save per‐run results. Default: "./result".
#' @param n_cores Integer or "max"; number of parallel worker processes. "max" uses all physical cores minus one. Default: "max".
#' @param cv_reps Integer number of repeated CV runs (used when \code{seeds} is \code{NULL}). Default: \code{NULL}.
#' @param blocksize Integer block size for CV fold creation. Default: 5.
#' @param maxdev Numeric or \code{NULL}; maximum deviation for LP constraints. Default: \code{NULL}.
#' @param n_restarts Integer number of random restarts for the initial GMR. Default: 1.
#'
#' @return A list with components:
#' \describe{
#'   \item{logs}{Character vector of log messages (one per CV run, plus summary).}
#'   \item{index_matrix}{Integer matrix with one row per run (fold × seed × penalties),
#'         showing \code{alpha_idx}, \code{theta_idx}, \code{seed_idx}, \code{fold_idx},
#'         and actual \code{lambda_alpha}, \code{lambda_theta} values.}
#' }
#'
#' @examples
#' \dontrun{
#' # Four folds, three repeats via explicit seeds
#' seeds <- c(101, 202, 303)
#' results <- cv_lcd_parallel(
#'   Y_bin               = Y_binned,
#'   X                   = X_mat,
#'   bin_mass            = counts_list,
#'   K                   = 2,
#'   cv_gridsize         = 4,
#'   nfold               = 4,
#'   seeds               = seeds,
#'   n_cores             = "max"
#' )
#'
#' # Or specify number of repeats instead of seeds
#' results2 <- cv_lcd_parallel(
#'   Y_bin               = Y_binned,
#'   X                   = X_mat,
#'   bin_mass            = counts_list,
#'   K                   = 2,
#'   cv_gridsize         = 4,
#'   nfold               = 4,
#'   n_cores             = 2,
#'   cv_reps             = 3
#' )
#' }
#' @export
cv_lcd_parallel <- function(
  Y_bin,
  X,
  bin_mass,
  K,
  lambda_alpha_range = c(1e-8, 1e-1),
  lambda_theta_range = c(1e-8, 1e-1),
  cv_gridsize        = 5,
  max_iter           = 30,
  iter_eta           = 1e-3,
  resp_threshold     = 1e-3,
  nfold              = 5,
  seeds              = NULL,
  trim_prob          = 0.05,
  save_dir           = "./result",
  n_cores            = "max",
  cv_reps            = NULL,
  blocksize          = 5,
  maxdev             = NULL,
  n_restarts         = 1
) {
  #— Create penalty grids —#
  alpha_lambdas <- flowmix::logspace(
    min    = lambda_alpha_range[1],
    max    = lambda_alpha_range[2],
    length = cv_gridsize
  ) |> sort()
  theta_lambdas <- flowmix::logspace(
    min    = lambda_theta_range[1],
    max    = lambda_theta_range[2],
    length = cv_gridsize
  ) |> sort()
  print(alpha_lambdas)
  print(theta_lambdas)
  
  #— Create CV folds —#
  folds <- flowmix::make_cv_folds(
    ylist     = Y_bin,
    nfold     = nfold,
    blocksize = blocksize
  )
  
  #— Build index matrix of all (fold × rep × penalties) combos —#
  
  if (is.null(seeds) & is.null(cv_reps)) {
    stop("`seeds` and `cv_reps` cannot be both NULL")
  }
  if (is.null(seeds)) {
    seeds = seq_len(cv_reps)
  }
  index_matrix <- make_cv_index_matrix(
    cv_gridsize    = cv_gridsize,
    nfold          = nfold,
    seeds          = seeds,
    alpha_lambdas  = alpha_lambdas,
    theta_lambdas  = theta_lambdas
  )
  print(index_matrix)
  save_path <- file.path(save_dir, "index_matrix.Rdata")
  save(index_matrix, file = save_path)
  
  #— Launch parallel cluster —#
  if (identical(n_cores, "max")) {
    workers <- parallel::detectCores(logical = FALSE) - 1
  } else {
    workers <- n_cores
  }
  cl <- parallel::makeCluster(workers)
  parallel::clusterExport(
    cl,
    varlist = c(
      "main",
      "binning",
      "initialize_model",
      "iteration",
      "compute_residuals",
      "pi_k",
      "mstep_update_alpha",
      "mstep_update_theta_log_concave",
      "mstep_update_intercepts",
      "mstep_estimate_log_concave_densities",
      "compute_surrogate_loglikelihood",
      "modified_logcondens",
      "e_step_log_concave",
      "weighted_quantile",
      "evaluate_lcd_model",
      "Y_bin", "X", "bin_mass", "K",
      "n_restarts", "max_iter", "iter_eta", "maxdev",
      "resp_threshold", "trim_prob", "save_dir",
      "folds", "index_matrix"
    ),
    envir = environment()
  )
  
  #— Run CV in parallel —#
  logs <- parallel::parLapply(
    cl,
    seq_len(nrow(index_matrix)),
    function(ii) {
      # Extract combination parameters
      alpha_idx   <- index_matrix[ii, "alpha_idx"]
      theta_idx   <- index_matrix[ii, "theta_idx"]
      seed_idx    <- index_matrix[ii, "seed_idx"]
      fold_idx    <- index_matrix[ii, "fold_idx"]
      lambda_alpha<- index_matrix[ii, "lambda_alpha"]
      lambda_theta<- index_matrix[ii, "lambda_theta"]
      
      log_msg <- paste0(
        "alpha=", lambda_alpha,
        ", theta=", lambda_theta,
        ", seed=", seed_idx,
        ", fold=", fold_idx, "\n"
      )
      
      # Split train/test
      test_i     <- folds[[fold_idx]]
      Y_tr       <- Y_bin[-test_i]
      bin_mass_tr<- bin_mass[-test_i]
      X_tr       <- X[-test_i, , drop = FALSE]
      
      # Fit model
      set.seed(seed_idx)
      out_log <- capture.output(
        res_ii <- tryCatch(
          main(
            Y               = Y_tr,
            X               = X_tr,
            biomass         = bin_mass_tr,
            binned          = TRUE,
            n_bins          = 0,
            K               = K,
            lambda_alpha    = lambda_alpha,
            lambda_theta    = lambda_theta,
            max_iter        = max_iter,
            iter_eta        = iter_eta,
            resp_threshold  = resp_threshold,
            debug           = TRUE,
            maxdev          = maxdev,
            n_restarts      = n_restarts
          ),
          error = function(e) {
            message("▶ Error fitting model: ", e$message)
            return(NA)
          }
        ),
        type = "message"      # capture both message() and cat()/print()
      )
      
      # append the captured log to your log_msg
      log_msg <- paste0(log_msg,
                        paste0(out_log, collapse = "\n"),
                        "\n")
      
      # now res_ii contains either the fit or NA
      if (!is.list(res_ii)) {
        prop_CV    <- NA
        trimmed_CV <- NA
        log_msg    <- paste0(log_msg, "Fit completely failed\n")
      } else if (!is.list(res_ii$iter)) {
        log_msg <- paste0("Error at EM iteration ", res_ii$iter_partial$failed_iter, 
                          ": ", res_ii$error, "\n")
        prop_CV    <- NA
        trimmed_CV <- NA
      }
      else {
        eval_res <- evaluate_lcd_model(
          model     = res_ii$iter,
          Y         = Y_bin[test_i],
          X         = X[test_i, , drop = FALSE],
          biomass   = bin_mass[test_i],
          trim_prob = trim_prob
        )
        prop_CV    <- eval_res$prop_infinite
        trimmed_CV <- eval_res$trimmed_avg_loglike
        log_msg    <- paste0(log_msg, "Fit succeeded\n")
      }
      
      # Save results
      save_path <- file.path(
        save_dir,
        sprintf("%d-%d-%d-%d.Rdata",
                alpha_idx, theta_idx, seed_idx, fold_idx)
      )
      save(prop_CV, trimmed_CV, log_msg, file = save_path)
      return(log_msg)
    }
  )
  
  parallel::stopCluster(cl)
  
  #— Summarize failures —#
  num_na <- sum(grepl("Error", logs))
  logs   <- c(
    logs,
    sprintf("Failures: %d/%d (%.1f%%)",
            num_na, length(logs), 100 * num_na / length(logs))
  )
  
  return(list(
    logs          = logs,
    index_matrix  = index_matrix
  ))
}
