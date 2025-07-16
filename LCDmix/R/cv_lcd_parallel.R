# Generated from create-LCDmix.Rmd: do not edit by hand

#' Cross‐validate the LCD mixture‐of‐experts model in parallel
#'
#' @description
#' Performs grid‐search cross‐validation over the L1 penalties \code{lambda_alpha}
#' and \code{lambda_theta} for the log‐concave mixture‐of‐experts model.  Creates
#' folds (with optional blocking), runs multiple repetitions per fold, and evaluates
#' models on held‐out data.  Computations are dispatched in parallel.
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
#' @param cv_grid_size Integer number of candidate lambdas per penalty. Default: 5.
#' @param n_restarts Integer number of random restarts for \code{flowmix::flowmix()}. Default: 1.
#' @param max_iter Integer maximum number of EM iterations per model fit. Default: 30.
#' @param iter_eta Numeric convergence threshold (relative change in Q). Default: 1e-3.
#' @param maxdev Numeric or \code{NULL}; max‐deviation constraint for LP updates. Default: \code{NULL}.
#' @param resp_threshold Numeric in \[0,1\]; responsibilities below this are zeroed. Default: 1e-3.
#' @param nfold Integer number of cross‐validation folds. Default: 5.
#' @param block_size Integer block size for folding. Default: 5.
#' @param cv_reps Integer number of repeated CV runs. Default: 5.
#' @param trim_prob Numeric in \[0,1\]; fraction of lowest‐likelihood points to trim when averaging. Default: 0.05.
#' @param save_dir Character; directory in which to save per‐run results. Default: "./result".
#' @param sparseMatrix Logical; if \code{TRUE}, build LP constraints as sparse matrices. Default: \code{TRUE}.
#' @param n_cores Integer or "max"; number of parallel worker processes. "max" uses all physical cores minus one. Default: "max".
#'
#' @return A list with components:
#' \describe{
#'   \item{logs}{Character vector of log messages (one per CV run, plus a summary).}
#'   \item{index_matrix}{Numeric matrix with one row per combination (fold × rep × penalties),
#'         showing indices and actual \code{lambda_alpha}, \code{lambda_theta} values.}
#' }
#'
#' @examples
#' \dontrun{
#' # 2‐component model, small grid, 3 folds, 2 reps
#' results <- cv_lcd_parallel(
#'   Y_bin               = Y_binned,
#'   X                   = X_mat,
#'   bin_mass            = bin_mass_list,
#'   K                   = 2,
#'   lambda_alpha_range  = c(1e-4, 1e-2),
#'   lambda_theta_range  = c(1e-4, 1e-2),
#'   cv_grid_size        = 4,
#'   nfold             = 3,
#'   cv_reps             = 2,
#'   n_cores             = "max"
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
  cv_grid_size       = 5,
  n_restarts         = 1,
  max_iter           = 30,
  iter_eta           = 1e-3,
  maxdev             = NULL,
  resp_threshold     = 1e-3,
  nfold            = 5,
  block_size         = 5,
  cv_reps            = 5,
  trim_prob          = 0.05,
  save_dir           = "./result",
  sparseMatrix      = TRUE,
  n_cores            = "max"
) {
  #— Create penalty grids —#
  alpha_lambdas <- flowmix::logspace(
    min    = lambda_alpha_range[1],
    max    = lambda_alpha_range[2],
    length = cv_grid_size
  ) |> sort()
  theta_lambdas <- flowmix::logspace(
    min    = lambda_theta_range[1],
    max    = lambda_theta_range[2],
    length = cv_grid_size
  ) |> sort()
  print(alpha_lambdas)
  print(theta_lambdas)
  
  #— Create CV folds —#
  folds <- flowmix::make_cv_folds(
    ylist     = Y_bin,
    nfold     = nfold,
    blocksize = block_size
  )
  
  #— Build index matrix of all (fold × rep × penalties) combos —#
  index_matrix <- make_cv_index_matrix(
    cv_gridsize   = cv_grid_size,
    nfold          = nfold,
    nrep           = cv_reps,
    alpha_lambdas  = alpha_lambdas,
    theta_lambdas  = theta_lambdas
  )
  print(index_matrix)
  
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
      "sparseMatrix", "folds", "index_matrix"
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
      rep_idx     <- index_matrix[ii, "rep_idx"]
      fold_idx    <- index_matrix[ii, "fold_idx"]
      lambda_alpha<- index_matrix[ii, "lambda_alpha"]
      lambda_theta<- index_matrix[ii, "lambda_theta"]
      
      log_msg <- paste0(
        "alpha=", lambda_alpha,
        ", theta=", lambda_theta,
        ", rep=", rep_idx,
        ", fold=", fold_idx, "\n"
      )
      
      # Split train/test
      test_i     <- folds[[fold_idx]]
      Y_tr       <- Y_bin[-test_i]
      bin_mass_tr<- bin_mass[-test_i]
      X_tr       <- X[-test_i, , drop = FALSE]
      
      # Fit model
      set.seed(rep_idx)
      fit_try <- tryCatch(
        main(
          Y               = Y_tr,
          X               = X_tr,
          biomass         = bin_mass_tr,
          binned          = TRUE,
          n_bins          = 0,
          K               = K,
          lambda_alpha    = lambda_alpha,
          lambda_theta    = lambda_theta,
          n_restarts      = n_restarts,
          max_iter        = max_iter,
          iter_eta        = iter_eta,
          maxdev          = maxdev,
          resp_threshold  = resp_threshold,
          sparseMatrix   = sparseMatrix
        ),
        error = function(e) e
      )
      
      if (inherits(fit_try, "error")) {
        log_msg <- paste0(log_msg, "Error: ", fit_try$message, "\n")
        prop_CV    <- NA
        trimmed_CV <- NA
      } else {
        # Evaluate on held‐out fold
        eval_res <- evaluate_lcd_model(
          model     = fit_try$iter,
          Y         = Y_bin[test_i],
          X         = X[test_i, , drop = FALSE],
          biomass   = bin_mass[test_i],
          trim_prob = trim_prob
        )
        prop_CV    <- eval_res$prop_infinite
        trimmed_CV <- eval_res$trimmed_avg_loglike
        log_msg <- paste0(log_msg, "Success\n")
      }
      
      # Save results
      save_path <- file.path(
        save_dir,
        sprintf("%d-%d-%d-%d.Rdata",
                alpha_idx, theta_idx, rep_idx, fold_idx)
      )
      save(prop_CV, trimmed_CV, file = save_path)
      return(log_msg)
    }
  )
  
  parallel::stopCluster(cl)
  
  #— Summarize failures —#
  num_na <- sum(grepl("Error:", logs))
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
