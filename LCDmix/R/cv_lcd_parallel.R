# Generated from create-LCDmix.Rmd: do not edit by hand

#' Cross-validate LCDmix over a penalty grid in parallel
#'
#' @description
#' Runs grid search cross-validation for the LCDmix model across the supplied
#' \code{alpha_lambdas} and \code{theta_lambdas}. For each combination of
#' (\code{lambda_alpha}, \code{lambda_theta}), CV fold, and seed, the function
#' fits the model on the training split, evaluates the held-out split with
#' \code{evaluate_lcd_model()}, and saves \code{prop_CV} and \code{trimmed_CV}
#' to \code{save_dir}. Workers load \pkg{flowmix} and \pkg{LCDmix} internally.
#'
#' @param Y_bin List of length \eqn{TT}; each element is an \eqn{M_t \times 1}
#'   matrix of binned responses at time \eqn{t}.
#' @param X Numeric \eqn{TT \times p} covariate matrix (rows = time points).
#' @param bin_mass List of length \eqn{TT}; each element is a numeric vector of
#'   biomass per bin at time \eqn{t}.
#' @param K Integer; number of mixture components.
#' @param alpha_lambdas Numeric vector of candidate \code{lambda_alpha} values.
#'   Default: \code{c(1e-4, 1e-3)}.
#' @param theta_lambdas Numeric vector of candidate \code{lambda_theta} values.
#'   Default: \code{c(1e-4, 1e-3)}.
#' @param max_iter Integer; maximum EM iterations per fit. Default: \code{30}.
#' @param iter_eta Numeric; convergence threshold on relative change in Q.
#'   Default: \code{1e-3}.
#' @param resp_threshold Numeric in \[0,1\]; responsibilities below this are set
#'   to zero for numerical stability. Default: \code{1e-3}.
#' @param nfold Integer; number of CV folds. Default: \code{5}.
#' @param seeds Integer vector of random seeds (one per repeat). If \code{NULL},
#'   you must supply \code{cv_reps}; seeds will be set to \code{1:cv_reps}.
#'   Default: \code{NULL}.
#' @param trim_prob Numeric in \[0,1\]; fraction of lowest-likelihood points to
#'   trim when averaging the held-out log-likelihood. Default: \code{0.05}.
#' @param save_dir Character; directory where the function writes the CV index
#'   matrix (\code{index_matrix.Rdata}) and per-run \code{*.Rdata} files of the
#'   form \code{"<alpha_idx>-<theta_idx>-<seed_idx>-<fold_idx>.Rdata"}.
#'   Default: \code{"./result"}.
#' @param n_cores Integer or \code{"max"}; number of parallel workers.
#'   \code{"max"} uses all physical cores minus one. Default: \code{"max"}.
#' @param cv_reps Integer; number of repeats used only when \code{seeds} is
#'   \code{NULL}. Default: \code{NULL}.
#' @param blocksize Integer; block size passed to \code{flowmix::make_cv_folds()}.
#'   Default: \code{5}.
#' @param maxdev Numeric or \code{NULL}; maximum-deviation constraint for the
#'   LP update in the theta-step. Default: \code{NULL}.
#' @param n_restarts Integer; number of random restarts for initialization.
#'   Default: \code{1}.
#'
#' @details
#' On each worker the function calls:
#' \enumerate{
#' \item \code{library(flowmix); library(LCDmix)} to load dependencies,
#' \item \code{main()} on the training split for each (\code{lambda_alpha},
#'       \code{lambda_theta}, seed, fold),
#' \item \code{evaluate_lcd_model()} on the hold-out split to compute
#'       \code{prop_CV} and \code{trimmed_CV}.
#' }
#' The function prints the sorted penalty vectors and writes the constructed
#' \code{index_matrix} to \code{save_dir/index_matrix.Rdata}.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{logs}}{Character vector of per-run log messages plus a final
#'         summary line reporting the number and percent of failures.}
#'   \item{\code{index_matrix}}{Numeric matrix with columns
#'         \code{alpha_idx}, \code{theta_idx}, \code{seed_idx}, \code{fold_idx},
#'         \code{lambda_alpha}, \code{lambda_theta}, corresponding to all CV jobs.}
#' }
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' # Suppose Y_bin, X, bin_mass already prepared
#' res <- cv_lcd_parallel(
#'   Y_bin          = Y_bin,
#'   X              = X,
#'   bin_mass       = bin_mass,
#'   K              = 2,
#'   alpha_lambdas  = c(1e-4, 1e-3, 1e-2),
#'   theta_lambdas  = c(1e-4, 1e-3, 1e-2),
#'   nfold          = 5,
#'   seeds          = 1:3,
#'   save_dir       = "cv_out",
#'   n_cores        = "max"
#' )
#' head(res$index_matrix)
#' tail(res$logs, 1)
#' }
#'
#' @export
cv_lcd_parallel <- function(
  Y_bin,
  X,
  bin_mass,
  K,
  alpha_lambdas      = c(1e-4, 1e-3),
  theta_lambdas      = c(1e-4, 1e-3),
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
  if (!dir.exists(save_dir)) dir.create(save_dir)
  
  alpha_lambdas = sort(alpha_lambdas)
  theta_lambdas = sort(theta_lambdas)
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
    workers <- as.integer(n_cores)
  }
  cl <- parallel::makeCluster(workers)
  parallel::clusterEvalQ(cl, {library(flowmix); library(LCDmix); NULL })
  parallel::clusterExport(
    cl,
    varlist = c(
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
      
      # Skip if already exists
      save_path <- file.path(
        save_dir,
        sprintf("%d-%d-%d-%d.Rdata",
                alpha_idx, theta_idx, seed_idx, fold_idx)
      )
      if (file.exists(save_path)) {
        return(paste0("Skip (cached): ", basename(save_path)))
      }
      
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
        type = "output"     # capture both message() and cat()/print()
      )

      # append the captured log to your log_msg
      log_msg <- paste0(log_msg,
                        paste0(out_log, collapse = "\n"),
                        "\n")
      
      # now res_ii contains either the fit or NA
      if (!is.list(res_ii)) {
        log_msg    <- paste0(log_msg, "Error: Fit completely failed\n")
        prop_CV    <- NA
        trimmed_CV <- NA
        
      } else if (is.null(res_ii$iter)) {
        log_msg <- paste0("Error at EM iteration ", res_ii$iter_partial$failed_iter, 
                          ": ", res_ii$iter_partial$error, "\n")
        prop_CV    <- NA
        trimmed_CV <- NA
        
      } else {
        
        log_msg    <- paste0(log_msg, "Fit succeeded\n")
        eval_res   <- evaluate_lcd_model(
          model     = res_ii$iter,
          Y         = Y_bin[test_i],
          X         = X[test_i, , drop = FALSE],
          biomass   = bin_mass[test_i],
          trim_prob = trim_prob
        )
        prop_CV    <- eval_res$prop_infinite
        trimmed_CV <- eval_res$trimmed_avg_loglike
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
