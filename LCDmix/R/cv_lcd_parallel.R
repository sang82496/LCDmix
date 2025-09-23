# Generated from create-LCDmix.Rmd: do not edit by hand

#' Cross-validate LCDmix over a penalty grid (parallel)
#'
#' Runs grid-search cross-validation across \code{alpha_lambdas} × \code{theta_lambdas}.
#' For each (\code{lambda_alpha}, \code{lambda_theta}) × (seed, fold), fits the training split
#' with \code{main()} and evaluates the held-out split via \code{evaluate_lcd_model()},
#' saving per-run results to disk. Existing result files are skipped on re-runs.
#'
#' @param Y_bin List of length \eqn{TT}; element \eqn{t} is a numeric vector of length \eqn{n_t}
#'   containing binned responses at time \eqn{t}.
#' @param X Numeric matrix \eqn{TT \times p}; covariates aligned by time.
#' @param bin_mass List of length \eqn{TT}; per-time weights for each observation/bin.
#' @param K Integer; number of mixture components.
#' @param alpha_lambdas Numeric vector of candidate \code{lambda_alpha} values. Default \code{c(1e-4, 1e-3)}.
#' @param theta_lambdas Numeric vector of candidate \code{lambda_theta} values. Default \code{c(1e-4, 1e-3)}.
#' @param max_iter Integer; maximum EM iterations per fit. Default \code{30}.
#' @param iter_eta Numeric; convergence threshold on relative change in the surrogate objective. Default \code{1e-3}.
#' @param resp_threshold Numeric in \eqn{[0,1]}; responsibilities below this are zeroed. Default \code{1e-3}.
#' @param nfold Integer; number of CV folds. Default \code{5}.
#' @param seeds Integer vector of random seeds (one per repeat). If \code{NULL}, use \code{cv_reps}. Default \code{NULL}.
#' @param trim_prob Numeric in \eqn{[0,1)}; trimming fraction used inside \code{evaluate_lcd_model()}. Default \code{0.01}.
#' @param save_dir Character; directory for outputs (\code{index_matrix.Rdata} and per-run \code{*.rds}). Default \code{"./result"}.
#' @param n_cores Integer or \code{"max"}; number of parallel workers. \code{"max"} uses all physical cores minus one. Default \code{"max"}.
#' @param cv_reps Integer; number of repeats, used only when \code{seeds} is \code{NULL}. Default \code{NULL}.
#' @param blocksize Integer; block size passed to \code{flowmix::make_cv_folds()}. Default \code{20}.
#'
#' @details
#' Folds are constructed with \code{flowmix::make_cv_folds()}, preserving temporal blocks.
#' Each worker loads \pkg{flowmix} and \pkg{LCDmix}. Results are cached per setting in
#' files named \code{"<alpha_idx>-<theta_idx>-<seed_idx>-<fold_idx>.rds"} under \code{save_dir}.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{logs}}{Character vector of per-run logs plus a final failure summary line.}
#'   \item{\code{index_matrix}}{Numeric matrix with columns \code{alpha_idx}, \code{theta_idx},
#'         \code{seed_idx}, \code{fold_idx}, \code{lambda_alpha}, \code{lambda_theta}.}
#' }
#'
#' @seealso \code{\link{make_cv_index_matrix}}, \code{\link{evaluate_lcd_model}},
#'   \code{\link[flowmix]{make_cv_folds}}
#'
#' @examples
#' \dontrun{
#' out <- cv_lcd_parallel(
#'   Y_bin          = Y_bin,
#'   X              = X,
#'   bin_mass       = bin_mass,
#'   K              = 2,
#'   alpha_lambdas  = 10^seq(-4, -2, length.out = 3),
#'   theta_lambdas  = 10^seq(-4, -2, length.out = 3),
#'   nfold          = 5,
#'   seeds          = 1:3,
#'   save_dir       = "cv_runs"
#' )
#' str(out$index_matrix)
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
  trim_prob          = 0.01,
  save_dir           = "./result",
  n_cores            = "max",
  cv_reps            = NULL,
  blocksize          = 20
) {
  # Ensure output directory exists (nested ok).
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  # Sort penalties for reproducible ordering in filenames/logs.
  alpha_lambdas <- sort(alpha_lambdas)
  theta_lambdas <- sort(theta_lambdas)

  # — Create CV folds —
  folds <- flowmix::make_cv_folds(
    ylist     = Y_bin,
    nfold     = nfold,
    blocksize = blocksize
  )

  # — Build index matrix of all (fold × rep × penalties) combos —
  if (is.null(seeds) && is.null(cv_reps)) {
    stop("`seeds` and `cv_reps` cannot be both NULL")
  }
  if (is.null(seeds)) {
    seeds <- seq_len(cv_reps)
  }
  index_matrix <- make_cv_index_matrix(
    nfold          = nfold,
    seeds          = seeds,
    alpha_lambdas  = alpha_lambdas,
    theta_lambdas  = theta_lambdas
  )

  save(index_matrix, file = file.path(save_dir, "index_matrix.Rdata"))

  # — Launch parallel cluster —
  n_workers <- if (identical(n_cores, "max")) parallel::detectCores(logical = FALSE) else as.integer(n_cores)
  cl <- parallel::makeCluster(n_workers)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

  parallel::clusterEvalQ(cl, { library(flowmix); library(LCDmix); NULL })
  parallel::clusterExport(
    cl,
    varlist = c(
      "Y_bin", "X", "bin_mass", "K", "max_iter", "iter_eta",
      "resp_threshold", "trim_prob", "save_dir", "folds", "index_matrix"
    ),
    envir = environment()
  )

  # — Run CV in parallel —
  logs <- parallel::parLapply(
    cl,
    seq_len(nrow(index_matrix)),
    function(ii) {
      # Extract combination parameters
      alpha_idx     <- index_matrix[ii, "alpha_idx"]
      theta_idx     <- index_matrix[ii, "theta_idx"]
      seed_idx      <- index_matrix[ii, "seed_idx"]
      fold_idx      <- index_matrix[ii, "fold_idx"]
      lambda_alpha  <- index_matrix[ii, "lambda_alpha"]
      lambda_theta  <- index_matrix[ii, "lambda_theta"]

      out_path <- file.path(
        save_dir,
        sprintf("%d-%d-%d-%d.rds", alpha_idx, theta_idx, seed_idx, fold_idx)
      )
      if (file.exists(out_path)) {
        return(paste0("Skipping existing: alpha=", lambda_alpha,
                      ", theta=", lambda_theta,
                      ", seed=",  seed_idx,
                      ", fold=",  fold_idx))
      }

      log_msg <- paste0(
        "alpha=", lambda_alpha,
        ", theta=", lambda_theta,
        ", seed=", seed_idx,
        ", fold=", fold_idx, "\n"
      )

      # Split train/test
      test_i      <- folds[[fold_idx]]
      Y_tr        <- Y_bin[-test_i]
      bin_mass_tr <- bin_mass[-test_i]
      X_tr        <- X[-test_i, , drop = FALSE]

      # Fit model
      set.seed(seed_idx)
      err_msg <- NULL
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
            trim_prob       = trim_prob,
            debug           = TRUE
          ),
          error = function(e) {
            err_msg <<- paste0("▶ Error fitting model: ", e$message)
            NULL
          }
        ),
        type = "output"
      )

      # Append captured logs
      log_msg <- paste0(log_msg, paste0(out_log, collapse = "\n"), "\n")

      if (!is.list(res_ii)) {
        log_msg       <- paste0(log_msg, "Error: Fit completely failed: ", err_msg, "\n")
        prop_inf      <- NA_real_
        trimmed_loglik<- NA_real_
        finite_loglik <- NA_real_
        L             <- NA_real_

      } else if (is.null(res_ii$iter)) {
        log_msg       <- paste0(log_msg, "Error at EM iteration ",
                                res_ii$iter_partial$failed_iter, ": ",
                                res_ii$iter_partial$error, "\n")
        prop_inf      <- NA_real_
        trimmed_loglik<- NA_real_
        finite_loglik <- NA_real_
        L             <- NA_real_

      } else {
        log_msg <- paste0(log_msg, "Fit succeeded\n")

        eval_res <- evaluate_lcd_model(
          model        = res_ii$iter,
          Y_test       = Y_bin[test_i],
          X_test       = X[test_i, , drop = FALSE],
          biomass_test = bin_mass[test_i],
          trim_prob    = trim_prob
        )
        prop_inf       <- eval_res$prop_inf
        trimmed_loglik <- eval_res$trimmed_loglik
        finite_loglik  <- eval_res$finite_loglik
        L              <- res_ii$L$loglik
      }

      # Save results
      saveRDS(
        list(
          prop_inf       = prop_inf,
          trimmed_loglik = trimmed_loglik,
          finite_loglik  = finite_loglik,
          L              = L,
          log_msg        = log_msg
        ),
        file = out_path
      )

      return(log_msg)
    }
  )

  # — Summarize failures —
  num_na <- sum(grepl("Error", logs, fixed = TRUE))
  logs   <- c(
    logs,
    sprintf("Failures: %d/%d (%.1f%%)", num_na, length(logs), 100 * num_na / length(logs))
  )

  return(list(
    logs         = logs,
    index_matrix = index_matrix
  ))
}
