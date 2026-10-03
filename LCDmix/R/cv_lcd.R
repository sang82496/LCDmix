# Generated from create-LCDmix.Rmd: do not edit by hand

#' Cross-validation of the penalty pair for one dataset
#'
#' @description
#' Runs \code{cv_lcd_onejob()} in parallel for every combination of fold,
#' seed, \code{lambda_alpha} and \code{lambda_theta}, and saves one file per
#' job in \code{save_dir}, together with \code{index_matrix.rds}. A job whose
#' file already exists is not run again, so an interrupted run can be resumed.
#' The folds are blocks of consecutive time points from
#' \code{flowmix::make_cv_folds()}. Use \code{cv_lcd_summary()} on
#' \code{save_dir} to choose the penalty pair, then \code{refit_lcd()}.
#'
#' The workers load the installed LCDmix package, so install the current
#' version before running.
#'
#' @inheritParams main
#' @inheritParams initialization
#' @param alpha_lambdas,theta_lambdas Candidate values of \code{lambda_alpha}
#'   and \code{lambda_theta}. They are sorted in increasing order first.
#' @param nfold Number of folds.
#' @param seeds Seeds (restarts) for every fold and penalty pair. If
#'   \code{NULL}, \code{1:cv_reps}.
#' @param save_dir Output directory; it is created if needed.
#' @param n_cores Number of worker processes, or \code{"max"} for all physical
#'   cores.
#' @param cv_reps Number of seeds, used only when \code{seeds} is \code{NULL}.
#' @param blocksize Block size passed to \code{flowmix::make_cv_folds()}.
#'
#' @return A list with \code{index_matrix} (from \code{cv_idx_mat()}) and
#'   \code{summary} (a text with the number of failed jobs).
#'
#' @seealso \code{\link{cv_lcd_simul}} for many datasets,
#'   \code{\link{cv_lcd_summary}}, \code{\link{refit_lcd}}
#'
#' @export
cv_lcd <- function(
  Y_bin,
  X,
  bin_mass,
  K,
  alpha_lambdas      = c(1e-4, 1e-3),
  theta_lambdas      = c(1e-4, 1e-3),
  max_iter           = 30,
  iter_eta           = 1e-4,
  resp_threshold     = 1e-3,
  nfold              = 5,
  seeds              = NULL,
  trim_prob          = 0.03,
  save_dir           = "./cv_saves",
  n_cores            = "max",
  cv_reps            = NULL,
  blocksize          = 20,
  lp_time_limit      = 3600,
  update             = c("lp", "optim"),    # NEW
  maxdev             = NULL                 # NEW (fixP) - appended last
) {
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  alpha_lambdas <- sort(alpha_lambdas)
  theta_lambdas <- sort(theta_lambdas)

  folds <- flowmix::make_cv_folds(ylist = Y_bin, nfold = nfold, blocksize = blocksize)

  if (is.null(seeds) && is.null(cv_reps)) stop("`seeds` and `cv_reps` cannot be both NULL")
  if (is.null(seeds)) seeds <- seq_len(cv_reps)

  index_matrix <- cv_idx_mat(
    nfold         = nfold,
    seeds         = seeds,
    alpha_lambdas = alpha_lambdas,
    theta_lambdas = theta_lambdas
  )
  saveRDS(index_matrix, file = file.path(save_dir, "index_matrix.rds"))

  n_workers <- if (identical(n_cores, "max")) parallel::detectCores(logical = FALSE) else as.integer(n_cores)
  cl <- parallel::makeCluster(n_workers)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

  parallel::clusterEvalQ(cl, { library(flowmix); library(LCDmix); NULL })
  parallel::clusterExport(
    cl,
    varlist = c("Y_bin","X","bin_mass","K","max_iter","iter_eta","resp_threshold",
                "trim_prob","save_dir","folds","index_matrix", "lp_time_limit", "update",
                "maxdev"),                                                   # NEW (fixP)
    envir = environment()
  )

  res <- parallel::parLapply(
    cl,
    seq_len(nrow(index_matrix)),
    function(ii) {
      job <- index_matrix[ii, , drop = FALSE][1, ]
      # ensure names present for worker
      job <- setNames(as.numeric(job), colnames(index_matrix))
      alpha_idx     <- job[["alpha_idx"]]
      theta_idx     <- job[["theta_idx"]]
      seed_idx      <- job[["seed_idx"]]
      fold_idx      <- job[["fold_idx"]]
      
      out_path <- file.path(
        save_dir,
        sprintf("%d-%d-%d-%d.rds", alpha_idx, theta_idx, seed_idx, fold_idx)
      )
      
      if (file.exists(out_path)) {
        res_ii = readRDS(out_path)
        return(!is.na(res_ii$fit_trimmed_loglik))
      }

      res_ii = cv_lcd_onejob(
        job            = job,
        Y_bin          = Y_bin,
        X              = X,
        bin_mass       = bin_mass,
        folds          = folds,
        K              = K,
        max_iter       = max_iter,
        iter_eta       = iter_eta,
        resp_threshold = resp_threshold,
        trim_prob      = trim_prob,
        save_dir       = save_dir,
        lp_time_limit  = lp_time_limit,
        update         = update,
        maxdev         = maxdev              # NEW (fixP)
      )
      return(res_ii)
    }
  )
  success   <- unlist(res, use.names = FALSE)
  summary   <- sprintf("Failures: %d/%d (%.1f%%)", sum(!success), 
                       length(success), 100 * sum(!success)/length(success))

  return(list(index_matrix = index_matrix, 
              summary = summary))
}
