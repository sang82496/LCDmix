# Generated from create-LCDmix.Rmd: do not edit by hand

#' Run cross-validation jobs across many simulations in parallel
#'
#' @description
#' Executes one CV job per row of \code{simul_idx_mat} in parallel. Each row
#' identifies a simulation file, a fold, a random seed, and a pair of
#' penalty indices/values. For each job, the function:
#' (1) loads the simulation from \code{sim_dir}, (2) rebuilds CV folds,
#' (3) fits the LCD mixture model on the train split with \code{main()},
#' (4) evaluates on the held-out split via \code{evaluate_lcd_model()},
#' and (5) saves \code{prop_CV}, \code{trimmed_CV}, and a text log to \code{save_dir}.
#'
#' @param simul_idx_mat Numeric matrix where each row is one job. Must contain
#'   columns \code{sim_idx}, \code{alpha_idx}, \code{theta_idx}, \code{seed_idx},
#'   \code{fold_idx}, and the corresponding \code{lambda_alpha}, \code{lambda_theta}
#'   values (e.g., as produced by \code{simul_idx_matrix()}).
#' @param sim_dir Character path to the directory with simulation files named
#'   \code{sim_<sim_idx>.rds}. Default: \code{"sim_data"}.
#' @param K Integer number of mixture components.
#' @param max_iter Integer; maximum EM iterations per fit. Default: 30.
#' @param iter_eta Numeric; convergence threshold on relative change in Q. Default: 1e-3.
#' @param resp_threshold Numeric in [0,1]; responsibilities below this are zeroed
#'   for numerical stability. Default: 1e-3.
#' @param nfold Integer number of CV folds (used when rebuilding folds for each simulation).
#'   Default: 5.
#' @param trim_prob Numeric in [0,0.5]; fraction of lowest-likelihood points to trim
#'   when averaging held-out log-likelihood. Default: 0.01.
#' @param save_dir Character path to write per-job results (\code{*.Rdata}).
#'   Files are saved as \code{<sim_idx>-<alpha_idx>-<theta_idx>-<seed_idx>-<fold_idx>.Rdata}.
#'   Default: \code{"cv_saves"}.
#' @param n_cores Integer or \code{"max"}; number of parallel workers. \code{"max"}
#'   uses all physical cores minus one. Default: \code{"max"}.
#' @param blocksize Integer block size passed to \code{flowmix::make_cv_folds()}.
#'   Default: 5.
#' @param maxdev Numeric or \code{NULL}; max-deviation constraint for LP updates
#'   in theta-step. Default: \code{NULL}.
#' @param n_restarts Integer; number of random restarts for initialization in \code{main()}.
#'   Default: 1.
#'
#' @return A character vector of per-job log messages, followed by a summary line
#'   reporting the number and fraction of failures (jobs where fitting or evaluation failed).
#'
#' @examples
#' \dontrun{
#' # Build an index matrix over sims × folds × seeds × penalties
#' idx <- simul_idx_matrix(
#'   n_sims        = 2,
#'   alpha_lambdas = c(1e-3, 1e-2),
#'   theta_lambdas = c(1e-3, 1e-2),
#'   seeds         = 1:2,
#'   nfold         = 3
#' )
#' # Run CV jobs and save results under "cv_saves/"
#' logs <- run_cv_simul(
#'   simul_idx_mat = idx,
#'   sim_dir       = "sim_data",
#'   K             = 2,
#'   save_dir      = "cv_saves",
#'   n_cores       = "max"
#' )
#' cat(tail(logs), sep = "\n")
#' }
#'
#' @export
run_cv_simul <- function(
  simul_idx_mat,
  sim_dir  = "sim_data",
  K,
  max_iter           = 30,
  iter_eta           = 1e-3,
  resp_threshold     = 1e-3,
  nfold              = 5,
  trim_prob          = 0.01,
  save_dir           = "cv_saves",
  n_cores            = "max",
  blocksize          = 5,
  maxdev             = NULL,
  n_restarts         = 1
) {
  # ensure output folder
  if (!dir.exists(save_dir)) dir.create(save_dir)
  
  # determine worker count
  n_workers <- if (identical(n_cores, "max")) parallel::detectCores(logical = FALSE) else as.integer(n_cores)
  cl <- parallel::makeCluster(n_workers)
  # export functions and data that workers need
  parallel::clusterEvalQ(cl, {library(flowmix); library(LCDmix); NULL })
  parallel::clusterExport(
    cl,
    varlist = c("simul_idx_mat", "sim_dir", "K", "n_restarts", "iter_eta", "maxdev", 
                "resp_threshold", "trim_prob", "nfold", "blocksize", "save_dir", "max_iter"),
    envir   = environment()
  )
  
  logs <- parallel::parLapply(cl, seq_len(nrow(simul_idx_mat)), function(i) {
    idx_row       <- simul_idx_mat[i, ]
    sim_idx       <- as.integer(idx_row["sim_idx"])
    seed_idx      <- as.integer(idx_row["seed_idx"])
    alpha_idx     <- as.integer(idx_row["alpha_idx"])
    theta_idx     <- as.integer(idx_row["theta_idx"])
    fold_idx      <- as.integer(idx_row["fold_idx"])
    lambda_alpha  <- as.numeric(idx_row["lambda_alpha"])
    lambda_theta  <- as.numeric(idx_row["lambda_theta"])
    
    log_msg  <- character(0)
    
    sim_file <- file.path(sim_dir, paste0("sim_", sim_idx, ".rds"))
    if (!file.exists(sim_file)) {
      log_msg <- paste0("Sim file missing: ", sim_file)
      return(log_msg)
    }
    
    # Skip if already exists
    save_path <- file.path(
      save_dir,
      sprintf("%d-%d-%d-%d-%d.rds",
              sim_idx, alpha_idx, theta_idx, seed_idx, fold_idx)
      )
    if (file.exists(save_path)) {
      log_msg <- paste0(log_msg, paste0("Skip (cached): ", basename(save_path)))
      return(log_msg)
    }

    
    # 1) load sim
    sim <- readRDS(sim_file)
    
    # 2) split into train / test
    Y_bin      <- sim$ylist
    X          <- sim$X
    bin_mass   <- sim$countslist
    
    # Create CV folds
    folds <- flowmix::make_cv_folds(
      ylist     = Y_bin,
      nfold     = nfold,
      blocksize = blocksize
    )
    
    # fold splits by time‐index
    test_i     <- folds[[fold_idx]]
    Y_tr       <- Y_bin[-test_i]
    bin_mass_tr<- bin_mass[-test_i]
    X_tr       <- X[-test_i, , drop = FALSE]
    
    log_msg <- paste0(
        "alpha=",   lambda_alpha,
        ", theta=", lambda_theta,
        ", seed=",  seed_idx,
        ", fold=",  fold_idx, "\n"
      )
    
    # 3) fit main()
    set.seed(seed_idx)
    err_msg <- NULL
    out_log <- capture.output(
      fit_try <- tryCatch({
        main(
          Y              = Y_tr,
          X              = X_tr,
          biomass        = bin_mass_tr,
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
        )
      }, error = function(e) {
        err_msg <<- paste0("Error fitting sim=", sim_idx,
                            ", alpha=", lambda_alpha,
                            ", theta=", lambda_theta, 
                              " seed=", seed_idx, 
                              " fold=", fold_idx, ": ",
                               e$message)
        return(NULL)
      }),
        type = "output"
      )
    
    log_msg <- paste0(log_msg,
                        paste0(out_log, collapse = "\n"),
                        "\n")
    
    # now fit_try contains either the fit or NA
      if (is.null(fit_try)) {
        log_msg    <- paste0(log_msg, "Error: Fit completely failed: ", err_msg, "\n")
        prop_CV    <- NA
        trimmed_CV <- NA
        Q_final    <- NA
        
      } else if (is.null(fit_try$iter)) {
        log_msg <- paste0(log_msg, "Error at EM iteration ", fit_try$iter_partial$failed_iter, 
                          ": ", fit_try$iter_partial$error, "\n")
        prop_CV    <- NA
        trimmed_CV <- NA
        Q_final    <- NA
        
      } else {
        ev <- evaluate_lcd_model(
          model     = fit_try$iter,
          Y         = Y_bin[test_i],
          X         = X[test_i, , drop = FALSE],
          biomass   = bin_mass[test_i],
          trim_prob = trim_prob
        )
        prop_CV    <- ev$prop_infinite
        trimmed_CV <- ev$trimmed_avg_loglike
        log_msg <- paste0(log_msg,
          "Success sim=", sim_idx,
          " fold=",       fold_idx,
          " seed=",       seed_idx,
          " alpha=",      round(lambda_alpha,5),
          " theta=",      round(lambda_theta,5),
          " -> prop=",    round(prop_CV,3),
          ", trimmed=",   round(trimmed_CV,3)
        )
        Q_final    <- tail(fit_try$iter$Q, 1)
      }
    
    # 5) save result
    saveRDS(list(prop_CV    = prop_CV,
                 trimmed_CV = trimmed_CV,
                 Q_final    = Q_final,
                 log_msg    = log_msg), file = save_path)
    return(log_msg)
  })
  parallel::stopCluster(cl)
  
  # 7) Summarize failures
  num_na <- sum(grepl("Error", logs))
  logs   <- c(
    logs,
    sprintf("Failures: %d/%d (%.1f%%)",
            num_na, length(logs), 100 * num_na / length(logs))
  )
  return(logs)
}
