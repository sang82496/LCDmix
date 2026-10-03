# Generated from create-LCDmix.Rmd: do not edit by hand

#' Cross-validation of the penalty pair for many simulated datasets
#'
#' @description
#' Runs every cross-validation job of every dataset in one parallel cluster
#' (no nested parallelism). For dataset \eqn{s}, the results go to
#' \code{base_dir/sim_<s>/}: one file per job from \code{cv_lcd_onejob()},
#' and \code{index_matrix.rds}. A job whose file already exists is not run
#' again, so an interrupted run can be resumed. Each worker reads its dataset
#' from \code{sim_files} and builds the folds with
#' \code{flowmix::make_cv_folds()}. Use \code{cv_lcd_summary_simul()} to
#' choose the penalty pair of each dataset, then \code{refit_lcd_simul()}.
#'
#' The workers load the installed LCDmix package, so install the current
#' version before running.
#'
#' @param sim_files Paths of \code{.rds} files, one per dataset; each holds a
#'   list with \code{Y_bin}, \code{X} and \code{bin_mass}, as written by
#'   \code{simulate_and_save()}.
#' @param base_dir Output directory.
#' @inheritParams cv_lcd
#' @inheritParams main
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{grand_jobs}}{Data frame with one row per job: \code{sim_idx}
#'     and the columns of \code{cv_idx_mat()}.}
#'   \item{\code{summary}}{A text with the number of failed jobs.}
#' }
#'
#' @seealso \code{\link{cv_lcd}}, \code{\link{cv_lcd_summary_simul}},
#'   \code{\link{refit_lcd_simul}}, \code{\link{eval_lcd}},
#'   \code{\link[flowmix]{make_cv_folds}}
#'
#' @examples
#' \dontrun{
#' # 25 datasets, 5 x 5 penalty grid, 10 seeds, 5 folds (31,250 jobs)
#' sims <- file.path("sim_data", sprintf("sim_%d.rds", 1:25))
#'
#' res <- cv_lcd_simul(
#'   sim_files     = sims,
#'   K             = 2,
#'   alpha_lambdas = 10^seq(-4, -2, length.out = 5),
#'   theta_lambdas = 10^seq(-4, -2, length.out = 5),
#'   nfold         = 5,
#'   seeds         = 1:10,
#'   base_dir      = "cv_runs",
#'   n_cores       = 64
#' )
#'
#' # Chosen penalty pair of each dataset (results in cv_runs/sim_<s>/)
#' summ <- cv_lcd_summary_simul(base_dir = "cv_runs", num_sims = 25)
#' opt_list <- lapply(summ, function(x) x$opt_lambdas)
#' }
#'
#' @export
cv_lcd_simul <- function(
  sim_files,
  K,
  alpha_lambdas,
  theta_lambdas,
  nfold          = 5,
  seeds          = NULL,
  cv_reps        = NULL,
  max_iter       = 30,
  iter_eta       = 1e-4,
  resp_threshold = 1e-3,
  trim_prob      = 0.03,
  blocksize      = 10,
  base_dir       = "./cv_saves",
  n_cores        = "max",
  lp_time_limit  = 600,
  calc_Q_every   = FALSE,                # NEW - appended
  update         = c("lp", "optim"),    # NEW
  maxdev         = NULL                 # NEW (fixP) - appended last
) {
  if (is.null(seeds) && is.null(cv_reps)) stop("`seeds` and `cv_reps` cannot be both NULL")
  if (is.null(seeds)) seeds <- seq_len(cv_reps)
  if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE)

  num_sims <- length(sim_files)
  alpha_lambdas <- sort(alpha_lambdas)
  theta_lambdas <- sort(theta_lambdas)

  # --- per-sim index matrices & subdirs ---
  per_sim_idx <- vector("list", num_sims)
  for (s in seq_len(num_sims)) {
    sim_dir <- file.path(base_dir, sprintf("sim_%d", s))
    if (!dir.exists(sim_dir)) dir.create(sim_dir, recursive = TRUE)

    # index matrix for this simulation (no data needed)
    idx <- cv_idx_mat(
      nfold         = nfold,
      seeds         = seeds,
      alpha_lambdas = alpha_lambdas,
      theta_lambdas = theta_lambdas
    )
    # persist for summaries
    saveRDS(idx, file = file.path(sim_dir, "index_matrix.rds"))
    per_sim_idx[[s]] <- cbind(sim_idx = s, idx)
  }

  # --- grand job matrix ---
  grand_jobs <- do.call(rbind, per_sim_idx)
  rownames(grand_jobs) <- NULL
  grand_jobs <- as.data.frame(grand_jobs)

  # One outer cluster
  n_workers <- if (identical(n_cores, "max")) parallel::detectCores(logical = FALSE) else as.integer(n_cores)
  cl <- parallel::makeCluster(n_workers)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

  parallel::clusterEvalQ(cl, { library(flowmix); library(LCDmix); NULL })
  parallel::clusterExport(
    cl,
    varlist = c("sim_files", "grand_jobs", "K", "max_iter", "iter_eta",
                "resp_threshold", "trim_prob", "blocksize", "base_dir", 
                "nfold", "lp_time_limit", "update", "calc_Q_every",
                "maxdev"),                                                   # NEW (fixP)
    envir = environment()
  )

   # --- worker: one job per task; return TRUE on success, FALSE on fail ---
  res <- parallel::parLapply(cl, seq_len(nrow(grand_jobs)), function(ii) {
    row <- grand_jobs[ii, ]
    s   <- as.integer(row[["sim_idx"]])
    a_i <- as.integer(row[["alpha_idx"]])
    t_i <- as.integer(row[["theta_idx"]])
    sd_i<- as.integer(row[["seed_idx"]])
    f_i <- as.integer(row[["fold_idx"]])
    la  <- as.numeric(row[["lambda_alpha"]])
    lt  <- as.numeric(row[["lambda_theta"]])
    
    # ensure sim_<s> directory exists
    sim_dir <- file.path(base_dir, sprintf("sim_%d", s))
    if (!dir.exists(sim_dir)) dir.create(sim_dir, recursive = TRUE)
    
    out_path <- file.path(
      sim_dir,
      sprintf("%d-%d-%d-%d.rds", a_i, t_i, sd_i, f_i)
    )
    if (file.exists(out_path)) {
      res_ii = readRDS(out_path)
      return(!is.na(res_ii$fit_trimmed_loglik))
    }

    # load the one simulation needed for this job
    sim = readRDS(sim_files[s])
    Y_bin   <- sim$Y_bin
    X       <- sim$X
    bin_mass<- sim$bin_mass

    # folds for this simulation (deterministic given inputs)
    folds <- flowmix::make_cv_folds(ylist = Y_bin, nfold = nfold, blocksize = blocksize)

    # build job descriptor expected by cv_lcd_onejob()
    job <- c(alpha_idx = a_i, theta_idx = t_i, seed_idx = sd_i, fold_idx = f_i, 
             lambda_alpha = la, lambda_theta = lt)

    # run one CV job; cv_lcd_onejob writes its own RDS and returns a log string
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
        save_dir       = sim_dir,
        lp_time_limit  = lp_time_limit,
        calc_Q_every   = calc_Q_every,       # NEW
        update         = update,
        maxdev         = maxdev              # NEW (fixP)
      )
      return(res_ii)
  })

  success   <- unlist(res, use.names = FALSE)
  summary   <- sprintf("Failures: %d/%d (%.1f%%)", sum(!success), length(success), 100 * sum(!success)/length(success))
  
  return(list(grand_jobs = grand_jobs,
              summary = summary))
}
