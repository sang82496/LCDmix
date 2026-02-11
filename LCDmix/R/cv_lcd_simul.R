# Generated from create-LCDmix.Rmd: do not edit by hand

#' Run LCDmix CV across many simulated datasets (outer-parallel, chunked)
#'
#' Processes a collection of simulated datasets using a single parallel cluster
#' (no nested parallel). For each simulation file, builds CV folds, constructs
#' an index over folds × seeds × (\code{lambda_alpha}, \code{lambda_theta}),
#' and evaluates all jobs. Results for each simulation are written into
#' \code{save_dir/sim_<s>/}, one \code{*.rds} per job, and an
#' \code{index_matrix.Rdata} for summary. Existing files are skipped, so runs
#' are resumable.
#'
#' @param sim_files Character vector of length \eqn{S}; paths to simulation files.
#'   Each file must provide objects named \code{Y_bin}, \code{X}, and \code{bin_mass}.
#'   Files may be \code{.rds} (containing a list) or \code{.RData}.
#' @param K Integer; number of mixture components.
#' @param alpha_lambdas Numeric vector of candidate \code{lambda_alpha} values.
#' @param theta_lambdas Numeric vector of candidate \code{lambda_theta} values.
#' @param nfold Integer; number of CV folds per simulation. Default \code{5}.
#' @param seeds Integer vector of seeds (one per repeat). If \code{NULL}, supply
#'   \code{cv_reps}. Default \code{NULL}.
#' @param cv_reps Integer; number of repeats used only when \code{seeds} is
#'   \code{NULL}. Default \code{NULL}.
#' @param max_iter Integer; maximum EM iterations per fit. Default \code{30}.
#' @param iter_eta Numeric; convergence threshold on relative change in the
#'   surrogate objective. Default \code{1e-3}.
#' @param resp_threshold Numeric in \eqn{[0,1]}; responsibilities below this are
#'   zeroed for stability. Default \code{1e-3}.
#' @param trim_prob Numeric in \eqn{[0,1)}; trimming fraction used inside
#'   \code{eval_lcd()}. Default \code{0.03}.
#' @param blocksize Integer; block size for \code{flowmix::make_cv_folds()}.
#'   Default \code{20}.
#' @param save_dir Character; base directory for outputs. Per-simulation results
#'   are written to \code{save_dir/sim_<s>/}. Default \code{"./result"}.
#' @param n_cores Integer or \code{"max"}; number of workers for the single
#'   outer cluster. \code{"max"} uses all physical cores minus one. Default \code{"max"}.
#'
#' @details
#' This function uses a single outer cluster to process a \emph{grand} index of
#' all simulation jobs in chunks, avoiding nested parallel. For reproducibility,
#' job seeds are  (\code{seed_idx})
#' before fitting. Each per-job RDS file contains:
#' \code{prop_inf}, \code{trimmed_loglik}, \code{finite_loglik}, \code{L},
#' and \code{log_msg}. After completion, run \code{\link{cv_lcd_summary}} on each
#' simulation's subdirectory to obtain per-sim CV scores and optimal penalties.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{grand_index}}{Data frame of all jobs with a \code{sim_idx} column.}
#'   \item{\code{logs}}{Character vector of per-chunk log strings.}
#' }
#'
#' @seealso \code{\link{cv_lcd}}, \code{\link{cv_lcd_summary}},
#'   \code{\link{eval_lcd}}, \code{\link[flowmix]{make_cv_folds}}
#'
#' @examples
#' \dontrun{
#' # 25 simulations, 5x5 lambda grid, 10 repeats, 5 folds (31,250 jobs total):
#' sims <- sprintf("sim_data_%02d.rds", 1:25)
#'
#' res <- cv_lcd_simul(
#'   sim_files     = sims,
#'   K             = 2,
#'   alpha_lambdas = 10^seq(-4, -2, length.out = 5),
#'   theta_lambdas = 10^seq(-4, -2, length.out = 5),
#'   nfold         = 5,
#'   seeds         = 1:10,
#'   save_dir      = "cv_runs",
#'   n_cores       = 64,          # USC CARC cap
#' )
#'
#' # Summarize each simulation (stored under cv_runs/sim_<s>/):
#' for (s in seq_along(sims)) {
#'   sim_dir <- file.path("cv_runs", sprintf("sim_%d", s))
#'   load(file.path(sim_dir, "index_matrix.Rdata"))  # loads index_matrix
#'   summ <- cv_lcd_summary(index_matrix, save_dir = sim_dir)
#'   print(list(sim = s, opt = summ$opt_lambdas, cv_score = max(summ$reduced_mat[, "cv_score"], na.rm = TRUE)))
#' }
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
  iter_eta       = 1e-3,
  resp_threshold = 1e-3,
  trim_prob      = 0.03,
  blocksize      = 10,
  base_dir       = "./cv_saves",
  n_cores        = "max",
  lp_time_limit  = 600
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
                "resp_threshold", "trim_prob", "blocksize", "base_dir", "nfold", "lp_time_limit"),
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
        lp_time_limit  = lp_time_limit
      )
      return(res_ii)
  })

  success   <- unlist(res, use.names = FALSE)
  summary   <- sprintf("Failures: %d/%d (%.1f%%)", sum(!success), length(success), 100 * sum(!success)/length(success))
  
  return(list(grand_jobs = grand_jobs,
              summary = summary))
}
