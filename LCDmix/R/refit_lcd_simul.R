# Generated from create-LCDmix.Rmd: do not edit by hand

#' Refit LCDmix across many simulations and record best files (no rereads)
#'
#' For each simulation, runs multiple refits across the given seeds using a
#' single outer parallel cluster (no nested parallel). Each refit is cached to
#' \code{base_dir/sim_<s>/refit/refit_<seed>.rds}. The function returns:
#' (i) a per-job table with refit metrics, (ii) logs, and (iii) a per-simulation
#' \emph{best table} giving the simulation index, its \code{lambda_alpha} and
#' \code{lambda_theta}, and the filename of the best refit chosen by the largest
#' training objective \eqn{L}. The best table is built directly from in-memory
#' results—no disk rereads.
#'
#' @param sim_files Character vector; paths to simulation files. Each file must
#'   contain \code{Y_bin}, \code{X}, and \code{bin_mass}. Files may be \code{.rds}
#'   (containing a list) or \code{.RData}.
#' @param opt_lambdas_list List of length \code{length(sim_files)}; each element
#'   is \code{c(lambda_alpha, lambda_theta)} to use for that simulation.
#' @param K Integer; number of mixture components.
#' @param seeds Integer vector of seeds to run. If \code{NULL}, supply \code{cv_reps}.
#' @param cv_reps Integer; number of repeats used only when \code{seeds} is
#'   \code{NULL} (seeds become \code{1:cv_reps}). Default \code{NULL}.
#' @param max_iter Integer; maximum EM iterations per refit. Default \code{30}.
#' @param iter_eta Numeric; convergence tolerance for the surrogate objective.
#'   Default \code{1e-3}.
#' @param resp_threshold Numeric in \eqn{[0,1]}; responsibilities below this are
#'   zeroed for stability. Default \code{1e-3}.
#' @param trim_prob Numeric in \eqn{[0,1)}; trimming fraction used during fitting.
#'   Default \code{0.01}.
#' @param base_dir Character; base directory for outputs. Refit caches are written
#'   under \code{base_dir/sim_<s>/refit/}. Default \code{"./cv_saves"}.
#' @param n_cores Integer or \code{"max"}; number of workers for the single outer
#'   cluster. Default \code{"max"}.
#' @param debug Logical; forwarded to \code{main()} for verbose diagnostics.
#'   Default \code{FALSE}.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{jobs}}{Data frame of per-job metrics with columns
#'         \code{sim_idx}, \code{seed}, \code{lambda_alpha}, \code{lambda_theta},
#'         \code{L}, and \code{file}.}
#'   \item{\code{logs}}{Character vector of chunk logs.}
#'   \item{\code{best_table}}{Data frame with columns
#'         \code{sim}, \code{lambda_alpha}, \code{lambda_theta}, \code{best_file};
#'         one row per simulation, where \code{best_file} is the cached filename
#'         of the refit with the largest \eqn{L} (or \code{NA} if none succeeded).}
#' }
#'
#' @seealso \code{\link{refit_onejob}}, \code{\link{cv_lcd_simul}}
#'
#' @examples
#' \dontrun{
#' sims <- sprintf("sim_%02d.rds", 1:10)
#' # suppose each sim has its own chosen penalties from CV:
#' opt_list <- replicate(10, c(1e-3, 1e-3), simplify = FALSE)
#' out <- refit_lcd_simul(
#'   sim_files        = sims,
#'   opt_lambdas_list = opt_list,
#'   K                = 2,
#'   seeds            = 1:8,
#'   base_dir         = "cv_saves",
#'   n_cores          = 32
#' )
#' out$best_table
#' }
#'
#' @export
refit_lcd_simul <- function(
  sim_files,
  opt_lambdas_list,
  K,
  seeds = NULL, 
  cv_reps = NULL,
  max_iter = 30, 
  iter_eta = 1e-3, 
  resp_threshold = 1e-3, 
  trim_prob = 0.01,
  base_dir = "./cv_saves",
  n_cores = "max"
) {
  if (is.null(seeds) && is.null(cv_reps)) stop("`seeds` or `cv_reps` required")
  if (is.null(seeds)) seeds <- seq_len(cv_reps)
  if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE)

  num_sims <- length(sim_files)
  if (length(opt_lambdas_list) != num_sims)
    stop("`opt_lambdas_list` length must equal `length(sim_files)`.")

  # Ensure per-simulation refit dirs exist
  for (s in seq_len(num_sims)) {
    sim_refit_dir <- file.path(base_dir, sprintf("sim_%d", s), "refit")
    if (!dir.exists(sim_refit_dir)) dir.create(sim_refit_dir, recursive = TRUE)
  }

  # Grand job list: one row per (sim, seed) with that sim's lambdas
  grand_jobs <- do.call(rbind, lapply(seq_len(num_sims), function(s) {
    data.frame(
      sim_idx      = s,
      seed         = seeds,
      lambda_alpha = opt_lambdas_list[[s]][1],
      lambda_theta = opt_lambdas_list[[s]][2]
    )
  }))
  rownames(grand_jobs) <- NULL

  # Cluster
  n_workers <- if (identical(n_cores, "max")) parallel::detectCores(logical = FALSE) else as.integer(n_cores)
  cl <- parallel::makeCluster(n_workers)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

  parallel::clusterEvalQ(cl, { library(LCDmix); NULL })
  parallel::clusterExport(
    cl,
    varlist = c("sim_files","grand_jobs","K","max_iter","iter_eta",
                "resp_threshold","trim_prob","base_dir"),
    envir = environment()
  )

  # Worker: run one refit job
  res <- parallel::parLapply(cl, seq_len(nrow(grand_jobs)), function(ii) {
    row <- grand_jobs[ii, ]
    s   <- as.integer(row[["sim_idx"]])
    sd  <- as.integer(row[["seed"]])
    la  <- as.numeric(row[["lambda_alpha"]])
    lt  <- as.numeric(row[["lambda_theta"]])
    
    sim_refit_dir <- file.path(base_dir, sprintf("sim_%d", s), "refit")
    out_path <- file.path(sim_refit_dir, sprintf("refit_%d.rds", sd))
    
    # cached?
    if (file.exists(out_path)) {
      obj <- readRDS(out_path)
      return(!is.null(obj$fit))
    }
    
    # load this simulation
    sim <- readRDS(sim_files[s])
    Y_bin    <- sim$Y_bin
    X        <- sim$X
    bin_mass <- sim$bin_mass
    
    # run one refit (writes cache under sim_refit_dir)
    res_ii <- refit_onejob(
      Y_bin = Y_bin, X = X, bin_mass = bin_mass, K = K,
      lambda_alpha = la, lambda_theta = lt,
      seed = sd, max_iter = max_iter, iter_eta = iter_eta,
      resp_threshold = resp_threshold, trim_prob = trim_prob,
      save_dir = sim_refit_dir)
    return(res_ii)
    })
  
  success   <- unlist(res, use.names = FALSE)
  summary   <- sprintf("Failures: %d/%d (%.1f%%)", sum(!success), length(success), 100 * sum(!success)/length(success))
  
  best_table <- data.frame(
    sim            = seq_len(num_sims),
    lambda_alpha   = vapply(opt_lambdas_list, function(x) as.numeric(x[1]), numeric(1)),
    lambda_theta   = vapply(opt_lambdas_list, function(x) as.numeric(x[2]), numeric(1)),
    best_idx       = NA_integer_,
    stringsAsFactors = FALSE
  )

  for (s in seq_len(num_sims)) {
    sim_refit_dir <- file.path(base_dir, sprintf("sim_%d", s), "refit")
    L_vec     <- rep(NA_real_, length(seeds))
    
    for (i in seq_along(seeds)) {
      sd <- as.integer(seeds[i])
      file_name  <- file.path(sim_refit_dir, sprintf("refit_%d.rds", sd))
      if (!file.exists(file_name)) next
      obj <- readRDS(file_name)
      if (!is.null(obj$fit)) {
        L_vec[i]  <- obj$fit$L$trimmed_loglik
      }
    }
    
    if (all(!is.finite(L_vec))) next
    best_table$best_idx[s] <- seeds[which.max(L_vec)]
  }

  return(list(grand_jobs = grand_jobs,
              summary    = summary,
              best_table = best_table))
}
