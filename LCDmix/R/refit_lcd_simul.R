# Generated from create-LCDmix.Rmd: do not edit by hand

#' Refit the model with several seeds for many simulated datasets
#'
#' @description
#' Runs \code{refit_onejob()} for every combination of dataset and seed in one
#' parallel cluster (no nested parallelism), each dataset with its own
#' penalty pair. Refit \code{seed} of dataset \eqn{s} is saved to
#' \code{base_dir/sim_<s>/refit/refit_<seed>.rds}; a file that already exists
#' is not computed again. After the run, the saved files are read and, for
#' each dataset, the seed with the largest \code{fit_L} (the penalized
#' training log-likelihood) is recorded.
#'
#' The workers load the installed LCDmix package, so install the current
#' version before running.
#'
#' @param sim_files Paths of \code{.rds} files, one per dataset; each holds a
#'   list with \code{Y_bin}, \code{X} and \code{bin_mass}.
#' @param opt_lambdas_list A list of the same length as \code{sim_files};
#'   element \eqn{s} is \code{c(lambda_alpha, lambda_theta)} for dataset
#'   \eqn{s}.
#' @param base_dir Output directory, usually the \code{base_dir} of
#'   \code{cv_lcd_simul()}.
#' @inheritParams refit_lcd
#' @inheritParams main
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{grand_jobs}}{Data frame with one row per refit:
#'     \code{sim_idx}, \code{seed}, \code{lambda_alpha}, \code{lambda_theta}.}
#'   \item{\code{summary}}{A text with the number of failed refits.}
#'   \item{\code{loglik_lst}}{List with one numeric vector per dataset:
#'     \code{fit_L} of each seed (\code{NA} for a failed seed).}
#'   \item{\code{best_table}}{Data frame with one row per dataset: \code{sim},
#'     \code{lambda_alpha}, \code{lambda_theta}, \code{best_idx} (seed of the
#'     best refit), \code{loglik} (its \code{fit_L}), \code{best_file} (path of
#'     its file) and \code{n_failed_seeds}. \code{best_idx}, \code{loglik}
#'     and \code{best_file} are \code{NA}, with a warning, for a dataset whose
#'     refits all failed.}
#' }
#'
#' @seealso \code{\link{refit_onejob}}, \code{\link{cv_lcd_simul}}
#'
#' @examples
#' \dontrun{
#' sims <- file.path("sim_data", sprintf("sim_%d.rds", 1:10))
#' summ <- cv_lcd_summary_simul(base_dir = "cv_saves", num_sims = 10)
#' out  <- refit_lcd_simul(
#'   sim_files        = sims,
#'   opt_lambdas_list = lapply(summ, function(x) x$opt_lambdas),
#'   K                = 2,
#'   seeds            = 1:10,
#'   base_dir         = "cv_saves",
#'   n_cores          = 32
#' )
#' out$best_table
#' best_fit_1 <- readRDS(out$best_table$best_file[1])$fit
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
  iter_eta = 1e-4, 
  resp_threshold = 1e-3, 
  trim_prob = 0.03,
  base_dir = "./cv_saves",
  n_cores = "max",
  lp_time_limit = 600,
  calc_Q_every = FALSE,
  update = c("lp", "optim"),    # NEW
  maxdev = NULL                 # NEW (fixP) - appended last
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
    varlist = c("sim_files","grand_jobs","K","max_iter","iter_eta","resp_threshold",
                "trim_prob","base_dir", "lp_time_limit", "update", "calc_Q_every",
                "maxdev"),   # NEW (fixP)
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
      save_dir = sim_refit_dir, lp_time_limit = lp_time_limit, 
      update = update, calc_Q_every = calc_Q_every,
      maxdev = maxdev)     # NEW (fixP)
    return(res_ii)
    })
  
  success   <- unlist(res, use.names = FALSE)
  summary   <- sprintf("Failures: %d/%d (%.1f%%)", sum(!success), length(success), 100 * sum(!success)/length(success))
  
   best_table <- data.frame(
     sim            = seq_len(num_sims),
     lambda_alpha   = vapply(opt_lambdas_list, function(x) as.numeric(x[1]), numeric(1)),
     lambda_theta   = vapply(opt_lambdas_list, function(x) as.numeric(x[2]), numeric(1)),
     best_idx       = NA_integer_,
     loglike        = NA_real_,
     best_file      = NA_character_,   # NEW: resolved path, or NA if none
     n_failed_seeds = NA_integer_,     # NEW: how many of `seeds` produced no fit
     stringsAsFactors = FALSE
   )
  loglik_lst = list()

  for (s in seq_len(num_sims)) {
    sim_refit_dir <- file.path(base_dir, sprintf("sim_%d", s), "refit")
    L_vec     <- rep(NA_real_, length(seeds))
    
    for (i in seq_along(seeds)) {
      sd <- as.integer(seeds[i])
      file_name  <- file.path(sim_refit_dir, sprintf("refit_%d.rds", sd))
      if (!file.exists(file_name)) next
      obj <- readRDS(file_name)
      if (!is.null(obj$fit_L)) {
        L_vec[i]  <- obj$fit_L
      }
    }
     loglik_lst[[s]] = L_vec
     best_table$n_failed_seeds[s] <- sum(!is.finite(L_vec))       # NEW
     if (all(!is.finite(L_vec))) {
       warning("refit_lcd_simul(): simulation ", s,
               " has no usable refit (all ", length(seeds)," seeds failed)")  # NEW
       next
     }
     L_vec[!is.finite(L_vec)] <- -Inf
     bi <- seeds[which.max(L_vec)]
     best_table$best_idx[s]  <- bi
     best_table$loglike[s]   <- max(L_vec)
     best_table$best_file[s] <- file.path(sim_refit_dir,             # NEW
                                          sprintf("refit_%d.rds", as.integer(bi)))
  }
  colnames(best_table)[5] = 'loglik'
  
   n_dead <- sum(is.na(best_table$best_idx))
   if (n_dead > 0)
     warning("refit_lcd_simul(): ", n_dead, " of ", num_sims,
             " simulations have no usable refit: ",
             paste(which(is.na(best_table$best_idx)), collapse = ", "))
   
  return(list(grand_jobs = grand_jobs,
              summary    = summary,
              loglik_lst = loglik_lst,
              best_table = best_table))
}
