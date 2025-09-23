# Generated from create-LCDmix.Rmd: do not edit by hand

#' Refit LCDmix per simulation (with best lambdas) and score via \code{mixture_metric}
#'
#' @description
#' For each simulation in \code{best_table}, this function refits the LCD
#' mixture-of-experts model using the provided optimal penalty pair
#' \code{(lambda_alpha, lambda_theta)} and multiple random seeds. Each fit is
#' evaluated with \code{mixture_metric()} against the simulation ground truth,
#' and a per-run result is written to \code{save_dir} as
#' \code{refit_sim<sim_idx>_seed<seed>.rds}. Results are computed in parallel.
#'
#' @param best_table A two- or three-column matrix/data.frame with at least
#'   columns \code{sim_idx}, \code{lambda_alpha}, \code{lambda_theta}. Typically
#'   the output of \code{summarize_simul_cv()}.
#' @param sim_dir Directory containing simulation files named
#'   \code{"sim_<sim_idx>.rds"} produced by \code{simulate_and_save()}.
#' @param save_dir Directory to write per-refit result files
#'   \code{"refit_sim<sim_idx>_seed<seed>.rds"}.
#' @param K Integer; number of mixture components.
#' @param max_iter Integer; maximum EM iterations in \code{main()}. Default \code{30}.
#' @param iter_eta Numeric; convergence tolerance (relative change in Q). Default \code{1e-3}.
#' @param resp_threshold Numeric in \[0,1\]; responsibilities below this are zeroed
#'   for numerical stability. Default \code{1e-3}.
#' @param refit_seeds Integer vector of RNG seeds for repeated refits per
#'   simulation. Default \code{1:1}.
#' @param n_cores Integer or \code{"max"}; number of parallel workers. \code{"max"}
#'   uses all physical cores minus one. Default \code{"max"}.
#'
#' @details
#' Each job:
#' \enumerate{
#'   \item Loads \code{sim_<sim_idx>.rds} (which contains \code{ylist}, \code{X},
#'         \code{countslist}, \code{prob}, and \code{dens_true}).
#'   \item Calls \code{main()} with the best \code{lambda_alpha/theta} and the
#'         given seed.
#'   \item Computes the weighted L1 mixture distance via \code{mixture_metric()}.
#'   \item Saves an \code{.rds} with \code{metric}, \code{L}, \code{fit_try},
#'         and log/error text.
#' }
#'
#' @return A matrix (one row per \code{sim_idx} × \code{seed_idx}) with columns:
#' \describe{
#'   \item{\code{sim_idx}}{Simulation ID.}
#'   \item{\code{lambda_alpha}}{Penalty used for mixture weights.}
#'   \item{\code{lambda_theta}}{Penalty used for component parameters.}
#'   \item{\code{seed_idx}}{Refit seed.}
#'   \item{\code{L}}{Final loglikelihood objective value (NA on failure).}
#'   \item{\code{metric}}{Weighted L1 distance to truth from \code{mixture_metric()} (NA on failure).}
#' }
#'
#' @seealso \code{\link{simulate_and_save}}, \code{\link{cv_lcd_simul}},
#'   \code{\link{summarize_simul_cv}}, \code{\link{mixture_metric}}
#'
#' @examples
#' \dontrun{
#' # Suppose 'best_tbl' came from summarize_simul_cv()
#' out <- simul_refit(
#'   best_table  = best_tbl,
#'   sim_dir     = "sim_data",
#'   save_dir    = "cv_saves",
#'   K           = 2,
#'   refit_seeds = 1:3,
#'   n_cores     = "max"
#' )
#' head(out)
#' }
#'
#' @export
simul_refit <- function(
  sim_files,
  opt_lambdas_list,      # list length = length(sim_files), each c(lambda_alpha, lambda_theta)
  K,
  seeds = NULL, cv_reps = NULL,
  max_iter = 30, iter_eta = 1e-3, resp_threshold = 1e-3, trim_prob = 0.01,
  base_dir = "./cv_saves",
  n_cores = "max",
  chunk_size = 20,
  debug = FALSE
) {
  if (is.null(seeds) && is.null(cv_reps)) stop("`seeds` or `cv_reps` required")
  if (is.null(seeds)) seeds <- seq_len(cv_reps)
  if (!dir.exists(base_dir)) dir.create(base_dir, recursive = TRUE)

  num_sims <- length(sim_files)
  if (length(opt_lambdas_list) != num_sims)
    stop("opt_lambdas_list must have length equal to length(sim_files).")

  # Load sims
  sims <- vector("list", num_sims)
  for (s in seq_len(num_sims)) {
    env <- new.env(parent = emptyenv())
    f <- sim_files[s]
    if (grepl("\\.rds$", f, ignore.case = TRUE)) {
      obj <- readRDS(f); if (is.list(obj)) base::list2env(obj, env) else stop("RDS must be a list with Y_bin, X, bin_mass.")
    } else {
      load(f, envir = env)
    }
    sims[[s]] <- env
    sim_refit_dir <- file.path(base_dir, sprintf("sim_%d", s), "refit")
    if (!dir.exists(sim_refit_dir)) dir.create(sim_refit_dir, recursive = TRUE)
  }

  # Build grand job index over (sim, seed)
  grand_jobs <- do.call(rbind, lapply(seq_len(num_sims), function(s) {
    data.frame(sim_idx = s, seed = seeds, lambda_alpha = opt_lambdas_list[[s]][1],
               lambda_theta = opt_lambdas_list[[s]][2])
  }))
  rownames(grand_jobs) <- NULL

  # Chunk rows
  rows <- seq_len(nrow(grand_jobs))
  chunks <- split(rows, ceiling(seq_along(rows) / chunk_size))

  # Cluster
  n_workers <- if (identical(n_cores, "max")) parallel::detectCores(logical = FALSE) else as.integer(n_cores)
  cl <- parallel::makeCluster(n_workers)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

  parallel::clusterEvalQ(cl, { library(LCDmix); NULL })
  parallel::clusterExport(
    cl,
    varlist = c("sims","grand_jobs","K","max_iter","iter_eta","resp_threshold",
                "trim_prob","base_dir","debug"),
    envir = environment()
  )

  logs <- parallel::parLapply(cl, chunks, function(rr) {
    out <- character(length(rr))
    for (i in seq_along(rr)) {
      gi <- rr[i]
      row <- grand_jobs[gi, ]
      s <- row[["sim_idx"]]
      sd <- row[["seed"]]
      la <- row[["lambda_alpha"]]
      lt <- row[["lambda_theta"]]

      env <- sims[[s]]
      sim_refit_dir <- file.path(base_dir, sprintf("sim_%d", s), "refit")

      res <- refit_onejob(
        Y_bin = env$Y_bin, X = env$X, bin_mass = env$bin_mass, K = K,
        lambda_alpha = la, lambda_theta = lt,
        seed = sd, max_iter = max_iter, iter_eta = iter_eta,
        resp_threshold = resp_threshold, trim_prob = trim_prob,
        save_dir = sim_refit_dir, debug = debug
      )
      out[i] <- paste0("[sim ", s, "] ", res$log_msg)
    }
    paste(out, collapse = "\n")
  })

  return(list(
    jobs = grand_jobs,
    logs = unlist(logs, use.names = FALSE)
  ))
}
