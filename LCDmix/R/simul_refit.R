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
#'   \item Saves an \code{.rds} with \code{metric}, \code{Q_final}, \code{fit_try},
#'         and log/error text.
#' }
#'
#' @return A matrix (one row per \code{sim_idx} × \code{seed_idx}) with columns:
#' \describe{
#'   \item{\code{sim_idx}}{Simulation ID.}
#'   \item{\code{lambda_alpha}}{Penalty used for mixture weights.}
#'   \item{\code{lambda_theta}}{Penalty used for component parameters.}
#'   \item{\code{seed_idx}}{Refit seed.}
#'   \item{\code{Q_final}}{Final surrogate objective value (NA on failure).}
#'   \item{\code{metric}}{Weighted L1 distance to truth from \code{mixture_metric()} (NA on failure).}
#' }
#'
#' @seealso \code{\link{simulate_and_save}}, \code{\link{run_cv_simul}},
#'   \code{\link{summarize_simul_cv}}, \code{\link{mixture_metric}}
#'
#' @examples
#' \dontrun{
#' # Suppose 'best_tbl' came from summarize_simul_cv()
#' out <- simul_refit(
#'   best_table  = best_tbl,
#'   sim_dir     = "sim_data",
#'   save_dir    = "cv_results",
#'   K           = 2,
#'   refit_seeds = 1:3,
#'   n_cores     = "max"
#' )
#' head(out)
#' }
#'
#' @export
simul_refit <- function(
    best_table, 
    sim_dir = "sim_data", 
    save_dir = "cv_results",
    K,
    max_iter           = 30,
    iter_eta           = 1e-3,
    resp_threshold     = 1e-3,
    refit_seeds        = 1:1,
    n_cores            = "max"
    ) {
  
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)
  
  # Use expand.grid to get every combination
  idx <- expand.grid(
    sim_idx        = best_table[, "sim_idx"],
    seed_idx       = refit_seeds
  )
  # join in the penalty columns
  best_table <- merge(
    idx,
    as.data.frame(best_table)[, c("sim_idx","lambda_alpha","lambda_theta")],
    by = "sim_idx",
    sort = FALSE
  )
  # ensure columns in desired order
  best_table <- as.matrix(best_table[, c("sim_idx","lambda_alpha","lambda_theta","seed_idx")])
  
  # determine worker count
  if (identical(n_cores, "max")) {
    n_workers <- parallel::detectCores(logical = FALSE) - 1
  } else {
    n_workers <- as.integer(n_cores)
  }
  cl <- parallel::makeCluster(n_workers)
  # export functions and data that workers need
  parallel::clusterEvalQ(cl, {library(flowmix); library(LCDmix); NULL })
  parallel::clusterExport(
    cl,
    varlist = c("best_table", "sim_dir", "save_dir",
                "K", "max_iter", "iter_eta", "resp_threshold", 
                "refit_seeds"),
    envir   = environment()
  )
  
  res_i <- parallel::parLapply(cl, seq_len(nrow(best_table)), function(i) {
    idx_row       <- best_table[i, ]
    sim_idx       <- as.integer(idx_row["sim_idx"])
    seed_idx      <- as.integer(idx_row["seed_idx"])
    lambda_alpha  <- as.numeric(idx_row["lambda_alpha"])
    lambda_theta  <- as.numeric(idx_row["lambda_theta"])

    sim     <- readRDS(file.path(sim_dir, paste0("sim_", sim_idx, ".rds")))
    log_msg <- paste0(
        "sim_idx=", sim_idx,
        ", alpha=", lambda_alpha,
        ", theta=", lambda_theta, 
        ", seed=",  seed_idx, "\n"
      )
    
    # Skip if already exists
    save_path <- file.path(
        save_dir,
        sprintf("refit-%d-%d.rds", sim_idx, seed_idx)
      )
    
    if (file.exists(save_path)) {
      print(paste0("Skip (cached): ", basename(save_path)))
      saved = readRDS(save_path)
      return(c(
        sim_idx      = sim_idx,
        lambda_alpha = lambda_alpha,
        lambda_theta = lambda_theta,
        seed_idx     = seed_idx,
        Q_final      = as.numeric(saved$Q_final),
        metric       = as.numeric(saved$metric)
      ))
    }
    
    set.seed(seed_idx)
    err_msg <- NULL
    out_log <- capture.output(
      fit_try <- tryCatch({
        main(
          Y               = sim$ylist,
          X               = sim$X,
          biomass         = sim$countslist,
          binned          = TRUE,
          n_bins          = 0,
          K               = K,
          lambda_alpha    = lambda_alpha,
          lambda_theta    = lambda_theta,
          max_iter        = max_iter,
          iter_eta        = iter_eta,
          resp_threshold  = resp_threshold,
          debug           = TRUE
        )
      }, error = function(e) {
        err_msg <<- paste0("Error fitting sim=", sim_idx,
                            ", alpha=", lambda_alpha,
                            ", theta=", lambda_theta, 
                             ", seed=", seed_idx, ": ",
                          e$message)
        return(NULL)
      }),
        type = "output"      #<<— capture only stdout (print/cat), not messages
      )

      # append the captured log to your log_msg
      log_msg <- paste0(log_msg,
                        paste0(out_log, collapse = "\n"),
                        "\n")
    
      if (is.null(fit_try) || !is.list(fit_try$iter)) {
        log_msg    <- paste0(log_msg, "✖ Error: ", err_msg, "\n")
        Q_final    <- NA
        metric     <- NA
      } else {
        log_msg    <- paste0(log_msg, "Fit succeeded\n")
        Q_final    <- tail(fit_try$iter$Q, 1)
        metric     <- mixture_metric(
                        sim$ylist, sim$X, sim$countslist,
                        est_res  = fit_try$iter,
                        true_res = list(prob=sim$prob, dens_true=sim$dens_true)
                      )$weighted
      }
    
    # save RDS
    saveRDS(
      list(metric  = metric,
           Q_final = Q_final,
           fit_try = fit_try,
           error   = err_msg,
           logs    = log_msg),
      file = file.path(
        save_dir,
        sprintf("refit-%d-%d.rds", sim_idx, seed_idx)
      )
    )
    
    # return summary
    return( c(sim_idx        = sim_idx,
              lambda_alpha   = lambda_alpha,
              lambda_theta   = lambda_theta,
              seed_idx       = seed_idx,
              Q_final        = Q_final,
              metric         = metric))
  })
  parallel::stopCluster(cl)
  
   #— Summarize failures —#
  res = vapply(res_i, function(x) {is.na(as.numeric(x["Q_final"]))}, logical(1))

  res_logs   <- sprintf("Failures: %d/%d (%.1f%%)",
            sum(res), length(res), 100 * sum(res) / length(res))
  
  # bind into a data.frame
  return(list(res_logs = res_logs,
              res_table = as.data.frame(do.call(rbind, res_i)))
        )
}
