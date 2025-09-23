# Generated from create-LCDmix.Rmd: do not edit by hand

#' Refit LCDmix multiple times and select the best training fit
#' @export
refit_lcd <- function(
  Y_bin, 
  X, 
  bin_mass, 
  K, 
  opt_lambdas,
  seeds = NULL, 
  cv_reps = NULL,
  max_iter = 30, 
  iter_eta = 1e-3, 
  resp_threshold = 1e-3, 
  trim_prob = 0.01,
  save_dir = "./refits", 
  n_cores = "max", 
  debug = FALSE
) {
  if (is.null(seeds) && is.null(cv_reps)) stop("`seeds` or `cv_reps` required")
  if (is.null(seeds)) seeds <- seq_len(cv_reps)
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  lambda_alpha <- opt_lambdas[1]
  lambda_theta <- opt_lambdas[2]

  n_workers <- if (identical(n_cores, "max")) parallel::detectCores(logical = FALSE) else as.integer(n_cores)
  cl <- parallel::makeCluster(n_workers)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

  parallel::clusterEvalQ(cl, {library(LCDmix); NULL })
  parallel::clusterExport(
    cl,
    varlist = c("Y_bin","X","bin_mass","K","lambda_alpha","lambda_theta",
                "max_iter","iter_eta","resp_threshold","trim_prob","save_dir","debug"),
    envir = environment()
  )

  res_list <- parallel::parLapply(cl, seeds, function(sd) {
    refit_onejob(
      Y_bin = Y_bin, X = X, bin_mass = bin_mass, K = K,
      lambda_alpha = lambda_alpha, lambda_theta = lambda_theta,
      seed = sd, max_iter = max_iter, iter_eta = iter_eta,
      resp_threshold = resp_threshold, trim_prob = trim_prob,
      save_dir = save_dir, debug = debug
    )
  })

  logs <- vapply(res_list, `[[`, character(1), "log_msg")
  scores <- vapply(res_list, function(x) x$L, numeric(1))

  if (all(is.na(scores))) {
    logs <- c(logs, "All refits failed.")
    return(list(logs = logs, refit_scores = scores, best_fit = NULL))
  }

  best_idx <- which.max(scores)
  logs <- c(logs, paste0(
    "Completed ", length(seeds), " repeats; ",
    sum(is.na(scores)), " failures; best L = ", round(scores[best_idx], 6),
    " (seed index ", best_idx, ").\n"
  ))

  return(list(logs = logs, refit_scores = scores, best_fit = res_list[[best_idx]]))
}
