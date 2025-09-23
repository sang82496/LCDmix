# Generated from create-LCDmix.Rmd: do not edit by hand

#' Refit LCDmix across many simulations (outer-parallel)
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
