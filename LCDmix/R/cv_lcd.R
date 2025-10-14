# Generated from create-LCDmix.Rmd: do not edit by hand

#' @export
cv_lcd <- function(
  Y_bin,
  X,
  bin_mass,
  K,
  alpha_lambdas      = c(1e-4, 1e-3),
  theta_lambdas      = c(1e-4, 1e-3),
  max_iter           = 30,
  iter_eta           = 1e-3,
  resp_threshold     = 1e-3,
  nfold              = 5,
  seeds              = NULL,
  trim_prob          = 0.01,
  save_dir           = "./cv_saves",
  n_cores            = "max",
  cv_reps            = NULL,
  blocksize          = 20
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
                "trim_prob","save_dir","folds","index_matrix"),
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
        job       = job,
        Y_bin     = Y_bin,
        X         = X,
        bin_mass  = bin_mass,
        folds     = folds,
        K         = K,
        max_iter  = max_iter,
        iter_eta  = iter_eta,
        resp_threshold = resp_threshold,
        trim_prob = trim_prob,
        save_dir  = save_dir
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
