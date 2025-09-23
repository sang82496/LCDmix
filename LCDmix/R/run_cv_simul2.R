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
#'   \code{evaluate_lcd_model()}. Default \code{0.01}.
#' @param blocksize Integer; block size for \code{flowmix::make_cv_folds()}.
#'   Default \code{20}.
#' @param save_dir Character; base directory for outputs. Per-simulation results
#'   are written to \code{save_dir/sim_<s>/}. Default \code{"./result"}.
#' @param n_cores Integer or \code{"max"}; number of workers for the single
#'   outer cluster. \code{"max"} uses all physical cores minus one. Default \code{"max"}.
#' @param chunk_size Integer; number of CV jobs each worker processes serially per
#'   task to reduce overhead. Default \code{25}.
#'
#' @details
#' This function uses a single outer cluster to process a \emph{grand} index of
#' all simulation jobs in chunks, avoiding nested parallel. For reproducibility,
#' job seeds are offset by the simulation index (\code{seed_idx + 1000 * sim_idx})
#' before fitting. Each per-job RDS file contains:
#' \code{prop_inf}, \code{trimmed_loglik}, \code{finite_loglik}, \code{L},
#' and \code{log_msg}. After completion, run \code{\link{LCD_cv_summary}} on each
#' simulation's subdirectory to obtain per-sim CV scores and optimal penalties.
#'
#' @return A list with:
#' \describe{
#'   \item{\code{grand_index}}{Data frame of all jobs with a \code{sim_idx} column.}
#'   \item{\code{logs}}{Character vector of per-chunk log strings.}
#' }
#'
#' @seealso \code{\link{cv_lcd_parallel2}}, \code{\link{LCD_cv_summary}},
#'   \code{\link{evaluate_lcd_model}}, \code{\link[flowmix]{make_cv_folds}}
#'
#' @examples
#' \dontrun{
#' # 25 simulations, 5x5 lambda grid, 10 repeats, 5 folds (31,250 jobs total):
#' sims <- sprintf("sim_data_%02d.rds", 1:25)
#'
#' res <- run_cv_simul2(
#'   sim_files     = sims,
#'   K             = 2,
#'   alpha_lambdas = 10^seq(-4, -2, length.out = 5),
#'   theta_lambdas = 10^seq(-4, -2, length.out = 5),
#'   nfold         = 5,
#'   seeds         = 1:10,
#'   save_dir      = "cv_runs",
#'   n_cores       = 64,          # USC CARC cap
#'   chunk_size    = 25           # tune 10–50 based on runtime
#' )
#'
#' # Summarize each simulation (stored under cv_runs/sim_<s>/):
#' for (s in seq_along(sims)) {
#'   sim_dir <- file.path("cv_runs", sprintf("sim_%d", s))
#'   load(file.path(sim_dir, "index_matrix.Rdata"))  # loads index_matrix
#'   summ <- LCD_cv_summary(index_matrix, save_dir = sim_dir)
#'   print(list(sim = s, opt = summ$opt_lambdas, cv_score = max(summ$reduced_mat[, "cv_score"], na.rm = TRUE)))
#' }
#' }
#'
#' @export
run_cv_simul2 <- function(
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
  trim_prob      = 0.01,
  blocksize      = 20,
  save_dir       = "./result",
  n_cores        = "max",
  chunk_size     = 25
) {
  if (is.null(seeds) && is.null(cv_reps)) stop("`seeds` and `cv_reps` cannot be both NULL")
  if (is.null(seeds)) seeds <- seq_len(cv_reps)
  if (!dir.exists(save_dir)) dir.create(save_dir, recursive = TRUE)

  S <- length(sim_files)
  sims <- vector("list", S)

  # Load sims + build per-sim folds and index matrices
  for (s in seq_len(S)) {
    env <- new.env(parent = emptyenv())
    f <- sim_files[s]
    if (grepl("\\.rds$", f, ignore.case = TRUE)) {
      obj <- readRDS(f); if (is.list(obj)) list2env(obj, env) else stop("RDS must be a list with Y_bin, X, bin_mass.")
    } else {
      load(f, envir = env)
    }
    env$folds <- flowmix::make_cv_folds(ylist = env$Y_bin, nfold = nfold, blocksize = blocksize)
    env$index_matrix <- make_cv_index_matrix(
      nfold          = nfold,
      seeds          = seeds,
      alpha_lambdas  = sort(alpha_lambdas),
      theta_lambdas  = sort(theta_lambdas)
    )
    # Persist per-sim index for later summary
    sim_dir <- file.path(save_dir, sprintf("sim_%d", s))
    if (!dir.exists(sim_dir)) dir.create(sim_dir, recursive = TRUE)
    saveRDS(env$index_matrix, file = file.path(sim_dir, "index_matrix.rds"))
    sims[[s]] <- env
  }

  # Grand index with sim_idx + seed offset to decorrelate RNG across sims
  grand_list <- lapply(seq_len(S), function(s) {
    im <- sims[[s]]$index_matrix
    cbind(sim_idx = s, im)
  })
  grand_index <- do.call(rbind, grand_list)
  rownames(grand_index) <- NULL
  grand_index <- as.data.frame(grand_index)

  # Chunk rows
  rows   <- seq_len(nrow(grand_index))
  chunks <- split(rows, ceiling(seq_along(rows) / chunk_size))

  # One outer cluster
  n_workers <- if (identical(n_cores, "max")) parallel::detectCores(logical = FALSE) else as.integer(n_cores)
  cl <- parallel::makeCluster(n_workers)
  on.exit(try(parallel::stopCluster(cl), silent = TRUE), add = TRUE)

  parallel::clusterEvalQ(cl, { library(flowmix); library(LCDmix); NULL })
  parallel::clusterExport(
    cl,
    varlist = c("sims","grand_index","K","max_iter","iter_eta","resp_threshold",
                "trim_prob","save_dir"),
    envir = environment()
  )

  logs <- parallel::parLapply(cl, chunks, function(chunk_rows) {
    out <- character(length(chunk_rows))
    for (i in seq_along(chunk_rows)) {
      gi  <- chunk_rows[i]
      row <- grand_index[gi, ]
      s   <- row[["sim_idx"]]

      # Build the 'job' vector for the worker
      job <- c(
        alpha_idx    = as.numeric(row[["alpha_idx"]]),
        theta_idx    = as.numeric(row[["theta_idx"]]),
        seed_idx     = as.numeric(row[["seed_idx"]]),
        fold_idx     = as.numeric(row[["fold_idx"]]),
        lambda_alpha = as.numeric(row[["lambda_alpha"]]),
        lambda_theta = as.numeric(row[["lambda_theta"]])
      )

      env      <- sims[[s]]
      sim_dir  <- file.path(save_dir, sprintf("sim_%d", s))

      out[i] <- cv_lcd_worker(
        job        = job,
        Y_bin      = env$Y_bin,
        X          = env$X,
        bin_mass   = env$bin_mass,
        folds      = env$folds,
        K          = K,
        max_iter   = max_iter,
        iter_eta   = iter_eta,
        resp_threshold = resp_threshold,
        trim_prob  = trim_prob,
        save_dir   = sim_dir
      )
    }
    paste(out, collapse = "\n")
  })

  return(list(
    grand_index = grand_index,
    logs        = unlist(logs, use.names = FALSE)
  ))
}
