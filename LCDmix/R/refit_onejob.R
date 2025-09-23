# Generated from create-LCDmix.Rmd: do not edit by hand

#' @title Internal: run one LCDmix refit (single seed) and cache
#' @keywords internal
#' @export
refit_onejob <- function(
  Y_bin,
  X,
  bin_mass,
  K,
  lambda_alpha,
  lambda_theta,
  seed,
  max_iter,
  iter_eta,
  resp_threshold,
  trim_prob,
  save_dir,
  debug = FALSE
) {
  seed_int <- as.integer(seed)
  out_path <- file.path(save_dir, sprintf("refit_%d.rds", seed_int))

  # Fast path: cached result
  if (file.exists(out_path)) {
    saved <- tryCatch(readRDS(out_path), error = function(e) NULL)
    if (is.list(saved)) {
      L_cached <- if (!is.null(saved$L)) as.numeric(saved$L) else NA_real_
      return(list(log_msg = saved$log_msg, L = L_cached, fit = saved$fit, file = out_path))
    }
    # if cache unreadable, fall through to re-fit
  }

  log_msg <- paste0(
    "▶ Refit seed ", seed_int,
    " | lambda_alpha=", lambda_alpha,
    ", lambda_theta=", lambda_theta, "\n"
  )

  set.seed(seed_int)
  err_msg <- NULL
  out_log <- utils::capture.output(
    fit <- tryCatch(
      main(
        Y               = Y_bin,
        X               = X,
        biomass         = bin_mass,
        binned          = TRUE,
        n_bins          = 0,
        K               = K,
        lambda_alpha    = lambda_alpha,
        lambda_theta    = lambda_theta,
        max_iter        = max_iter,
        iter_eta        = iter_eta,
        resp_threshold  = resp_threshold,
        trim_prob       = trim_prob,
        debug           = debug
      ),
      error = function(e) { err_msg <<- e$message; NULL }
    ),
    type = "output"
  )
  log_msg <- paste0(log_msg, paste(out_log, collapse = "\n"), "\n")

  if (!is.list(fit) || is.null(fit$L) || is.null(fit$L$loglik)) {
    log_msg <- paste0(log_msg, "✖ Refit failed: ", err_msg)
    saveRDS(list(log_msg = log_msg, L = NA_real_, fit = NULL), file = out_path)
    return(list(log_msg = log_msg, L = NA_real_, fit = NULL, file = out_path))
  }

  L <- as.numeric(fit$L$loglik)
  log_msg <- paste0(
    log_msg,
    "✔ Completed seed ", seed_int,
    "; final L = ", round(L, 6), "\n"
  )

  saveRDS(list(log_msg = log_msg, L = L, fit = fit), file = out_path)

  return(list(log_msg = log_msg, L = L, fit = fit, file = out_path))
}
