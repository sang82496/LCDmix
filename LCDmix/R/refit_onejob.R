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

  # failure path
    if (!is.list(fit)) {
      log_msg <- paste0(log_msg, "✖ Refit failed: ", err_msg)
      saveRDS(list(
        fit      = NULL,
        log_msg  = log_msg
      ), file = out_path)
      return(FALSE)
    }
  
  fit_L <- as.numeric(fit$L$trimmed_loglik)
  log_msg <- paste0(log_msg, "✔ Completed seed ", seed_int, "; final trimmed loglikelihood = ", 
                    round(fit_L, 6), "\n")
  saveRDS(list(
      fit     = fit,
      log_msg = log_msg
    ), file = out_path)
  return(TRUE)
}
