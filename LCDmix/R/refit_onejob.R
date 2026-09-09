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
  lp_time_limit,
  update = c("lp", "optim")     # NEW - must be LAST
) {
  update <- match.arg(update)           # NEW
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
        debug           = TRUE,
        lp_time_limit   = lp_time_limit,
        update          = update
      ),
      error = function(e) { err_msg <<- e$message; NULL }
    ),
    type = "output"
  )
  log_msg <- paste0(log_msg, paste(out_log, collapse = "\n"), "\n")

  # failure path
    if (!is.list(fit) | is.null(fit$iter)) {
      err_txt <- .fit_error_text(fit, err_msg)                     # NEW
      log_msg <- paste0(log_msg, "✖ Refit failed: ", err_txt)
      saveRDS(list(
        fit         = NULL,
        fit_L       = NULL,
        err_msg     = err_txt,                                     # NEW
        failed_iter = if (is.list(fit) && !is.null(fit$iter_partial))  # NEW
                        fit$iter_partial$failed_iter else NA_integer_, # NEW
        log_msg     = log_msg
      ), file = out_path)
      return(FALSE)
    }
  fit_L = as.numeric(fit$L$loglik + fit$L$penalty) # unpenalized loglik
  log_msg <- paste0(log_msg, "✔ Completed seed ", seed_int, "; final loglikelihood = ", 
                      round(fit_L, 6), "\n")
  
  saveRDS(list(
      fit     = fit,
      fit_L   = fit_L,
      log_msg = log_msg
    ), file = out_path)
  return(TRUE)
}
