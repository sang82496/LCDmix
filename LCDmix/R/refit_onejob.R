# Generated from create-LCDmix.Rmd: do not edit by hand

#' Fit the model once with one seed and save the result
#'
#' @description
#' Calls \code{set.seed(seed)} and \code{main()} on the full data with the
#' given penalty pair, and saves the result to
#' \code{save_dir/refit_<seed>.rds}. Called by \code{refit_lcd()} and
#' \code{refit_lcd_simul()}, which run several seeds (restarts) and keep the
#' best.
#'
#' @inheritParams main
#' @inheritParams initialization
#' @param seed Seed (restart index); it selects the random flowmix start.
#' @param save_dir Directory for the result file.
#'
#' @return \code{TRUE} if the fit succeeded, \code{FALSE} otherwise. The saved
#'   file is a list with components:
#' \describe{
#'   \item{fit}{The \code{main()} result; \code{NULL} on failure.}
#'   \item{fit_L}{\code{fit$L$loglik}, the penalized training log-likelihood
#'     (bins scored by probability), used to choose among seeds; \code{NULL}
#'     on failure.}
#'   \item{fit_L_unpen}{The same without the penalty, for reporting;
#'     \code{NULL} on failure.}
#'   \item{maxdev_diag}{The result of \code{maxdev_summary()} (success only).}
#'   \item{err_msg, failed_iter}{The text from \code{fit_error_text()} and the
#'     iteration at which the fit failed (failure only).}
#'   \item{log_msg}{Printed output of the fit.}
#' }
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
  calc_Q_every = FALSE,         # NEW - appended
  update = c("lp", "optim"),    # NEW
  maxdev = NULL                 # NEW (fixP) - appended last
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
        calc_Q_every    = calc_Q_every,   # NEW
        lp_time_limit   = lp_time_limit,
        update          = update,
        maxdev          = maxdev          # NEW (fixP)
      ),
      error = function(e) { err_msg <<- e$message; NULL }
    ),
    type = "output"
  )
  log_msg <- paste0(log_msg, paste(out_log, collapse = "\n"), "\n")

  # failure path
    if (!is.list(fit) | is.null(fit$iter)) {
      err_txt <- fit_error_text(fit, err_msg)                     # NEW
      log_msg <- paste0(log_msg, "✖ Refit failed: ", err_txt)
      saveRDS(list(
        fit         = NULL,
        fit_L       = NULL,
        fit_L_unpen = NULL,
        err_msg     = err_txt,                                     # NEW
        failed_iter = if (is.list(fit) && !is.null(fit$iter_partial))  # NEW
                        fit$iter_partial$failed_iter else NA_integer_, # NEW
        log_msg     = log_msg
      ), file = out_path)
      return(FALSE)
    }
   fit_L       = as.numeric(fit$L$loglik)                  # penalized: what the EM maximizes; picks the restart
   fit_L_unpen = as.numeric(fit$L$loglik + fit$L$penalty)  # unpenalized, kept for reporting only
   log_msg <- paste0(log_msg, "✔ Completed seed ", seed_int, "; final penalized objective = ",
                       round(fit_L, 6), "\n")
  
  saveRDS(list(
      fit         = fit,
      fit_L       = fit_L,
      fit_L_unpen = fit_L_unpen,
      maxdev_diag = maxdev_summary(fit$iter),   # NEW (fixP)
      log_msg     = log_msg
    ), file = out_path)
  return(TRUE)
}
