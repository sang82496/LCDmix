# Generated from create-LCDmix.Rmd: do not edit by hand

#' Run one cross-validation job and save its result
#'
#' @description
#' Fits the model with one penalty pair and one seed on all time points except
#' one fold, evaluates it on the held-out fold with \code{eval_lcd()}, and
#' saves the result to
#' \code{save_dir/<alpha_idx>-<theta_idx>-<seed_idx>-<fold_idx>.rds}. It calls
#' \code{set.seed(seed_idx)} before \code{main()}, so the seed selects the
#' random flowmix start. Called by \code{cv_lcd()} and \code{cv_lcd_simul()}.
#'
#' @param job A named numeric vector with \code{alpha_idx}, \code{theta_idx},
#'   \code{seed_idx}, \code{fold_idx}, \code{lambda_alpha} and
#'   \code{lambda_theta}: one row of \code{cv_idx_mat()}.
#' @param folds A list with one element per fold, the indices of the held-out
#'   time points, as returned by \code{flowmix::make_cv_folds()}.
#' @param save_dir Directory for the result file.
#' @inheritParams main
#' @inheritParams initialization
#'
#' @return \code{TRUE} if the fit succeeded, \code{FALSE} otherwise. The saved
#'   file is a list with components:
#' \describe{
#'   \item{eval_trimmed_loglik, eval_finite_loglik, eval_med_loglik}{Held-out
#'     log-likelihoods from \code{eval_lcd()} (bins scored by probability)
#'     with the penalty added back, that is, unpenalized.}
#'   \item{eval_prop_inf, eval_penalty, eval_w, eval_trimmed_w}{Held-out share
#'     of weight with log-likelihood \code{-Inf}, penalty, total held-out
#'     weight, and held-out weight kept after trimming.}
#'   \item{fit_loglik, fit_trimmed_loglik, fit_med_loglik}{Training
#'     log-likelihoods from \code{fit$L}, penalized. \code{cv_lcd_summary()}
#'     uses \code{fit_trimmed_loglik} to choose among seeds.}
#'   \item{iter_num}{Number of EM iterations.}
#'   \item{n_outside_total}{Sum of \code{n_outside_every}.}
#'   \item{lp_max_over, lp_n_out_total}{Largest distance outside the support
#'     and total count from \code{lp_check_every}; \code{NA} unless
#'     \code{calc_Q_every = TRUE}.}
#'   \item{n_ascent_violations}{Number of iterations with
#'     \eqn{Q_E < Q_{ref}}; \code{NA} unless \code{calc_Q_every = TRUE}.}
#'   \item{Q_every, n_outside_every, lp_check_every}{Traces from
#'     \code{iteration()}; \code{NULL} unless \code{calc_Q_every = TRUE}.}
#'   \item{maxdev_diag}{The result of \code{maxdev_summary()}.}
#'   \item{err_msg, failed_iter}{\code{NA} on success; on failure, the text
#'     from \code{fit_error_text()} and the iteration at which the fit
#'     failed.}
#'   \item{log_msg}{Printed output of the fit.}
#' }
#' On failure all numeric fields are \code{NA} and the traces are \code{NULL}.
#' @keywords internal
#'
#' @export
cv_lcd_onejob <- function(
  job,
  Y_bin, 
  X, 
  bin_mass, 
  folds,
  K,
  max_iter, 
  iter_eta, 
  resp_threshold, 
  trim_prob,
  save_dir,
  lp_time_limit = 3600,
  calc_Q_every  = FALSE,                 # NEW
  update         = c("lp", "optim"),    # NEW
  maxdev         = NULL                 # NEW (fixP) - appended last
) {
  update <- match.arg(update)           # NEW
  alpha_idx     <- job[["alpha_idx"]]
  theta_idx     <- job[["theta_idx"]]
  seed_idx      <- job[["seed_idx"]]
  fold_idx      <- job[["fold_idx"]]
  lambda_alpha  <- job[["lambda_alpha"]]
  lambda_theta  <- job[["lambda_theta"]]

  out_path <- file.path(
    save_dir,
    sprintf("%d-%d-%d-%d.rds", alpha_idx, theta_idx, seed_idx, fold_idx)
  )

  log_msg <- paste0(
    "alpha=", lambda_alpha,
    ", theta=", lambda_theta,
    ", seed=", seed_idx,
    ", fold=", fold_idx, "\n"
  )

  # Split train/test
  test_i      <- folds[[fold_idx]]
  Y_tr        <- Y_bin[-test_i]
  bin_mass_tr <- bin_mass[-test_i]
  X_tr        <- X[-test_i, , drop = FALSE]

  # Fit
  set.seed(seed_idx)
  err_msg <- NULL
  out_log <- utils::capture.output(
    fit <- tryCatch(
      main(
        Y               = Y_tr,
        X               = X_tr,
        biomass         = bin_mass_tr,
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
        calc_Q_every    = calc_Q_every,      # NEW
        lp_time_limit   = lp_time_limit,
        update          = update,
        maxdev          = maxdev             # NEW (fixP)
      ),
      error = function(e) { err_msg <<- e$message; NULL }
    ),
    type = "output"
  )

  log_msg <- paste0(log_msg, paste0(out_log, collapse = "\n"), "\n")

  if (!is.list(fit) || is.null(fit$iter)) {
    err_txt <- fit_error_text(fit, err_msg)                       # NEW
    failed_iter <- if (is.list(fit) && !is.null(fit$iter_partial)) # NEW
                     fit$iter_partial$failed_iter else NA_integer_ # NEW
    log_msg <- paste0(log_msg, "✖ Fit failed: ", err_txt)
    saveRDS(
      list(
        eval_prop_inf       = NA_real_,
        eval_trimmed_loglik = NA_real_,
        eval_finite_loglik  = NA_real_,
        eval_med_loglik     = NA_real_,
        eval_penalty        = NA_real_,
        eval_w              = NA_real_,
        eval_trimmed_w      = NA_real_,
        fit_loglik          = NA_real_,
        fit_trimmed_loglik  = NA_real_,
        fit_med_loglik      = NA_real_,
        err_msg             = err_txt,                             # NEW
        failed_iter         = failed_iter,                         # NEW
       iter_num = NA_integer_, 
       n_outside_total = NA_real_,
       lp_max_over = NA_real_, 
       lp_n_out_total = NA_real_,
       n_ascent_violations = NA_integer_,
       Q_every = NULL, 
       n_outside_every = NULL, 
       lp_check_every = NULL,
       maxdev_diag    = NULL,
        log_msg             = log_msg
      ),
      file = out_path
    )
    return(FALSE)
  }

  # Evaluate on hold-out
  eval_res <- eval_lcd(
    model        = fit$iter,
    Y_test       = Y_bin[test_i],
    X_test       = X[test_i, , drop = FALSE],
    biomass_test = bin_mass[test_i],
    trim_prob    = trim_prob
  )
  
  log_msg <- paste0(log_msg, "✔ Saved: ", basename(out_path))

  saveRDS(
    list(
      eval_prop_inf       = eval_res$prop_inf,
      eval_trimmed_loglik = eval_res$trimmed_loglik + eval_res$penalty,
      eval_finite_loglik  = eval_res$finite_loglik + eval_res$penalty,
      eval_med_loglik     = eval_res$med_loglik + eval_res$penalty,
      eval_penalty        = eval_res$penalty,
      eval_w              = eval_res$sum_w,
      eval_trimmed_w      = eval_res$sum_trimmed_w,
      fit_loglik          = fit$L$loglik,
      fit_trimmed_loglik  = fit$L$trimmed_loglik,
      fit_med_loglik      = fit$L$med_loglik,
       # --- NEW: diagnostics -------------------------------------------------
       iter_num            = fit$iter$iter_num,
       n_outside_total     = sum(fit$iter$n_outside_every),
       lp_max_over         = if (length(fit$iter$lp_check_every))
                               max(vapply(fit$iter$lp_check_every,
                                          function(m) max(m[, "max_over"]), numeric(1))) else NA_real_,
       lp_n_out_total      = if (length(fit$iter$lp_check_every))
                               sum(vapply(fit$iter$lp_check_every,
                                          function(m) sum(m[, "n_out"]), numeric(1))) else NA_real_,
       n_ascent_violations = if (isTRUE(calc_Q_every)) {
                               q <- fit$iter$Q_every; it <- fit$iter$iter_num
                               sum(q[5 * seq_len(it) + 1] - q[5 * seq_len(it) - 3] < 0)
                             } else NA_integer_,
       Q_every             = if (isTRUE(calc_Q_every)) fit$iter$Q_every         else NULL,
       n_outside_every     = if (isTRUE(calc_Q_every)) fit$iter$n_outside_every else NULL,
       lp_check_every      = if (isTRUE(calc_Q_every)) fit$iter$lp_check_every  else NULL,
       maxdev_diag         = maxdev_summary(fit$iter),              # NEW (fixP)
       # ----------------------------------------------------------------------
       err_msg             = NA_character_,
       failed_iter         = NA_integer_,
       log_msg             = log_msg
    ),
    file = out_path
  )
  return(TRUE)
}
