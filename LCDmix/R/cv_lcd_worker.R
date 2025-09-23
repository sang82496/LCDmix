# Generated from create-LCDmix.Rmd: do not edit by hand

#' @title Internal: run a single CV job (one fold/seed/λ-pair)
#' @description Executes one CV job: fit on train split, evaluate on hold-out,
#'   and write a per-job RDS file. Returns a concise log string.
#' @param job named numeric vector/data.frame row with fields:
#'   alpha_idx, theta_idx, seed_idx, fold_idx, lambda_alpha, lambda_theta
#' @param Y_bin list of responses; @param X matrix; @param bin_mass list of weights
#' @param folds list of fold indices as from flowmix::make_cv_folds()
#' @param K,max_iter,iter_eta,resp_threshold,trim_prob standard LCDmix args
#' @param save_dir directory to write "<alpha>-<theta>-<seed>-<fold>.rds"
#' @return character (one log line); writes result file as side-effect
#' @keywords internal
#' 
#' @export
cv_lcd_worker <- function(
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
  save_dir
) {
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
  if (file.exists(out_path)) {
    return(paste0("Skipping existing: alpha=", lambda_alpha,
                  ", theta=", lambda_theta,
                  ", seed=",  seed_idx,
                  ", fold=",  fold_idx))
  }

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
    fit_try <- tryCatch(
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
        debug           = TRUE
      ),
      error = function(e) { err_msg <<- e$message; NULL }
    ),
    type = "output"
  )

  log_msg <- paste0(log_msg, paste0(out_log, collapse = "\n"), "\n")

  if (!is.list(fit_try) || is.null(fit_try$iter)) {
    log_msg <- paste0(log_msg, "✖ Fit failed: ", err_msg)
    saveRDS(
      list(
        prop_inf       = NA_real_,
        trimmed_loglik = NA_real_,
        finite_loglik  = NA_real_,
        L              = NA_real_,
        log_msg        = log_msg
      ),
      file = out_path
    )
    return(log_msg)
  }

  # Evaluate on hold-out
  eval_res <- evaluate_lcd_model(
    model        = fit_try$iter,
    Y_test       = Y_bin[test_i],
    X_test       = X[test_i, , drop = FALSE],
    biomass_test = bin_mass[test_i],
    trim_prob    = trim_prob
  )

  prop_inf       <- eval_res$prop_inf
  trimmed_loglik <- eval_res$trimmed_loglik
  finite_loglik  <- eval_res$finite_loglik
  L              <- fit_try$L$loglik

  log_msg <- paste0(log_msg, "✔ Saved: ", basename(out_path))

  saveRDS(
    list(
      prop_inf       = prop_inf,
      trimmed_loglik = trimmed_loglik,
      finite_loglik  = finite_loglik,
      L              = L,
      log_msg        = log_msg
    ),
    file = out_path
  )

  return(log_msg)
}
