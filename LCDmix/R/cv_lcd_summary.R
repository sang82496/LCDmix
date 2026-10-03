# Generated from create-LCDmix.Rmd: do not edit by hand

#' Choose the penalty pair from saved cross-validation results
#'
#' @description
#' Reads the files that \code{cv_lcd_onejob()} saved in \code{save_dir} and
#' computes a cross-validation score for every penalty pair:
#'
#' 1. For each pair and fold, keeps the seed (restart) with the largest
#'    training log-likelihood: \code{fit_trimmed_loglik} (penalized) if
#'    \code{cv_by_trimmed = TRUE}, otherwise \code{fit_med_loglik}.
#' 2. Averages the held-out score of the kept fits over the folds:
#'    \code{eval_trimmed_loglik} (or \code{eval_med_loglik}), which is
#'    unpenalized and scores bins by their probability.
#' 3. Chooses the pair with the largest score. \code{NA} scores count as
#'    \code{-Inf}.
#'
#' A missing file gives a warning and is treated as a failed fit.
#'
#' The two stages use different criteria on purpose, as in flowmix
#' (\code{flowmix::one_job()} scores held-out data with both penalties set to
#' 0, and \code{flowmix::cv_aggregate()} keeps the restart with the best
#' penalized training objective). Step 1 compares restarts of the same
#' penalized problem, so it uses a penalized criterion. Step 2
#' estimates how well each fit predicts new data; the penalty is part of the
#' fitting objective, not of the predictive density, so it is left out.
#'
#' @param index_matrix The matrix from \code{cv_idx_mat()}, saved as
#'   \code{index_matrix.rds} by \code{cv_lcd()} and \code{cv_lcd_simul()}.
#' @param save_dir Directory with the per-job files.
#' @param cv_by_trimmed \code{TRUE} to use trimmed log-likelihoods,
#'   \code{FALSE} to use weighted medians.
#' @param equal_w \code{TRUE} to give each fold the same weight in step 2;
#'   \code{FALSE} to weight each fold by its held-out weight (the weight kept
#'   after trimming when \code{cv_by_trimmed = TRUE}).
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{CVmat}}{Data frame: \code{index_matrix} with the
#'     \code{fit_*} and \code{eval_*} values of each job.}
#'   \item{\code{reduced_mat}}{Matrix with one row per pair: \code{alpha_idx},
#'     \code{theta_idx}, \code{lambda_alpha}, \code{lambda_theta},
#'     \code{cv_score}, and \code{counts} (jobs, over all seeds and folds,
#'     with finite training and held-out values).}
#'   \item{\code{opt_lambdas}}{\code{c(lambda_alpha, lambda_theta)} of the
#'     chosen pair.}
#'   \item{\code{max_prop_inf}}{Largest \code{eval_prop_inf} over all jobs.}
#'   \item{\code{cv_by_trimmed}, \code{equal_w}}{The arguments, returned
#'     unchanged.}
#' }
#'
#' @seealso \code{\link{cv_lcd}}, \code{\link{cv_lcd_summary_simul}},
#'   \code{\link{eval_lcd}}, \code{\link{cv_idx_mat}}
#'
#' @examples
#' \dontrun{
#' # After cv_lcd(..., save_dir = "cv_runs"):
#' index_matrix <- readRDS(file.path("cv_runs", "index_matrix.rds"))
#' summ <- cv_lcd_summary(index_matrix, save_dir = "cv_runs")
#' summ$opt_lambdas
#' }
#'
#' @export
cv_lcd_summary <- function(
  index_matrix,
  save_dir = "./cv_saves",
  cv_by_trimmed = T,
  equal_w = T
) {
  n_runs <- nrow(index_matrix)
  # Initialize CV matrix with placeholders
  CVmat <- cbind(
    index_matrix        = index_matrix,
    fit_trimmed_loglik  = rep(NA_real_, n_runs),
    fit_med_loglik      = rep(NA_real_, n_runs),
    fit_loglik          = rep(NA_real_, n_runs),
    eval_prop_inf       = rep(NA_real_, n_runs),
    eval_trimmed_loglik = rep(NA_real_, n_runs),
    eval_med_loglik     = rep(NA_real_, n_runs),
    eval_penalty        = rep(NA_real_, n_runs),
    eval_w              = rep(NA_real_, n_runs),
    eval_trimmed_w      = rep(NA_real_, n_runs)
  )
  CVmat <- as.data.frame(CVmat)
  
  # Load each run's results
  for (i in seq_len(n_runs)) {
    alpha_idx  <- index_matrix[i, "alpha_idx"]
    theta_idx  <- index_matrix[i, "theta_idx"]
    seed_idx   <- index_matrix[i, "seed_idx"]
    fold_idx   <- index_matrix[i, "fold_idx"]
    
    file_name  <- file.path(
    save_dir,
    sprintf("%d-%d-%d-%d.rds",
            alpha_idx, theta_idx, seed_idx, fold_idx))
    
    # skip runs whose file never got written
    if (!file.exists(file_name)) {
      warning("CV file missing: ", file_name, "; leaving NA")
      next
    }
    # load
    mat = readRDS(file_name)
    CVmat[i, "fit_trimmed_loglik"]  <- mat$fit_trimmed_loglik
    CVmat[i, "fit_med_loglik"]      <- mat$fit_med_loglik
    CVmat[i, "fit_loglik"]          <- mat$fit_loglik
    CVmat[i, "eval_prop_inf"]       <- mat$eval_prop_inf 
    CVmat[i, "eval_trimmed_loglik"] <- mat$eval_trimmed_loglik
    CVmat[i, "eval_med_loglik"]     <- mat$eval_med_loglik
    CVmat[i, "eval_penalty"]        <- mat$eval_penalty
    CVmat[i, "eval_w"]              <- mat$eval_w
    CVmat[i, "eval_trimmed_w"]      <- mat$eval_trimmed_w
  }
  ## 4) For each (alpha_idx, theta_idx, fold_idx), keep the row with largest fit_loglik
  gid <- interaction(CVmat$alpha_idx, CVmat$theta_idx, CVmat$fold_idx, drop = TRUE)
  row_groups <- split(seq_len(nrow(CVmat)), gid)
  
  if (cv_by_trimmed) {
    pick <- CVmat$fit_trimmed_loglik
  } else {
    pick <- CVmat$fit_med_loglik
  }

  pick_idx <- vapply(row_groups, function(idx) {
    v <- pick[idx]
    if (all(!is.finite(v))) return(NA_integer_)  # no usable value in this group
    idx[which.max(v)]                 # pick the row index with max loglik
  }, integer(1L))
  
  selected <- CVmat[pick_idx[!is.na(pick_idx)], , drop = FALSE]
  
  # Group selected rows by (alpha_idx, theta_idx)
  gid2 <- interaction(selected$alpha_idx, selected$theta_idx, drop = TRUE)
  row_groups2 <- split(seq_len(nrow(selected)), gid2)

  ## 5) Average those per-fold maxima across folds → cv_score per (alpha,theta)
  cv_score <- vapply(row_groups2, function(idx) {
    if (equal_w) {
      if (cv_by_trimmed) {
        mean(selected$eval_trimmed_loglik[idx], na.rm = TRUE)
      } else {
        mean(selected$eval_med_loglik[idx], na.rm = TRUE)
      }
    } else {
      if (cv_by_trimmed) {
        weighted.mean(selected$eval_trimmed_loglik[idx], selected$eval_trimmed_w[idx], na.rm = TRUE)
      } else {
        weighted.mean(selected$eval_med_loglik[idx], selected$eval_w[idx], na.rm = TRUE)
      }
    }
  }, double(1L))
  
  mat <- cbind(
    do.call(rbind, lapply(strsplit(names(row_groups2), "[.]"), as.numeric)),
    cv_score
  )
  colnames(mat)[1:2] <- c("alpha_idx", "theta_idx")
  
  ## 6) Unique lambda values per combo (from the original matrix)
  unique_lams <- unique(CVmat[, c("alpha_idx","theta_idx","lambda_alpha","lambda_theta")])

  ## 7) Count rows in CVmat with non-NA and finite loglik per (alpha_idx, theta_idx)
  if (cv_by_trimmed) {
    counts_df <- aggregate(
      I(is.finite(eval_trimmed_loglik) & is.finite(fit_trimmed_loglik)) ~ alpha_idx + theta_idx,
      data = CVmat, FUN  = sum)
  } else {
    counts_df <- aggregate(
      I(is.finite(eval_med_loglik) & is.finite(fit_med_loglik)) ~ alpha_idx + theta_idx,
      data = CVmat, FUN  = sum)
  }
  names(counts_df)[3] <- "counts"
  
  ## 8) cbind
  reduced_mat <- merge(unique_lams, mat, by = c("alpha_idx","theta_idx"), all.x = TRUE)
  reduced_mat <- merge(reduced_mat, counts_df, by = c("alpha_idx","theta_idx"), all.x = TRUE)
  
  ## 9) Pick optimal lambdas
  score_vec <- reduced_mat$cv_score
  score_vec[!is.finite(score_vec)] <- -Inf   # handles NA and NaN
  opt_lambdas <- as.numeric(reduced_mat[which.max(score_vec), c("lambda_alpha", "lambda_theta")])

  return(list(
    CVmat        = CVmat,
    reduced_mat  = as.matrix(reduced_mat),
    opt_lambdas  = opt_lambdas,
    max_prop_inf = max(CVmat[, "eval_prop_inf"], na.rm = TRUE),
    cv_by_trimmed = cv_by_trimmed,
    equal_w       = equal_w
  ))
}
