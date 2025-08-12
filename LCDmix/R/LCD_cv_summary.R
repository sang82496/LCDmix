# Generated from create-LCDmix.Rmd: do not edit by hand

#' Summarize CV results and select optimal penalty values
#'
#' @description
#' Loads saved per‐run CV results (from \code{save_dir}) according to the
#' supplied \code{index_matrix}, warns on any missing or failed loads, and then
#' aggregates:
#' 1. Per‐seed average of \code{trimmed_CV} across folds,
#' 2. The maximum of those seed‐averages per (\code{alpha_idx},\code{theta_idx}),
#' 3. The count of non‐missing runs per combination,
#' and merges these with the corresponding penalty values.  Finally selects the
#' optimal \code{lambda_alpha} and \code{lambda_theta} that maximize the CV score.
#'
#' @param index_matrix Integer matrix with one row per CV run, containing columns
#'   \code{alpha_idx}, \code{theta_idx}, \code{seed_idx}, \code{fold_idx},
#'   \code{lambda_alpha}, and \code{lambda_theta}.
#' @param save_dir     Character; directory where the CV run \*.Rdata files are stored.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{CVmat}}{Numeric matrix of dimension \code{n_runs × 8}, combining
#'     the original \code{index_matrix} columns with two new columns:
#'     \code{prop_CV} and \code{trimmed_CV}, one row per run.}
#'   \item{\code{reduced_mat}}{Numeric matrix with one row per unique
#'     (\code{alpha_idx},\code{theta_idx}) combination, and columns:
#'     \code{alpha_idx}, \code{theta_idx}, \code{lambda_alpha},
#'     \code{lambda_theta}, \code{cv_score} (max seed‐average),
#'     and \code{count} (number of non‐NA runs).}
#'   \item{\code{opt_lambdas}}{Numeric vector of length 2 giving the
#'     optimal \code{lambda_alpha} and \code{lambda_theta}.}
#'   \item{\code{max_NA_prop}}{Numeric; the maximum \code{prop_CV} value
#'     across all runs, ignoring \code{NA}.}
#' }
#'
#' @examples
#' \dontrun{
#' # Given `index_matrix` and CV .Rdata files in "./result":
#' summary <- LCD_cv_summary(index_matrix, save_dir = "./result")
#' print(summary$opt_lambdas)
#' head(summary$reduced_mat)
#' }
#'
#' @export
LCD_cv_summary <- function(
  index_matrix,
  save_dir = "./result",
  simul    = FALSE,
  sim_idx  = NA
) {
  n_runs <- nrow(index_matrix)
  # Initialize CV matrix with placeholders
  CVmat <- cbind(
    index_matrix,
    prop_CV     = rep(NA_real_, n_runs),
    trimmed_CV  = rep(NA_real_, n_runs)
  )
  
  # Load each run's results
  for (i in seq_len(n_runs)) {
    alpha_idx  <- index_matrix[i, "alpha_idx"]
    theta_idx  <- index_matrix[i, "theta_idx"]
    seed_idx   <- index_matrix[i, "seed_idx"]
    fold_idx   <- index_matrix[i, "fold_idx"]
    
    if (!simul){
      file_name  <- file.path(
      save_dir,
      sprintf("%d-%d-%d-%d.Rdata",
              alpha_idx, theta_idx, seed_idx, fold_idx))
      } else {
      file_name  <- file.path(
      save_dir,
      sprintf("%d-%d-%d-%d-%d.Rdata",
              sim_idx, alpha_idx, theta_idx, seed_idx, fold_idx))
    }
    
    # skip runs whose file never got written
    if (!file.exists(file_name)) {
      warning("CV file missing: ", file_name, "; leaving NA")
      next
    }
    # load
    load(file_name, envir = environment())
    CVmat[i, "prop_CV"]    <- prop_CV
    CVmat[i, "trimmed_CV"] <- trimmed_CV
  }
  good <- !is.na(CVmat[,"trimmed_CV"])
  if (!any(good)) {
    stop("No successful CV runs found: all trimmed_CV are NA")
  }
  
  # 4) Compute per‐seed mean over folds, *only* on those good rows
  rep_means <- aggregate(
    trimmed_CV ~ alpha_idx + theta_idx + seed_idx,
    data    = CVmat,
    subset  = good,
    FUN     = mean  # no need for na.rm=TRUE since we filtered out NA
  )
  
  # 5) Take max across seeds → cv_score
  best_by_combo <- aggregate(
    trimmed_CV ~ alpha_idx + theta_idx,
    data   = rep_means,
    FUN    = max
  )
  names(best_by_combo)[3] <- "cv_score"
  
  # 6) Unique lambda values per combo
  unique_lams <- unique(
    CVmat[, c("alpha_idx","theta_idx","lambda_alpha","lambda_theta")]
  )
  
  # 7) Count non‐NA original runs per combo
  counts_df <- aggregate(
    I(!is.na(trimmed_CV)) ~ alpha_idx + theta_idx,
    data = CVmat,
    FUN  = sum
  )
  names(counts_df)[3] <- "count"
  
  # 8) Merge & sort
  reduced_mat <- merge(unique_lams,  best_by_combo, by = c("alpha_idx","theta_idx"))
  reduced_mat <- merge(reduced_mat, counts_df,     by = c("alpha_idx","theta_idx"))
  reduced_mat <- reduced_mat[order(reduced_mat$alpha_idx, reduced_mat$theta_idx), ]
  
  # 9) Pick optimal lambdas
  opt_row     <- which.max(reduced_mat$cv_score)
  opt_lambdas <- as.numeric(reduced_mat[opt_row, c("lambda_alpha","lambda_theta")])
  
  list(
    CVmat        = CVmat,
    reduced_mat  = as.matrix(reduced_mat),
    opt_lambdas  = opt_lambdas,
    max_NA_prop  = max(CVmat[, "prop_CV"], na.rm = TRUE)
  )
}
