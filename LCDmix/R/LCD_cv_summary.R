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
    L           = rep(NA_real_, n_runs),
    prop_CV     = rep(NA_real_, n_runs),
    trimmed_CV  = rep(NA_real_, n_runs)
  )
  CVmat <- as.data.frame(CVmat)
  
  # Load each run's results
  for (i in seq_len(n_runs)) {
    alpha_idx  <- index_matrix[i, "alpha_idx"]
    theta_idx  <- index_matrix[i, "theta_idx"]
    seed_idx   <- index_matrix[i, "seed_idx"]
    fold_idx   <- index_matrix[i, "fold_idx"]
    
    if (!simul){
      file_name  <- file.path(
      save_dir,
      sprintf("%d-%d-%d-%d.rds",
              alpha_idx, theta_idx, seed_idx, fold_idx))
      } else {
      file_name  <- file.path(
      save_dir,
      sprintf("%d-%d-%d-%d-%d.rds",
              sim_idx, alpha_idx, theta_idx, seed_idx, fold_idx))
    }
    
    # skip runs whose file never got written
    if (!file.exists(file_name)) {
      warning("CV file missing: ", file_name, "; leaving NA")
      next
    }
    # load
    mat = readRDS(file_name)
    CVmat[i, "L"]          <- mat$L
    CVmat[i, "prop_CV"]    <- mat$prop_CV
    CVmat[i, "trimmed_CV"] <- mat$trimmed_CV
  }
  
  ## 4) For each (alpha_idx, theta_idx, fold_idx), keep the row with largest L
  # Order so that within each (alpha,theta,fold) the first row has max L (NAs go last)
  ord <- with(CVmat,
    order(alpha_idx, theta_idx, fold_idx, -L, na.last = TRUE)  # NA Q go last
  )
  selected <- CVmat[ord, ]
  selected <- selected[!duplicated(selected[c("alpha_idx","theta_idx","fold_idx")]), ]
  
  ## Warn if any selected rows have NA trimmed_CV
  if (any(is.na(selected$trimmed_CV))) {
    warning("Selected rows contain NA trimmed_CV values; they will be ignored when averaging.")
  }
  
  ## 5) Average the selected trimmed_CV across folds -> cv_score per (alpha,theta)
  best_by_combo <- aggregate(
    trimmed_CV ~ alpha_idx + theta_idx,
    data = selected,
    FUN  = function(x) mean(x, na.rm = TRUE)
  )
  names(best_by_combo)[3] <- "cv_score"
  
  ## 6) Unique lambda values per combo (from the original matrix)
  unique_lams <- unique(CVmat[, c("alpha_idx","theta_idx","lambda_alpha","lambda_theta")])

  ## 7) Count rows in CVmat with non-NA L per (alpha_idx, theta_idx)
  counts_df <- aggregate(
    I(!is.na(L)) ~ alpha_idx + theta_idx,
    data = CVmat,
    FUN  = sum
  )
  names(counts_df)[3] <- "counts"
  
  ## 8) Merge & sort  (unchanged except for the new column name)
  reduced_mat <- merge(unique_lams, best_by_combo, by = c("alpha_idx","theta_idx"), all.x = TRUE)
  reduced_mat <- merge(reduced_mat, counts_df, by = c("alpha_idx","theta_idx"), all.x = TRUE)
  reduced_mat <- reduced_mat[order(reduced_mat$alpha_idx, reduced_mat$theta_idx), ]
  
  ## 9) Pick optimal lambdas (same as before)
  opt_row     <- which.max(reduced_mat$cv_score)
  opt_lambdas <- as.numeric(reduced_mat[opt_row, c("lambda_alpha","lambda_theta")])
  
  return(list(
    CVmat        = CVmat,
    reduced_mat  = as.matrix(reduced_mat),
    opt_lambdas  = opt_lambdas,
    max_NA_prop  = max(CVmat[, "prop_CV"], na.rm = TRUE)
  ))
}
