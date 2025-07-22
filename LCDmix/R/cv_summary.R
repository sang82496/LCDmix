# Generated from create-LCDmix.Rmd: do not edit by hand

#' Summarize cross‐validation results for the LCD mixture model
#'
#' @description
#' Loads per‐run CV results saved as \code{.Rdata} files and assembles them into a
#' full CV matrix.  Computes a reduced summary of the best cross‐validation score
#' for each combination of \code{lambda_alpha} and \code{lambda_theta}, and selects
#' the optimal penalties.
#'
#' @param index_matrix Numeric matrix with columns:
#'   \describe{
#'     \item{\code{alpha_idx}}{Index of \code{lambda_alpha} in the grid.}
#'     \item{\code{theta_idx}}{Index of \code{lambda_theta} in the grid.}
#'     \item{\code{rep_idx}}{Cross‐validation repetition index.}
#'     \item{\code{fold_idx}}{Cross‐validation fold index.}
#'     \item{\code{lambda_alpha}}{Value of the \code{lambda_alpha} penalty.}
#'     \item{\code{lambda_theta}}{Value of the \code{lambda_theta} penalty.}
#'   }
#' @param save_dir Character; directory where per‐run \code{.Rdata} files are stored.
#'   Files should be named “\code{alpha_idx-theta_idx-rep_idx-fold_idx.Rdata}” and
#'   contain objects \code{prop_CV} and \code{trimmed_CV}.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{CVmat}}{Numeric matrix of dimension \code{n_runs × 8}, containing
#'     the original \code{index_matrix} plus two new columns:
#'     \code{prop_CV} and \code{trimmed_CV}.}
#'   \item{\code{reduced_mat}}{Numeric \code{n_alpha*n_theta × 3} matrix, with columns
#'     \code{lambda_alpha}, \code{lambda_theta}, and the best (max) mean
#'     \code{trimmed_CV} across repetitions for each penalty pair.}
#'   \item{\code{opt_lambdas}} Numeric vector of length 2: the optimal
#'     \code{lambda_alpha} and \code{lambda_theta} (row of \code{reduced_mat} with
#'     largest CV score).}
#'   \item{\code{max_NA_prop}} Numeric; the maximum proportion of infinite (NA)
#'     CV log‐likelihoods across all runs.}
#'
#' @examples
#' \dontrun{
#' # After running `cv_lcd_parallel()` which saved .Rdata files in "./result":
#' idx_mat <- make_cv_index_matrix(5, 5, 3, alpha_lambdas, theta_lambdas)
#' summary <- cv_summary(idx_mat, save_dir = "./result")
#' str(summary)
#' }
#' @export
cv_summary <- function(
  index_matrix,
  save_dir = "./result"
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
    rep_idx    <- index_matrix[i, "rep_idx"]
    fold_idx   <- index_matrix[i, "fold_idx"]
    file_name  <- file.path(
      save_dir,
      sprintf("%d-%d-%d-%d.Rdata",
              alpha_idx, theta_idx, rep_idx, fold_idx)
    )
    # skip runs whose file never got written
    if (!file.exists(file_name)) {
      warning("CV file missing: ", file_name, "; leaving NA")
      next
    }
    
    # try to load, but if it errors just warn & continue
    tryCatch({
      load(file_name, envir = environment())
      CVmat[i, "prop_CV"]    <- prop_CV
      CVmat[i, "trimmed_CV"] <- trimmed_CV
    }, error = function(e) {
      warning("Failed to load or assign from ", file_name, ": ", e$message)
      # CVmat row stays NA
    })
  }
  
  # Compute max proportion of infinite/NA log‐likelihoods
  max_NA_prop <- max(CVmat[, "prop_CV"], na.rm = TRUE)
  
  # Determine grid dimensions
  n_rep   <- length(unique(index_matrix[, "rep_idx"]))
  n_fold  <- length(unique(index_matrix[, "fold_idx"]))
  n_alpha <- length(unique(index_matrix[, "alpha_idx"]))
  n_theta <- length(unique(index_matrix[, "theta_idx"]))
  n_chunk <- n_alpha * n_theta
  block_sz <- n_rep * n_fold
  
  # Build reduced summary: one row per (alpha, theta) pair
  reduced_mat <- matrix(
    NA_real_,
    nrow = n_chunk,
    ncol = 3,
    dimnames = list(NULL, c("lambda_alpha", "lambda_theta", "cv_score"))
  )
  
  for (j in seq_len(n_chunk)) {
    rows_j <- ((j - 1) * block_sz + 1):(j * block_sz)
    submat <- CVmat[rows_j, , drop = FALSE]
    # Compute mean trimmed_CV per repetition, then take the maximum across reps
    mean_by_rep <- tapply(
      submat[, "trimmed_CV"],
      submat[, "rep_idx"],
      mean,
      na.rm = TRUE
    )
    best_score <- max(mean_by_rep, na.rm = TRUE)
    reduced_mat[j, ] <- c(
      submat[1, "lambda_alpha"],
      submat[1, "lambda_theta"],
      best_score
    )
  }
  
  # Select optimal penalties (max CV score)
  opt_row     <- which.max(reduced_mat[, "cv_score"])
  opt_lambdas <- reduced_mat[opt_row, c("lambda_alpha", "lambda_theta")]
  
  return(list(
    CVmat        = CVmat,
    reduced_mat  = reduced_mat,
    opt_lambdas  = opt_lambdas,
    max_NA_prop  = max_NA_prop
  ))
}
