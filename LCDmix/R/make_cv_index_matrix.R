# Generated from create-LCDmix.Rmd: do not edit by hand

#' Generate cross‐validation index matrix for penalty grid search
#'
#' @description
#' Constructs a matrix of all combinations of cross‐validation fold, random seeds (repeats),
#' and penalty indices for tuning \code{lambda_alpha} and \code{lambda_theta}.
#' 
#'
#' @param cv_gridsize Integer; number of candidate values for each penalty (\code{lambda_alpha} and \code{lambda_theta}).
#' @param nfold Integer; number of folds in cross‐validation.
#' @param seeds Integer vector. Random seeds (or repeat identifiers) to use for each repetition.
#' @param alpha_lambdas Numeric vector of length \code{cv_gridsize}; the candidate values of \code{lambda_alpha}.
#' @param theta_lambdas Numeric vector of length \code{cv_gridsize}; the candidate values of \code{lambda_theta}.
#'
#' @return
#' A numeric matrix with one row per combination (total rows = \code{cv_gridsize^2 * nfold * nrep})
#' and six columns:
#' \describe{
#'   \item{\code{alpha_idx}}{Index in \code{alpha_lambdas} (1..\code{cv_gridsize}).}
#'   \item{\code{theta_idx}}{Index in \code{theta_lambdas} (1..\code{cv_gridsize}).}
#'   \item{\code{seed_idx}}{The seed or repeat identifier for this run.}
#'   \item{\code{fold_idx}}{Fold index (1..\code{nfold}).}
#'   \item{\code{lambda_alpha}}{Value of \code{alpha_lambdas[alpha_idx]}.}
#'   \item{\code{lambda_theta}}{Value of \code{theta_lambdas[theta_idx]}.}
#' }
#'
#' @examples
#' \dontrun{
#' # 5 candidate lambdas for each penalty, 4 folds, 3 repetitions
#' alpha_vals <- seq(0.01, 0.1, length.out = 5)
#' theta_vals <- seq(0.001, 0.01, length.out = 5)
#' seeds      <- c(101, 202, 303)  # three repeats
#' idx_mat <- make_cv_index_matrix(
#'   cv_gridsize   = 5,
#'   nfold         = 4,
#'   seeds         = seeds,
#'   alpha_lambdas = alpha_vals,
#'   theta_lambdas = theta_vals
#' )
#' head(idx_mat)
#' }
#' @export
make_cv_index_matrix <- function(
  cv_gridsize,
  nfold,
  seeds,
  alpha_lambdas,
  theta_lambdas
) {
  # Create all combinations of fold, repetition, theta‐index, and alpha‐index
  grid <- expand.grid(
    fold_idx   = seq_len(nfold),
    seed_idx   = seeds,
    theta_idx  = seq_len(cv_gridsize),
    alpha_idx  = seq_len(cv_gridsize)
  )
  
  # Assemble the index matrix with corresponding lambda values
  iimat <- cbind(
    alpha_idx    = grid$alpha_idx,
    theta_idx    = grid$theta_idx,
    seed_idx     = grid$seed_idx,
    fold_idx     = grid$fold_idx,
    lambda_alpha = alpha_lambdas[grid$alpha_idx],
    lambda_theta = theta_lambdas[grid$theta_idx]
  )
  return(iimat)
}
