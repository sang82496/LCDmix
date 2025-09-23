# Generated from create-LCDmix.Rmd: do not edit by hand

#' Create a cross-validation index matrix (folds × seeds × penalty grids)
#'
#' @description
#' Builds the full factorial grid over CV folds, random seeds, and the indices
#' of \code{alpha_lambdas} and \code{theta_lambdas}, and attaches the actual
#' penalty values for each row. Each row corresponds to one CV job.
#'
#' @param nfold Integer. Number of CV folds.
#' @param seeds Integer vector of random seeds (or repeat identifiers). One
#'   row will be created for each seed × fold × penalty-pair combination.
#' @param alpha_lambdas Numeric vector of candidate \code{lambda_alpha} values.
#' @param theta_lambdas Numeric vector of candidate \code{lambda_theta} values.
#'
#' @return A numeric matrix with columns:
#' \describe{
#'   \item{\code{alpha_idx}}{Index into \code{alpha_lambdas} (1..length).}
#'   \item{\code{theta_idx}}{Index into \code{theta_lambdas} (1..length).}
#'   \item{\code{seed_idx}}{Seed used for this run (from \code{seeds}).}
#'   \item{\code{fold_idx}}{Fold index (1..\code{nfold}).}
#'   \item{\code{lambda_alpha}}{Actual \code{alpha_lambdas[alpha_idx]}.}
#'   \item{\code{lambda_theta}}{Actual \code{theta_lambdas[theta_idx]}.}
#' }
#'
#' @examples
#' \dontrun{
#' idx <- cv_idx_mat(
#'   nfold         = 5,
#'   seeds         = c(101, 202, 303),
#'   alpha_lambdas = c(1e-4, 1e-3, 1e-2),
#'   theta_lambdas = c(1e-4, 1e-3, 1e-2)
#' )
#' head(idx)
#' }
#'
#' @export
cv_idx_mat <- function(
  nfold,
  seeds,
  alpha_lambdas,
  theta_lambdas
) {
  alpha_size = length(alpha_lambdas)
  theta_size = length(theta_lambdas)
  # Create all combinations of fold, repetition, theta‐index, and alpha‐index
  grid <- expand.grid(
    fold_idx   = seq_len(nfold),
    seed_idx   = seeds,
    theta_idx  = seq_len(theta_size),
    alpha_idx  = seq_len(alpha_size)
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
