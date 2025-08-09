# Generated from create-LCDmix.Rmd: do not edit by hand

#' Build a CV job index matrix across simulations, folds, seeds, and penalties
#'
#' @description
#' Creates the full factorial grid over simulation IDs, CV folds, random seeds,
#' and the indices of \code{alpha_lambdas} and \code{theta_lambdas}. It also
#' attaches the actual \code{lambda_alpha} and \code{lambda_theta} values for
#' each row. Each row corresponds to one cross-validation job.
#'
#' @param n_sims Integer number of simulated datasets (files \code{sim_<i>.rds}).
#' @param alpha_lambdas Numeric vector of candidate \code{lambda_alpha} values.
#' @param theta_lambdas Numeric vector of candidate \code{lambda_theta} values.
#' @param seeds Vector of random seeds to use for repeated CV runs.
#' @param nfold Integer number of CV folds per simulation.
#'
#' @return A numeric matrix with columns:
#' \describe{
#'   \item{sim_idx}{Simulation index (1..\code{n_sims}).}
#'   \item{alpha_idx}{Index into \code{alpha_lambdas}.}
#'   \item{theta_idx}{Index into \code{theta_lambdas}.}
#'   \item{seed_idx}{Random seed used for that job.}
#'   \item{fold_idx}{CV fold index (1..\code{nfold}).}
#'   \item{lambda_alpha}{Actual \code{alpha_lambdas[alpha_idx]}.}
#'   \item{lambda_theta}{Actual \code{theta_lambdas[theta_idx]}.}
#' }
#'
#' @examples
#' \dontrun{
#' idx_mat <- simul_idx_matrix(
#'   n_sims        = 3,
#'   alpha_lambdas = c(1e-3, 1e-2),
#'   theta_lambdas = c(1e-3, 1e-2, 1e-1),
#'   seeds         = 1:2,
#'   nfold         = 5
#' )
#' head(idx_mat)
#' }
#'
#' @export
simul_idx_matrix <- function(
  n_sims,
  alpha_lambdas,
  theta_lambdas,
  seeds,
  nfold
) {
  grid <- expand.grid(
    fold_idx   = seq_len(nfold),
    seed_idx   = seeds,
    theta_idx  = seq_len(length(theta_lambdas)),
    alpha_idx  = seq_len(length(alpha_lambdas)),
    sim_idx    = seq_len(n_sims),
    stringsAsFactors = FALSE
  )
  grid <- grid[, c("sim_idx","alpha_idx","theta_idx","seed_idx","fold_idx")]
  grid$lambda_alpha <- alpha_lambdas[grid$alpha_idx]
  grid$lambda_theta <- theta_lambdas[grid$theta_idx]
  return(as.matrix(grid))
}
