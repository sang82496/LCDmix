# Generated from create-LCDmix.Rmd: do not edit by hand

#' Summarize CV results across simulations and pick best lambdas
#'
#' @description
#' For each unique simulation ID in \code{simul_idx_mat}, this helper subsets the
#' corresponding rows, calls \code{LCD_cv_summary()} to aggregate the CV results for
#' that simulation, and extracts the optimal penalty pair
#' \code{(lambda_alpha, lambda_theta)}. It returns a compact matrix with one row
#' per simulation.
#'
#' @param simul_idx_mat Numeric matrix of CV jobs (e.g., from
#'   \code{simul_idx_matrix()}). Must contain a \code{sim_idx} column; all other
#'   columns required by your \code{LCD_cv_summary()} are passed along after removing
#'   \code{sim_idx}.
#' @param save_dir Character path to the directory containing the per–CV-run
#'   result files produced by \code{run_cv_simul()}.
#'
#' @return A numeric matrix with three columns:
#' \describe{
#'   \item{\code{sim_idx}}{Simulation identifier.}
#'   \item{\code{lambda_alpha}}{Selected \code{lambda_alpha} for that simulation.}
#'   \item{\code{lambda_theta}}{Selected \code{lambda_theta} for that simulation.}
#' }
#'
#' @examples
#' \dontrun{
#' # Suppose you already ran:
#' # logs <- run_cv_simul(idx_mat, sim_dir="sim_data", K=2, save_dir="cv_saves")
#'
#' best_tbl <- summarize_simul_cv(
#'   simul_idx_mat = idx_mat,
#'   save_dir      = "cv_saves"
#' )
#' head(best_tbl)
#' }
#'
#' @export
summarize_simul_cv <- function(simul_idx_mat, save_dir) {
  sim_idxs <- unique(simul_idx_mat[,"sim_idx"])
  bests   <- matrix(0, nrow = length(sim_idxs), ncol = 3)
  colnames(bests) = c("sim_idx", "lambda_alpha", "lambda_theta")
  for (sim in sim_idxs) {
    idx_mat <- simul_idx_mat[simul_idx_mat[,"sim_idx"] == sim, -1]
    out     <- LCD_cv_summary(idx_mat, save_dir, simul = TRUE, sim_idx = sim)
    bests[sim,] <- c(sim, out$opt_lambdas[1], out$opt_lambdas[2])
  }
  return(bests)
}
