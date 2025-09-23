# Generated from create-LCDmix.Rmd: do not edit by hand

#' Compute L1 distance between true and estimated parameter matrices with label‐switching correction
#'
#' @description
#' For two‐component mixture models, calculates the L1 norm of the difference
#' between a true parameter matrix and an estimated one, taking into account
#' the possibility that the component labels may be swapped. Returns the smaller
#' of the direct L1 difference and the L1 difference under a swap of the two columns.
#'
#' @param true_params Numeric matrix of dimension \eqn{p \times 2}, where each column
#'   corresponds to the true parameters of one component.
#' @param est_params Numeric matrix of the same dimension \eqn{p \times 2}, where each
#'   column corresponds to the estimated parameters of one component.
#'
#' @return A single numeric: the minimal L1 distance between \code{true_params} and
#'   \code{est_params} under the two possible labelings.
#'
#' @examples
#' \dontrun{
#' # True slopes for two components (p = 3)
#' true_mat <- matrix(c(1,2,3,  -1,-2,-3), ncol = 2)
#' # Estimated with potential label swap
#' est_mat1 <- true_mat + 0.1  # small noise, same ordering
#' est_mat2 <- true_mat[, 2:1] + 0.2  # swapped + noise
#' d1 <- comp_l1_dist(true_mat, est_mat1)
#' d2 <- comp_l1_dist(true_mat, est_mat2)
#' stopifnot(d1 < d2)
#' }
#' @export
comp_l1_dist <- function(
  true_params,
  est_params
) {
  # Direct L1 difference
  diff_direct <- sum(abs(true_params - est_params))
  # L1 difference under swapping the two component columns
  diff_swapped <- sum(abs(true_params - est_params[, c(2, 1)]))
  # Return the minimal distance
  return(min(diff_direct, diff_swapped))
}
