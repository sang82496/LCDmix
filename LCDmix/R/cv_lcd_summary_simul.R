# Generated from create-LCDmix.Rmd: do not edit by hand

#' Choose the penalty pair for each simulated dataset
#'
#' @description
#' Calls \code{cv_lcd_summary()} on \code{base_dir/sim_<s>/} for
#' \eqn{s = 1, \ldots,} \code{num_sims}, the layout that
#' \code{cv_lcd_simul()} writes. Folds get equal weight
#' (\code{equal_w = TRUE}).
#'
#' @param base_dir The \code{base_dir} given to \code{cv_lcd_simul()}.
#' @param num_sims Number of datasets.
#' @param cv_by_trimmed Passed to \code{cv_lcd_summary()}.
#'
#' @return A list of length \code{num_sims}; element \eqn{s} is the
#'   \code{cv_lcd_summary()} result of dataset \eqn{s}, so
#'   \code{out[[s]]$opt_lambdas} is its chosen pair.
#'
#' @export
cv_lcd_summary_simul <- function(base_dir, num_sims, cv_by_trimmed = T) {
  out <- vector("list", num_sims)
  for (s in seq_len(num_sims)) {
    sim_dir <- file.path(base_dir, sprintf("sim_%d", s))
    idx <- readRDS(file.path(sim_dir, "index_matrix.rds"))
    out[[s]] <- cv_lcd_summary(idx, save_dir = sim_dir, cv_by_trimmed = cv_by_trimmed)
  }
  return(out)
}
