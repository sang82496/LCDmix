# Generated from create-LCDmix.Rmd: do not edit by hand

#' Simulate multiple skew‐normal mixture datasets and save to disk
#'
#' @description
#' For each combination of \code{seed}, \code{intercept_gap}, and
#' \code{skew_alpha}, this function calls \code{generate_skewed_data()}
#' with the supplied simulation parameters and writes the result out
#' (one list per file) as \code{sim_<i>.rds} in \code{sim_dir}.  
#' It returns the data‐frame of all parameter combinations (in the same
#' row‐order as the saved files).
#'
#' @param seeds          Integer vector of random seeds.
#' @param intercept_gaps Numeric vector of change‐point intercept gaps.
#' @param skew_alphas    Numeric vector of skew‐normal shape parameters.
#' @param sim_dir       Character; directory to hold \code{sim_<i>.rds} files.
#'                       Will be created if it does not exist.  Default: \code{"sim_data"}.
#' @param n_per_time     Number of observations per time‐point in the second half.
#'                       Passed to \code{generate_skewed_data()}. Default: 1000.
#' @param beta_par       The slope magnitude for the “baseline” covariate.
#'                       Default: 0.5.
#' @param n_covariates   Number of covariates (including baseline & change‐point).
#'                       Default: 10.
#' @param grid_size      Number of histogram bins for each time‐point.
#'                       Default: 30.
#' @param n_time         Total number of time‐points to simulate. Default: 100.
#'
#' @return A data.frame with columns \code{seed}, \code{intercept_gap},
#'   \code{skew_alpha}, where each row \code{i} corresponds to file
#'   \code{sim_dir/sim_i.rds}.
#'
#' @export
simulate_and_save <- function(
  seeds,
  intercept_gaps,
  skew_alphas,
  sim_dir    = "sim_data",
  n_per_time  = 1000,
  beta_par    = 0.5,
  n_covariates= 10,
  grid_size   = 30,
  n_time      = 100
) {
  combos <- expand.grid(
    seed          = seeds,
    intercept_gap = intercept_gaps,
    skew_alpha    = skew_alphas,
    stringsAsFactors = FALSE
  )
  if (!dir.exists(sim_dir)) dir.create(sim_dir, recursive = TRUE)

  for (i in seq_len(nrow(combos))) {
    sc  <- combos[i, ]
    sim <- generate_skewed_data(
      seed          = sc$seed,
      intercept_gap = sc$intercept_gap,
      skew_alpha    = sc$skew_alpha,
      n_per_time    = n_per_time,
      beta_par      = beta_par,
      n_covariates  = n_covariates,
      grid_size     = grid_size,
      n_time        = n_time
    )
    saveRDS(sim, file = file.path(sim_dir, paste0("sim_", i, ".rds")))
  }
  return(combos)
}
