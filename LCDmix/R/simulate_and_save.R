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
#' @param sim_seeds      Integer vector of random seeds.
#' @param gaps           Numeric vector of change‐point intercept gaps.
#' @param skew_alphas    Numeric vector of skew‐normal shape parameters.
#' @param sim_dir        Character; directory to hold \code{sim_<i>.rds} files.
#'                       Will be created if it does not exist.  Default: \code{"sim_data"}.
#' @param nt             Number of observations per time‐point in the second half.
#'                       Passed to \code{generate_skewed_data()}. Default: 1000.
#' @param TT             Total number of time‐points to simulate. Default: 100.               
#' @param theta_par      The slope magnitude for the “baseline” covariate.
#'                       Default: 0.5.
#' @param p              Number of covariates (including baseline & change‐point).
#'                       Default: 10.
#' @param B              Number of histogram bins for each time‐point.
#'                       Default: 30.
#' 
#' @return A data.frame with columns \code{seed}, \code{intercept_gap},
#'   \code{skew_alpha}, where each row \code{i} corresponds to file
#'   \code{sim_dir/sim_i.rds}.
#'
#' @export
simulate_and_save <- function(
  sim_seeds,
  gaps,
  is_heavytail = FALSE,
  df = NULL,
  skew_alphas = NULL,
  sim_dir    = "sim_data",
  nt         = 1000,
  TT         = 100,
  theta_par  = 0.5,
  p          = 10,
  B          = 30,
  sim_helper_dir = '.'
) {
  if (!dir.exists(sim_dir)) dir.create(sim_dir, recursive = TRUE)
  
  if (is_heavytail){
    assertthat::assert_that(!is.null(df))
    jobs <- expand.grid(sim_seed = sim_seeds, gap = gaps, 
                          df = df, stringsAsFactors = F)
    for (i in seq_len(nrow(jobs))) {
      row  <- jobs[i, ]
      sim <- gen_simul_data(sim_seed = row$sim_seed, nt, TT, theta_par, p, B, is_heavytail = T, 
                            df = row$df, gap = row$gap, sim_helper_dir = sim_helper_dir)
      saveRDS(sim, file = file.path(sim_dir, paste0("sim_", i, ".rds")))
    }
    
  } else { # if skewed
    
    assertthat::assert_that(!is.null(skew_alphas))
    jobs <- expand.grid(sim_seed = sim_seeds, gap = gaps, 
                          skew_alpha = skew_alphas, stringsAsFactors = F)
    for (i in seq_len(nrow(jobs))) {
      row  <- jobs[i, ]
      sim <- gen_simul_data(sim_seed = row$sim_seed, nt, TT, theta_par, p, B, is_heavytail = F, 
                            skew_alpha = row$skew_alpha, gap = row$gap, sim_helper_dir = sim_helper_dir)
      saveRDS(sim, file = file.path(sim_dir, paste0("sim_", i, ".rds")))
    }
  }
  return(jobs)
}
