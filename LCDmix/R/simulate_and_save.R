# Generated from create-LCDmix.Rmd: do not edit by hand

#' Simulate many datasets and save them to disk
#'
#' @description
#' Calls \code{gen_simul_data()} once for every combination of
#' \code{sim_seeds}, \code{gaps} and, depending on \code{noisetype},
#' \code{skew_alphas} (\code{"skewed"}) or \code{df} (\code{"heavytail"}), and
#' saves dataset \eqn{i} as \code{sim_dir/sim_<i>.rds}.
#'
#' @param sim_seeds Seeds, one dataset each.
#' @param gaps Values of \code{gap}, the intercept of component 2.
#' @param df Degrees of freedom for \code{"heavytail"}; a vector gives one
#'   dataset per value.
#' @param skew_alphas Shape parameters for \code{"skewed"}; a vector gives one
#'   dataset per value.
#' @param sim_dir Output directory; it is created if needed.
#' @inheritParams gen_simul_data
#'
#' @return The data frame of combinations from \code{expand.grid()}, with
#'   columns \code{sim_seed} and \code{gap}, and \code{skew_alpha}
#'   (\code{"skewed"}) or \code{df} (\code{"heavytail"}). Row \eqn{i} belongs
#'   to the file \code{sim_<i>.rds}; \code{sim_seed} varies fastest.
#'
#' @export
simulate_and_save <- function(
  sim_seeds,
  gaps,
  noisetype  = 'gaussian',
  df         = NULL,
  skew_alphas = NULL,
  sim_dir    = "sim_data",
  nt         = 1000,
  TT         = 100,
  theta_par  = 0.5,
  p          = 10,
  B          = 30,
  sim_helper_dir = '.'
) {
  ## Setup and basic checks
  assertthat::assert_that(nt %% 5 == 0)
  assertthat::assert_that(noisetype %in% c('heavytail', 'skewed', 'laplace', 'exponential', 'gaussian'))

  if (!dir.exists(sim_dir)) dir.create(sim_dir, recursive = TRUE)
  
  if (noisetype == 'heavytail'){
    assertthat::assert_that(!is.null(df))
    jobs <- expand.grid(sim_seed = sim_seeds, gap = gaps, df = df, stringsAsFactors = F)
    for (i in seq_len(nrow(jobs))) {
      row  <- jobs[i, ]
      sim <- gen_simul_data(sim_seed = row$sim_seed, nt, TT, theta_par, p, B, noisetype = noisetype, 
                            df = row$df, gap = row$gap, sim_helper_dir = sim_helper_dir)
      saveRDS(sim, file = file.path(sim_dir, paste0("sim_", i, ".rds")))
    }
    
  } else if (noisetype == 'skewed'){
    assertthat::assert_that(!is.null(skew_alphas))
    jobs <- expand.grid(sim_seed = sim_seeds, gap = gaps, 
                          skew_alpha = skew_alphas, stringsAsFactors = F)
    for (i in seq_len(nrow(jobs))) {
      row  <- jobs[i, ]
      sim <- gen_simul_data(sim_seed = row$sim_seed, nt, TT, theta_par, p, B, noisetype = noisetype, 
                            skew_alpha = row$skew_alpha, gap = row$gap, sim_helper_dir = sim_helper_dir)
      saveRDS(sim, file = file.path(sim_dir, paste0("sim_", i, ".rds")))
    }
    
  } else if (noisetype == 'laplace'){
    jobs <- expand.grid(sim_seed = sim_seeds, gap = gaps, stringsAsFactors = F)
    for (i in seq_len(nrow(jobs))) {
      row  <- jobs[i, ]
      sim <- gen_simul_data(sim_seed = row$sim_seed, nt, TT, theta_par, p, B, noisetype = noisetype, 
                            gap = row$gap, sim_helper_dir = sim_helper_dir)
      saveRDS(sim, file = file.path(sim_dir, paste0("sim_", i, ".rds")))
    }
  
  } else if (noisetype == 'exponential'){
    jobs <- expand.grid(sim_seed = sim_seeds, gap = gaps, stringsAsFactors = F)
    for (i in seq_len(nrow(jobs))) {
      row  <- jobs[i, ]
      sim <- gen_simul_data(sim_seed = row$sim_seed, nt, TT, theta_par, p, B, noisetype = noisetype, 
                            gap = row$gap, sim_helper_dir = sim_helper_dir)
      saveRDS(sim, file = file.path(sim_dir, paste0("sim_", i, ".rds")))
    }
    
  } else { #gaussian
    jobs <- expand.grid(sim_seed = sim_seeds, gap = gaps, stringsAsFactors = F)
    for (i in seq_len(nrow(jobs))) {
      row  <- jobs[i, ]
      sim <- gen_simul_data(sim_seed = row$sim_seed, nt, TT, theta_par, p, B, noisetype = noisetype, 
                            gap = row$gap, sim_helper_dir = sim_helper_dir)
      saveRDS(sim, file = file.path(sim_dir, paste0("sim_", i, ".rds")))
    }
  }
  return(jobs)
}
