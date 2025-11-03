# Generated from create-LCDmix.Rmd: do not edit by hand

#' @export
simulate_and_save_laplace <- function(
  sim_seeds,
  gaps,
  scale      = NULL,
  sim_dir    = "sim_data",
  nt         = 1000,
  TT         = 100,
  theta_par  = 0.5,
  p          = 10,
  B          = 30,
  sim_helper_dir = '.'
) {
  if (!dir.exists(sim_dir)) dir.create(sim_dir, recursive = TRUE)
  
  assertthat::assert_that(!is.null(scale))
  jobs <- expand.grid(sim_seed = sim_seeds, gap = gaps, 
                        scale = scale, stringsAsFactors = F)
  for (i in seq_len(nrow(jobs))) {
    row  <- jobs[i, ]
    sim <- gen_simul_laplace(sim_seed = row$sim_seed, nt, TT, theta_par, p, B, 
                             scale = row$scale, gap = row$gap, sim_helper_dir = sim_helper_dir)
    saveRDS(sim, file = file.path(sim_dir, paste0("sim_", i, ".rds")))
  }

  return(jobs)
}
