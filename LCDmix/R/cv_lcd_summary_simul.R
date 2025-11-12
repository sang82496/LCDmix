# Generated from create-LCDmix.Rmd: do not edit by hand

#' @export
cv_lcd_summary_simul <- function(base_dir, num_sims, trimmed = T) {
  out <- vector("list", num_sims)
  for (s in seq_len(num_sims)) {
    sim_dir <- file.path(base_dir, sprintf("sim_%d", s))
    idx <- readRDS(file.path(sim_dir, "index_matrix.rds"))
    out[[s]] <- cv_lcd_summary(idx, save_dir = sim_dir, trimmed = trimmed)
  }
  return(out)
}
