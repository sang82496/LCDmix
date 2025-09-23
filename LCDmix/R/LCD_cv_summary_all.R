# Generated from create-LCDmix.Rmd: do not edit by hand

#' @export
LCD_cv_summary_all <- function(base_dir, S) {
  out <- vector("list", S)
  
  for (s in seq_len(S)) {
    sim_dir <- file.path(base_dir, sprintf("sim_%d", s))
    idx <- readRDS(file.path(sim_dir, "index_matrix.rds"))
    out[[s]] <- LCD_cv_summary(idx, save_dir = sim_dir, simul = FALSE)
  }
  return(out)
}
