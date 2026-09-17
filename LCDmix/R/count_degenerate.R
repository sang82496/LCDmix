# Generated from create-LCDmix.Rmd: do not edit by hand

#' @export
count_degenerate <- function(base_dir = "cv_saves") {
  fs  <- list.files(base_dir, pattern = "^[0-9]+-[0-9]+-[0-9]+-[0-9]+\\.rds$",
                    recursive = TRUE, full.names = TRUE)
  msg <- vapply(fs, function(f) {
    m <- readRDS(f)$err_msg
    if (is.null(m) || is.na(m)) NA_character_ else m
  }, character(1))
  deg <- grepl("^degenerate component", msg)
  cat(sprintf("degenerate-component terminations: %d of %d cells (%.2f%%)\n",
              sum(deg), length(fs), 100 * mean(deg)))
  if (any(deg))
    print(table(stage = ifelse(grepl("failed at iteration", msg[deg]),
                               "iteration", "initialization")))
  invisible(data.frame(file = fs[deg], err_msg = msg[deg], stringsAsFactors = FALSE))
}
