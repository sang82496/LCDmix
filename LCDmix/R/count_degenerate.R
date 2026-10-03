# Generated from create-LCDmix.Rmd: do not edit by hand

#' Count cross-validation fits that stopped at a degenerate component
#'
#' @description
#' Reads every cross-validation result file under \code{base_dir}
#' (recursively; files named \code{<a>-<t>-<s>-<f>.rds}, as written by
#' \code{cv_lcd_onejob()}) and counts the fits whose error message starts
#' with "degenerate component". \code{mstep_g()} gives this error when a
#' component has fewer than two distinct residuals. Prints the count and, if
#' there are any, how many happened at initialization and how many during the
#' EM iterations.
#'
#' @param base_dir Directory to search.
#'
#' @return Invisibly, a data frame with the file name and the error message of
#'   each degenerate fit.
#'
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
