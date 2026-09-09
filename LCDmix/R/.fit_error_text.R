# Generated from create-LCDmix.Rmd: do not edit by hand

#' @keywords internal
#' 
#' @export
.fit_error_text <- function(fit, err_msg) {
  ## Error thrown out of main() -- tryCatch already captured it.
  if (!is.null(err_msg)) return(err_msg)
  ## main() returned normally with iter = NULL: iteration() caught it inside.
  if (is.list(fit) && !is.null(fit$iter_partial)) {
    e <- fit$iter_partial$error
    i <- fit$iter_partial$failed_iter
    if (!is.null(e))
      return(paste0(e, if (!is.null(i)) paste0(" [failed at iteration ", i, "]") else ""))
    return("main() returned iter = NULL but iter_partial$error was empty")
  }
  "unknown failure (fit was not a list, or had no iter_partial)"
}
