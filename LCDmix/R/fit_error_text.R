# Generated from create-LCDmix.Rmd: do not edit by hand

#' Error message of a failed fit
#'
#' @description
#' Builds the error text that \code{cv_lcd_onejob()} and
#' \code{refit_onejob()} save when a fit fails. A fit fails in one of two
#' ways: \code{main()} stops with an error, which the caller catches
#' (\code{err_msg}); or, with \code{debug = TRUE}, \code{iteration()} catches
#' the error and \code{main()} returns normally with \code{iter = NULL}
#' (\code{fit$iter_partial$error}).
#'
#' @param fit The value returned by \code{main()}, or \code{NULL} if
#'   \code{main()} stopped with an error.
#' @param err_msg The caught error message, or \code{NULL}.
#'
#' @return A character string: \code{err_msg} if it is not \code{NULL};
#'   otherwise the error stored in \code{fit$iter_partial}, followed by
#'   " [failed at iteration i]"; otherwise a fixed text saying that the cause
#'   is unknown.
#' @keywords internal
#'
#' @export
fit_error_text <- function(fit, err_msg) {
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
