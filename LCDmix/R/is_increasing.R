# Generated from create-LCDmix.Rmd: do not edit by hand

#' Test whether a numeric vector is non‐decreasing
#'
#' @description
#' Checks if each element of the input vector is greater than or equal to the
#' previous element, i.e., verifies that the sequence is non‐decreasing.
#'
#' @param v A numeric vector.
#'
#' @return Logical scalar: \code{TRUE} if \code{v} is non‐decreasing (or of length ≤ 1), 
#'   otherwise \code{FALSE}.
#'
#' @examples
#' \dontrun{
#' is_increasing(c(1, 2, 2, 5))  # TRUE
#' is_increasing(c(3, 1, 4))     # FALSE
#' is_increasing(42)             # TRUE
#' }
#' @export
is_increasing <- function(v) {
  n <- length(v)
  if (n <= 1) {
    return(TRUE)
  }
  for (i in seq(2, n)) {
    if (v[i] < v[i - 1]) {
      return(FALSE)
    }
  }
  return(TRUE)
}
