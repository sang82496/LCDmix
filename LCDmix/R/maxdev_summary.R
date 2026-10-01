# Generated from create-LCDmix.Rmd: do not edit by hand

#' Summarize how often the maxdev guard was active in one fit
#'
#' @description
#' Reads the per-LP-call records that \code{mstep_theta()} stores in
#' \code{iter$theta_diag}. The guard is active in an LP call when the current
#' slopes exceed \code{maxdev} at one or more time points (by more than 1e-8).
#'
#' @param iter The \code{iter} element of a \code{main()} fit.
#'
#' @return A named numeric vector: \code{n_lp_calls} (LP calls, one per
#'   iteration and component), \code{n_calls_active} (calls with the guard
#'   active), \code{n_relaxed_total} (time points relaxed, summed over calls),
#'   and \code{max_excess} (largest excess over \code{maxdev}). All but
#'   \code{n_lp_calls} are \code{NA} when \code{maxdev} was \code{NULL}.
#' @export
maxdev_summary <- function(iter) {
  d <- unlist(iter$theta_diag, recursive = FALSE)   # one entry per (iteration, component)
  n_rel  <- vapply(d, function(z) if (is.null(z$maxdev_n_relaxed)) NA_real_ else as.numeric(z$maxdev_n_relaxed), numeric(1))
  excess <- vapply(d, function(z) if (is.null(z$maxdev_max_excess)) NA_real_ else as.numeric(z$maxdev_max_excess), numeric(1))
  if (!length(d) || all(is.na(n_rel))) {
    return(c(n_lp_calls = length(d), n_calls_active = NA, n_relaxed_total = NA, max_excess = NA))
  }
  return(c(n_lp_calls      = length(d),
           n_calls_active  = sum(n_rel > 0, na.rm = TRUE),
           n_relaxed_total = sum(n_rel, na.rm = TRUE),
           max_excess      = max(excess, na.rm = TRUE)))
}
