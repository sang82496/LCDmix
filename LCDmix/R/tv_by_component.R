# Generated from create-LCDmix.Rmd: do not edit by hand

#' Total variation distance between estimated and true error densities
#'
#' @description
#' For each component, in the order of \code{lcd_component_order()} or
#' \code{flow_component_order()}, computes \code{tv_density()} between the
#' estimated error density and the true error density of \code{sim}. The
#' integration range comes from \code{tv_grid()} and covers the estimated
#' support (for flowmix, 6 standard deviations on each side of 0). A
#' component whose computation fails gets \code{NA}.
#'
#' @param sim A dataset from \code{gen_simul_data()}.
#' @param fit The \code{iter} element of an LCDmix fit, or a flowmix fit.
#' @param n Number of grid points for \code{tv_density()}.
#'
#' @return A numeric vector of length \eqn{K}, in the order of the true
#'   components.
#'
#' @export
tv_by_component <- function(sim, fit, n = 2001) {
  d_true <- function(v) LCDmix::err_true_fun(sim, v)
  q_true <- function(p) LCDmix::err_true_quantile(sim, p)
  is_lcd <- !is.null(fit$g_new)
  K <- if (is_lcd) length(fit$g_new) else dim(fit$mn)[3]

  ord <- if (is_lcd) lcd_component_order(fit) else flow_component_order(fit)

  return(vapply(seq_len(K), function(k) {
    kk <- ord[k]
    if (is_lcd) {
      g <- fit$g_new[[kk]]
      fhat <- lcd_err_density(g)
      sup  <- range(g$x)
    } else {
      fhat <- flow_err_density(as.numeric(fit$sigma)[kk])
      sd_k <- sqrt(as.numeric(fit$sigma)[kk])
      sup  <- c(-6 * sd_k, 6 * sd_k)
    }
    rg <- tv_grid(q_true, est_support = sup)
    return(tryCatch(tv_density(fhat, d_true, rg[1], rg[2], n = n),
                    error = function(e) NA_real_))
  }, numeric(1)))
}
