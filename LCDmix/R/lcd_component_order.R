# Generated from create-LCDmix.Rmd: do not edit by hand

#' Order of the components of an LCDmix fit
#'
#' @description
#' Orders the components by intercept, and breaks ties by the slopes in
#' turn. \code{flow_component_order()} applies the same rule to a flowmix
#' fit, and \code{res_LCDmix_simul()} uses it too, so that the components of
#' different fits are matched to the true components. In
#' \code{gen_simul_data()} component 1 has the smaller intercept.
#'
#' @param res The \code{iter} element of a \code{main()} fit (fields
#'   \code{theta0_new} and \code{theta_new}).
#'
#' @return An integer vector of length \eqn{K}: element 1 is the index of the
#'   component with the smallest intercept, and so on.
#'
#' @export
lcd_component_order <- function(res) {
  theta_mat <- rbind(unlist(res$theta0_new), do.call(cbind, res$theta_new))
  return(do.call(order, as.data.frame(t(theta_mat))))
}
