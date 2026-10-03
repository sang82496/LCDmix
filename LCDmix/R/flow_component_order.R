# Generated from create-LCDmix.Rmd: do not edit by hand

#' Order of the components of a flowmix fit
#'
#' @description
#' The rule of \code{lcd_component_order()} applied to a flowmix fit: orders
#' the components by intercept, and breaks ties by the slopes in turn.
#'
#' @param flow A flowmix fit (field \code{beta}: a list of \eqn{K} vectors,
#'   intercept first).
#'
#' @return An integer vector of length \eqn{K}: element 1 is the index of the
#'   component with the smallest intercept, and so on.
#'
#' @export
flow_component_order <- function(flow) {
  theta_mat <- do.call(cbind, flow$beta)
  return(do.call(order, as.data.frame(t(theta_mat))))
}
