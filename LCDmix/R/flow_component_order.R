# Generated from create-LCDmix.Rmd: do not edit by hand

#' @export
flow_component_order <- function(flow) {
  theta_mat <- do.call(cbind, flow$beta)
  return(do.call(order, as.data.frame(t(theta_mat))))
}
