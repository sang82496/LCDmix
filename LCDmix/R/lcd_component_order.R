# Generated from create-LCDmix.Rmd: do not edit by hand

#' Component order used throughout the pipeline: order by intercept, then by
#' slopes. Copied from the analysis function in create-LCDmix.Rmd so that the
#' two never drift apart.
#' @export
lcd_component_order <- function(res) {
  theta_mat <- rbind(unlist(res$theta0_new), do.call(cbind, res$theta_new))
  return(do.call(order, as.data.frame(t(theta_mat))))
}
