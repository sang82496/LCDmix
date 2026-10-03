# Generated from create-LCDmix.Rmd: do not edit by hand

#' Bin responses into equal-width bins and sum their weights
#'
#' @description
#' Pools the responses of all time points, cuts their range into
#' \code{n_bins} equal-width bins, and returns for each time point the centers
#' of the occupied bins and the total weight in each. All time points share
#' the same grid. With \code{n_bins = 0} there is no binning: each distinct
#' response value is its own bin, and the weights of equal values are summed.
#'
#' @param Y A list of length \code{TT}; \code{Y[[t]]} is a numeric vector or a
#'   one-column matrix of responses at time \eqn{t}.
#' @param biomass A list of length \code{TT}; \code{biomass[[t]]} holds one
#'   nonnegative weight per element of \code{Y[[t]]}.
#' @param n_bins Number of equal-width bins, or \code{0} for no binning.
#'
#' @return A list with components:
#' \describe{
#'   \item{Y_bin}{List of length \code{TT}; element \eqn{t} is an
#'     \eqn{M_t \times 1} matrix of the centers of the occupied bins, in
#'     increasing order.}
#'   \item{bin_mass}{List of length \code{TT}; element \eqn{t} is a numeric
#'     vector of length \eqn{M_t}, the total weight in each of those bins.}
#' }
#' The names of \code{Y} are kept when \code{n_bins > 0}.
#'
#' @examples
#' \dontrun{
#' Y_list       <- list(rnorm(100), rnorm(150))
#' biomass_list <- lapply(Y_list, function(y) runif(length(y), 0.5, 2))
#' # 20 equal-width bins
#' binned <- binning(Y_list, biomass_list, n_bins = 20)
#' str(binned)
#' # No binning: each distinct Y value is its own bin
#' binned2 <- binning(Y_list, biomass_list, n_bins = 0)
#' }
#' @export
binning <- function(
  Y,
  biomass,
  n_bins = 40
) {
  TT <- length(Y)
  
  # If no binning: each unique Y value becomes its own bin
  if (n_bins == 0) {
    Y_bin    <- vector("list", TT)
    bin_mass <- vector("list", TT)
    
    for (t in seq_len(TT)) {
      # Sum weights by unique response
      w_by_y <- tapply(biomass[[t]], factor(Y[[t]]), sum)
      centers   <- as.numeric(names(w_by_y))
      masses    <- as.numeric(w_by_y)
      
      Y_bin[[t]]    <- matrix(centers, ncol = 1)
      bin_mass[[t]] <- masses
    }
    
  } else {
    # Determine global response range and equal-width cutpoints
    pooled_range <- range(unlist(Y))
    cuts         <- seq(pooled_range[1], pooled_range[2], length = n_bins + 1)
    
    # Assign each Y to a bin index (1..n_bins)
    binned_idx <- lapply(
      Y,
      findInterval,
      vec = cuts,
      rightmost.closed = TRUE
    )
    binned_idx <- lapply(binned_idx, factor, levels = seq_len(n_bins))
    
    # Precompute bin centers
    bin_centers <- (cuts[-1] + cuts[-length(cuts)]) / 2
    
    Y_bin    <- vector("list", TT)
    bin_mass <- vector("list", TT)
    
    for (t in seq_len(TT)) {
      # Sum weights within each bin
      masses_t <- tapply(biomass[[t]], binned_idx[[t]], sum)
      
      # Keep only bins with non-NA (nonzero) mass
      occupied <- !is.na(masses_t)
      
      Y_bin[[t]]    <- matrix(bin_centers[occupied], ncol = 1)
      bin_mass[[t]] <- as.numeric(masses_t[occupied])
    }
    
    # Preserve any names on the input list
    names(Y_bin)    <- names(Y)
    names(bin_mass) <- names(Y)
  }
  
  # Return binned responses and corresponding biomass sums
  return(list(
    Y_bin    = Y_bin,
    bin_mass = bin_mass
  ))
}
