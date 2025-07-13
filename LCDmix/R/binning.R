# Generated from create-LCDmix.Rmd: do not edit by hand

#' Bin responses by biomass weights into fixed or data-driven bins
#'
#' @description
#' Given a list of response vectors/matrices \code{Y} and corresponding biomass weights,
#' this function either treats them as already binned (when \code{n_bins = 0}) or
#' partitions the range of pooled responses into \code{n_bins} equal-width bins,
#' computing for each time point \(t\):
#' - \code{Y_bin[[t]]}: the bin centers with nonzero total weight  
#' - \code{bin_mass[[t]]}: the sum of biomass weights in each occupied bin  
#'
#' @param Y A list of length \code{TT}, where \code{Y[[t]]} is a numeric vector or single-column matrix of responses at time \(t\).
#' @param biomass A list of the same length \code{TT}, where \code{biomass[[t]]} is a numeric vector of weights for each observation in \code{Y[[t]]}.
#' @param n_bins Integer number of equal-width bins to create across the pooled range of \code{Y}.  
#'   If \code{n_bins = 0}, no binning is performed and each unique response value is treated as its own “bin.”  
#'   Default: \code{40}.
#'
#' @return A list with components:
#' \describe{
#'   \item{Y_bin}{List of length \code{TT}; each element is a single-column matrix of bin centers (numeric) with positive total weight.}
#'   \item{bin_mass}{List of length \code{TT}; each element is a numeric vector of total biomass weights per bin, aligned with \code{Y_bin[[t]]}.}
#' }
#'
#' @examples
#' \dontrun{
#' Y_list     <- list(rnorm(100), rnorm(150))
#' biomass_list <- lapply(Y_list, function(y) runif(length(y), 0.5, 2))
#' # Data-driven bins
#' binned <- binning(Y_list, biomass_list, n_bins = 20)
#' str(binned)
#' # No binning: treat each unique Y as its own bin
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
  list(
    Y_bin    = Y_bin,
    bin_mass = bin_mass
  )
}
