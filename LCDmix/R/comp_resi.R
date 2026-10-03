# Generated from create-LCDmix.Rmd: do not edit by hand

#' Compute residuals for binned responses given the expert coefficients
#'
#' @description
#' For each time point \eqn{t}, bin \eqn{i} and component \eqn{k}, computes
#' \deqn{u_{tik} = y_{ti} - \theta_{0k} - X_t^\top \theta_k.}
#'
#' @param Y_bin A list of length \code{TT}; each element is an
#'   \eqn{M_t \times 1} matrix of bin centers at time \eqn{t}.
#' @param X A numeric \eqn{TT \times p} matrix; row \code{X[t, ]} is the
#'   covariate vector at time \eqn{t}.
#' @param intercepts A list of length \code{K}; \code{intercepts[[k]]} is the
#'   intercept \eqn{\theta_{0k}}.
#' @param slopes A list of length \code{K}; \code{slopes[[k]]} is the slope
#'   vector \eqn{\theta_k} of length \eqn{p}.
#'
#' @return A list of length \code{TT}; element \eqn{t} is an
#'   \eqn{M_t \times K} matrix whose entry \code{[i, k]} is \eqn{u_{tik}}.
#'
#' @examples
#' \dontrun{
#' TT         <- 3; K <- 2; p <- 2; n_bins <- 5
#' Y_bin      <- lapply(1:TT, function(i) matrix(rnorm(n_bins), ncol = 1))
#' X          <- matrix(rnorm(TT * p), nrow = TT, ncol = p)
#' intercepts <- list(0.5, -1.2)
#' slopes     <- list(c(1.0, 0.5), c(-0.5, 2.0))
#' res_list   <- comp_resi(Y_bin, X, intercepts, slopes)
#' str(res_list)
#' }
#' @export
comp_resi <- function(
  Y_bin,
  X,
  intercepts,
  slopes
) {
  TT <- length(Y_bin)
  K  <- length(intercepts)
  residuals <- vector("list", TT)

  for (t in seq_len(TT)) {
    M_t <- nrow(Y_bin[[t]])
    # initialize M_t x K matrix for residuals at time t
    res_t <- matrix(0, nrow = M_t, ncol = K)

    for (k in seq_len(K)) {
      # predicted value for component k at time t
      pred_tk <- intercepts[[k]] + sum(X[t, ] * slopes[[k]])
      # residuals: observed minus predicted
      res_t[, k] <- as.numeric(Y_bin[[t]][, 1]) - pred_tk
    }

    residuals[[t]] <- res_t
  }

  return(residuals)
}
