# Generated from create-LCDmix.Rmd: do not edit by hand

#' Update mixture‐weight parameters via penalized multinomial regression
#'
#' @description
#' Fits a weighted multinomial logistic regression of component‐assignment weights
#' on covariates with an L1 penalty to update the mixture‐weight parameters \code{alpha}.
#'
#' @param X A numeric \eqn{TT \times p} matrix of covariates (rows = time points).
#' @param weights A list of length \eqn{TT}, each element an \eqn{M_t \times K}
#'   matrix of posterior weights (e.g.\ responsibilities × biomass) for each bin and component.
#' @param idx A list of length \eqn{TT}, each element an \eqn{M_t \times K}
#'   logical matrix indicating which bins have effectively nonzero posterior weight.
#' @param lambda_alpha Positive numeric L1 penalty on the non‐intercept mixture‐weight coefficients.
#'
#' @return A numeric \eqn{K \times (p+1)} matrix \code{alpha}, where each row corresponds to
#'   one mixture component, the first column is the intercept, and the remaining \eqn{p}
#'   columns are the slope coefficients.
#'
#' @examples
#' \dontrun{
#' TT <- 5; p <- 3; K <- 2
#' X <- matrix(rnorm(TT * p), nrow = TT, ncol = p)
#' # Simulate posterior weights and threshold mask
#' post_w <- lapply(1:TT, function(i) matrix(runif(4 * K), ncol = K))
#' mask  <- lapply(post_w, function(w) w > 0.1)
#' alpha <- mstep_alpha(X, post_w, mask, lambda_alpha = 1e-3)
#' }
#' @export
mstep_alpha <- function(
  X,
  weights,
  idx,
  lambda_alpha
) {
  TT <- nrow(X)
  K  <- ncol(idx[[1]])
  
  # Sequence of candidate lambda values (from large to target)
  lambda_max <- lambda_alpha * 100
  lambda_seq <- exp(seq(log(lambda_max), log(lambda_alpha), length.out = 30))
  
  # Sum posterior weights per time point and component
  weight_sum <- matrix(0, nrow = TT, ncol = K)
  for (t in seq_len(TT)) {
    for (k in seq_len(K)) {
      mask_tk         <- idx[[t]][, k]
      weight_sum[t, k] <- sum(weights[[t]][mask_tk, k])
    }
  }
  
  # Fit penalized multinomial regression
  fit <- glmnet::glmnet(
    x         = X,
    y         = weight_sum,
    lambda    = lambda_seq,
    family    = "multinomial",
    intercept = TRUE
  )
  
  # Extract coefficients at the requested lambda
  coefs_list <- glmnet::coef.glmnet(fit, s = lambda_alpha)
  alpha_mat  <- t(as.matrix(do.call(cbind, coefs_list)))
  alpha_mat  <- sweep(alpha_mat, 2, alpha_mat[1, ], FUN = "-")
  
  # Return K x (p+1) matrix: rows = components, cols = intercept + slopes
  return(alpha_mat)
}
