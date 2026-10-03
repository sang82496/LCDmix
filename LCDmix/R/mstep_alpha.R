# Generated from create-LCDmix.Rmd: do not edit by hand

#' Update the gate coefficients by penalized multinomial regression
#'
#' @description
#' Fits an L1-penalized multinomial logistic regression with \code{glmnet}.
#' The response at time \eqn{t} is the vector of summed posterior weights of
#' the active bins, \eqn{s_{tk} = \sum_{i \in I_{tk}} w_{tik}}; glmnet treats
#' each row as proportions with weight \eqn{\sum_k s_{tk}}. glmnet runs on a
#' path of 30 values from \code{100 * lambda_alpha} down to
#' \code{lambda_alpha}, and the coefficients at \code{lambda_alpha} are
#' returned, centered so that component 1 is the reference (row 1 is 0).
#'
#' If the proportions \eqn{s_{tk} / \sum_l s_{tl}} are the same at every time
#' point (to 1e-10), or if glmnet fails, the exact intercept-only solution is
#' returned: the log of the total weight of each component, with all slopes 0.
#'
#' glmnet uses its default \code{standardize = TRUE}, so the penalty acts on
#' the coefficients of the standardized covariates. See the Details of
#' \code{main()} on covariate scale.
#'
#' @param X A numeric \eqn{TT \times p} covariate matrix.
#' @param weights A list of length \code{TT}; each element is an
#'   \eqn{M_t \times K} matrix of posterior weights.
#' @param idx A list of length \code{TT}; each element is an
#'   \eqn{M_t \times K} logical matrix of active bins.
#' @param lambda_alpha Positive L1 penalty on the non-intercept gate
#'   coefficients.
#'
#' @return A numeric \eqn{K \times (p+1)} matrix: row \eqn{k} holds the
#'   intercept and the \eqn{p} slopes of component \eqn{k}; row 1 is 0.
#'
#' @examples
#' \dontrun{
#' TT <- 5; p <- 3; K <- 2
#' X <- matrix(rnorm(TT * p), nrow = TT, ncol = p)
#' post_w <- lapply(1:TT, function(i) matrix(runif(4 * K), ncol = K))
#' mask   <- lapply(post_w, function(w) w > 0.1)
#' alpha  <- mstep_alpha(X, post_w, mask, lambda_alpha = 1e-3)
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
  lambda_max <- lambda_alpha * 100
  lambda_seq <- exp(seq(log(lambda_max), log(lambda_alpha), length.out = 30))

  weight_sum <- matrix(0, nrow = TT, ncol = K)
  for (t in seq_len(TT)) {
    for (k in seq_len(K)) {
      mask_tk <- idx[[t]][, k]
      weight_sum[t, k] <- sum(weights[[t]][mask_tk, k])
    }
  }

  # NEW: exact MLE when class proportions do not vary across time points
  intercept_only <- function() {
    tot <- colSums(weight_sum)
    a0  <- log(pmax(tot, .Machine$double.xmin) / sum(tot))
    a   <- cbind(a0, matrix(0, K, ncol(X)))
    sweep(a, 2, a[1, ], FUN = "-")
  }
  rs <- rowSums(weight_sum)
  pr <- weight_sum[rs > 0, , drop = FALSE] / rs[rs > 0]
  if (!nrow(pr) || max(abs(sweep(pr, 2, colMeans(pr)))) < 1e-10) return(intercept_only())

  fit <- tryCatch(                                             # NEW
    glmnet::glmnet(x = X, y = weight_sum, lambda = lambda_seq,
                   family = "multinomial", intercept = TRUE),
    error = function(e) NULL)
  if (is.null(fit)) return(intercept_only())                  # NEW

  coefs_list <- glmnet::coef.glmnet(fit, s = lambda_alpha)
  alpha_mat  <- t(as.matrix(do.call(cbind, coefs_list)))
  alpha_mat  <- sweep(alpha_mat, 2, alpha_mat[1, ], FUN = "-")
  return(alpha_mat)
}
