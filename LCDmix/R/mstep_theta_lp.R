# Generated from create-LCDmix.Rmd: do not edit by hand

#' M‐step update of intercept and slopes via linear programming under log-concavity
#'
#' @description
#' Updates the intercept \eqn{\theta_{0k}} and slope vector \eqn{\theta_k} for mixture
#' component \eqn{k} by solving a linear program that enforces the log‐concave density
#' constraints given by \code{density_k}, maximizes the weighted log‐likelihood, and
#' applies an L1 penalty on the slopes.
#'
#' @param Y_bin List of length \code{TT}; each element is an \eqn{M_t \times 1} matrix of binned response values at time \eqn{t}.
#' @param X Numeric \eqn{TT \times p} matrix of covariates (rows = time points).
#' @param weights List of length \code{TT}; each element is an \eqn{M_t \times K} matrix of posterior weights for each bin and component.
#' @param residuals List of length \code{TT}; each element is an \eqn{M_t \times K} matrix of residuals for each bin and component.
#' @param density_k An object returned by \code{modified_logcondens()} for component \eqn{k}, containing fields \code{x}, \code{phi}, and \code{IsKnot}.
#' @param idx List of length \code{TT}; each element is an \eqn{M_t \times K} logical matrix where \code{TRUE} indicates bins contributing to component \eqn{k}.
#' @param intercept_k Numeric scalar; current intercept parameter \eqn{\theta_{0k}}.
#' @param slopes_k Numeric vector of length \eqn{p}; current slope parameters \eqn{\theta_k}.
#' @param lambda_theta Nonnegative numeric; L1 penalty on the slope parameters.
#' @param component Integer in \(\{1,\dots,K\}\); index of the component to update.
#'
#' @return A list with components:
#' \describe{
#'   \item{theta0_k}{Numeric; updated intercept for component \eqn{k}.}
#'   \item{theta_k}{Numeric vector of length \eqn{p}; updated slopes for component \eqn{k}.}
#' }
#'
#' @export
mstep_theta_lp <- function(
  Y_bin,
  X,
  weights,
  residuals,
  density_k,
  idx,
  intercept_k,
  slopes_k,
  lambda_theta,
  component,
  lp_time_limit = 3600
) {
  TT <- length(Y_bin)
  p  <- ncol(X)

  # Collect residuals, weights, responses, and covariates for component
  res_k <- numeric(0)
  w_k   <- numeric(0)
  Y_k   <- numeric(0)
  X_k   <- matrix(nrow = 0, ncol = p)
  skip_ts <- integer(0)
  
  for (t in seq_len(TT)) {
    idx_tk <- idx[[t]][, component]
    if (any(idx_tk)) {
      res_k <- c(res_k, residuals[[t]][idx_tk, component])
      w_k   <- c(w_k, weights[[t]][idx_tk, component])
      Y_k   <- c(Y_k, Y_bin[[t]][idx_tk, 1])
      X_k   <- rbind(
                X_k,
                matrix(rep(X[t, ], sum(idx_tk)), nrow = sum(idx_tk), byrow = TRUE)
              )
    } else {
      skip_ts <- c(skip_ts, t)
    }
  }
  
  n <- length(Y_k)
  
  # Extract knots and slopes of piecewise linear density
  x_m     <- density_k$x[as.logical(density_k$IsKnot)]
  phi_m   <- density_k$phi[as.logical(density_k$IsKnot)]
  J       <- length(x_m) - 1
  b       <- diff(phi_m) / diff(x_m)
  beta0   <- b * x_m[-length(x_m)] - phi_m[-length(phi_m)]
  
  # Build constraint matrix and RHS
  const_mat <- NULL
  const_vec <- numeric(0)
  
  # Epigraph constraints for log‐concavity
  for (j in seq_len(J)) {
    tmp   <- cbind(Matrix::Diagonal(n), b[j], b[j] * X_k)
    block <- cbind(tmp, -tmp)
    const_mat <- rbind(const_mat, block)
    const_vec <- c(const_vec, b[j] * Y_k - beta0[j])
  }
  
  # Feasibility constraints: ensure predictions lie within range of residuals
  if (length(skip_ts) == 0) {
    TT_new <- TT
    X_new  <- X
  } else {
    TT_new <- TT - length(skip_ts)
    X_new  <- X[-skip_ts, , drop = FALSE]
  }
  tmp1 <- cbind(matrix(0, nrow = TT_new, ncol = n), 1, X_new)
  block1 <- cbind(tmp1, -tmp1)
  block2 <- cbind(-tmp1, tmp1)
  const_mat <- rbind(const_mat, block1, block2)
  
  # RHS for feasibility: bounding by min/max of Y_k
  L <- min(res_k); U <- max(res_k)
  tmp_vec <- numeric(2 * TT_new)
  cnt <- 1
  for (t in seq_len(TT)) {
    idx_tk <- idx[[t]][, component]
    if (any(idx_tk)) {
      tmp_vec[cnt]           <- min(Y_bin[[t]][idx_tk]) - L
      tmp_vec[TT_new + cnt]  <- U - max(Y_bin[[t]][idx_tk])
      cnt <- cnt + 1
    }
  }
  const_vec <- c(const_vec, tmp_vec)
  
  #–– Debugging: print size and memory usage of constraint matrix ––#
#  print(dim(const_mat))
#  print(format(object.size(as.matrix(const_mat)), "Gb"))
#  print(format(object.size(const_mat), "Mb"))
  
  # Objective: maximize w_k^T z - N*lambda_theta * |theta_k| - w_k^T z'
  N_total <- sum(unlist(weights))
  obj_coef <- c(
    w_k,
    0,
    rep(-N_total * lambda_theta, p),
    -w_k,
    0,
    rep(-N_total * lambda_theta, p)
  )
  
  const_dir <- rep("<=", nrow(const_mat))
  
  # Solve the linear program
  lp_res <- Rsymphony::Rsymphony_solve_LP(
    obj = obj_coef, 
    mat = const_mat, 
    dir = const_dir, 
    rhs = const_vec,  
    max = TRUE,
    time_limit = lp_time_limit)
  
  if (lp_res$status != 0) {
    print("No solution has been stored by Rsymphony. Change the LP solver to lpSolve")
    lp_res <- lpSolve::lp(
    direction    = "max",
    objective.in = obj_coef,
    const.mat    = const_mat,
    const.dir    = const_dir,
    const.rhs    = const_vec
  )
  }
  
  if (lp_res$status != 0) {
    stop("LP did not find an optimal solution (status = ", lp_res$status, ")")
  }
  
  sol <- lp_res$solution
  theta0_new <- sol[n + 1] - sol[2 * n + p + 2]
  theta_new  <- sol[(n + 2):(n + p + 1)] -
                sol[(2 * n + p + 3):(2 * (n + p + 1))]
  
  return(list(
    theta0_k = theta0_new,
    theta_k  = theta_new
  ))
}
