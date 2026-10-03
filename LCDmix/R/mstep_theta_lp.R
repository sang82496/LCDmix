# Generated from create-LCDmix.Rmd: do not edit by hand

#' Update the coefficients of one expert by linear programming
#'
#' @description
#' The fitted log-density \eqn{\hat g_k} of component \eqn{k} is concave and
#' piecewise linear, so maximizing
#' \deqn{\sum_i w_i \hat g_k(y_i - \theta_{0k} - x_i^\top \theta_k) - N \lambda_\theta \|\theta_k\|_1}
#' over \eqn{(\theta_{0k}, \theta_k)} is a linear program. The sum runs over the
#' active bins of component \eqn{k}, and \eqn{N} is the sum of all posterior
#' weights, so the objective is \eqn{N} times the part of \eqn{Q} that depends
#' on \eqn{\theta_k}. The linear program has
#'
#' - one variable \eqn{z_i} per active bin, with
#'   \eqn{z_i \le \phi_j + b_j (u_i - x_j)} for every linear piece \eqn{j}
#'   of \eqn{\hat g_k} (knot \eqn{x_j}, value \eqn{\phi_j}, slope \eqn{b_j});
#' - support constraints: at every time point with active bins, all residuals
#'   stay in \eqn{[L_k, U_k]}, the range of the current residuals, which is
#'   the support of \eqn{\hat g_k};
#' - optionally, the \code{maxdev} constraints;
#' - the objective \eqn{\sum_i w_i z_i - N \lambda_\theta \|\theta_k\|_1}.
#'
#' Every variable is split into a positive and a negative part. The program is
#' solved with \code{Rsymphony::Rsymphony_solve_LP()}; if that does not return
#' an optimal solution, with \code{lpSolve::lp()}. The function stops with an
#' error if neither does.
#'
#' @param Y_bin A list of length \code{TT}; each element is an
#'   \eqn{M_t \times 1} matrix of bin centers at time \eqn{t}.
#' @param X A numeric \eqn{TT \times p} covariate matrix.
#' @param weights A list of length \code{TT}; each element is an
#'   \eqn{M_t \times K} matrix of posterior weights.
#' @param residuals A list of length \code{TT}; each element is an
#'   \eqn{M_t \times K} matrix of current residuals. Only their range over the
#'   active bins of the component is used, as \eqn{[L_k, U_k]}.
#' @param density_k The \code{modified_logcondens()} fit of the component
#'   (fields \code{x}, \code{phi} and \code{IsKnot} are used).
#' @param idx A list of length \code{TT}; each element is an
#'   \eqn{M_t \times K} logical matrix of active bins. \code{iteration()}
#'   passes the active bins of the previous E-step, the bins that
#'   \code{density_k} was fitted to.
#' @param intercept_k Current intercept. Not used by the linear program; it is
#'   an argument so that both updates take the same arguments.
#' @param slopes_k Current slope vector of length \eqn{p}. Used only by the
#'   \code{maxdev} guard.
#' @param lambda_theta Nonnegative L1 penalty on the slopes.
#' @param component Index \eqn{k} of the component to update.
#' @param lp_time_limit Time limit in seconds for Rsymphony. If it is reached,
#'   lpSolve is tried without a time limit.
#' @param maxdev \code{NULL} (default, no constraint) or a positive number.
#'   Bounds the deviation of the component mean from its intercept,
#'   \eqn{|X_t^\top \theta_k| \le} \code{maxdev} at every time point \eqn{t},
#'   as in flowmix. If the current slopes already exceed \code{maxdev} at some
#'   \eqn{t}, the bound there is relaxed to the current deviation (the guard),
#'   so the linear program stays feasible and the update cannot decrease
#'   \eqn{Q}. The return value reports how often that happened.
#'
#' @return A list with components:
#' \describe{
#'   \item{theta0_k}{Updated intercept.}
#'   \item{theta_k}{Updated slope vector of length \eqn{p}.}
#'   \item{maxdev_n_relaxed}{Number of time points where the guard relaxed the
#'     bound (current deviation above \code{maxdev} by more than 1e-8);
#'     \code{NA} when \code{maxdev} is \code{NULL}.}
#'   \item{maxdev_max_excess}{Largest amount by which the current deviation
#'     exceeds \code{maxdev} (0 if never); \code{NA} when \code{maxdev} is
#'     \code{NULL}.}
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
  lp_time_limit,
  maxdev = NULL                     # NEW (fixP) - appended last
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
  
  # NEW (fixP): maximum-deviation constraints (flowmix), |X_t' theta_k| <= maxdev
  # for all t. Only the slopes enter, so mstep_shift() cannot break them later.
  # Guard: if the current slopes exceed maxdev at some t, the bound there is
  # the current deviation instead, so slopes_k stays feasible, the LP cannot
  # become infeasible, and the theta-block cannot decrease Q. The guard is
  # counted as active only above 1e-8, so rounding error is not counted.
  maxdev_n_relaxed  <- NA_integer_
  maxdev_max_excess <- NA_real_
  if (!is.null(maxdev)) {
    if (!is.numeric(maxdev) || length(maxdev) != 1L || !(maxdev > 0))
      stop("mstep_theta_lp(): maxdev must be NULL or a single positive number.")
    dev_now <- abs(as.vector(X %*% slopes_k))
    maxdev_n_relaxed  <- sum(dev_now > maxdev + 1e-8)
    maxdev_max_excess <- max(0, dev_now - maxdev)
    dev_rhs <- pmax(maxdev, dev_now)
    tmp2 <- cbind(matrix(0, nrow = TT, ncol = n + 1), X)
    const_mat <- rbind(const_mat, cbind(tmp2, -tmp2), cbind(-tmp2, tmp2))
    const_vec <- c(const_vec, dev_rhs, dev_rhs)
  }
  
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
    # fixS: lpSolve::lp() does not accept the sparse Matrix built above
    # ("argument is not a matrix"), and a dense copy can need gigabytes, so
    # pass the constraints as (row, column, value) triplets instead.
    trip <- Matrix::summary(Matrix::Matrix(const_mat, sparse = TRUE))
    lp_res <- lpSolve::lp(
      direction    = "max",
      objective.in = obj_coef,
      dense.const  = cbind(trip$i, trip$j, trip$x),
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
    theta_k  = theta_new,
    maxdev_n_relaxed  = maxdev_n_relaxed,     # NEW (fixP)
    maxdev_max_excess = maxdev_max_excess     # NEW (fixP)
  ))
}
