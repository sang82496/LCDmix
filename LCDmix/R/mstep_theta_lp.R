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
#' - one variable \eqn{s_t} per time point with active bins, the shift of the
#'   component at that time, tied to the coefficients by the equality
#'   \eqn{s_t = \theta_{0k} + X_t^\top \theta_k};
#' - one variable \eqn{z_i} per active bin, with
#'   \eqn{z_i \le \phi_j + b_j (u_i - x_j)} for every linear piece \eqn{j}
#'   of \eqn{\hat g_k} (knot \eqn{x_j}, value \eqn{\phi_j}, slope \eqn{b_j}),
#'   where \eqn{u_i = y_i - s_t} is the residual of bin \eqn{i} and \eqn{t}
#'   is its time point;
#' - support constraints: at every time point with active bins, all residuals
#'   stay in \eqn{[L_k, U_k]}, the range of the current residuals, which is
#'   the support of \eqn{\hat g_k}; these are a lower and an upper bound on
#'   \eqn{s_t};
#' - optionally, the \code{maxdev} constraints;
#' - the objective \eqn{\sum_i w_i z_i - N \lambda_\theta \|\theta_k\|_1}.
#'
#' The variables \eqn{s_t} are not needed to state the program: the
#' constraints of the bins could be written with
#' \eqn{(\theta_{0k}, \theta_k)} directly. With \eqn{s_t} each of those
#' constraints involves only \eqn{z_i} and \eqn{s_t}, not \eqn{z_i},
#' \eqn{\theta_{0k}} and all \eqn{p} slopes, and the program has the same
#' solutions. With many covariates this makes it much faster to solve.
#'
#' The slopes are split into a positive and a negative part, which makes the
#' L1 penalty linear; \eqn{z_i}, \eqn{s_t} and \eqn{\theta_{0k}} are free
#' variables. The program is solved with
#' \code{Rsymphony::Rsymphony_solve_LP()}. The function stops with an error
#' if the solver does not return an optimal solution, for example when it
#' reaches its time limit.
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
#'   the function stops with an error.
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

  # Time points with active bins of the component (the others are skipped),
  # and the residuals, weights and responses of those bins
  keep_t <- which(vapply(idx, function(m) any(m[, component]), logical(1)))
  if (length(keep_t) == 0L) {
    stop("mstep_theta_lp(): component ", component, " has no active bins")
  }
  n_t   <- vapply(idx[keep_t], function(m) sum(m[, component]), numeric(1))
  res_k <- unlist(lapply(keep_t, function(t) residuals[[t]][idx[[t]][, component], component]))
  w_k   <- unlist(lapply(keep_t, function(t) weights[[t]][idx[[t]][, component], component]))
  Y_k   <- unlist(lapply(keep_t, function(t) Y_bin[[t]][idx[[t]][, component], 1]))
  Y_min <- vapply(keep_t, function(t) min(Y_bin[[t]][idx[[t]][, component], 1]), numeric(1))
  Y_max <- vapply(keep_t, function(t) max(Y_bin[[t]][idx[[t]][, component], 1]), numeric(1))
  X_keep <- X[keep_t, , drop = FALSE]

  n      <- length(Y_k)               # number of active bins
  T_k    <- length(keep_t)            # number of time points with active bins
  time_i <- rep(seq_len(T_k), n_t)    # for bin i, the position of its time point in keep_t

  # Extract knots and slopes of piecewise linear density
  x_m     <- density_k$x[as.logical(density_k$IsKnot)]
  phi_m   <- density_k$phi[as.logical(density_k$IsKnot)]
  J       <- length(x_m) - 1
  b       <- diff(phi_m) / diff(x_m)
  beta0   <- b * x_m[-length(x_m)] - phi_m[-length(phi_m)]

  # fixAA: the shift s_t = theta0 + X_t' theta of each time point is a variable.
  # Columns:  z (n) | s (T_k) | theta0 (1) | theta+ (p) | theta- (p)
  # z, s and theta0 are free; theta+ and theta- are >= 0.
  n_free <- n + T_k + 1
  col_s  <- n + seq_len(T_k)
  col_t0 <- n + T_k + 1
  col_tp <- n_free + seq_len(p)
  col_tm <- n_free + p + seq_len(p)
  n_col  <- n_free + 2 * p

  # The constraints are collected as (row, column, value) triplets.
  # Epigraph constraints for log-concavity, one per (active bin, linear piece):
  #   z_i + b_j * s_t <= b_j * y_i - beta0_j
  i_e     <- rep(seq_len(n), times = J)
  j_e     <- rep(seq_len(J), each = n)
  row_epi <- seq_len(n * J)
  # Link constraints, one per time point: s_t - theta0 - X_t' (theta+ - theta-) = 0
  row_link <- n * J + seq_len(T_k)

  tri_i <- c(row_epi, row_epi,
             row_link, row_link, rep(row_link, times = p), rep(row_link, times = p))
  tri_j <- c(i_e, col_s[time_i[i_e]],
             col_s, rep(col_t0, T_k), rep(col_tp, each = T_k), rep(col_tm, each = T_k))
  tri_x <- c(rep(1, n * J), b[j_e],
             rep(1, T_k), rep(-1, T_k), -as.vector(X_keep), as.vector(X_keep))
  const_vec <- c(b[j_e] * Y_k[i_e] - beta0[j_e], rep(0, T_k))
  const_dir <- c(rep("<=", n * J), rep("==", T_k))
  n_row <- n * J + T_k

  # Feasibility constraints: the residuals of each time point stay within the
  # range [L, U] of the current residuals. These are bounds on the shift:
  #   max(Y_tk) - U <= s_t <= min(Y_tk) - L
  L <- min(res_k); U <- max(res_k)
  s_lo <- Y_max - U
  s_hi <- Y_min - L
  # The current shift satisfies both bounds, because L and U are the smallest
  # and the largest current residual. When the active bins of one time point
  # span the whole support the two bounds are equal, and rounding can put them
  # in the wrong order by about 1e-16. Rsymphony rejects crossed bounds, so
  # bounds that cross by rounding only are set to their midpoint.
  crossed <- s_lo > s_hi
  if (any(crossed)) {
    if (max(s_lo[crossed] - s_hi[crossed]) > 1e-8) {
      stop("mstep_theta_lp(): the support constraints of component ", component, " cannot be satisfied")
    }
    s_mid <- (s_lo[crossed] + s_hi[crossed]) / 2
    s_lo[crossed] <- s_mid
    s_hi[crossed] <- s_mid
  }

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
    row_dev <- n_row + seq_len(TT)
    tri_i <- c(tri_i, rep(row_dev, times = p), rep(row_dev, times = p),
               rep(TT + row_dev, times = p), rep(TT + row_dev, times = p))
    tri_j <- c(tri_j, rep(col_tp, each = TT), rep(col_tm, each = TT),
               rep(col_tp, each = TT), rep(col_tm, each = TT))
    tri_x <- c(tri_x, as.vector(X), -as.vector(X), -as.vector(X), as.vector(X))
    const_vec <- c(const_vec, dev_rhs, dev_rhs)
    const_dir <- c(const_dir, rep("<=", 2 * TT))
    n_row <- n_row + 2 * TT
  }

  keep  <- tri_x != 0
  tri_i <- tri_i[keep]; tri_j <- tri_j[keep]; tri_x <- tri_x[keep]
  const_mat <- Matrix::sparseMatrix(i = tri_i, j = tri_j, x = tri_x, dims = c(n_row, n_col))

  # Objective: maximize w_k^T z - N*lambda_theta * (theta+ + theta-)
  N_total <- sum(unlist(weights))
  obj_coef <- numeric(n_col)
  obj_coef[seq_len(n)] <- w_k
  obj_coef[c(col_tp, col_tm)] <- -N_total * lambda_theta

  # Solve the linear program. z and theta0 have no lower bound; s_t lies
  # between s_lo and s_hi; theta+ and theta- keep the default bound, >= 0.
  col_free <- c(seq_len(n), col_t0)
  lp_res <- Rsymphony::Rsymphony_solve_LP(
    obj = obj_coef,
    mat = const_mat,
    dir = const_dir,
    rhs = const_vec,
    bounds = list(lower = list(ind = c(col_free, col_s), val = c(rep(-Inf, length(col_free)), s_lo)),
                  upper = list(ind = col_s, val = s_hi)),
    max = TRUE,
    time_limit = lp_time_limit)

  # No second solver: a fit whose program is not solved to optimality stops here.
  if (lp_res$status != 0) {
    stop("LP did not find an optimal solution (status = ", lp_res$status, ")")
  }

  sol <- lp_res$solution
  theta0_new <- sol[col_t0]
  theta_new  <- sol[col_tp] - sol[col_tm]

  return(list(
    theta0_k = theta0_new,
    theta_k  = theta_new,
    maxdev_n_relaxed  = maxdev_n_relaxed,     # NEW (fixP)
    maxdev_max_excess = maxdev_max_excess     # NEW (fixP)
  ))
}
