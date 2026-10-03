# Generated from create-LCDmix.Rmd: do not edit by hand

#' Compute starting values from a Gaussian mixture of regressions
#'
#' @description
#' Computes the starting point of the EM-type algorithm:
#'
#' 1. Fits a Gaussian mixture of regressions with \code{flowmix::flowmix()}
#'    (one random start, \code{nrep = 1}).
#' 2. Takes the gate coefficients from that fit, centered so that component 1
#'    is the reference (row 1 is 0), and the intercepts and slopes.
#' 3. Computes the residuals and the starting responsibilities from the
#'    Gaussian component densities \eqn{N(0, \sigma_k^2)}, evaluated at the
#'    residuals and multiplied by the gate probabilities.
#' 4. Re-centers the intercepts with \code{mstep_shift()}.
#' 5. Fits the starting log-concave densities with \code{mstep_g()}.
#' 6. Computes the starting value of the surrogate objective \eqn{Q} with
#'    \code{comp_Q()}.
#'
#' The flowmix start is random. Call \code{set.seed()} first for reproducible
#' values.
#'
#' @param Y_bin A list of length \code{TT}; each element is an
#'   \eqn{M_t \times 1} matrix of bin centers at time \eqn{t}.
#' @param X A numeric \eqn{TT \times p} covariate matrix; row \code{X[t, ]}
#'   belongs to time \eqn{t}.
#' @param bin_mass A list of length \code{TT}; each element is a numeric vector
#'   of length \eqn{M_t}, the total weight (biomass) in each bin.
#' @param K Number of mixture components.
#' @param lambda_alpha L1 penalty on the gate coefficients. Passed to flowmix
#'   as \code{prob_lambda} and used in \eqn{Q}.
#' @param lambda_theta L1 penalty on the expert slopes. Passed to flowmix as
#'   \code{mean_lambda} and used in \eqn{Q}.
#' @param resp_threshold Bins whose responsibility for component \eqn{k} is
#'   above this value are active for \eqn{k} (\code{idx_init}).
#' @param maxdev \code{NULL} or a positive number; passed to
#'   \code{flowmix::flowmix()}. See \code{main()}.
#'
#' @return A list with components:
#' \describe{
#'   \item{flow}{The fit returned by \code{flowmix::flowmix()}.}
#'   \item{idx_init}{List of \code{TT} logical \eqn{M_t \times K} matrices;
#'     \code{TRUE} where \code{resp_init} is above \code{resp_threshold}.}
#'   \item{resp_init}{List of \code{TT} \eqn{M_t \times K} matrices of
#'     responsibilities. They are not thresholded; each row sums to 1.}
#'   \item{weight_init}{List of \code{TT} \eqn{M_t \times K} matrices of
#'     posterior weights, \code{resp_init} times \code{bin_mass}.}
#'   \item{theta0_init}{List of \code{K} intercepts, after \code{mstep_shift()}.}
#'   \item{theta_init}{List of \code{K} slope vectors, from flowmix.}
#'   \item{alpha_init}{\eqn{K \times (p+1)} matrix of gate coefficients; row 1
#'     is 0.}
#'   \item{resi_init}{List of \code{TT} \eqn{M_t \times K} matrices of
#'     residuals for \code{theta0_init} and \code{theta_init}.}
#'   \item{g_init}{List of \code{K} log-concave fits from \code{mstep_g()}.}
#'   \item{Q}{Starting value of \eqn{Q}: a number named \code{"Q"} that
#'     carries the attributes of \code{comp_Q()}.}
#'   \item{Q_every}{The same value; \code{iteration()} appends to it.}
#' }
#'
#' @examples
#' \dontrun{
#' TT       <- 4; K <- 2; p <- 3; n_bins <- 5
#' Y_bin    <- lapply(1:TT, function(t) matrix(rnorm(n_bins), ncol = 1))
#' bin_mass <- lapply(Y_bin, function(y) runif(nrow(y)))
#' X        <- matrix(rnorm(TT * p), nrow = TT, ncol = p)
#' init     <- initialization(
#'   Y_bin, X, bin_mass, K,
#'   lambda_alpha   = 1e-3,
#'   lambda_theta   = 1e-3,
#'   resp_threshold = 1e-3
#' )
#' str(init)
#' }
#' @export
initialization <- function(
  Y_bin,
  X,
  bin_mass,
  K,
  lambda_alpha,
  lambda_theta,
  resp_threshold,
  maxdev = NULL                         # NEW (fixP) - appended
) {
  # 1) Fit Gaussian mixture regression via flowmix
  flow_res <- flowmix::flowmix(
    ylist        = Y_bin,
    X            = X,
    countslist   = bin_mass,
    numclust     = K,
    prob_lambda  = lambda_alpha,
    mean_lambda  = lambda_theta,
    maxdev       = maxdev,              # fixP: was NULL
    nrep         = 1
  )
  
  message("✔ flowmix initialization complete")
  
  # 2) Extract initial parameters
  alpha_init  <- flow_res$alpha
  alpha_init  <- sweep(alpha_init, 2, alpha_init[1, ], FUN = "-")
  
  theta0_init <- lapply(flow_res$beta, `[[`, 1)
  theta_init  <- lapply(flow_res$beta, function(b) b[-1])
  
  # 3) Compute initial residuals
  resi_init <- comp_resi(Y_bin, X, theta0_init, theta_init)
  
  # 4) Compute initial responsibilities and weights
  TT           <- length(Y_bin)
  resp_init    <- vector("list", TT)
  idx_init     <- vector("list", TT)
  weight_init  <- vector("list", TT)
  pi_mat       <- pi_k(X, alpha_init)  # TT x K matrix
  
  for (t in seq_len(TT)) {
    M_t      <- nrow(Y_bin[[t]])
    likeli_t <- matrix(0, nrow = M_t, ncol = K)
    for (k in seq_len(K)) {
      likeli_t[, k] <- dnorm(
        x    = resi_init[[t]][, k],
        mean = 0,    # fixR: resi_init is already Y - mn; was mean = flow_res$mn[t, 1, k]
        sd   = sqrt(flow_res$sigma[k])
      ) * pi_mat[t, k]
    }
    resp_t           <- likeli_t / rowSums(likeli_t)
    resp_init[[t]]   <- resp_t
    idx_init[[t]]    <- resp_t > resp_threshold
    weight_init[[t]] <- resp_t * bin_mass[[t]]
  }
  
  # 5) One M-step update for intercepts
  theta0_init <- mstep_shift(Y_bin, X, weight_init, idx_init, theta_init)
  resi_init   <- comp_resi(Y_bin, X, theta0_init, theta_init)
  
  # 6) Initial log-concave density estimates
  g_init <- mstep_g(resi_init, weight_init, idx_init)
  
  # 7) Compute initial surrogate log-likelihood
  Q_init <- comp_Q(
    X, g_init, resi_init, theta_init,
    alpha_init, idx_init, weight_init,
    lambda_alpha, lambda_theta,
    Y_bin = Y_bin, intercepts = theta0_init          # NEW
    )
  names(Q_init) <- "Q"       # note: <- not =, and attributes survive naming
  
  # Return initialization results
  return(list(
    flow        = flow_res,
    idx_init    = idx_init,
    resp_init   = resp_init,
    weight_init = weight_init,
    theta0_init = theta0_init,
    theta_init  = theta_init,
    alpha_init  = alpha_init,
    resi_init   = resi_init,
    g_init      = g_init,
    Q           = Q_init,
    Q_every     = Q_init
  ))
}
