# Generated from create-LCDmix.Rmd: do not edit by hand

#' Generate synthetic time‐series data from a 2‐component skewed mixture model
#'
#' @description
#' Simulates a univariate response observed over \code{n_time} time points from a
#' two‐component mixture‐of‐experts model with skew‐normal noise.  At each time
#' point, a changing‐probability mixture of two regression lines generates
#' responses; the mixture weights vary with a binary “change‐point” covariate
#' and an optional smooth baseline trend.  The resulting continuous observations
#' are then binned into histograms for use with the LCDmix pipeline.
#'
#' @param seed Optional integer; random seed for reproducibility. Default: \code{NULL}.
#' @param n_per_time Integer; base number of observations per time point in the
#'   second half of the series.  The first half uses \code{0.8 * n_per_time}.
#'   Must be a multiple of 5. Default: \code{200}.
#' @param beta_par Numeric; slope coefficient on the baseline covariate for component 1
#'   (component 2 has intercept offset given by \code{intercept_gap}). Default: \code{0.5}.
#' @param n_covariates Integer ≥ 3; total number of covariates (including baseline and change‐point). Default: \code{3}.
#' @param grid_size Integer; number of grid cells per axis for binning the data. Default: \code{30}.
#' @param skew_alpha Numeric; shape parameter \(\alpha\) for the skew‐normal noise. Default: \code{10}.
#' @param intercept_gap Numeric; difference in intercepts between component 2 and 1. Default: \code{3}.
#' @param n_time Integer; total number of time points (must be even). Default: \code{100}.
#'
#' @return A list with components:
#' \describe{
#'   \item{\code{ylist}}{List of length \code{n_time}, each an \eqn{M_t \times 1}
#'     matrix of binned counts per time point.}
#'   \item{\code{X}}{\eqn{n_time \times n_covariates} covariate matrix:  
#'     column 1 = scaled smooth baseline trend,  
#'     column 2 = binary change‐point indicator (0 for first half, 1 for second),  
#'     columns 3:\(n_{\text{covariates}}\) = i.i.d.\ noise covariates.}
#'   \item{\code{countslist}}{List of length \code{n_time} of numeric vectors of raw
#'     histogram counts (before normalization).}
#'   \item{\code{mnmat}}{\eqn{n_time \times 2} matrix of true component means at each time.}
#'   \item{\code{prob}}{\eqn{n_time \times 2} matrix of true mixture probabilities.}
#'   \item{\code{mn}}{\eqn{n_time \times 1 \times 2} array of true means (for plotting consistency).}
#'   \item{\code{numclust}}{Integer number of mixture components (always 2).}
#'   \item{\code{alpha}}{\eqn{(n_{\text{covariates}}+1)\times 2} matrix of true mixture‐weight coefficients.}
#'   \item{\code{beta}}{\eqn{(n_{\text{covariates}}+1)\times 2} matrix of true regression coefficients.}
#'   \item{\code{sigma}}{\eqn{2\times1\times1} array of true noise variances for each component.}
#' }
#'
#' @examples
#' \dontrun{
#' # Simulate a small dataset
#' dat <- generate_skewed_data(
#'   seed          = 42,
#'   n_per_time    = 50,
#'   beta_par      = 1.0,
#'   n_covariates  = 4,
#'   grid_size     = 20,
#'   skew_alpha    = 5,
#'   intercept_gap = 2,
#'   n_time        = 20
#' )
#' str(dat)
#' }
#' @export
generate_skewed_data <- function(
  seed           = NULL,
  n_per_time     = 200,
  beta_par       = 0.5,
  n_covariates   = 3,
  grid_size      = 30,
  skew_alpha     = 10,
  intercept_gap  = 3,
  n_time         = 100
) {
  #— Argument checks —#
  assertthat::assert_that(n_per_time %% 5 == 0, msg = "`n_per_time` must be a multiple of 5")
  assertthat::assert_that(n_covariates >= 3, msg = "`n_covariates` must be at least 3")
  assertthat::assert_that(n_time %% 2 == 0, msg = "`n_time` must be even")
  
  # Set seed if requested
  if (!is.null(seed)) set.seed(seed)
  
  # Number of clusters
  numclust <- 2
  
  # Determine sample sizes per time point
  n_first_half  <- rep(0.8 * n_per_time, n_time / 2)
  n_second_half <- rep(n_per_time,     n_time / 2)
  n_list        <- c(n_first_half, n_second_half)
  
  #— Simulate covariates —#
  # Baseline trend: smooth version of a fixed sequence
  par_raw <- c(
  1807.884, 168.2681, 0.0006315789, -0.0139, -0.01336364, -0.014, -0.013125, -0.014,
  -0.0128, 0.218, 0.3867778, 23.893, 3320.086, 3278.974, 3584.604, 3769.884,
  3124.176, 3210.607, 2222.561, 2619.597, 1457.061, 87.64753, 0.03627778, -0.009875,
  -0.002125, 0.372, 0.526, 32.36786, 317.9361, 784.968, 1448.081, 1624.513,
  2499.522, 2059.033, 2005.622, 1719.437, 1811.938, 2868.744, 1851.804, 1101.775,
  646.0992, 27.92433, 0.6711875, -0.013, -0.01711765, -0.014125, -0.0134, -0.01175,
  0.01726316, 0.3675625, 0.4200588, 35.718, 272.8423, 655.5318, 872.1162, 1685.27,
  2138.353, 2799.593, 3552.781, 3003.885, 2739.476, 2525.467, 1645.232, 1257.148,
  322.8541, 54.98407, 1.811, -0.006733333, -0.009428571, -0.0069, -0.005692308, -0.0107,
  -0.009923077, 0.2817778, 1.216714, 76.83171, 211.934, 373.5088, 518.5455, 792.2462,
  931.5338, 1421.289, 913.1166, 841.0782, 922.5017, 788.3413, 685.737, 291.6708,
  122.616, 34.555, 0.422, -0.009333333, -0.0095, -0.01375, -0.01181818, -0.01192308,
  -0.009, -0.01041667, 0.6727857, 81.47508, 279.9459, 870.0117, 1304.603, 1713.278,
  2235.767, 1578.55, 1455.561, 1173.019, 386.7722, 57.29713, 0.8726, -0.006333333,
  -0.008833333, -0.008857143, -0.009375, -0.00925, -0.015, 0.3832, 0.4872308, 30.92275,
  213.5052, 462.9413, 692.3243, 1132.996, 2123.787, 2373.636, 2605.431, 2730.907,
  2340.096, 1557.706, 860.8426, 743.79, 315.8718, 83.95737, 0.2305, -0.005583333,
  -0.002230769, -0.003636364, -0.006181818, -0.001235294, 0.1075625, 0.3570769, 0.4778889, 49.719,
  354.4997, 817.8346, 1259.956, 2148.988, 2673.856, 2670.551, 2690.914, 2648.071,
  2759.417, 2783.573, 2852.975, 2533.883, 1810.079, 422.1133, 0.9876316, -0.0025,
  0.00006666667, -0.005777778, -0.004642857, -0.005333333, -0.0072, 0.2985, 0.4086667, 75.63286,
  353.0761, 724.7677, 1230.961, 1749.818, 2318.679, 2721.202, 3232.953, 3081.986,
  2675.911, 2377.742, 1600.023, 999.8954, 0.061, -0.01078571, -0.0075, -0.008666667,
  -0.008777778, -0.007166667, -0.008823529, -0.009875, 0.01844444, 26.37, 179.0059, 551.9136,
  1019.848, 1489.251, 1818.534, 2081.818, 2362.124, 2111.173, 1787.832, 1234.41,
  603.6168, 219.7619, 40.71359, 0.04358824, -0.00525, -0.005, -0.004, -0.00575,
  -0.008266667, -0.003777778, 0.3596875, 0.3949444, 42.80317, 262.8376, 800.8294, 1031.592,
  1414.766, 2902.38, 3331.281, 3102.55, 3314.92, 3258.176, 2502.952, 1913.286,
  908.4718, 317.1317, 44.07631, 0.01391667, -0.010875, -0.008, -0.005352941, -0.003705882,
  -0.003263158, 0.01389474, 0.3741875, 0.40085, 97.6352, 526.4492, 1717.746, 2481.61,
  2660.157, 3574.407, 3236.346, 3144.131, 3512.318, 3191.223, 3414.783, 2315.886,
  1544.328, 845.7028, 48.32487, 0.0098, -0.004, -0.0035, -0.006, -0.0053125,
  -0.00525, -0.0007272727, 0.251125, 0.3925625, 31.10895, 349.4825, 989.3238, 1673.217,
  3008.452, 3327.474, 3668.375, 3520.383, 3809.708, 3307.647, 3284.979, 3025.367,
  2197.985, 1805.234, 36.2584, -0.0096, -0.01033333, -0.006142857, -0.007416667, -0.0052,
  -0.007666667, -0.006375, -0.009923077, -0.005222222, 31.77713, 260.5398, 854.8976, 2440.643,
  3061.953, 3348.637, 3232.407, 3626.341, 4077.914, 2544.177, 1044.943, 464.6145,
  21.11, -0.01592308, -0.0135, -0.01288889, -0.01535714, -0.01515789, -0.01415, -0.01206667,
  -0.01313333, -0.01109091, 1818.87, 1880.722
  )
  par_smoothed <- stats::ksmooth(
    x        = seq_along(par_raw),
    y        = par_raw,
    bandwidth = 5,
    x.points  = seq_along(par_raw)
  )$y
  baseline <- scale(par_smoothed[1:n_time])
  
  # Change‐point indicator (0/1)
  cp <- c(rep(0, n_time / 2), rep(1, n_time / 2))
  
  # Additional noise covariates
  noise_covs <- replicate(
    n_covariates - 2,
    stats::rnorm(n_time),
    simplify = FALSE
  )
  X <- cbind(
    baseline,
    cp,
    do.call(cbind, noise_covs)
  )
  colnames(X) <- c("baseline", "cp", paste0("noise", seq_len(n_covariates - 2)))
  
  #— True model parameters —#
  # Regression coefficients for each component (intercept + covariates)
  beta <- matrix(0, nrow = n_covariates + 1, ncol = numclust)
  # Component 1
  beta[1, 1] <- 0
  beta[2, 1] <- beta_par
  # Component 2
  beta[1, 2] <- intercept_gap
  beta[2, 2] <- -beta_par
  # Remaining slopes left at zero
  rownames(beta) <- c("intercept", colnames(X))
  colnames(beta) <- paste0("clust", 1:numclust)
  
  # Mixture‐weight (alpha) coefficients
  alpha <- matrix(0, nrow = n_covariates + 1, ncol = numclust)
  # Shift mixture in second component
  alpha[1, 2] <- -10
  alpha[3, 2] <-  10 + log(1/4)
  rownames(alpha) <- c("intercept", colnames(X))
  colnames(alpha) <- paste0("clust", 1:numclust)
  
  # Compute true means and mixture probabilities over time
  design_mat <- cbind(1, X)
  mnmat      <- design_mat %*% beta
  pi_mat     <- exp(design_mat %*% alpha)
  pi_mat     <- pi_mat / rowSums(pi_mat)
  
  #— Simulate raw data and add skew‐normal noise —#
  ylist <- vector("list", n_time)
  for (t in seq_len(n_time)) {
    # Sample cluster memberships
    clust_ids <- sample(
      seq_len(numclust),
      size      = n_list[t],
      replace   = TRUE,
      prob      = pi_mat[t, ]
    )
    # True component means
    means_t <- mnmat[t, clust_ids]
    # Skew‐normal noise centered at zero
    omega   <- sqrt(1 / (1 - 2 * (1/pi) * skew_alpha^2 / (1 + skew_alpha^2)))
    shift   <- omega * skew_alpha / sqrt(1 + skew_alpha^2) * sqrt(2/pi)
    noise_t <- sn::rsn(
      n      = n_list[t],
      xi     = 0,
      omega  = omega,
      alpha  = skew_alpha
    ) - shift
    # Observations
    ylist[[t]] <- matrix(means_t + noise_t, ncol = 1)
  }
  
  #— Bin into histograms for LCDmix —#
  dat_grid <- flowmix::make_grid(ylist, gridsize = grid_size)
  binned   <- flowmix::bin_many_cytograms(
    ylist       = ylist,
    manual.grid = dat_grid,
    mc.cores    = 4,
    verbose     = TRUE
  )
  ylist         <- binned$ybin_list
  countslist    <- binned$counts_list
  
  #— Pack return object —#
  mn_array <- array(
    dim = c(n_time, 1, numclust),
    data = as.vector(mnmat)
  )
  sigma_array <- array(1, dim = c(numclust, 1, 1))
  
  return(list(
    ylist       = ylist,
    X           = X,
    countslist  = countslist,
    mnmat       = mnmat,
    prob        = pi_mat,
    mn          = mn_array,
    numclust    = numclust,
    alpha       = alpha,
    beta        = beta,
    sigma       = sigma_array
  ))
}
