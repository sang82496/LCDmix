# Generated from create-LCDmix.Rmd: do not edit by hand

#' @export
gen_simul_laplace <- function(
  sim_seed       = NULL,
  nt             = 1000,
  TT             = 100,
  theta_par      = 0.5,
  p              = 10,
  B              = 30,
  scale          = 1,
  gap            = 4,
  sim_helper_dir = "."
) {

  ## Setup and basic checks
  assertthat::assert_that(nt %% 5 ==0)
  K = 2
  stopifnot(p >= 3)
  if(!is.null(sim_seed)) set.seed(sim_seed)
  ntlist = c(rep(0.8 * nt, TT/2), rep(nt, TT/2))
  noisetype = 'laplace'

  ## Generate covariate
  par = readRDS(file.path(sim_helper_dir, "simul_helper.rds"))

  Xrest = do.call(cbind, lapply(1:(p-2), function(ii) rnorm(TT)) )
  X = cbind(scale(par[1:TT]), c(rep(0, TT/2), rep(1, TT/2)), Xrest)
  colnames(X) = c("par", "cp", paste0("noise", 1:(p-2)))

  ## theta coefficients
  theta = matrix(0, ncol = K, nrow = p+1)
  theta[0+1,1] = 0
  theta[1+1,1] = theta_par
  theta[0+1,2] = gap
  theta[1+1,2] = -theta_par
  colnames(theta) = paste0("clust", 1:K)
  rownames(theta) = c("intercept", "par", "cp", paste0("noise", 1:(p-2)))

  ## alpha coefficients
  alpha = matrix(0, ncol = K, nrow = p+1)
  alpha[0+1, 2] = -10
  alpha[2+1, 2] = 10 + log(1/4)

  colnames(alpha) = paste0("clust", 1:K)
  rownames(alpha) = c("intercept", "par", "cp", paste0("noise", 1:(p-2)))

  ## Generate means and probabilities
  mnmat = cbind(1, X) %*% theta
  prob = pi_k(X, t(alpha))
  
  ## Samples |nt| memberships out of (1:K) according to the probs in prob.
  ## Data is a probabilistic mixture from these two means, over time.
  
  ylist = lapply(1:TT, function(tt){
     draws = sample(1:K, size = ntlist[tt], replace = TRUE,
                    prob = c(prob[tt,1], prob[tt,2]))
     mns = mnmat[tt,]
     means = mns[draws]
     ## Add noise to obtain data points.
     noise = rlaplace(ntlist[tt], 0, scale = scale)
     datapoints = means + noise
     cbind(datapoints)
   })

  # Binning
  biomass = vector("list", TT)
  biomass = lapply(1:TT, function(t){biomass[[t]] = rep(1, length(ylist[[t]]))})
  binned = LCDmix::binning(ylist, biomass, n_bins = B)

  return(list(Y_bin = binned$Y_bin, 
              X = X,
              bin_mass = binned$bin_mass,
              ## The true generating model:
              noisetype = noisetype,
              mnmat = mnmat,
              prob = prob,
              alpha = alpha,
              theta = theta,
              scale = scale
              ))
}
