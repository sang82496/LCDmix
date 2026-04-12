# Generated from create-LCDmix.Rmd: do not edit by hand

#' @export
res_LCDmix_simul <- function(s){
  mat = summ[[s]]$reduced_mat
  print(summ[[s]]$opt_lambdas)
  sim = readRDS(paste0('sim_data/sim_', s, '.rds'))
  flow_bests = readRDS('flow_res_bests.rds')
  flow = flow_bests[[s]]
  res = refit_bests[[s]]
  X = res$X
  Y_bin = res$Y_bin
  bin_mass = res$bin_mass
  K = 2
  TT = length(Y_bin)
  
  
  
  ## 1. Lambda heatmap ##
  cv_grid <- as.data.frame(mat)
  best <- cv_grid %>% slice_max(cv_score, n = 1)
  
  ggplot(cv_grid, aes(x = lambda_alpha, y = lambda_theta, fill = cv_score)) +
    geom_tile() +
    geom_point(data = best, shape = 21, size = 3, color = "white", fill = "red") +
    geom_text(data = best, aes(label = sprintf("max=%.3f", cv_score)),
              nudge_x = 0, nudge_y = 0, vjust = -1, size = 3.2, color = "black") +
    scale_x_log10(labels = label_scientific(digits = 2)) +
    scale_y_log10(labels = label_scientific(digits = 2)) +
    scale_fill_gradient2(
      name = "cv_score",
      low = "yellow", mid = "white", high = "white",
      midpoint = max(cv_grid$cv_score)
    ) +
    labs(x = expression(lambda[alpha]), y = expression(lambda[theta]),
         title = "CV score over (lambda_alpha, lambda_theta)") +
    coord_fixed() +
    theme_minimal(base_size = 12) +
    theme(panel.grid.minor = element_blank())
  
  
  
  
  
  
  ## 2. Surrogate Loglikelihood ##
  is_increasing(res$iter$Q) 
  # Should be non‐decreasing if the EM update behaved correctly
  
  par(mfrow = c(1,2), pty = "s", cex = 0.7)
  plot(
    res$iter$Q, type = "b", xlab = "Iteration", ylab = "Surrogate Q",
    main = "Convergence of Q"
  )
  print("Surrogate Q values:")
  print(res$iter$Q)
  print(paste("The number of iterations:", res$iter$iter_num))
  
  
  
  
  
  ## 3. Parameter estimates ##
  #   theta : intercept + slopes;   alpha : mixture‐weight coefficients
  
  # true coefficients
  alpha_true = sim$alpha
  theta_true = sim$theta
  
  alpha_est_flow <- t(sweep(flow$alpha, 2, flow$alpha[1, ], FUN = "-")) # (p+1) × K
  alpha_est_lcd  <- t(res$iter$alpha_new)
  
  theta_est_flow <- do.call(cbind, flow$beta)
  theta_est_lcd  <- rbind(
    do.call(cbind, res$iter$theta0_new),  # intercepts
    do.call(cbind, res$iter$theta_new)    # slopes
  )
  
  rownames(theta_est_lcd) = rownames(alpha_true)
  colnames(theta_est_lcd) = colnames(alpha_true)
  
  rownames(theta_est_flow) = rownames(alpha_true)
  colnames(theta_est_flow) = colnames(alpha_true)
  
  rownames(theta_true) = rownames(alpha_true)
  colnames(theta_true) = colnames(alpha_true)
  
  rownames(alpha_est_lcd) = rownames(alpha_true)
  colnames(alpha_est_lcd) = colnames(alpha_true)
  
  rownames(alpha_est_flow) = rownames(alpha_true)
  colnames(alpha_est_flow) = colnames(alpha_true)
  
  rownames(alpha_true) = rownames(alpha_true)
  colnames(alpha_true) = colnames(alpha_true)
  
  
  # label matching
  ord1 = do.call(order, as.data.frame(t(theta_est_lcd)))
  theta_est_lcd = theta_est_lcd[, ord1]
  alpha_est_lcd = alpha_est_lcd[, ord1]
  
  ord2 = do.call(order, as.data.frame(t(theta_est_flow)))
  theta_est_flow = theta_est_flow[, ord2]
  alpha_est_flow = alpha_est_flow[, ord2]
  
  g_est = list()
  flow_mn = list()
  flow_sigma = rep(0, K)
  for (k in 1:K) {
    g_est[[k]] = res$iter$g_new[[ord1[k]]]
    flow_mn[[k]] = flow$mn[,,ord2[k]]
    flow_sigma[k] = flow$sigma[ord2[k],,1]
  }
  for (k in 1:K) {
    res$iter$theta0_new[[k]] = theta_est_lcd[1,k]
    res$iter$theta_new[[k]] = theta_est_lcd[-1,k]
    res$iter$alpha_new[k,] = alpha_est_lcd[,k]
    res$iter$g_new[[k]] = g_est[[k]]
    flow$alpha[k,] = alpha_est_flow[,k]
    flow$beta[[k]] = theta_est_flow[[k]]
    flow$mn[,,k] = flow_mn[[k]]
    flow$sigma[k,,1] = flow_sigma[k]
  }
  
  print("Theta estimates (rows=intercept+slopes, cols=components):")
  print(theta_est_flow)
  print(theta_est_lcd)
  print(theta_true)
  
  print("Theta sparsity (proportion of < 1e-6):")
  print(mean(abs(theta_est_flow[-1,]) < 1e-6)) # 0.35
  print(mean(abs(theta_est_lcd[-1,])  < 1e-6)) # 0.95
  print(mean(abs(theta_true[-1,])  < 1e-6))    # 0.9
  
  print("Alpha estimates (rows=components, cols=intercept+slopes):")
  print(alpha_est_flow)
  print(alpha_est_lcd)
  print(alpha_true)
  
  print("Alpha sparsity (proportion of < 1e-6):")
  print(mean(abs(as.numeric(alpha_est_flow[-1, -1, drop = FALSE])) < 1e-6, na.rm = TRUE)) # 0
  print(mean(abs(as.numeric(alpha_est_lcd[-1, -1, drop = FALSE])) < 1e-6, na.rm = TRUE))  # 1
  print(mean(abs(as.numeric(alpha_true[-1, -1, drop = FALSE])) < 1e-6, na.rm = TRUE))     # 0.9
  
  
  mean_curves_true = cbind(1, X) %*% theta_true
  mean_curves_flow <- cbind(1, X) %*% theta_est_flow
  mean_curves_lcd  <- cbind(1, X) %*% theta_est_lcd
  
  
  
  ## 4. Binned Y variables over time ##
  #  Heatmap of binned mass over time
  mids   <- sort(unique(unlist(res$Y_bin)))
  heat_mat <- matrix(0, nrow = length(mids), ncol = TT)
  for (t in seq_len(TT)) {
    rows <- match(res$Y_bin[[t]], mids)
    heat_mat[rows, t] <- res$bin_mass[[t]]
  }
  zmat <- t(heat_mat)  # transpose for image()
  layout(matrix(1:2, nrow = 1), widths = c(4, 1.5))
  image(
    1:TT, mids, zmat,
    col       = rev(heat.colors(100)),
    xlab      = "Time (t)",
    ylab      = "Response Y (bin midpoints)",
    main      = "Heatmap of bin_mass over time"
  #  useRaster = TRUE
  )
  par(mar = c(5, 2, 4, 4))    # bottom, left, top, right
  breaks <- seq(min(heat_mat), max(heat_mat), length.out = 100)
  image(
    x    = 1,               # dummy x position
    y    = breaks,          # vertical scale = your data range
    z    = matrix(breaks, nrow = 1, ncol = length(breaks)),
    col  = rev(heat.colors(100)),
    axes = FALSE, xlab = "", ylab = ""
  )
  axis(
    side = 4, at = pretty(breaks), labels = pretty(breaks),
    las = 1,cex.axis = 0.7)
  mtext("bin_mass", side = 4, line = 2.5, cex  = 0.7)
  layout(1)
  
  
  # Bubble size ~ bin_mass
  par(mfrow = c(1,2), pty = "s", cex = 0.7)
  plot(
    NA, xlim = c(1, TT), ylim = range(res$Y_bin),
    xlab = "t", ylab = "y",
    main = "Observed data (bubble ~ mass)"
  )
  for (t in seq_len(TT)) {
    points(
      rep(t, length(res$Y_bin[[t]])),
      res$Y_bin[[t]],
      cex = 10 * res$bin_mass[[t]] / sum(res$bin_mass[[t]])
    )
  }
  
  # only when K is correctly specified
  lines(1:TT, mean_curves_true[,1], col = 3, lwd = 2)
  lines((TT/2):TT, mean_curves_true[(TT/2):TT,2], col = 3, lwd = 2)
  
  
  
  
  ## 5. Mean‐curves ##
  
  ## flowmix
  par(mfrow = c(1,2), pty = "s", cex = 0.7)
  plot(
    NA, xlim = c(1, TT), ylim = range(res$Y_bin),
    xlab = "t", ylab = "y",
    main = "Flowmix mean curves"
  )
  for (t in seq_len(TT)) {
      points(
        rep(t, length(res$Y_bin[[t]])),
        res$Y_bin[[t]],
        cex = 10 * res$bin_mass[[t]] / sum(res$bin_mass[[t]])
      )
  }
  # true
  lines(1:TT, mean_curves_true[,1], col = 3, lwd = 2)
  lines((TT/2):TT, mean_curves_true[(TT/2):TT,2], col = 3, lwd = 2)
  # estimation
  for (k in seq_len(K)) {
    lines(seq_len(TT), mean_curves_flow[, k], lwd = 2, col = 2)
  }
  
  ## LCDmix
  plot(
    NA, xlim = c(1, TT), ylim = range(res$Y_bin),
    xlab = "t", ylab = "y",
    main = "LCDmix mean curves"
  )
  for (t in seq_len(TT)) {
      points(
        rep(t, length(res$Y_bin[[t]])),
        res$Y_bin[[t]],
        cex = 10 * res$bin_mass[[t]] / sum(res$bin_mass[[t]])
      )
  }
  # true
  lines(1:TT, mean_curves_true[,1], col = 3, lwd = 2)
  lines((TT/2):TT, mean_curves_true[(TT/2):TT,2], col = 3, lwd = 2)
  # estimation
  for (k in seq_len(K)) {
    lines(seq_len(TT), mean_curves_lcd[, k], lwd = 2, col = 2)
  }
  
  
  
  
  
  
  
  ## 6. Median curves ##
  median_curves_lcd  <- matrix(0, nrow = TT, ncol = K)
  
  for (k in 1:K) {
    for (t in 1:TT){
      resp = res$iter$resp_new[[t]][,ord1[k]]
      median_curves_lcd[t,k] = weighted_quantile(res$Y_bin[[t]] * resp, res$bin_mass[[t]] * resp, prob = 0.5)
    }
  }
  
  ## flowmix
  par(mfrow = c(1,2), pty = "s", cex = 0.7)
  plot(
    NA, xlim = c(1, TT), ylim = range(res$Y_bin),
    xlab = "t", ylab = "y",
    main = "Flowmix median curves"
  )
  for (t in seq_len(TT)) {
      points(
        rep(t, length(res$Y_bin[[t]])),
        res$Y_bin[[t]],
        cex = 10 * res$bin_mass[[t]] / sum(res$bin_mass[[t]])
      )
  }
  # estimation
  for (k in seq_len(K)) {
    lines(seq_len(TT), mean_curves_flow[, k], lwd = 2, col = 2)
  }
  
  ## LCDmix
  plot(
    NA, xlim = c(1, TT), ylim = range(res$Y_bin),
    xlab = "t", ylab = "y",
    main = "LCDmix median curves"
  )
  for (t in seq_len(TT)) {
      points(
        rep(t, length(res$Y_bin[[t]])),
        res$Y_bin[[t]],
        cex = 10 * res$bin_mass[[t]] / sum(res$bin_mass[[t]])
      )
  }
  # estimation
  for (k in seq_len(K)) {
    lines(seq_len(TT), median_curves_lcd[,k], lwd = 2, col = 2)
  }
  
  
  
  
  
  ## 7. Pi‐curves ##
  pi_flow  <- pi_k(X, t(alpha_est_flow))
  pi_lcd   <- pi_k(X, t(alpha_est_lcd))
  pi_true  <- pi_k(X, t(alpha_true))
  
  # flowmix
  par(mfrow = c(1,2), pty = "s", cex = 0.7)
  plot(
    1:TT, pi_flow[,1], type = "l", col = 2, ylim = c(0,1),
    xlab = "t", ylab = expression(pi[t,k]),
    main = "Flowmix pi vs t"
  )
  for (k in 1:K) lines(1:TT, pi_true[,k], col = 3)
  for (k in 2:K) lines(1:TT, pi_flow[,k], col = 2)
  
  # LCDmix
  plot(
    1:TT, pi_lcd[,1], type = "l", col = 2, ylim = c(0,1),
    xlab = "t", ylab = expression(pi[t,k]),
    main = "LCDmix pi vs t"
  )
  for (k in 1:K) lines(1:TT, pi_true[,k], col = 3)
  for (k in 2:K) lines(1:TT, pi_lcd[,k], col = 2)
  
  
  
  
  ## 8. experts (errors) estimates ##
  rg = range(g_est[[1]]$xn, g_est[[2]]$xn)
  gridk    <- seq(
    rg[1], rg[2], length = 200
  )
  err_true_den = err_true_fun(sim, gridk)
  
  
  par(mfrow = c(1,2), pty = "s", cex = 0.7)
  for (k in seq_len(K)) {
    flow_den = dnorm(gridk, 0, sqrt(flow_sigma[k]))
    LCD_den = logcondens::evaluateLogConDens(gridk, g_est[[k]])[,3]
    ymax = max(flow_den, LCD_den, err_true_den)
    
    par(mfrow = c(1,2), pty = "s", cex = 0.7)
    # flowmix
    plot(
      gridk, flow_den,
      type = "l", main = paste("flowmix expert", k),
      xlab = "residual", ylab = "density",
      ylim = c(0, ymax), lwd = 2
    )
    lines(gridk, err_true_den, col = 3, lwd = 2)
    legend("topright", legend = c("Flowmix", "Truth"),
    col = c(1, 3), lty = 1, lwd = 2, bty = "n")
  
    # LCDmix
    plot(
      gridk, LCD_den,
      type = "l", main = paste("LCDmix expert", k),
      xlab = "residual", ylab = "density",
      ylim = c(0, ymax), lwd = 2
    )
    lines(gridk, err_true_den, col = 3, lwd = 2)
    legend("topright", legend = c("LCDmix", "Truth"),
    col = c(1, 3), lty = 1, lwd = 2, bty = "n")
  }
  
  
  
  
  
  ## 9. Mixture density estimates ##
  times <- c(1, 51)
  par(mfrow = c(1, 2), pty = "s", cex = 0.7)
  for (t in times) {
    dens_flow  <- dens_est_fun(flow, t, Y_bin[[t]], X)
    mix_flow   <- dens_flow %*% pi_flow[t, ]
    dens_lcd  <- dens_est_fun(res$iter, t, Y_bin[[t]], X)
    mix_lcd   <- dens_lcd %*% pi_lcd[t, ]
    dens_true <- dens_true_fun(sim, t, Y_bin[[t]])
    mix_true  <- dens_true %*% pi_true[t, ]
    
    freq      <- res$bin_mass[[t]] / sum(res$bin_mass[[t]])
    mix_flow  <- mix_flow / sum(mix_flow)
    mix_lcd   <- mix_lcd  / sum(mix_lcd)
    mix_true <- mix_true / sum(mix_true)
    
    ymax = max(freq, mix_flow, mix_lcd)
    plot(
      res$Y_bin[[t]], freq, type = "b", main = paste("T =", t), 
      xlab = "Y", ylab = "density", ylim = c(0, ymax), col = 'white')
    lines(res$Y_bin[[t]], mix_flow, col = 4, lwd = 2)
    points(res$Y_bin[[t]], mix_true, col = 3)
    lines(res$Y_bin[[t]], mix_lcd,  col = 2,  lwd = 2)
    points(res$Y_bin[[t]], mix_true, col = 3)
    
    # place legend in upper‐right inside the plot
    legend("topright", legend = c("Flowmix", "LCDmix", "Truth"),
      col = c(4, 2, 3), lty = 1, lwd = 2, bty = "n")
  }
}
