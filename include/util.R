
# - 
mixture_normal <- function(N, mu, sigma, p)
{
  K = length(mu)
  out = numeric(N)
  
  for (nn in 1:N) {
    z = sample(1:K, 1, prob = p, replace = FALSE)
    out[nn] = rnorm(1,  mu[z], sigma[z])
  }
  return(out)
}


# circadian parametric function
circadian <- function(tt, R, theta)
{
  m = theta$m
  a = theta$a
  alpha = theta$alpha
  beta = theta$beta
  phi = theta$phi
  
  out = m + a*expit(beta*(cos((tt/R - phi) * (2*pi)/24) - alpha))
  
  return(out)
}



getPosteriorPredictive_covariates <- function(parms, sbj, T_all, R,
                                              n_draws = 50, theta = FALSE)
{
  m_sample = parms$m
  beta_sample = parms$beta
  sigma_sample = parms$sigma
  alpha_sample = parms$alpha
  
  if (theta == TRUE) {
    a_sample = parms$theta_a
    phi_sample = parms$theta_phi
  } else {
    a_sample = parms$a 
    phi_sample = parms$phi
  }

  
  N = T_all[sbj]
  
  n_MCMC = dim(a_sample)[1]
  y_hat = matrix(NA, n_draws, N)
  idxs_MCMC = sample(1:n_MCMC, n_draws)
  
  for (ii in 1:n_draws) {
    
    jj = idxs_MCMC[ii]
    theta = list()
    
    theta$m = m_sample[jj, sbj]
    theta$alpha = alpha_sample[jj, sbj]
    theta$beta = beta_sample[jj, sbj]
    theta$a = a_sample[jj, sbj]
    theta$phi = phi_sample[jj, sbj]
    sigma = sigma_sample[jj, sbj]
    
    mu = sapply(1:N, circadian, R, theta)
    y_hat[ii, ] = rnorm(N, mean = mu, sd = sigma)
    
  }
  
  return(y_hat)
  
}

getPosteriorPredictive <- function(parms, N, R, n_draws = 50)
{
  m_sample = parms$m
  beta_sample = parms$beta
  phi_sample = parms$phi
  sigma_sample = parms$sigma
  a_sample = parms$a 
  alpha_sample = parms$alpha
  
  n_MCMC = length(m_sample)
  y_hat = matrix(NA, n_draws, N)
  idxs_MCMC = sample(1:n_MCMC, n_draws)
  
  for (ii in 1:n_draws) {
    
    jj = idxs_MCMC[ii]
    
    theta = list()
    theta$m = m_sample[jj]
    theta$a = a_sample[jj]
    theta$alpha = alpha_sample[jj]
    theta$beta = beta_sample[jj]
    theta$phi = phi_sample[jj]
    
    mu = sapply(1:N, circadian, R, theta)
    y_hat[ii, ] = rnorm(N, mean = mu, sd = sigma_sample[jj])
    
  }
  return(y_hat)
}


plotPosteriorPredictive <- function(obs_single, y_hat, plt_main = "")
{
  cols = c("gray95", "gray90", "gray85", "gray80", "gray75")
  plt_pch = 20
  plt_cex = 1.0
  
  N = length(obs_single)
  
  probs <- seq(from=0.1, to=0.9, by=0.1)
  cred <- sapply(1:N, function(t) quantile(y_hat[, t], probs = probs))
  
  plot(1:N, type = "n",  
       ylab = "log(activity+1)", 
       main = plt_main, 
       xlim = c(1, N),
       ylim = c(0, 8),
       cex.lab = 1.1, xlab = "")
  polygon(c(1:N, rev(1:N)), c(cred[1,], rev(cred[9,])),
          col = scales::alpha(cols[1], 0.9), border = NA)
  polygon(c(1:N, rev(1:N)), c(cred[2,], rev(cred[8,])),
          col = scales::alpha(cols[2], 0.9), border = NA)
  polygon(c(1:N, rev(1:N)), c(cred[3,], rev(cred[7,])),
          col = scales::alpha(cols[3], 0.9), border = NA)
  polygon(c(1:N, rev(1:N)), c(cred[4,], rev(cred[6,])),
          col = scales::alpha(cols[4], 0.9), border = NA)
  
  points(1:N, obs_single, col=scales::alpha("white", .9),
         pch=plt_pch, cex= (plt_cex + 0.2))
  points(1:N, obs_single, pch = plt_pch,
         cex = plt_cex, col =  scales::alpha("black", .7))
  
  
  lines(1:N, cred[5,], col= scales::alpha("red", 0.8), lwd=3)
}

# = 
DExp_upperQuantile_MC <- function(X, N_MC, Q, sparsity_level, sigma_x) 
{
  theta <- rep(NA, Q)
  out <- rep(NA, N_MC)
  N_sbj <- nrow(X)
  
  magnitude_all <- array(NA, c(N_sbj, N_MC))
  
  for (ii in 1:N_MC) {
    zeros <- sample(c(0, 1), Q, replace = TRUE, 
                    prob = c(1 - sparsity_level, sparsity_level))
    eta_tmp <- (1-zeros)*rlaplace(Q, 0, sigma_x)
    
    for (nn in 1:N_sbj) {
      magnitude_all[nn, ii] = as.numeric(exp(X[nn, ] %*% eta_tmp))
    }
    out[ii] <- quantile(magnitude_all[, ii], probs = c(0.95))
  }
  return(out)
}


# - 
l1_ball_projection <- function(beta, r){
  p <- length(beta)
  abs_beta <- abs(beta)
  abs_beta_sort_index <- order(abs_beta, decreasing = TRUE)
  mu <- pmax(cumsum(abs_beta[abs_beta_sort_index]) - r, 0)
  c <- max(which((abs_beta[abs_beta_sort_index] > mu/(1:p)) == TRUE))
  theta <- sign(beta)*pmax(abs(beta) - mu[c]/c, 0)
  return(theta)
}


# 
generate_theta_l1ball <- function(p, r_alpha, beta_gen, N_MC = 5000){
  # r ~ Exp(r_alpha)
  theta <- matrix(NA, nrow = N_MC, ncol = p)
  for(j in 1:N_MC){
    r <- rexp(1, r_alpha)
    beta <- beta_gen(p)
    theta[j, ] <- l1_ball_projection(beta, r)
  }
  return(theta)
}

# - 
setHyperametersCALC <- function(X, settingsHyperParmsSearch)
{
  Q = dim(X)[2]
  ampl_upper <- settingsHyperParmsSearch$ampl_upper
  phase_upper <- settingsHyperParmsSearch$phase_upper 
  sparsity_level <- settingsHyperParmsSearch$sparsity_level
  N_MC <- settingsHyperParmsSearch$N_MC
  n_grid <- settingsHyperParmsSearch$n_grid
  grid_upper_sigma <- settingsHyperParmsSearch$grid_upper_sigma
  grid_upper_lambda <- settingsHyperParmsSearch$grid_upper_lambda
  
  sigma_a_grid <- seq(from = 0.1, to = grid_upper_sigma, length.out = n_grid)
  sigma_phi_grid <- seq(from = 0.1, to = grid_upper_sigma, length.out = n_grid)
  lambda_a_grid <- seq(from = 0.1, to = grid_upper_lambda, length.out = 2*n_grid)
  lambda_phi_grid <- seq(from = 0.1, to = grid_upper_lambda, length.out = 2*n_grid)
  
  a_tmp <- rep(NA, n_grid)
  phi_tmp <- rep(NA, n_grid)
  sparsity_a_tmp <- rep(NA, n_grid)
  sparsity_phi_tmp <- rep(NA, n_grid)
  
  
  # - setting \sigma_x ## (DE(0, \sigma_x)), x in {a, \phi}
  cat("..setting sigma_x...")
  for (jj in 1:n_grid) {
    cat(jj, "/", n_grid, "\n")
    magnitude_a_sim<- DExp_upperQuantile_MC(X, N_MC, Q, sparsity_level, sigma_a_grid[jj])
    magnitude_phi_sim <- DExp_upperQuantile_MC(X, N_MC, Q, sparsity_level, sigma_phi_grid[jj])
    a_tmp[jj] <- mean(magnitude_a_sim)
    phi_tmp[jj] <- mean(magnitude_phi_sim)
  }
  sigma_a_set <- sigma_a_grid[which.min(abs(a_tmp - ampl_upper))]
  sigma_phi_set <- sigma_phi_grid[which.min(abs(phi_tmp - phase_upper))]
  
  # setting r_x | \sigma_a ## (DE(0, r_x) x in {a, \phi}
  cat("..setting r_x...")
  for (jj in 1:(2*n_grid)) {
    cat(jj, "/", (2*n_grid), "\n")
    a_l1ball <- generate_theta_l1ball(Q, r_alpha = lambda_a_grid[jj], 
                                      beta_gen = function(p){rlaplace(p, 0, sigma_a_set)},
                                      N_MC)
    phi_l1ball <- generate_theta_l1ball(Q, r_alpha = lambda_phi_grid[jj], 
                                        beta_gen = function(p){rlaplace(p, 0, sigma_phi_set)},
                                        N_MC)
    sparsity_a_tmp[jj] <- mean(a_l1ball == 0)
    sparsity_phi_tmp[jj] <- mean(phi_l1ball == 0)
  }
  lambda_a_set <- lambda_a_grid[which.min(abs(sparsity_a_tmp - sparsity_level))]
  lambda_phi_set <- lambda_phi_grid[which.min(abs(sparsity_phi_tmp - sparsity_level))]
  
  out = list()
  out$sigma_a = sigma_a_set 
  out$sigma_phi = sigma_phi_set
  out$lambda_a = lambda_a_set 
  out$lambda_phi = lambda_phi_set
  
  return(out)
}


# - 


getMetrics <- function(predicted, true_values) 
{
  
  out <- list()
  # Create confusion matrix
  confusion_matrix <- confusionMatrix(factor(predicted), factor(true_values))
  
  # Extract metrics
  out$accuracy <- confusion_matrix$overall['Accuracy']
  #out$specificity <- confusion_matrix$byClass['Spec']
  out$precision <- confusion_matrix$byClass['Pos Pred Value']
  #out$sensitivity <- confusion_matrix$byClass['Sens']
  out$MCC <- mcc(predicted, true_values)
  out$F1 <- F1_Score(predicted, true_values)
  
  return(out)
}


# - 
addNightBlocks <- function(sbj, T_all, lower_rct = -0.05, upper_rct = -0.01) 
{
  rect(0, lower_rct, 96, upper_rct, col="black") 
  i = 96 + 12*12
  rect(0 + i, lower_rct, i + 12*12, upper_rct, col="black") 
  
  for (k in 1:(ceil(T_all[sbj]/288)-2)) {
    i = i + 2*12*12
    rect(0 + i, lower_rct, i + 12*12, upper_rct, col="black") 
  }
  i = i + 2*12*12
  rect(0 + i, lower_rct, N, upper_rct, col="black") 
  
}


