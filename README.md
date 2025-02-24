# CALC
Bayesian covariate-dependent anti-logistic circadian (CALC) model for activity data. 

`CALC` is an R package for Bayesian modeling of circadian rest-activity rhythms using a covariate-dependent anti-logistic framework. Designed to integrate wearable actigraphy data with demographic and clinical predictors, `CALC` enables cohort-level analysis with enhanced flexibility and interpretability. The model incorporates covariate-driven amplitude and phase parameters while enforcing sparsity via an l1-ball projection prior, facilitating efficient estimation of significant predictors. It supports robust Bayesian inference using Hamiltonian Monte Carlo and is implemented in Stan. More details can be found in Hadj-Amar, B., Krishnan, V., & Vannucci, M. (2025) "Bayesian Covariate-Dependent Circadian Modeling of Rest-Activity Rhythms", published in Data Science in Science.

This software allows for the following four Bayesian modeling option 

1. CALC $l_1$-ball
2. CALC Lasso
3. RAR
4. Cosinor
  
## Example 

We provide a snapshot of `tutorial.R`, which contains a tutorial for using our software in R

```r
# stan model
model <- stan_model(file = "stan/CALC_l1ball.stan")

# stan data
data_stan <- list(N_sbj = N_sbj,
                  T_all = T_all,
                  T_max = T_max,
                  Q = dim(X)[2],
                  R = R,
                  Y = obs_all,
                  X = X,
                  sd_sigma = 1,
                  mean_eta = 0,
                  sd_eta_a = hyperParms$sigma_a,
                  sd_eta_phi =  hyperParms$sigma_phi,
                  r_lambda_a = hyperParms$lambda_a,
                  r_lambda_phi = hyperParms$lambda_phi)
# stan init
init_fun <- function(...) list(eta_a = c(1, rep(0, Q-1)), 
                               eta_phi = c(1, rep(0, Q-1)))

# stan fit
fit_sparseRAR <- sampling(object = model,
                          data = data_stan,
                          seed = 108,
                          chains = 1,
                          iter = N_warmup + N_MCMC,
                          warmup = N_warmup,
                          init = init_fun)
```

```r
# - Representative Traceplots

par(mfrow = c(4, 4))
par(mai=rep(0.4, 4))
for (ii in 1:Q) {
  plot(theta_a_sims[, ii], type = "l", xlab = "Iterations", col = "goldenrod3", 
       main = bquote(eta[.(paste(" a,", ii))]), cex.main = 1.5)
}
```

<p align="center">
<img src="https://github.com/Beniamino92/sparseVARHSMM/blob/main/figures/postpred_training.png" width="700" heigth="100"/> 
</p>

```
# - PPI and Magnitude
PPI_a = numeric(Q)
PPI_phi = numeric(Q)
for (jj in 1:Q) {
  PPI_a[jj] = sum(theta_a_sims[, jj] != 0)/N_MCMC
  PPI_phi[jj] = sum(theta_phi_sims[, jj] != 0)/N_MCMC
}
```

<p align="center">
<img src="https://github.com/Beniamino92/sparseVARHSMM/blob/main/figures/postpred_training.png" width="700" heigth="100"/> 
</p>

``` r 
# - Posterior Predictive 

par(mfrow = c(2, 2))
par(mai=rep(0.4, 4))
for (sbj in sbj_plots) {
  N = T_all[sbj]
  y_hat = getPosteriorPredictive_covariates(parms, sbj, T_all, R, n_draws = 100)
  plot(1:N, obs_all[sbj, 1:T_all[sbj]],  
       ylab = "y", xlab = "", 
       xlim = c(1, N), cex.lab = 1.1, 
       ylim = c(0, 7),
       cex = 0.6,
       main = paste("Subject ", sbj, sep = ''))
  
  points(1:N, obs_all[sbj, 1:T_all[sbj]],
         cex = 0.6, col = "blue")
  points(1:N, obs_all[sbj, 1:T_all[sbj]],
         pch = 20, cex = 0.7, col = "lightblue")
  
  probs <- seq(from=0.1, to=0.9, by=0.1)
  cred <- sapply(1:N, function(t) quantile(y_hat[, t], probs = probs))
  polygon(c(1:N, rev(1:N)), c(cred[1,], rev(cred[9,])),
          col = scales::alpha(c_light, 0.5), border = NA)
  polygon(c(1:N, rev(1:N)), c(cred[2,], rev(cred[8,])),
          col = scales::alpha(c_light_highlight, 0.5), border = NA)
  polygon(c(1:N, rev(1:N)), c(cred[3,], rev(cred[7,])),
          col = scales::alpha(c_mid, 0.5), border = NA)
  polygon(c(1:N, rev(1:N)), c(cred[4,], rev(cred[6,])),
          col = scales::alpha(c_mid_highlight, 0.5), border = NA)
}
```
<p align="center">
<img src="https://github.com/Beniamino92/sparseVARHSMM/blob/main/figures/postpred_training.png" width="700" heigth="100"/> 
</p>
