# CALC
Bayesian covariate-dependent anti-logistic circadian (CALC) model for activity data. 

`CALC` is an R package for Bayesian modeling of circadian rest-activity rhythms using a covariate-dependent anti-logistic framework. Designed to integrate wearable actigraphy data with demographic and clinical predictors, `CALC` enables cohort-level analysis with enhanced flexibility and interpretability. The model incorporates covariate-driven amplitude and phase parameters while enforcing sparsity via an l1-ball projection prior, facilitating efficient estimation of significant predictors. It supports robust Bayesian inference using Hamiltonian Monte Carlo and is implemented in Stan. More details can be found in Hadj-Amar, B., Krishnan, V., & Vannucci, M. (2025) "Bayesian Covariate-Dependent Circadian Modeling of Rest-Activity Rhythms", published in Data Science in Science.

This software allows for the following four Bayesian modeling option 

1. CALC $l_1$-ball
2. CALC Lasso
3. RAR
4. Cosinor
  
## Example 

We provide a snapshot of `tutorial.Rmd`, which contains a tutorial for using our software in R

* Run RAPGP sampler
  ```R
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
<p align="center">
<img src="https://github.com/Beniamino92/BayesRPAGP/blob/main/plots/example.png" width="500" heigth="240"/> 
</p>

