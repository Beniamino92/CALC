# CALC
Bayesian covariate-dependent anti-logistic circadian (CALC) model for activity data. 

BayesCALC is an R package for Bayesian modeling of circadian rest-activity rhythms using a covariate-dependent anti-logistic framework. Designed to integrate wearable actigraphy data with demographic and clinical predictors, BayesCALC enables cohort-level analysis with enhanced flexibility and interpretability. The model incorporates covariate-driven amplitude and phase parameters while enforcing sparsity via an l1-ball projection prior, facilitating efficient estimation of significant predictors. It supports robust Bayesian inference using Hamiltonian Monte Carlo and is implemented in Stan. More details can be found in Hadj-Amar, B., Krishnan, V., & Vannucci, M. (2025) "Bayesian Covariate-Dependent Circadian Modeling of Rest-Activity Rhythms", published in Data Science in Science.


## Example 

We provide a snapshot of `tutorial.Rmd`, which contains a tutorial for using our software in R

* Run RAPGP sampler
  ```R
  results <- fit_rpagp(y = dat$y, n_iter = n_iter,
                         theta0 = theta0, hyperparam = hyperparam,
                         pinned_point = pinned_point,
                         pinned_value = pinned_value)
  ```
<p align="center">
<img src="https://github.com/Beniamino92/BayesRPAGP/blob/main/plots/example.png" width="500" heigth="240"/> 
</p>

