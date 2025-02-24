functions{
  
  real circadian(int tt, real R, real m, real a, real phi) 
  {
    real out; 
    out = m + a*cos((tt/R - phi) * (2*pi())/24);
    return out; 
  }
  
  real llk_lp(int N, real R, vector y, real m, real a, 
              real phi, real sigma)
  {
    real out;
    vector[N] mu;
    
    for (tt in 1:N) {
      mu[tt] = circadian(tt, R, m, a, phi);
    }
    
    out = normal_lpdf(y | mu, sigma);
    return out;
  }
}


data {
  int<lower=0> N;
  real<lower=1> R; 
  vector[N] y;
  real<lower=0> mean_m;
  real<lower=0> sd_m;
  real<lower=0> mean_a;
  real<lower=0> sd_a;
  real<lower=0> mean_sigma;
  real<lower=0> sd_sigma;
}


parameters {
  real<lower=0> m;
  real<lower=0> a;
  real<lower=0,upper=24> phi;
  real<lower=0> sigma;
}


model {
  // priors
  target += normal_lpdf(m | mean_m, sd_m);
  target += normal_lpdf(a | mean_a, sd_a);
  target += cauchy_lpdf(sigma | mean_sigma, sd_sigma);
  // llk 
  target += llk_lp(N, R, y, m, a, phi, sigma);
}

