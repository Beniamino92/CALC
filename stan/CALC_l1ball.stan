
functions {
  real expit(real x) 
  {
    real out;
    out = exp(x)/(1+exp(x)); 
    return out; 
  }
  
  real circadian(int tt, real R, real m, real a, 
                 real alpha, real beta, real phi) 
  {
    real out; 
    out = m + a*expit(beta*(cos((tt/R - phi) * (2*pi())/24) - alpha));
    return out; 
  }
  
  vector l1_ball_projection(vector eta, real r) {
    
     if(sum(abs(eta)) <= r){return eta;}
     else{
      int p = num_elements(eta);
      vector[p] abs_eta_sort = sort_desc(abs(eta));
      vector[p] cumsum_abs_eta_sort = cumulative_sum(abs_eta_sort);
      vector[p] mu;
      vector[p] theta;
      int c = p;
      for(j in 1:p){
         mu[j] = max({cumsum_abs_eta_sort[j] - r, 0});
         if(abs_eta_sort[j] > mu[j]/j){
            c = j;
         }
      }
      for(i in 1:p){
         theta[i] = abs(eta[i])/eta[i]*max({abs(eta[i]) - mu[c]/c, 0});
      }      
      return theta;
     }
    }
  
  real llk_lp(array[] int T_all, int N_sbj, real R, array[] vector Y,
              vector m, vector a, vector alpha, vector beta, vector phi, 
              vector sigma)
  {
    real out = 0.0;
    for (nn in 1:N_sbj) {
      int T = T_all[nn];
      vector[T] mu;
      for (tt in 1:T) {
        mu[tt] = circadian(tt, R, m[nn], a[nn], alpha[nn], beta[nn], phi[nn]);
      }
      out += normal_lpdf(Y[nn, 1:T] | mu, sigma[nn]);
    }
    return out;
  }
  
}


data {
  int<lower=0> N_sbj; // number of subjects
  array[N_sbj] int<lower=1> T_all; // length of time series (for each subject)
  int<lower=1> T_max; // maximum length of time series
  int<lower=1> Q; // number of covariates (including intercept)
  real<lower=1> R; // sampling rate (interval per hour)
  array[N_sbj] vector[T_max] Y; // time series for all sbj
  array[N_sbj] vector[Q] X; // covariates
  
  // hyperparms
  // real<lower=0> mean_beta;
  // real<lower=0> sd_beta;
  // real<lower=0> mean_sigma;
  real<lower=0> sd_sigma;
  real mean_eta;
  real<lower=0> sd_eta_a;
  real<lower=0> sd_eta_phi;
  real<lower=0> r_lambda_a;
  real<lower=0> r_lambda_phi;
}

parameters {
  vector[N_sbj] m_aux;
  vector[N_sbj] beta_aux;
  vector<lower=-1,upper=1>[N_sbj] alpha;
  //vector<lower=0>[N_sbj] beta;
  vector[Q] eta_a; 
  vector[Q] eta_phi;
  real<lower=0> ray_a;
  real<lower=0> ray_phi;
  real mean_m; 
  real<lower=0> sd_m; 
  real mean_beta; 
  real<lower=0> sd_beta; 
  vector<lower=0>[N_sbj] sigma;
}

transformed parameters {
  vector<lower=0>[N_sbj] m; 
  vector<lower=0>[N_sbj] beta; 
  vector<lower=0, upper=24>[N_sbj] phi;
  vector<lower=0>[N_sbj] a;
  vector[Q] theta_a; 
  vector[Q] theta_phi; 
  
  theta_a = l1_ball_projection(eta_a, ray_a);
  theta_phi = l1_ball_projection(eta_phi, ray_phi);
  
  for (n in 1:N_sbj) {
    beta[n] = exp(beta_aux[n]); 
    m[n] = exp(m_aux[n]);
    a[n] = exp(X[n]'*theta_a);
    phi[n] = exp(X[n]'*theta_phi);
  }
  // probably I can write this as matrix multiplication
}



model {
  // priors
  target += exponential_lpdf(ray_a | r_lambda_a);
  target += exponential_lpdf(ray_phi | r_lambda_phi);
  target += double_exponential_lpdf(eta_a | mean_eta, sd_eta_a);
  target += double_exponential_lpdf(eta_phi | mean_eta, sd_eta_phi);
  target += normal_lpdf(m_aux | mean_m, sd_m);
  target += normal_lpdf(mean_m |0, sqrt(10)); 
  target += cauchy_lpdf(sd_m |0, 1); 
  target += normal_lpdf(beta_aux | mean_beta, sd_beta);
  target += normal_lpdf(mean_beta |5, 3); 
  target += cauchy_lpdf(sd_m |0, 1); 
  target += cauchy_lpdf(sigma | 0, sd_sigma);
  // llk
  target += llk_lp(T_all, N_sbj, R, Y, m, a, alpha, beta, phi, sigma);
}

// generated quantities {
//   vector[Q] theta_a; 
//   vector[Q] theta_phi; 
//   theta_a = l1_ball_projection(eta_a, ray_a);
//   theta_phi = l1_ball_projection(eta_phi, ray_phi);
// }
