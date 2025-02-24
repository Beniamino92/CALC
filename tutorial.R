my_dir <- "/Users/beniamino/Desktop_New/Research/My_Papers/Methodological/DSS_2025+/Untitled/CALC/"
setwd(my_dir)

c_light = "gray95"
c_light_highlight = "gray90"
c_mid = "gray85"
c_mid_highlight = "gray80"
c_dark = "gray75"

library("farff")
library("rstan")
library("bayesplot")
library("bridgesampling")
library("locfit")
library("RColorBrewer")
library("ggplot2")
library("dplyr")
library("tidyr")
library("LaplacesDemon")
library("scales")

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())

source("include/util.R")
source("include/util_data.R")
source("include/util_stan.R")


# generating parameter
N_sbj <- 30
N <- 1000
R <- 12.0
T_all = rep(N, N_sbj)
T_max = N

# generating covariate 
set.seed(108)
M = 11
obs_all <- matrix(NA, nrow = N_sbj, ncol = N)
X <- generate_covariate_mix(N_sbj, M)
Q <- ncol(X)
head(X)
colnames(X) <- 1:15
head(X)


# true effects
eta_phi = c(c(1.0, 0.9, 0.01, 0), rep(0, M))
eta_a = c(c(0.6, 0.0, 0.0, 0.5, -0.6, 0, 0.1), rep(0, M-3))
indicator_phi = as.integer(eta_phi != 0)
indicator_a = as.integer(eta_a != 0)


theta_all <- list()
theta <- list()
theta$m = 3 # > 0
theta$a = NULL # > 0
theta$beta =  10 # > 0
theta$alpha = -0.3 # in (-1,1)
theta$phi = NULL # in [0, 24]


# -- generating data all subjects 
for (nn in 1:N_sbj) {
  
  # sbj specifc acrophase phi
  theta_all[[nn]] <- theta
  theta_all[[nn]]$phi <- as.numeric(exp(X[nn, ] %*% eta_phi))
  theta_all[[nn]]$a <- as.numeric(exp(X[nn, ] %*% eta_a))
  
  # RAR + error 
  signal <- sapply(1:N, circadian , R, theta_all[[nn]])
  error <- rnorm(N, sd = 0.5)
  obs_all[nn, ] <- signal + error
}


# - 
phi_all = numeric(N_sbj)
a_all = numeric(N_sbj)
for (nn in 1:N_sbj) {
  phi_all[nn] = theta_all[[nn]]$phi
  a_all[nn] = theta_all[[nn]]$a
}


# ------------------  plot data paper
sbj_plots <- c(6, 15, 20, 25)
n_draw = length(sbj_plots)
my_pallete <- brewer.pal(n_draw, "Dark2")

# - prep df
df_obs = as.data.frame(obs_all)
colnames(df_obs) <- 1:N
df_obs$sbj = as.character(1:N_sbj)
df_tmp <- df_obs %>%  as_tibble %>% 
  pivot_longer(-1001, names_to = "time") %>% 
  filter(sbj %in% sbj_plots) %>%  
  mutate(time = as.integer(time))

# 
df_tmp %>% 
  ggplot(aes(x = time, y = value, color = sbj)) +  
  geom_line(alpha = 0.7) + 
  geom_point(size = 0.5, alpha = 0.7) +
  scale_color_manual(name = "Subject", values=my_pallete) + 
  theme_bw() + 
  labs(fill = "", y = "y", x = "Time") + 
  theme(panel.grid = element_blank(),
        axis.text = element_text(size = 20, color = "black"),
        axis.title = element_text(size = 20),
        legend.title = element_text(size=15),
        legend.text = element_text(size = 15), 
        panel.border = element_blank(),
        axis.line.x = element_line(),
        axis.line.y = element_line(), 
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) 



# - MCMC info
N_MCMC = 2000
N_warmup = 1000

# - settings hyperparameters with/without prior information. 
priorInfoCALC <- FALSE
if (priorInfoCALC) {
  settingsHyperParmsSearch <- list()
  settingsHyperParmsSearch$ampl_upper <- 6
  settingsHyperParmsSearch$phase_upper <- 10
  settingsHyperParmsSearch$sparsity_level <- 0.75
  settingsHyperParmsSearch$N_MC <- 5000
  settingsHyperParmsSearch$n_grid <- 30
  settingsHyperParmsSearch$grid_upper_sigma <- 0.5
  settingsHyperParmsSearch$grid_upper_lambda <- 5
  hyperParms <- setHyperametersCALC(X, settingsHyperParmsSearch)
} else {
  hyperParms <- list()
  hyperParms$sigma_a <- 1;
  hyperParms$sigma_phi <- 1;
  hyperParms$lambda_a <- 1;
  hyperParms$lambda_phi <- 1;
}


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
print(fit_sparseRAR, pars = c("theta_a", "theta_phi"))

load("PAPER_illustrative_example.RData")

# extract output
parms <- rstan::extract(fit_sparseRAR) 
eta_a_sims <- parms$eta_a
eta_phi_sims <- parms$eta_phi
theta_a_sims <- parms$theta_a
theta_phi_sims <- parms$theta_phi

# trace plot 
par(mfrow = c(4, 4))
par(mai=rep(0.4, 4))
for (ii in 1:Q) {
  plot(theta_a_sims[, ii], type = "l", xlab = "Iterations", col = "goldenrod3", 
       main = bquote(eta[.(paste(" a,", ii))]), cex.main = 1.5)
}


# trace plots 
traceplot(fit_sparseRAR, 
          pars = c("theta_a"))
traceplot(fit_sparseRAR, 
          pars = c("theta_phi"))
traceplot(fit_sparseRAR, 
          pars = c("eta_a"))
traceplot(fit_sparseRAR, 
          pars = c("eta_phi"))
traceplot(fit_sparseRAR, 
          pars = c("ray_a"))
traceplot(fit_sparseRAR, 
          pars = c("ray_phi"))
traceplot(fit_sparseRAR, 
          pars = c("sigma"))
traceplot(fit_sparseRAR, 
          pars = c("m"))
traceplot(fit_sparseRAR, 
          pars = c("alpha"))
traceplot(fit_sparseRAR, 
          pars = c("lp__"))


# ---- plot posterior predictive
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


# - PPI 
PPI_a = numeric(Q)
PPI_phi = numeric(Q)
for (jj in 1:Q) {
  PPI_a[jj] = sum(theta_a_sims[, jj] != 0)/N_MCMC
  PPI_phi[jj] = sum(theta_phi_sims[, jj] != 0)/N_MCMC
}


# df a and phi
df_a <- reshape2::melt(data.frame(theta_a_sims))
df_phi <- reshape2::melt(data.frame(theta_phi_sims))

# ------ plot PPI + boxplots 
col_a = rep("grey", Q); col_a[PPI_a >= 0.5] = "red"
col_phi = rep("grey", Q); col_phi[PPI_phi >= 0.5] = "red"


# ------------ PPI and Magnitude: 
par(mfrow = c(2, 2))
labels <- 1:15
# Plot PPI (amplitude)
plot(PPI_a, type = "h", ylim = c(0, 1),
     cex = 1.5,  pch = 20, xaxt = "n", xlab = "Covariates",
     col = col_a, ylab = "PPI", lwd = 2, 
     main = "PPI - Amplitude")
points(indicator_a, cex = 1.5, col = col_a, pch = 20)
abline(h = 0.5, lty = "dotted", col = "black")
ticks <- seq(from = 1, to = Q)
#labels <- colnames(X)
axis(side = 1, ticks, labels, cex.axis = 0.9, tck=-0.02, las = 2)

# Magnitude (amplitude) - boxplot
boxplot(value ~ variable, 
        data = df_a,
        col = col_a,
        ylab = expression(eta[a]), 
        xlab = "Covariates",
        xaxt = "n", main = "Magnitude Effect - Amplitude")
axis(side = 1, ticks, labels, cex.axis = 0.9, tck=-0.02, las = 2)
x0s <- which(indicator_a == 1) - 0.4
x1s <- which(indicator_a == 1) + 0.4
y0s = eta_a[which(indicator_a == 1) ]
segments(x0 = x0s, x1 = x1s, y0 = y0s, col = "blue", lwd = 3)

# Plot PPI (phase)
plot(PPI_phi, type = "h", ylim = c(0, 1),
     cex = 1.5,  pch = 20, xaxt = "n", xlab = "Covariates",
     col = col_phi, ylab = "PPI", lwd = 2, 
     main = "PPI - Acrophase")
points(indicator_phi, cex = 1.5, col = col_phi, pch = 20)
abline(h = 0.5, lty = "dotted", col = "black")
ticks <- seq(from = 1, to = Q)
axis(side = 1, ticks, labels, cex.axis = 0.9, tck=-0.02, las = 2)

# Magnitude (Acrophase) - boxplot
boxplot(value ~ variable, 
        data = df_phi,
        col = col_phi,
        ylab = expression(eta[phi]),
        xlab = "Covariates",
        main = "Magnitude Effect - Acrophase",
        xaxt = "n")
axis(side = 1, ticks, labels, cex.axis = 0.9, tck=-0.02, las = 2)
x0s <- which(indicator_phi == 1) - 0.4
x1s <- which(indicator_phi == 1) + 0.4
y0s = eta_phi[which(indicator_phi == 1) ]
segments(x0 = x0s, x1 = x1s, y0 = y0s, col = "blue", lwd = 3)


# --------- scatter plot (true vs estimated) for each subject
amplitude_est <- numeric(N_sbj)
phase_est <- numeric(N_sbj)
m_est <- numeric(N_sbj)
alpha_est <- numeric(N_sbj)
beta_est <- numeric(N_sbj)
for (nn in 1:N_sbj) {
  amplitude_est[nn] = exp(apply(theta_a_sims, 2, mean)  %*% X[nn, ])
  phase_est[nn] = exp(apply(theta_phi_sims, 2, mean)  %*% X[nn, ])
}


# plot 
par(mfrow = c(1, 2))
plot(a_all, amplitude_est, 
     xlab = "True", ylab = "Estimate", 
     main = "Amplitude", 
     pch = 20, col = "blue", cex = 1.1)
abline(coef = c(0,1))
plot(phi_all, phase_est,
     xlab = "True", 
     ylab = "Estimate", main = "Phase", 
     pch = 20, col = "blue", cex = 1.1)
abline(coef = c(0,1))

# - 
par(mfrow = c(2, 2))
par(mai=rep(0.4, 4))
hist(parms$m, col = "green2", freq = F, breaks = 20, 
     main = "m")
abline(v = theta$m, col = "red", lty = "dotted", lwd = 3)
hist(parms$alpha, col = "green2", freq = F, 
     main = expression(alpha))
abline(v = theta$alpha, col = "red", lty = "dotted", lwd = 3)
hist(parms$beta,  col = "green2", freq = F, xlim = c(2, 14), 
     main = expression(beta))
abline(v = theta$beta, col = "red", lty = "dotted", lwd = 3)
hist(parms$sigma, col = "green2", freq = F, main = expression(sigma))
abline(v = 0.5, col = "red", lty = "dotted", lwd = 3)


# estimate subject specific parameters
mean(apply(parms$m, 2, mean))
mean(apply(parms$alpha, 2, mean))
mean(apply(parms$beta, 2, mean))
mean(apply(parms$sigma, 2, mean))



