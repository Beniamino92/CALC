
# - 
downSampleAllSubjects <- function(obs_all, n_skip)
{
  N_sbj = dim(obs_all)[2]
  T_sbj = dim(obs_all)[1]
  T_sub = length(seq(from = 1, to = T_sbj, by = n_skip))-1
  out = matrix(NA, N_sbj, T_sub)
  
  for (ii in 1:N_sbj) {
    out[ii, ] = downSample(obs_all[, ii], n_skip) 
  }
  
  return(out)
}



getActigraphData <- function(day_min, window_smoothing, n_skip, 
                             intercept = TRUE, interaction = TRUE, 
                             standardize_covariate = TRUE)
{
  df <- read_excel("data_Vaishnav/data.xlsx", range = cell_cols("C:BA"))
  N_sbj = dim(df)[2]
  # - actigraph data, all subjects, all time points
  obs_all_raw = df[3:dim(df)[1], ]
  # - covariates (age & sex)
  X_aux = t(as.matrix(df[1:2, ]))
  
  # - length time series, for each subjet
  T_all = numeric(N_sbj)
  for (ii in 1:N_sbj) {
    first_NA = which(is.na(obs_all_raw[, ii]))[1]
    T_all[ii] = first_NA-1
  }
  T_all[which(is.na(T_all))] = dim(obs_all_raw)[1] 
  
  # - index subjects analyzed
  T_min = day_min*1440
  idxs_sbj = which(T_all >= T_min)
  T_max = max(T_all[idxs_sbj])
  N_sbj_sub = length(idxs_sbj)
  
  
  # list: sub set data analyzed + log transformation
  obs_list = list()
  for (ii in 1:N_sbj_sub) {
    sbj = idxs_sbj[ii]
    obs_list[[ii]] = log(as.numeric(as.matrix(obs_all_raw[1:T_all[sbj], sbj])) + 1)
  }
  
  # list: smoothed data
  obs_smooth_list = list()
  for (ii in 1:N_sbj_sub) {
    if (window_smoothing > 1) {
      temp = movavg(obs_list[[ii]], window_smoothing, "s")
      obs_smooth_list[[ii]] = temp[window_smoothing:length(temp)]
    } else {
      obs_smooth_list[[ii]] = obs_list[[ii]]
    }
  }
  
  # down sample 
  obs_downsampled_list = list()
  for (ii in 1:N_sbj_sub) {
    obs_downsampled_list[[ii]] = downSample(obs_smooth_list[[ii]], n_skip)
  }
  
  # -- -from list let's get matrix and length of each time series
  T_all_out = numeric(N_sbj_sub)
  for (ii in 1:N_sbj_sub) {
    T_all_out[ii] = length(obs_downsampled_list[[ii]])
  }
  
  # -output 
  obs_out = matrix(-1, N_sbj_sub, max(T_all_out))
  for (ii in 1:N_sbj_sub) {
    obs_out[ii, 1:T_all_out[ii]] = obs_downsampled_list[[ii]]
  }
  
  # - creating matrix of covariates.
  X = matrix(NA, nrow = N_sbj, ncol = 2)
  X[which(X_aux[, 1] == "F"), 1] = 1
  X[which(X_aux[, 1] == "M"), 1] = 0
  X[, 2] = as.integer(X_aux[, 2])
  X_age = X[, 2]
  
  # - include intercept
  if (intercept) {
    X = cbind(rep(1, N_sbj), X)
    # - standardizing
    if (standardize_covariate) {
      X[, 3] = (X[, 3] - mean(X[, 3]))/sd(X[, 3])
    }
  } else {
    # - standardizing
    if (standardize_covariate) {
      X[, 2] = (X[, 2] - mean(X[, 2]))/sd(X[, 2])
    }
  }
  # interaction term (gender x age)
  interaction_gender_age = X[, 2]*X[, 3]
  if (interaction) {
    X = cbind(X, interaction_gender_age)
  }
  X = X[idxs_sbj, ]
  
  
  
  covariates = list()
  covariates$X = X
  covariates$age = X_age[idxs_sbj]
  
  return(list(Y = obs_out, T_all = T_all_out, X = covariates))
}


getActigraphDataFullCovariates <- function(day_min, window_smoothing, n_skip)
{
  df_old <- read_excel("data_Vaishnav/data.xlsx", range = cell_cols("C:BA"))
  df_new <- read_excel("data_Vaishnav/data_with_covariates.xlsx", 
                       range = cell_cols("A:BH"))
  IDs_old <- names(df_old)
  IDs_new <- names(df_new)[c(-1)]
  IDs_different <- setdiff(IDs_new, IDs_old)
  IDs_common <- intersect(IDs_new, IDs_old)
  
  # - extracting data 
  temp = df_new[16:dim(df_new)[1], 2:60]
  obs_all_raw = temp[, IDs_common]
  N_sbj = dim(obs_all_raw)[2]
  
  # - length time series, for each subjet
  T_all = numeric(N_sbj)
  for (ii in 1:N_sbj) {
    first_NA = which(is.na(obs_all_raw[, ii]))[1]
    T_all[ii] = first_NA-1
  }
  T_all[which(is.na(T_all))] = dim(obs_all_raw)[1] 
  
  # - index subjects analyzed
  T_min = day_min*1440
  idxs_sbj = which(T_all >= T_min)
  T_max = max(T_all[idxs_sbj])
  N_sbj_sub = length(idxs_sbj)
  
  
  # list: sub set data analyzed + log transformation
  obs_list = list()
  for (ii in 1:N_sbj_sub) {
    sbj = idxs_sbj[ii]
    obs_list[[ii]] = log(as.numeric(as.matrix(obs_all_raw[1:T_all[sbj], sbj])) + 1)
  }
  
  # list: smoothed data
  obs_smooth_list = list()
  for (ii in 1:N_sbj_sub) {
    if (window_smoothing > 1) {
      temp = movavg(obs_list[[ii]], window_smoothing, "s")
      obs_smooth_list[[ii]] = temp[window_smoothing:length(temp)]
    } else {
      obs_smooth_list[[ii]] = obs_list[[ii]]
    }
  }
  
  # - time stamps (every minute, for each subject)
  timeStamps_all <- list()
  temp <- seq(
    from=as.POSIXct("2012-1-1 0:00", tz="UTC"),
    to=as.POSIXct("2012-2-3 23:00", tz="UTC"),
    by="1 min")  
  for (sbj in 1:N_sbj) {
    timeStamps_all[[sbj]] <- as_hms(temp[window_smoothing:T_all[sbj]])
  }
  for (sbj in 1:N_sbj) {
    aux = seq(from = 1, to = length(timeStamps_all[[sbj]]), by = n_skip)
    timeStamps_all[[sbj]] = timeStamps_all[[sbj]][aux]
  }
  timeStamps_out = list()
  for (ii in 1:N_sbj_sub){
    sbj = idxs_sbj[ii]
    timeStamps_out[[ii]] = timeStamps_all[[sbj]]
  }
  
  # down sample 
  obs_downsampled_list = list()
  for (ii in 1:N_sbj_sub) {
    obs_downsampled_list[[ii]] = downSample(obs_smooth_list[[ii]], n_skip)
  }
  
  # -- -from list let's get matrix and length of each time series
  T_all_out = numeric(N_sbj_sub)
  for (ii in 1:N_sbj_sub) {
    T_all_out[ii] = length(obs_downsampled_list[[ii]])
  }
  
  # -output 
  obs_out = matrix(-1, N_sbj_sub, max(T_all_out))
  for (ii in 1:N_sbj_sub) {
    obs_out[ii, 1:T_all_out[ii]] = obs_downsampled_list[[ii]]
  }
  
  
  # - extracting covariates
  temp = t(as.matrix(df_new[1:14, ]))
  covariate_names <- gsub(" ", "", temp[1, ])
  X_temp <- temp[2:60, ]
  X_covariates <- X_temp[IDs_common, ]
  # recoding variable
  X = matrix(NA, nrow = N_sbj, ncol = dim(X_covariates)[2])
  for (qq in 1:dim(X_covariates)[2]) {
    if (X_covariates[1, qq] %in% c("F", "M", "NO", "YES")) {
      X[, qq] = recodeCovariates(X_covariates[, qq])
    } else {
      X[, qq] = as.integer(X_covariates[, qq])
    }
  }
  X <- matrix(as.numeric(X), ncol = ncol(X))
  colnames(X) <- covariate_names
  # standardizing
  X_std = X
  for (qq in 1:dim(X)[2]) {
    if(!setequal(sort(unique(X[, qq])), c(0, 1))) {
      X_std[, qq] = (X[, qq] - mean(X[, qq], na.rm = T))/sd(X[which(!is.na(X[, qq])), qq], na.rm = T)
    }
  }
  # missing values 
  X_std_imp <- missForest(X_std)$ximp
  # intercept
  X = cbind(rep(1, N_sbj), X)
  X_std = cbind(rep(1, N_sbj), X_std_imp)
  colnames(X)[1] <- "Intercept"
  colnames(X_std)[1] <- "Intercept"
  
  # output 
  covariates = list()
  covariates$X = X[idxs_sbj, ]
  covariates$X_std = X_std[idxs_sbj, ]
  covariates$IDs = IDs_common[idxs_sbj]
    timeStamps <- timeStamps_out
  
  return(list(Y = obs_out, T_all = T_all_out, covariates = covariates, 
              timeStamps =  timeStamps_out, IDs_sbj = IDs_common[idxs_sbj]))
}


# function generating covariates

generate_covariate <- function(N, intercept = FALSE)
{
  X_aux <- data.frame(gender = c(rep(1, as.integer(N/2)), rep(0, as.integer(N/2))), 
                      height = mixture_normal(N, mu = c(180, 160), sigma = c(2, 2), p = c(0.5, 0.5)), 
                      weight = mixture_normal(N, mu = c(75, 55), sigma = c(2, 2), p = c(0.5, 0.5)))
  if (intercept) {
    X <- cbind(cbind(rep(1, N), X_aux[, 1]),
               standardize(X_aux[, 2:3]))
  } else {
    X <- cbind(X_aux[, 1], standardize(X_aux[, 2:3]))
  }
  
  X <- as.data.frame(X)
  names(X) <- c("intercept", "gender", "heigth", "weigth")
  X <- as.matrix(X)
  
  return(X)
}


# - 
standardize <- function(obs, only_sd = FALSE)
{
  N <- nrow(obs)
  D <- ncol(obs)
  out <- matrix(NA, nrow = N, ncol = D)
  for (dd in 1:D) {
    if (only_sd) {
      out[, dd] <- (obs[, dd])/sd(obs[, dd])
    } else {
      out[, dd] <- (obs[, dd] - mean(obs[, dd]))/sd(obs[, dd])
    }
  }
  return(out)
}

# - 
filter_causal <- function(y, win_s) {
  
  N <- length(y)
  N_smooth <- N - win_s + 1
  y_smooth <- c()
  for (t in 1:N_smooth) {
    y_smooth <- c(y_smooth, mean(y[t:(t+win_s-1)], na.rm = T))
  }
  y_smooth <- c(rep(NA, win_s  - 1), y_smooth)
  
  return(y_smooth)
}

# - 
downSample <- function(obs, n_skip)
{
  N = length(obs)[1]
  aux = seq(from = 1, to = N, by = n_skip)
  obs_temp = numeric(length(aux)-1)
  for (jj in 1:(length(aux)-1)) {
    a = aux[jj]
    if (jj == (length(aux)-1)) {
      b = aux[jj+1]
    } else {
      b = aux[jj+1]-1
    }
    obs_temp[jj] = mean(obs[a:b])
  }
  return(obs_temp)
}


# - 
get_legend_quantiles <- function(x, n_levels)
{
  quantile_aux <- quantile(x, probs = seq(0, 1, by = 1/n_levels))
  quantile_aux <- as.integer(round(quantile_aux, 1))
  out <- character(n_levels)
  for (ii in 1:n_levels) {
    temp <- as.character(c(quantile_aux[ii], quantile_aux[ii+1]))
    out[ii] <-   paste(temp[1], temp[2],  sep = "-")
  }
  return(out)
}

# - 
recodeCovariates <- function(X)
{
  X[X == "F"] = 1
  X[X == "M"] = 0
  X[X == "NO"] = 0
  X[X == "YES"] = 1
  
  return(X)
}

# generate_covariate_mix <- function(N_sbj, Q, n_binary)
# {
#   X = matrix(NA, nrow = N_sbj, ncol = Q)
#   # first covariate from mixture
#   X[, 1] = mixture_normal(N_sbj, mu = c(-10, 10),
#                           sigma = c(1, 1),
#                           p = c(0.5, 0.5))
#   # other from normal
#   for (qq in 2:(Q-n_binary)) {
#     X[, qq] <- rnorm(N_sbj, 0, 1)
#   }
# 
#   X[, 1:(Q-n_binary)] = standardize(X[, 1:(Q-n_binary)])
#   for (qq in (Q-n_binary+1):Q) {
#     X[, qq] <- rbinom(N_sbj, 1, prob = 0.5)
#   }
#   return(X)
# }

# 
generate_covariate_mix <- function(N_sbj, M = 11)
{
  # meaningful covariates
  X_meaningful <- data.frame(gender = c(rep(1, as.integer(N_sbj/2)), rep(0, as.integer(N_sbj/2))),
                             height = mixture_normal(N_sbj, mu = c(180, 160), sigma = c(2, 2), p = c(0.5, 0.5)),
                             weight = mixture_normal(N_sbj, mu = c(75, 55), sigma = c(2, 2), p = c(0.5, 0.5)))

  # non meaningful covariates
  X_nonsense <- matrix(rnorm(N_sbj*M), nrow = N_sbj, ncol = M)

  # bind + standardzie
  X = cbind(X_meaningful, X_nonsense)
  X = cbind(X[, 1], standardize(X[, 2:ncol(X)]))

  # binary covariates
  X[, 7] = rbinom(N_sbj, 1, prob = 0.5)
  X[, 10] = rbinom(N_sbj, 1, prob = 0.5)

  # intercept
  X <- cbind(rep(1, N_sbj), X)
  head(X)

  # df and name col
  X <- as.data.frame(X)
  names(X) <- c(c("intercept", "gender", "heigth", "weigth"), LETTERS[1:M])
  X <- as.matrix(X)

  return(X)
}
