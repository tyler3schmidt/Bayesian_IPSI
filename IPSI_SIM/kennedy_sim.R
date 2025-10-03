library(SuperLearner)

# simulation from Kennedy (2018)
SL.library <- c("SL.gam", "SL.earth", "SL.ksvm", 
                "SL.ranger", "SL.glm", "SL.glm.interaction", "SL.stepAIC")
#******************************************************************************#
#******************************************************************************#
# functions 

##########################################################################
# function: generate_data
# takes in: sample size n and seed 
# returns: simulated data frame sim_data
generate_data <- function(n, seed=NULL){
  if(!is.null(seed)) set.seed(seed) # if given seed
  
  # creating covariates
  X <- matrix(rnorm(n*4), nrow = n, ncol = 4)
  
  # naming columns
  colnames(X) <- paste0("X", 1:4) 
  
  # linear combination of covariates
  logit_p <- -X[,1] + 0.5 * X[,2] - 0.25 * X[,3] -0.1 * X[,4] 
  
  # propensity score
  p <- 1 / (1 + exp(-logit_p))
  
  # Treatment simulation
  a <- rbinom(n, 1, prob = p)
  
  # outcome simulation 
  mu <- 200 + a* (10 + 13.7 * (2 * X[,1] + X[,2] + X[,3] + X[,4]))
  y <- rnorm(n,  mean= mu, sd = 1)
  
  sim_data <- data.frame(a = a, y = y, X1 = X[,1], X2 = X[,2], X3 = X[,3], X4 = X[,4])
  return(sim_data)
}




########################################################################
# function: data_transformation
# takes in: data frame dat
# returns: transformed data frame transformed_data
# note: function for transforming covariates based on kang and schafer (2007)
data_transformation <- function(dat){
  # creating covariates
  n <- nrow(dat)
  Z <- matrix(rnorm(n*4), nrow = n, ncol = 4)
  
  # naming columns
  colnames(Z) <- paste0("Z", 1:4) 
  
  #transformations
  Z[,1] = exp(dat$X1 / 2)
  Z[,2] = dat$X2 / (1 + exp(dat$X1)) + 10
  Z[,3] = ((dat$X1 * dat$X3) / 25 + 0.6)^3 
  Z[,4] = (dat$X2 + dat$X4 + 20)^2
  
  # placing transformed covariates in a data frame
  transformed_data <- data.frame(a = dat$a, y = dat$y, 
                                 X1 = Z[,1], X2 = Z[,2], X3 = Z[,3], X4 = Z[,4])
  return(transformed_data)
}



#################################################################
# function: param_nuisance
# takes in: a data frame dat
# returns: pi_hat, mu1_hat, mu0_hat
param_nuisance <- function(dat) {
  pi_model <- suppressWarnings(glm(a ~ X1 + X2 + X3 + X4, data = dat, family = binomial()))
  pi_hat<- suppressWarnings(predict(pi_model, type="response"))
  
  mu1_model <- lm(y ~  X1 + X2 + X3 + X4, data = dat, subset = (a==1))
  mu1_hat <- predict(mu1_model, newdata = dat)
  
  mu0_model <- lm(y ~ X1 + X2 + X3 + X4, data = dat, subset = (a == 0))
  mu0_hat <- predict(mu0_model, newdata = dat)
  
  return(list(pi_hat = pi_hat, mu1_hat = mu1_hat, mu0_hat = mu0_hat))
}




################################################################
# function: nonparam_nuisance
# takes in: a data frame dat and number of sample splits K
# returns: pi_hat, mu1_hat, mu0_hat, and fold assignment folds
nonparam_nuisance <- function(dat, K) {
  n <- nrow(dat)
  
  # randomly splits sample into k folds
  good_folds <- FALSE
  fold_try_index = 0
  while (!good_folds) {
    # error handling 
    if(fold_try_index == 100000) {
      return(NA)
    }
    good_folds <- TRUE
    folds <- sample(rep(1:K, length.out = n))
    
    for (k in 1:K) {
      train <- which(folds != k) 
      # Check if treated and untreated units exist in training set
      if (sum(dat$a[train] == 1) == 0 | sum(dat$a[train] == 0) == 0) {
        good_folds <- FALSE
        break  # no need to check further if one fold fails
      }
    }
    fold_try_index <- fold_try_index + 1
  }
  
  
  # storage vectors 
  pi_hat <- numeric(n)
  mu1_hat <- numeric(n)
  mu0_hat <- numeric(n)
  
  for (k in 1:K) {
    train <- which(folds != k) 
    test <- which(folds ==k)
    
    # propensity model 
    sl_pi <- SuperLearner(Y = dat$a[train],
                          X = dat[train, c("X1", "X2", "X3", "X4")],
                          family = binomial(),
                          SL.library = SL.library)
    pi_hat[test] <- predict(sl_pi, newdata = dat[test, c("X1", "X2", "X3", "X4")])$pred
    
    # Outcome regression with treatment indicator included
    sl_mu <- SuperLearner(
      Y = dat$y[train],
      X = dat[train, c("X1", "X2", "X3", "X4", "a")],
      family = gaussian(),
      SL.library = SL.library
    )
    
    # To get mu1(x): set a=1 for test observations
    newdata1 <- dat[test, c("X1", "X2", "X3", "X4")]
    newdata1$a <- 1
    mu1_hat[test] <- predict(sl_mu, newdata = newdata1)$pred
    
    # To get mu0(x): set a=0 for test observations
    newdata0 <- dat[test, c("X1", "X2", "X3", "X4")]
    newdata0$a <- 0
    mu0_hat[test] <- predict(sl_mu, newdata = newdata0)$pred
  }
  
  return(list(pi_hat = pi_hat, mu1_hat = mu1_hat, mu0_hat = mu0_hat, folds = folds))
}




###################################################
# function: plug_in_est
# takes in a data frame, delta,pi_hat, mu1_hat, mu0_hat. 
# returns psi_hat
plug_in_est <- function(a, y, delta, pi_hat, mu1_hat, mu0_hat) {
  
  # plug in estimator calculation
  numerator <- delta * pi_hat * mu1_hat + (1 - pi_hat) * mu0_hat
  denominator <- delta * pi_hat + 1 - pi_hat
  
  psi_hat <- mean(numerator / denominator)
  return(psi_hat)
}




###########################################
# function: ipw_est
# takes in: a data frame, delta, pi_hat
# returns: psi_hat
ipw_est <- function(a, y, delta, pi_hat) {
  # ipw estimator calculation
  numerator <- (delta * a + 1 - a) * y
  denominator <- delta * pi_hat + 1 - pi_hat
  
  psi_hat <- mean(numerator / denominator)
  return(psi_hat)
}





###############################################################
# function: proposed_est
# takes in a data frame, delta,binary variable for parametric/ nonparametric,
# number of sample splits, pi_hat, mu1_hat, mu0_hat, folds assignment.
# returns: psi_hat
proposed_est <- function(a, y, delta, parametric, K, pi_hat, mu1_hat, mu0_hat, folds) {
  
  n <- length(a)
  phi_vals <- numeric(n)
  if (parametric == 1) {
    #Parametric nuisance estimation
    
    # same denominator used for time-dependent weights and pseudo-outcome
    denom <- delta * pi_hat + (1 - pi_hat)
    
    # W_t
    w_t = (delta * a + 1 - a) / denom
    
    # first term of phi_val
    ipw_term <- w_t * y
    
    # pseudo-outcome R_t
    r_t <-(delta * pi_hat * mu1_hat + (1-pi_hat)*mu0_hat) / denom
    
    V_t_num <- (a * (1- pi_hat) - (1-a)*delta * pi_hat) 
    V_t_den <- delta /(1-delta)
    
    V_t <- V_t_num / V_t_den
    
    correction_term <- w_t * V_t * r_t
    
    phi_vals <- ipw_term + correction_term
    psi_hat <- mean(phi_vals)
    
  } else {
    psi_k <- numeric(K)
    # we iterate over k-folds
    for (k in 1:K) {
      train <- which(folds != k)
      test <- which(folds == k)
      
      a_test <- a[test]
      y_test <- y[test]
      pi_hat_test <- pi_hat[test]
      mu1_hat_test <- mu1_hat[test]
      mu0_hat_test<- mu0_hat[test]
      
      
      
      
      # same denominator used for time-dependent weights and pseudo-outcome
      denom <- delta * pi_hat_test + (1 - pi_hat_test)
      
      # time dependent weights Wt in D1
      w_t <- (delta * a_test + 1 - a_test) / denom
      
      # first term of phi_val
      ipw_term <- w_t * y_test
      
      # pseudo-outcome R_t
      r_t <-(delta * pi_hat_test * mu1_hat_test + (1-pi_hat_test)*mu0_hat_test) / denom
      
      V_t_num <- (a_test * (1- pi_hat_test) - (1-a_test)*delta * pi_hat_test) 
      V_t_den <- delta /(1-delta)
      
      V_t <- V_t_num / V_t_den
      
      correction_term <- w_t * V_t * r_t
      
      phi_vals_test <- ipw_term + correction_term
      psi_k[k] <- mean(phi_vals_test)
      
      phi_vals[test] = phi_vals_test
    }
    psi_hat = mean(psi_k)
  }
  
  
  
  
  return(list(psi_hat = psi_hat, phi_vals = phi_vals))
}





######################################################################
# function: delta_sequence
# takes in:  length of the sequence, I
# returns: vector of length I equally spaced on the log scale

delta_sequence <- function(I) {
  log_lower <- -2.3
  log_upper <- 2.3
  
  # equally spaced on log scale
  log_delta_seq <- seq(log_lower, log_upper, length.out = I)  
  
  # back-transform to original scale
  delta_seq <- exp(log_delta_seq)  
  
  return(delta_seq)
}




#######################################################################
# function: psi_true
# takes in: a delta sequence value 
# returns: true incremental effect true_psi
# note: for a given delta finds the true expectation by using very large n
psi_true <- function(Delta_seq) {
  I <- length(Delta_seq)
  true_psi <- numeric(I)
  
  n_true <- 1e7
  X_big <- matrix(rnorm(n_true * 4), ncol = 4)
  
  lin_pred <- -X_big[,1] + 0.5*X_big[,2] - 0.25*X_big[,3] - 0.1*X_big[,4]
  pi_true_big <- plogis(lin_pred)
  
  mu0 <- 200
  mu1 <- mu0 + 10 + 13.7 * (2*X_big[,1] + X_big[,2] + X_big[,3] + X_big[,4])
  
  # to itearate over delta sequence
  for (i in 1:I) {
    
    numerator <- Delta_seq[i] * pi_true_big * mu1 + (1 - pi_true_big)* mu0
    denominator <- Delta_seq[i]* pi_true_big + (1 - pi_true_big)
    
    true_psi[i] <- mean(numerator / denominator)
    
  }
  return(true_psi)
}





######################################################################
# function: ci_coverage
# takes in: phi_vals, psi_hat, true_psi, alpha, number of bootstrap replications B, n
# returns: boolean coverage
# note: phi_vals matrix must be n x I
ci_coverage <- function(phi_vals, psi_hat, true_psi, alpha, B, n) {
  I <- length(true_psi)
  
  # sigma_hat estimate
  sigma_hat <- numeric(I)
  
  # average product
  avg_prod <- numeric(I)
  
  # rademacher random variables
  rad_rvs <- numeric(n)
  
  max_vals <- numeric(B)
  
  
  
  for (i in 1:I) {
    diff_vec <- phi_vals[,i] - psi_hat[i]
    
    # for debugging 
    #cat("Mean of diff_vec:", mean(diff_vec), "\n")
    
    sigma_sqr <- mean(diff_vec^2)
    sigma_hat[i] <- sqrt(sigma_sqr)
  }
  
  for (b in 1:B) {
    rad_rvs <- 2 * rbinom(n, 1, 0.5) -1
    for (i in 1:I) {
      avg_prod[i] <- mean(rad_rvs * ((phi_vals[,i] - psi_hat[i]) / sigma_hat[i]))
    }
    max_vals[b] <- max(abs(sqrt(n) * avg_prod))
  }
  crit_value <- quantile(max_vals, 1-alpha)
  eff_ll <- psi_hat - (crit_value * sigma_hat/ sqrt(n))
  eff_ul <- psi_hat + (crit_value * sigma_hat/ sqrt(n))
  
  
  if (all(true_psi >= eff_ll & true_psi <= eff_ul)) {
    return(1)
  } else {
    return(0)
  }
}



################################################################
# function: compute simulations
# takes in: sample size n, number of delta values I, number of simulations J, 
# number of bootstrap replications B
# returns: results matrix, ci_matrix
# note: results matrix has the following order
# columns 1-4 propesd estimator, columns 5-8 plug-in, columns 9-12 ipw. 
# 1 is reg param, 2 is trans param, 3 reg nonparam, 4 trans nonparam, same pattern 5-12
# 13 delta_index, 14 sim_index, 15 k, 16 delta_val, 17 true_val. phi_val_matrix is used
# to find unifrom confidence bands
compute_simulation <- function(n, I, J) {
  sim_seed = 70 + J
  set.seed(sim_seed)
  
  delta_seq <- delta_sequence(I=I)
  
  
  # computes the true psi for each delta
  truth <- psi_true(Delta_seq = delta_seq)
  
  results_matrix <- matrix(NA, nrow = I , ncol = 17)
  coverage_vec <- numeric(4)
  
  # each column is an influence function vector from a delta_value. 
  reg_param_influence_matrix <- matrix(NA, nrow = n, ncol = I)
  trans_param_influence_matrix <- matrix(NA, nrow = n, ncol = I)
  reg_nonparam_influence_matrix <- matrix(NA, nrow = n, ncol = I)
  trans_nonparam_influence_matrix <- matrix(NA, nrow = n, ncol = I)
  
  # proposed
  proposed_reg_param_psi_hat <- numeric(I)
  proposed_trans_param_psi_hat <- numeric(I)
  proposed_reg_nonparam_psi_hat <- numeric(I)
  proposed_trans_nonparam_psi_hat <- numeric(I)
  
  #plug-in
  plug_reg_param_psi_hat <- numeric(I)
  plug_trans_param_psi_hat <- numeric(I)
  plug_reg_nonparam_psi_hat <- numeric(I)
  plug_trans_nonparam_psi_hat <- numeric(I)
  
  # ipw 
  ipw_reg_param_psi_hat <- numeric(I)
  ipw_trans_param_psi_hat <- numeric(I)
  ipw_reg_nonparam_psi_hat <- numeric(I)
  ipw_trans_nonparam_psi_hat <- numeric(I)
  
  
  
  # data generation 
  regData <-  generate_data(n=n)
  transData <- data_transformation(dat = regData)
  
  # regular parametric nuisance
  reg_param_nuis <- param_nuisance(dat = regData)
  reg_param_pi_hat <- reg_param_nuis$pi_hat
  reg_param_mu1_hat <- reg_param_nuis$mu1_hat
  reg_param_mu0_hat <- reg_param_nuis$mu0_hat
  
  # transformed parametric nuisance
  trans_param_nuis <- param_nuisance(dat = transData)
  trans_param_pi_hat <- trans_param_nuis$pi_hat
  trans_param_mu1_hat <- trans_param_nuis$mu1_hat
  trans_param_mu0_hat <- trans_param_nuis$mu0_hat
  
  # regular nonparametric nuisance
  reg_nonparam_nuis <- nonparam_nuisance(dat = regData, K=2)
  reg_nonparam_pi_hat <- reg_nonparam_nuis$pi_hat
  reg_nonparam_mu1_hat <- reg_nonparam_nuis$mu1_hat
  reg_nonparam_mu0_hat <- reg_nonparam_nuis$mu0_hat
  reg_folds <- reg_nonparam_nuis$folds
  
  # transformed nonparametric nuisance
  trans_nonparam_nuis <- nonparam_nuisance(dat = transData, K=2)
  trans_nonparam_pi_hat <- trans_nonparam_nuis$pi_hat
  trans_nonparam_mu1_hat <- trans_nonparam_nuis$mu1_hat
  trans_nonparam_mu0_hat <- trans_nonparam_nuis$mu0_hat
  trans_folds <- trans_nonparam_nuis$folds
  
  # treatment and outcome
  a <- regData$a
  y <- regData$y
  
  for (i in 1:I) {
    delta <- delta_seq[i]
    true_psi <- truth[i]
    
    
    
    # adding estimates
    
    # proposed reg param
    result1 <- proposed_est(a = a, y = y, delta = delta, parametric = 1, 
                            K = 2, pi_hat = reg_param_pi_hat, mu1_hat = reg_param_mu1_hat,
                            mu0_hat= reg_param_mu0_hat, folds = reg_folds)
    reg_param_influence_matrix[,i] <- result1$phi_vals
    proposed_reg_param_psi_hat[i] <- result1$psi_hat
    
    # proposed trans param 
    result2 <- proposed_est(a = a, y = y, delta = delta, parametric = 1, 
                            K = 2, pi_hat = trans_param_pi_hat, mu1_hat = trans_param_mu1_hat,
                            mu0_hat= trans_param_mu0_hat, folds = trans_folds)
    trans_param_influence_matrix[, i] <- result2$phi_vals
    proposed_trans_param_psi_hat[i] <- result2$psi_hat
    
    # proposed reg nonparam 
    result3 <- proposed_est(a = a, y = y, delta = delta, parametric = 0, 
                            K = 2, pi_hat = reg_nonparam_pi_hat, mu1_hat = reg_nonparam_mu1_hat,
                            mu0_hat= reg_nonparam_mu0_hat, folds = reg_folds)
    reg_nonparam_influence_matrix[,i] <- result3$phi_vals
    proposed_reg_nonparam_psi_hat[i] <- result3$psi_hat
    
    # proposed trans nonparam 
    result4 <- proposed_est(a = a, y = y, delta = delta, parametric = 0, 
                            K = 2, pi_hat = trans_nonparam_pi_hat, mu1_hat = trans_nonparam_mu1_hat,
                            mu0_hat= trans_nonparam_mu0_hat, folds = trans_folds)
    trans_nonparam_influence_matrix[,i] <- result4$phi_vals
    proposed_trans_nonparam_psi_hat[i] <- result4$psi_hat
    
    # plug-in reg param 
    plug_reg_param_psi_hat[i] <- plug_in_est(a = a, y = y, delta = delta, pi_hat = reg_param_pi_hat, 
                                        mu1_hat = reg_param_mu1_hat, mu0_hat = reg_param_mu0_hat)
    
    # plug-in trans param 
    plug_trans_param_psi_hat[i] <- plug_in_est(a = a, y = y, delta = delta, pi_hat = trans_param_pi_hat, 
                                        mu1_hat = trans_param_mu1_hat, mu0_hat = trans_param_mu0_hat)
    
    # plug-in reg nonparam 
    plug_reg_nonparam_psi_hat[i] <- plug_in_est(a = a, y = y, delta = delta, pi_hat = reg_nonparam_pi_hat, 
                                        mu1_hat = reg_nonparam_mu1_hat, mu0_hat = reg_nonparam_mu0_hat)
    
    # plug-in trans nonparam 
    plug_trans_nonparam_psi_hat[i] <- plug_in_est(a = a, y = y, delta = delta, pi_hat = trans_nonparam_pi_hat, 
                                        mu1_hat = trans_nonparam_mu1_hat, mu0_hat = trans_nonparam_mu0_hat)
    
    # ipw reg param
    ipw_reg_param_psi_hat[i] <- ipw_est(a = a, y = y, delta = delta, pi_hat = reg_param_pi_hat)
    
    # ipw trans param
    ipw_trans_param_psi_hat[i] <- ipw_est(a = a, y = y, delta = delta, pi_hat = trans_param_pi_hat)
    
    # ipw reg nonparam
    ipw_reg_nonparam_psi_hat[i] <- ipw_est(a = a, y = y, delta = delta, pi_hat = reg_nonparam_pi_hat)
    
    # ipw trans nonparam
    ipw_trans_nonparam_psi_hat[i] <- ipw_est(a = a, y = y, delta = delta, pi_hat = trans_nonparam_pi_hat)
    
 
  }
  
  # multiplier bootstrap code for each simulation 
  # influence function matrices and psi hats are obtained earlier in the inner loop 
  proposed_reg_param_coverage <- ci_coverage(phi_vals = reg_param_influence_matrix, 
                                             psi_hat = proposed_reg_param_psi_hat, 
                                 true_psi = truth, alpha = 0.05, B=2000, n=n)
  
  proposed_trans_param_coverage <- ci_coverage(phi_vals = trans_param_influence_matrix, 
                                               psi_hat = proposed_trans_param_psi_hat, 
                                 true_psi = truth, alpha = 0.05, B=2000, n=n)
  
  proposed_reg_nonparam_coverage <- ci_coverage(phi_vals = reg_nonparam_influence_matrix, 
                                                psi_hat = proposed_reg_nonparam_psi_hat, 
                                 true_psi = truth, alpha = 0.05, B=2000, n=n)
  
  proposed_trans_nonparam_coverage <- ci_coverage(phi_vals = trans_nonparam_influence_matrix,
                                                  psi_hat = proposed_trans_nonparam_psi_hat, 
                                 true_psi = truth, alpha = 0.05, B=2000, n=n)
  
  
  
  
  
  return(list(
    proposed_reg_param_psi_hat = proposed_reg_param_psi_hat,
    proposed_trans_param_psi_hat = proposed_trans_param_psi_hat,
    proposed_reg_nonparam_psi_hat = proposed_reg_nonparam_psi_hat,
    proposed_trans_nonparam_psi_hat = proposed_trans_nonparam_psi_hat,
    plug_reg_param_psi_hat = plug_reg_param_psi_hat,
    plug_trans_param_psi_hat = plug_trans_param_psi_hat,
    plug_reg_nonparam_psi_hat = plug_reg_nonparam_psi_hat,
    plug_trans_nonparam_psi_hat = plug_trans_nonparam_psi_hat,
    ipw_reg_param_psi_hat = ipw_reg_param_psi_hat,
    ipw_trans_param_psi_hat = ipw_trans_param_psi_hat,
    ipw_reg_nonparam_psi_hat = ipw_reg_nonparam_psi_hat,
    ipw_trans_nonparam_psi_hat = ipw_trans_nonparam_psi_hat,
    proposed_reg_param_coverage = proposed_reg_param_coverage,
    proposed_trans_param_coverage = proposed_trans_param_coverage,
    proposed_reg_nonparam_coverage = proposed_reg_nonparam_coverage,
    proposed_trans_nonparam_coverage = proposed_trans_nonparam_coverage))
  
}



sim_id <- as.numeric(commandArgs(TRUE))
data <- compute_simulation(n = 500, I=100, J = sim_id)

filename <- paste0("sim", sim_id, ".rds")
saveRDS(data, file=filename)
