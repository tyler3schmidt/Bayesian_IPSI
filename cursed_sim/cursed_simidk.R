library(SuperLearner)

# simulation from Kennedy (2018)
SL.library <- c("SL.gam", "SL.earth", "SL.ksvm", 
                "SL.ranger", "SL.glm", "SL.glm.interaction", "SL.stepAIC")



###################################################################################
# ********************************************************************************#
###################################################################################
# Start of Functions




##########################################################################
# function: generate_data
# takes in: sample size n and seed 
# returns: simulated matrix, first column treatment, second outcome, and the final four are covariates.
generate_data <- function(n, seed=NULL){
  sim_matrix <- matrix(NA, nrow = n, ncol = 6)
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
  
  # updating columns of sim_matrix
  sim_matrix[,1] <- a
  sim_matrix[,2] <- y
  sim_matrix[,3] <- X[,1]
  sim_matrix[,4] <- X[,2]
  sim_matrix[,5] <- X[,3]
  sim_matrix[,6] <- X[,4]
  return(sim_matrix)
}




####################################################################
# function: data_transformation
# takes in: a matrix dat
# returns: transformed matrix based on kang & schafer (2007)
data_transformation <- function(dat){
  # creating covariates
  n <- nrow(dat)
  Z <- matrix(NA, nrow = n, ncol = 6)
  
  
  
  #transformations
  Z[,1] = dat[,1]
  Z[,2] = dat[,2]
  Z[,3] = exp(dat[,3] / 2)
  Z[,4] = dat[,4] / (1 + exp(dat[,3])) + 10
  Z[,5] = ((dat[,3] * dat[,5]) / 25 + 0.6)^3 
  Z[,6] = (dat[,4] + dat[,6] + 20)^2
  
  
  return(Z)
}



##################################################################################
# function: ipsi_efficient_influence 
# takes in: Z, y, pi_hat_draw, delta, m_t, t
# returns efficient influence function estimate psi_delta
# note: this algebraic form is from the longitudinal form definition simplified 
#       to T = 1.
ipsi_efficient_influence <- function(Z, y, pi_hat, delta, m_t_1, m_t_0, ratio) {
  score_num <- Z *(1- pi_hat) - (1-Z) * delta * pi_hat
  score_denom <- delta / (1-delta)
  score_term <- score_num / score_denom
  
  shift_num <- delta * pi_hat * m_t_1 + (1-pi_hat) * m_t_0
  shift_denom <- delta * pi_hat + 1 - pi_hat
  shift_term <- shift_num / shift_denom
  
  cum_weight_num <- delta * Z + 1 - Z 
  cum_weight_denom <- delta * pi_hat + 1 - pi_hat
  cum_weight_term <- cum_weight_num/ cum_weight_denom
  
  y_num <- (delta * Z + 1 -Z) * y
  y_denom <- delta *pi_hat + 1 - pi_hat
  y_term <- y_num / y_denom
  
  phi_hat <- score_term * shift_term * cum_weight_term + y_term - mean(ratio)
  
  #dQ_num <- Z * delta * pi_hat + (1-Z) * (1-pi_hat)
  #dQ_denom <- delta * pi_hat + 1 - pi_hat
  #dQ <- dQ_num / dQ_denom
  
  #dP <- pi_hat^Z * (1 - pi_hat)^(1 - Z)
  
  # debug
  #weights <- dQ / dP
  
  # Plot the weights
  #plot(weights, type = "h", 
  #     main = "Weights dQ / dP",
  #     xlab = "Index", ylab = "Weight",
  #     col = "blue", lwd = 2)
  #abline(h = 1, col = "red", lty = 2)  # reference line at weight = 1
  
  return(phi_hat)
}




###################################################################################
# function: freq_nonparam_nuisance
# takes in: data frame dat, number of sample splits K
# returns: nuisance function estimates pi_hat, mu1_hat, mu0_hat, and fold assignment folds
freq_nonparam_nuisance <- function(dat, K) {
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
  
  # Plot overlap
  # pi_plot <- ggplot(data.frame(pi_hat = pi_hat, A = factor(dat$a)), aes(x = pi_hat, fill = A)) +
  #  geom_density(alpha = 0.5) +
  #  scale_x_continuous(limits = c(0, 1)) +
  # labs(title = "Propensity Score Overlap by Treatment",
  #     x = "pi_hat",
  #    y = "Density",
  #   fill = "Treatment") +
  #  theme_minimal()
  #print(pi_plot)
  return(list(pi_hat = pi_hat, mu1_hat = mu1_hat, mu0_hat = mu0_hat, folds = folds))
}










freq_nonparam_nuisance2 <- function(dat, K = NULL) {
  n <- nrow(dat)
  
  # propensity model
  sl_pi <- SuperLearner(Y = dat$a,
                        X = dat[, c("X1", "X2", "X3", "X4")],
                        family = binomial(),
                        SL.library = SL.library)
  pi_hat <- predict(sl_pi, newdata = dat[, c("X1", "X2", "X3", "X4")])$pred
  
  # outcome regression
  sl_mu <- SuperLearner(Y = dat$y,
                        X = dat[, c("X1", "X2", "X3", "X4", "a")],
                        family = gaussian(),
                        SL.library = SL.library)
  
  newdata1 <- dat[, c("X1", "X2", "X3", "X4")]; newdata1$a <- 1
  newdata0 <- dat[, c("X1", "X2", "X3", "X4")]; newdata0$a <- 0
  
  mu1_hat <- predict(sl_mu, newdata = newdata1)$pred
  mu0_hat <- predict(sl_mu, newdata = newdata0)$pred
  
  return(list(pi_hat = pi_hat, mu1_hat = mu1_hat, mu0_hat = mu0_hat, folds = rep(1, n)))
}



###################################################
# function: plug_in_est_ipsi
# takes in a data frame, delta_seq, pi_hat, mu1_hat, mu0_hat. 
# returns psi_hat
plug_in_est_ipsi <- function(a, y, delta_seq, nuisance) {
  pi_hat <- nuisance$pi_hat
  mu1_hat <- nuisance$mu1_hat
  mu0_hat <- nuisance$mu0_hat
  
  I <- length(delta_seq)
  psi_hat <- numeric(I)
  # plug in estimator calculation
  for (i in 1:I) {
    delta <- delta_seq[i]
    
    numerator <- delta * pi_hat * mu1_hat + (1 - pi_hat) * mu0_hat
    denominator <- delta * pi_hat + 1 - pi_hat
    
    psi_hat[i] <- mean(numerator / denominator)
  }
  
  return(psi_hat)
}




###############################################################
# function: proposed_est_ipsi
# takes in a data frame, delta,
# number of sample splits, pi_hat, mu1_hat, mu0_hat, folds assignment.
# returns: psi_hat
proposed_est_ipsi <- function(a, y, delta_seq, K, nuisance) {
  pi_hat <- nuisance$pi_hat
  mu1_hat <- nuisance$mu1_hat
  mu0_hat <- nuisance$mu0_hat
  folds <- nuisance$folds
  
  n <- length(a)
  I <- length(delta_seq)
  
  phi_vals <- matrix(NA, nrow = n, ncol = I)
  psi_hat <- numeric(I)
  
  for (i in 1:I) {
    delta <- delta_seq[i]
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
      
      phi_vals[test, i] = phi_vals_test
    }
    psi_hat[i] = mean(psi_k)
    
    
  }
  
  
  return(list(psi_hat = psi_hat, phi_vals = phi_vals))
}






######################################################################
# function: ci_coverage
# takes in: phi_vals, psi_hat, true_psi, alpha, number of bootstrap replications B, n
# returns: boolean coverage
# note: phi_vals matrix must be n x I
frequentist_coverage <- function(phi_vals, psi_hat, true_vec, B) {
  I <- length(true_vec)
  n <- nrow(phi_vals)
  # sigma_hat estimate
  sigma_hat <- numeric(I)
  
  # average product
  avg_prod <- numeric(I)
  
  # rademacher random variables
  rad_rvs <- numeric(n)
  
  max_vals <- numeric(B)
  
  
  pointwise_coverage <- numeric(I)
  pointwise_length <- numeric(I)
  
  for (i in 1:I) {
    diff_vec <- phi_vals[,i] - psi_hat[i]
    
    # for debugging 
    #cat("Mean of diff_vec:", mean(diff_vec), "\n")
    
    sigma_sqr <- mean(diff_vec^2)
    sigma_hat[i] <- sqrt(sigma_sqr)
    
    pointwise_ll <- psi_hat[i] - 1.96 * sigma_hat[i] / sqrt(n) 
    pointwise_ul <- psi_hat[i] + 1.96 * sigma_hat[i] / sqrt(n) 
    
    pointwise_length[i] <- pointwise_ul - pointwise_ll
    
    if(true_vec[i] >= pointwise_ll & true_vec[i] <= pointwise_ul) {
      pointwise_coverage[i] <- 1
    } else {
      pointwise_coverage[i] <- 0
    }
  }
  
  for (b in 1:B) {
    rad_rvs <- 2 * rbinom(n, 1, 0.5) -1
    for (i in 1:I) {
      avg_prod[i] <- mean(rad_rvs * ((phi_vals[,i] - psi_hat[i]) / sigma_hat[i]))
    }
    max_vals[b] <- max(abs(sqrt(n) * avg_prod))
  }
  crit_value <- quantile(max_vals, 0.95)
  uniform_ll <- psi_hat - (crit_value * sigma_hat/ sqrt(n))
  uniform_ul <- psi_hat + (crit_value * sigma_hat/ sqrt(n))
  
  uniform_length = uniform_ul - uniform_ll
  if(all(true_vec >= uniform_ll & true_vec <= uniform_ul)) {
    uniform_coverage <- 1
  } else {
    uniform_coverage <- 0
  }
  
  return(list( pointwise_coverage = pointwise_coverage, uniform_coverage = uniform_coverage, 
               pointwise_length = pointwise_length, uniform_length = uniform_length))
}







######################################################################
# function: delta_sequence
# takes in:  length of the sequence, I
# returns: vector of length I equally spaced on the log scale
# note: for the ipsi case
delta_sequence_ipsi <- function(I) {
  log_lower <- -2.3
  log_upper <- 2.3
  
  # equally spaced on log scale
  log_delta_seq <- seq(log_lower, log_upper, length.out = I)  
  
  # back-transform to original scale
  delta_seq <- exp(log_delta_seq)  
  
  return(delta_seq)
}





#######################################################################
# function: psi_true_ipsi
# takes in: a delta sequence value 
# returns: true incremental effect true_psi
# note: for a given delta finds the true expectation by using very large n
#       uses seed 42
psi_true_ipsi <- function(Delta_seq) {
  set.seed(42)
  I <- length(Delta_seq)
  true_psi <- numeric(I)
  
  n_true <- 1e6
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
# function: coverage
# takes in: (B x I) 
# returns if 95% credible interval contains true psi
bayesian_coverage <- function(est_matrix, true_vec) {
  psi_hat <- colMeans(est_matrix)  # posterior mean at each delta
  
  B <- nrow(est_matrix)
  I <- length(true_vec)
  
  # Posterior standard deviation at each delta (manual)
  sigma_hat <- numeric(I)
  for (i in 1:I) {
    sigma_hat[i] <- sqrt(sum((est_matrix[, i] - psi_hat[i])^2) / (B - 1))
  }
  
  max_vals <- numeric(B)
  
  # Pointwise coverage and length
  pointwise_coverage <- numeric(I)
  pointwise_length <- numeric(I)
  
  for (i in 1:I) {
    eff_ll <- quantile(est_matrix[, i], 0.025)
    eff_ul <- quantile(est_matrix[, i], 0.975)
    pointwise_length[i] <- eff_ul - eff_ll
    pointwise_coverage[i] <- as.numeric(true_vec[i] >= eff_ll & true_vec[i] <= eff_ul)
  }
  
  # Uniform coverage (studentized)
  for (b in 1:B) {
    max_vals[b] <- max(abs(est_matrix[b, ] - psi_hat) / sigma_hat)
  }
  
  crit_val <- quantile(max_vals, 0.95)
  
  uniform_ll <- psi_hat - crit_val * sigma_hat
  uniform_ul <- psi_hat + crit_val * sigma_hat
  uniform_length <- uniform_ul - uniform_ll
  
  uniform_coverage <- as.numeric(all(true_vec >= uniform_ll & true_vec <= uniform_ul))
  
  return(list(
    pointwise_coverage = pointwise_coverage,
    uniform_coverage = uniform_coverage,
    pointwise_length = pointwise_length,
    uniform_length = uniform_length
  ))
}




######################################################################
# function: IF_ipsi
# takes in: delta, pi_hat, mu1_hat, mu0_hat
# returns: phi_vals
# note: plug_in influence function values for multiplier bootstrap
IF_ipsi <- function(delta_seq, nuisance) {
  pi_hat <- nuisance$pi_hat
  mu1_hat <- nuisance$mu1_hat
  mu0_hat <- nuisance$mu0_hat
  
  n <- length(pi_hat)
  I <- length(delta_seq)
  
  phi_vals <- matrix(NA, nrow = n, ncol = I)
  for(i in 1:I) {
    delta <- delta_seq[i]
    
    numerator <- delta * pi_hat * mu1_hat + (1-pi_hat) * mu0_hat
    denominator <- delta * pi_hat + 1 - pi_hat
    phi_vals[, i] <- numerator / denominator 
    
  }
  
  return(phi_vals)
}




###################################################################################
# function: bayes_boot
# takes in: covar X, treatment Z, outcome y, number of bootstrap iterations B,
#           list from nuisance fit function, delta sequence,
#.           intervention (either "ipsi" or "static"), reg_X (non transformed)
# returns: (B x I) matrices psi_delta_matrix and efficient_influence_matrix 
#.          for the plug-in and one-step estimators respectively.
# note:     
bayes_boot <- function(X, Z, y, B, nuisance_fit, delta_seq) {
  
  n <- length(Z) 
  I <- length(delta_seq) 
  psi_delta_matrix <- matrix(NA, nrow = B, ncol = I)
  efficient_influence_matrix <- matrix(NA, nrow = B, ncol = I)
  
  
  
  pi_hat <- nuisance_fit$pi_hat
  mu0_hat <- nuisance_fit$mu0_hat
  mu1_hat <- nuisance_fit$mu1_hat
  
  
  
  
  
  
  for (b in 1:B) {
    # dirchlet distribution is the same as the standardized exponential distribution
    random_exp <- rexp(n)
    random_dirch <- random_exp / sum(random_exp)
    
    
    
    for (i in 1:I) {
      delta <- delta_seq[i]
      # denominator
      denom <- delta * pi_hat + (1 - pi_hat)
      ratio <- (delta * pi_hat * mu1_hat + (1 - pi_hat) * mu0_hat) / denom
      
      psi_delta_matrix[b, i] <- sum(ratio * random_dirch)
      
      efficient_influence_value <- ipsi_efficient_influence(Z=Z, y=y, pi_hat=pi_hat,
                                                            delta=delta, m_t_1=mu1_hat,m_t_0=mu0_hat, ratio = ratio)
      
      efficient_influence_matrix[b, i] <- mean(ratio) + sum(efficient_influence_value * random_dirch)                                                  
    }
    
    ## progress update every 100 iterations
    if (b %% 100 == 0) {
      cat("Finished bootstrap", b, "of", B, "\n")
      flush.console() 
    }
    
  }
  return(list(psi_delta_matrix = psi_delta_matrix, efficient_influence_matrix = efficient_influence_matrix))
}










#####################################################################################
# function: Compute Simulation
# takes in: sample size n, delta_index I, sim number J, Version, intervention
# returns: list of estimate and coverage for plug in and one-step estimators
#Note: The following versions are: "frequentist", "bart", "softbart", "softbcf"
#       static breaks when I isnt 100. can be easily fixed in the function 
# each psi_hat is a vector of length I
compute_simulation <- function(n, I, J, Version) {
  
  ipsi_delta_seq <-delta_sequence_ipsi((I=I))
  
  
  
  
  
  # always uses seed 42 internally
  ipsi_psi_true <- psi_true_ipsi(Delta_seq = ipsi_delta_seq) 
  
  
  # seed setting
  sim_seed = 70 + J
  set.seed(sim_seed)
  
  # data generating
  sim <- generate_data(n=n)
  Z <- sim[,1]
  y <- sim[,2]
  X <- sim[,3:6]
  
  ############# transformed data
  trans_sim <- data_transformation(dat = sim)
  trans_Z <- trans_sim[,1]
  trans_y <- trans_sim[,2]
  trans_X <- trans_sim[,3:6]
  
  reg_dat = data.frame(a = Z, y = y, X1 = X[,1], X2 = X[,2], X3 = X[,3], X4 = X[,4])
  trans_dat = data.frame(a = Z, y = y, X1 = trans_X[,1], X2 = trans_X[,2], 
                         X3 = trans_X[,3], X4 = trans_X[,4])
  
  
  
  if(Version=="frequentist"){
  
  # frequentist uses a different series of functions compared to bayesian 
  
  
  
  reg_freq_nuisance <- freq_nonparam_nuisance(dat = reg_dat, K=2)
  trans_freq_nuisance <- freq_nonparam_nuisance(dat = trans_dat, K=2)
  
  
  
  # psi_hats
  ipsi_reg_psi_hat <- plug_in_est_ipsi(a = Z, y = y, delta_seq = ipsi_delta_seq, nuisance = reg_freq_nuisance)
  ipsi_trans_psi_hat <- plug_in_est_ipsi(a = Z, y = y, delta_seq = ipsi_delta_seq, nuisance = trans_freq_nuisance)
  
  
  ipsi_reg_efficient_psi_hat_container <- proposed_est_ipsi(a = Z, y = y, delta_seq = ipsi_delta_seq, K = 2, nuisance = reg_freq_nuisance)
  ipsi_trans_efficient_psi_hat_container <- proposed_est_ipsi(a = Z, y = y, delta_seq = ipsi_delta_seq, K = 2, nuisance = trans_freq_nuisance)
  
  ipsi_reg_efficient_psi_hat <- ipsi_reg_efficient_psi_hat_container$psi_hat
  ipsi_trans_efficient_psi_hat <- ipsi_trans_efficient_psi_hat_container$psi_hat
  
  
  
  # coverage 
  ipsi_reg_if <- IF_ipsi(delta_seq = ipsi_delta_seq, nuisance = reg_freq_nuisance)
  ipsi_trans_if <- IF_ipsi(delta_seq = ipsi_delta_seq, nuisance = trans_freq_nuisance)
  
  
  
  # ipsi
  ipsi_reg_coverage <- frequentist_coverage(phi_vals = ipsi_reg_if,
                                            psi_hat = ipsi_reg_psi_hat, true_vec = ipsi_psi_true, B = 10000)
  ipsi_trans_coverage <- frequentist_coverage(phi_vals = ipsi_trans_if,
                                              psi_hat = ipsi_trans_psi_hat, true_vec = ipsi_psi_true, B = 10000)
  
  
  
  # efficient ipsi
  ipsi_reg_efficient_coverage <- frequentist_coverage(phi_vals = ipsi_reg_efficient_psi_hat_container$phi_vals,
                                                      psi_hat = ipsi_reg_efficient_psi_hat, true_vec = ipsi_psi_true, B = 10000)
  ipsi_trans_efficient_coverage <- frequentist_coverage(phi_vals = ipsi_trans_efficient_psi_hat_container$phi_vals,
                                                        psi_hat = ipsi_trans_efficient_psi_hat, true_vec = ipsi_psi_true, B = 10000)
  } else if (Version == "bayesian"){
  
  # bayesian fuckery
    reg_freq_nuisance2 <- freq_nonparam_nuisance2(dat = reg_dat, K=2)
    trans_freq_nuisance2 <- freq_nonparam_nuisance2(dat = trans_dat, K=2)
  
  # bayesian bootstrap
  
  # ipsi
  ipsi_reg_psi_matrix <- bayes_boot(X = X, Z = Z, y = y, B = 10000, nuisance_fit = reg_freq_nuisance2,
                                    delta_seq = ipsi_delta_seq)
  ipsi_trans_psi_matrix <- bayes_boot(X = trans_X, Z = Z, y = y, B = 10000, nuisance_fit = trans_freq_nuisance2,
                                      delta_seq = ipsi_delta_seq)
  
  ipsi_reg_psi_hat <- colMeans(ipsi_reg_psi_matrix$psi_delta_matrix)
  ipsi_trans_psi_hat <- colMeans(ipsi_trans_psi_matrix$psi_delta_matrix)
  
  ipsi_reg_efficient_psi_hat <- colMeans(ipsi_reg_psi_matrix$efficient_influence_matrix)
  ipsi_trans_efficient_psi_hat <- colMeans(ipsi_trans_psi_matrix$efficient_influence_matrix)
  
  ipsi_reg_coverage <- bayesian_coverage(est_matrix = ipsi_reg_psi_matrix$psi_delta_matrix, true_vec = ipsi_psi_true)
  ipsi_trans_coverage <- bayesian_coverage(est_matrix = ipsi_trans_psi_matrix$psi_delta_matrix, true_vec = ipsi_psi_true)
  
  ipsi_reg_efficient_coverage <- bayesian_coverage(est_matrix = ipsi_reg_psi_matrix$efficient_influence_matrix, true_vec = ipsi_psi_true)  
  ipsi_trans_efficient_coverage <- bayesian_coverage(est_matrix = ipsi_trans_psi_matrix$efficient_influence_matrix, true_vec = ipsi_psi_true)
  } else {
    cat("not correct version see function header")
  }
  
  return(list(
    # frequentist
    ipsi_reg_psi_hat = ipsi_reg_psi_hat,
    ipsi_reg_efficient_psi_hat = ipsi_reg_efficient_psi_hat,
    ipsi_reg_coverage = ipsi_reg_coverage,
    ipsi_reg_efficient_coverage = ipsi_reg_efficient_coverage,
    ipsi_trans_psi_hat = ipsi_trans_psi_hat,
    ipsi_trans_efficient_psi_hat = ipsi_trans_efficient_psi_hat,
    ipsi_trans_coverage = ipsi_trans_coverage,
    ipsi_trans_efficient_coverage = ipsi_trans_efficient_coverage
  ))
  
}
  # different estimator runs
sim_id <- as.numeric(commandArgs(TRUE))
# use for local test
sim_id <- 42

# manual change for each run
Size = "n1"
n = 500
# Version = frequentist
Version = "frequentist"
data <- compute_simulation(n = n, I=100, J = sim_id, Version = Version)
#filename <- paste0(Version, Size, "_sim", sim_id, ".rds")
#saveRDS(data, file=filename)

# Version = bayesian
Version = "bayesian"
data <- compute_simulation(n = n, I=100, J = sim_id, Version = Version)
#filename <- paste0(Version, Size, "_sim", sim_id, ".rds")
#saveRDS(data, file=filename)