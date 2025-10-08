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




#########################################################################
# function: bias_est
# takes in: est_matrix, true_vec, I, J
# returns: bias
bias_est <- function(est_matrix, true_vec, I, J) {
  
  # holds the absolute value of internal loop for averaging
  ans_vec = numeric(I)
  
  
  
  for (i in 1:I) {
    est_vec <- est_matrix[, i]
    diff_vec = est_vec - true_vec[i]
    sum = 0
    for (j in 1:J) {
      sum = sum + diff_vec[j]
    }
    ans_vec[i] <- abs(sum / J)
  }
  bias = mean(ans_vec)
  return(bias)
}




##########################################################################
# function: RMSE_est
# takes in: est_vec, true_vec, I, J, n
# returns: RMSE
RMSE_est <- function(est_matrix, true_vec, I, J, n) {
  
  # holds the absolute value of internal loop for averaging
  ans_vec = numeric(I)
  
  
  
  
  for (i in 1:I) {
    est_vec <- est_matrix[, i]
    diff_vec = est_vec - true_vec[i]
    sum = 0
    for (j in 1:J) {
      sum = sum + diff_vec[j]^2
    }
    ans_vec[i] <- sqrt(sum / J)
  }
  RMSE = sqrt(n) * mean(ans_vec)
  return(RMSE)
}

################### extract sims
data_cleaner <- function(n, I, J) {
  delta <- delta_sequence(I=I)
  true_psi <- psi_true(Delta_seq = delta)
  
  all_sims <- vector("list", J)  # preallocate a list
  
  for (i in 1:J) {
    filename <- paste0("sim", i, ".rds")
    all_sims[[i]] <- readRDS(filename)
  }
  
  # data storage setup
  psi_hat <- matrix(NA, nrow = J, ncol = I)
  efficient_psi_hat <- matrix(NA, nrow = J, ncol = I)
  reg_uniform_coverage <- numeric(J)
  reg_pointwise_coverage <- matrix(NA, nrow = J, ncol = I)
  
  reg_efficient_uniform_coverage <- numeric(J)
  reg_efficient_pointwise_coverage <- matrix(NA, nrow = J, ncol = I)
  
  # trans
  trans_psi_hat <-  matrix(NA, nrow = J, ncol = I)
  trans_efficient_psi_hat <- matrix(NA, nrow = J, ncol = I)
  
  trans_uniform_coverage <- numeric(J)
  trans_pointwise_coverage <- matrix(NA, nrow = J, ncol = I)
  
  trans_efficient_uniform_coverage <- numeric(J)
  trans_efficient_pointwise_coverage <- matrix(NA, nrow = J, ncol = I)
  
  for (j in 1:J) {
    psi_hat[j,] <- all_sims[[j]]$psi_hat
    efficient_psi_hat[j,] <- all_sims[[j]]$efficient_psi_hat
    reg_uniform_coverage[j] <- all_sims[[j]]$reg_coverage$uniform_coverage
    reg_pointwise_coverage[j, ] <- all_sims[[j]]$reg_coverage$pointwise_coverage
    
    reg_efficient_uniform_coverage[j] <- all_sims[[j]]$reg_efficient_coverage$uniform_coverage
    reg_efficient_pointwise_coverage[j, ] <- all_sims[[j]]$reg_efficient_coverage$pointwise_coverage
    
    trans_psi_hat[j,] <- all_sims[[j]]$trans_psi_hat
    trans_efficient_psi_hat[j,] <- all_sims[[j]]$trans_efficient_psi_hat
    trans_uniform_coverage[j] <- all_sims[[j]]$trans_coverage$uniform_coverage
    trans_pointwise_coverage[j, ] <- all_sims[[j]]$trans_coverage$pointwise_coverage
    
    trans_efficient_uniform_coverage[j] <- all_sims[[j]]$trans_efficient_coverage$uniform_coverage
    trans_efficient_pointwise_coverage[j, ] <- all_sims[[j]]$trans_efficient_coverage$pointwise_coverage
  }
  
  # bias
  bias <- bias_est(est_matrix = psi_hat, true_vec = true_psi, I= I, J=J)
  efficient_bias <- bias_est(est_matrix = efficient_psi_hat, true_vec = true_psi, I= I, J=J)
  
  # rmse
  rmse <- RMSE_est(est_matrix = psi_hat, true_vec = true_psi, I= I, J=J, n=n)
  efficient_rmse <- RMSE_est(est_matrix = efficient_psi_hat, true_vec = true_psi, I= I, J=J, n=n)
  
  # trans bias
  trans_bias <- bias_est(est_matrix = trans_psi_hat, true_vec = true_psi, I= I, J=J)
  trans_efficient_bias <- bias_est(est_matrix = trans_efficient_psi_hat, true_vec = true_psi, I= I, J=J)
  
  # trans rmse 
  trans_rmse <- RMSE_est(est_matrix = trans_psi_hat, true_vec = true_psi, I= I, J=J, n=n)
  trans_efficient_rmse <- RMSE_est(est_matrix = trans_efficient_psi_hat, true_vec = true_psi, I= I, J=J, n=n)
  
  return(list(
    reg_uniform_coverage = reg_uniform_coverage,
    reg_pointwise_coverage = reg_pointwise_coverage,
    reg_efficient_uniform_coverage = reg_efficient_uniform_coverage,
    reg_efficient_pointwise_coverage = reg_efficient_pointwise_coverage,
    reg_bias = bias, 
    reg_efficient_bias = efficient_bias, 
    reg_rmse = rmse, 
    reg_efficient_rmse = efficient_rmse,
    trans_uniform_coverage = trans_uniform_coverage,
    trans_pointwise_coverage = trans_pointwise_coverage,
    trans_efficient_uniform_coverage = trans_efficient_uniform_coverage,
    trans_efficient_pointwise_coverage = trans_efficient_pointwise_coverage,
    trans_bias = trans_bias, 
    trans_efficient_bias = trans_efficient_bias,
    trans_rmse = trans_rmse,
    trans_efficient_rmse = trans_efficient_rmse
  ))
}

bart_n3_data <- data_cleaner(n = 5000, I = 100, J = 500)

filename <-"bart_n3_data.rds"
saveRDS(bart_n3_data, file=filename)
