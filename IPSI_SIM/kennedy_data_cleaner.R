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
  proposed_reg_param_psi_hat <- matrix(NA, nrow = J, ncol = I)
  proposed_trans_param_psi_hat <- matrix(NA, nrow = J, ncol = I)
  proposed_reg_nonparam_psi_hat <- matrix(NA, nrow = J, ncol = I)
  proposed_trans_nonparam_psi_hat <- matrix(NA, nrow = J, ncol = I)
  
  plug_reg_param_psi_hat <- matrix(NA, nrow = J, ncol = I)
  plug_trans_param_psi_hat <- matrix(NA, nrow = J, ncol = I)
  plug_reg_nonparam_psi_hat <- matrix(NA, nrow = J, ncol = I)
  plug_trans_nonparam_psi_hat <- matrix(NA, nrow = J, ncol = I)
  
  ipw_reg_param_psi_hat <- matrix(NA, nrow = J, ncol = I)
  ipw_trans_param_psi_hat <- matrix(NA, nrow = J, ncol = I)
  ipw_reg_nonparam_psi_hat <- matrix(NA, nrow = J, ncol = I)
  ipw_trans_nonparam_psi_hat <- matrix(NA, nrow = J, ncol = I)
  
  proposed_reg_param_coverage <- numeric(J)
  proposed_trans_param_coverage <- numeric(J)
  proposed_reg_nonparam_coverage <- numeric(J)
  proposed_trans_nonparam_coverage <-numeric(J)
  
  for (j in 1:J) {
    proposed_reg_param_psi_hat[j, ] <- all_sims[[j]]$proposed_reg_param_psi_hat
    proposed_trans_param_psi_hat[j, ] <- all_sims[[j]]$proposed_trans_param_psi_hat
    proposed_reg_nonparam_psi_hat[j, ] <- all_sims[[j]]$proposed_reg_nonparam_psi_hat
    proposed_trans_nonparam_psi_hat[j, ] <- all_sims[[j]]$proposed_trans_nonparam_psi_hat
    
    plug_reg_param_psi_hat[j, ] <- all_sims[[j]]$plug_reg_param_psi_hat
    plug_trans_param_psi_hat[j, ] <- all_sims[[j]]$plug_trans_param_psi_hat
    plug_reg_nonparam_psi_hat[j, ] <- all_sims[[j]]$plug_reg_nonparam_psi_hat
    plug_trans_nonparam_psi_hat[j, ] <- all_sims[[j]]$plug_trans_nonparam_psi_hat
    
    ipw_reg_param_psi_hat[j, ] <- all_sims[[j]]$ipw_reg_param_psi_hat
    ipw_trans_param_psi_hat[j, ] <- all_sims[[j]]$ipw_trans_param_psi_hat
    ipw_reg_nonparam_psi_hat[j, ] <- all_sims[[j]]$ipw_reg_nonparam_psi_hat
    ipw_trans_nonparam_psi_hat[j, ] <- all_sims[[j]]$ipw_trans_nonparam_psi_hat
    
    proposed_reg_param_coverage[j] <- all_sims[[j]]$proposed_reg_param_coverage
    proposed_trans_param_coverage[j] <- all_sims[[j]]$proposed_trans_param_coverage
    proposed_reg_nonparam_coverage[j] <- all_sims[[j]]$proposed_reg_nonparam_coverage
    proposed_trans_nonparam_coverage[j] <- all_sims[[j]]$proposed_trans_nonparam_coverage
  }
  
 # bias 
  proposed_reg_param_bias <- bias_est(est_matrix = proposed_reg_param_psi_hat, true_vec = true_psi, I = I, J = J)
  proposed_trans_param_bias <- bias_est(est_matrix = proposed_trans_param_psi_hat, true_vec = true_psi, I = I, J = J)
  proposed_reg_nonparam_bias <- bias_est(est_matrix = proposed_reg_nonparam_psi_hat, true_vec = true_psi, I = I, J = J)
  proposed_trans_nonparam_bias <- bias_est(est_matrix = proposed_trans_nonparam_psi_hat, true_vec = true_psi, I = I, J = J)
  
  plug_reg_param_bias <- bias_est(est_matrix = plug_reg_param_psi_hat, true_vec = true_psi, I = I, J = J)
  plug_trans_param_bias <- bias_est(est_matrix = plug_trans_param_psi_hat, true_vec = true_psi, I = I, J = J)
  plug_reg_nonparam_bias <- bias_est(est_matrix = plug_reg_nonparam_psi_hat, true_vec = true_psi, I = I, J = J)
  plug_trans_nonparam_bias <- bias_est(est_matrix = plug_trans_nonparam_psi_hat, true_vec = true_psi, I = I, J = J)
  
  ipw_reg_param_bias <- bias_est(est_matrix = ipw_reg_param_psi_hat, true_vec = true_psi, I = I, J = J)
  ipw_trans_param_bias <- bias_est(est_matrix = ipw_trans_param_psi_hat, true_vec = true_psi, I = I, J = J)
  ipw_reg_nonparam_bias <- bias_est(est_matrix = ipw_reg_nonparam_psi_hat, true_vec = true_psi, I = I, J = J)
  ipw_trans_nonparam_bias <- bias_est(est_matrix = ipw_trans_nonparam_psi_hat, true_vec = true_psi, I = I, J = J)
  
  ###########################################################################################################
  # rmse 
  proposed_reg_param_RMSE <- RMSE_est(est_matrix = proposed_reg_param_psi_hat, true_vec = true_psi, I = I, J = J)
  proposed_trans_param_RMSE <- RMSE_est(est_matrix = proposed_trans_param_psi_hat, true_vec = true_psi, I = I, J = J)
  proposed_reg_nonparam_RMSE <- RMSE_est(est_matrix = proposed_reg_nonparam_psi_hat, true_vec = true_psi, I = I, J = J)
  proposed_trans_nonparam_RMSE <- RMSE_est(est_matrix = proposed_trans_nonparam_psi_hat, true_vec = true_psi, I = I, J = J)
  
  plug_reg_param_RMSE <- RMSE_est(est_matrix = plug_reg_param_psi_hat, true_vec = true_psi, I = I, J = J)
  plug_trans_param_RMSE <- RMSE_est(est_matrix = plug_trans_param_psi_hat, true_vec = true_psi, I = I, J = J)
  plug_reg_nonparam_RMSE <- RMSE_est(est_matrix = plug_reg_nonparam_psi_hat, true_vec = true_psi, I = I, J = J)
  plug_trans_nonparam_RMSE <- RMSE_est(est_matrix = plug_trans_nonparam_psi_hat, true_vec = true_psi, I = I, J = J)
  
  ipw_reg_param_RMSE <- RMSE_est(est_matrix = ipw_reg_param_psi_hat, true_vec = true_psi, I = I, J = J)
  ipw_trans_param_RMSE <- RMSE_est(est_matrix = ipw_trans_param_psi_hat, true_vec = true_psi, I = I, J = J)
  ipw_reg_nonparam_RMSE <- RMSE_est(est_matrix = ipw_reg_nonparam_psi_hat, true_vec = true_psi, I = I, J = J)
  ipw_trans_nonparam_RMSE <- RMSE_est(est_matrix = ipw_trans_nonparam_psi_hat, true_vec = true_psi, I = I, J = J)
  
  return(list(
    proposed_reg_param_coverage = proposed_reg_param_coverage,
    proposed_trans_param_coverage = proposed_trans_param_coverage,
    proposed_reg_nonparam_coverage = proposed_reg_nonparam_coverage,
    proposed_trans_nonparam_coverage = proposed_trans_nonparam_coverage,
    proposed_reg_param_bias = proposed_reg_param_bias,
    proposed_trans_param_bias = proposed_trans_param_bias,
    proposed_reg_nonparam_bias = proposed_reg_nonparam_bias,
    proposed_trans_nonparam_bias = proposed_trans_nonparam_bias,
    plug_reg_param_bias = plug_reg_param_bias,
    plug_trans_param_bias = plug_trans_param_bias,
    plug_reg_nonparam_bias = plug_reg_nonparam_bias,
    plug_trans_nonparam_bias = plug_trans_nonparam_bias,
    ipw_reg_param_bias = ipw_reg_param_bias,
    ipw_trans_param_bias = ipw_trans_param_bias,
    ipw_reg_nonparam_bias = ipw_reg_nonparam_bias,
    ipw_trans_nonparam_bias = ipw_trans_nonparam_bias,
    proposed_reg_param_RMSE = proposed_reg_param_RMSE,
    proposed_trans_param_RMSE = proposed_trans_param_RMSE,
    proposed_reg_nonparam_RMSE = proposed_reg_nonparam_RMSE,
    proposed_trans_nonparam_RMSE = proposed_trans_nonparam_RMSE,
    plug_reg_param_RMSE = plug_reg_param_RMSE,
    plug_trans_param_RMSE = plug_trans_param_RMSE,
    plug_reg_nonparam_RMSE = plug_reg_nonparam_RMSE,
    plug_trans_nonparam_RMSE = plug_trans_nonparam_RMSE,
    ipw_reg_param_RMSE = ipw_reg_param_RMSE,
    ipw_trans_param_RMSE = ipw_trans_param_RMSE,
    ipw_reg_nonparam_RMSE = ipw_reg_nonparam_RMSE,
    ipw_trans_nonparam_RMSE = ipw_trans_nonparam_RMSE
  ))
  
}

n1_data <- data_cleaner(n = 500, I = 100, J = 500)

filename <-"n1_data.rds"
saveRDS(n1_data, file=filename)