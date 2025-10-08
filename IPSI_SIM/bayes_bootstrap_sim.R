library(BART)
library(SoftBart)


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
# returns: transformed matrix
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




#######################################################################
fit_probit <- function(X, Y, num_tree, num_iter) {
  hypers <- Hypers(X, Y, k = 1/6, num_tree = num_tree, sigma_hat = 1)
  opts   <- Opts(update_sigma = FALSE)
  probit_forest <- MakeForest(hypers, opts)
  
  # store posterior draws of pi(x) = Pr(Y=1|X)
  pi_train <- matrix(nrow = num_iter, ncol = nrow(X))
  
  # initialize
  r <- probit_forest$do_predict(X)
  upper <- ifelse(Y == 0, 0, Inf)
  lower <- ifelse(Y == 0, -Inf, 0)
  Z <- truncnorm::rtruncnorm(n = length(Y), a = lower, b = upper, mean = r, sd = 1)
  
  for(i in 1:num_iter) {
    r <- probit_forest$do_gibbs(X, Z, X, 1)
    Z <- truncnorm::rtruncnorm(n = length(Y), a = lower, b = upper, mean = r, sd = 1)
    pi_train[i, ] <- pnorm(r)  # posterior draw of nuisance function
  }
  pi_hat <- colMeans(pi_train[(num_iter/2 +1): num_iter, ])
  return(list(pi_train = pi_train, pi_hat = pi_hat))
}



#######################################################################
# function: fit_vc_bart
# takes in: prognositc forest alpha_forest, treatment effect forest beta_forst,
# outcome vector y, prognostic covariate matrix prognostic_X, treatment_X, treatment vector Z, num_iter
# returns alpha, beta, sigma, mu
fit_vc_bart <- function(alpha_forest, beta_forest, y, prognostic_X, treatment_X,
                        Z, num_iter) {
  ## Variables to save
  
  n <- nrow(prognostic_X)
  
  alpha_out <- matrix(NA, nrow = num_iter, ncol = n)
  beta_out <- matrix(NA, nrow = num_iter, ncol = n)
  sigma_out <- numeric(num_iter)
  
  
  ## Initializing alpha vector
  alpha <- alpha_forest$do_predict(prognostic_X)
  for(i in 1:num_iter) {
    R <- (y - alpha) / Z
    beta <- beta_forest$do_gibbs_weighted(treatment_X, R, Z^2, treatment_X, 1)
    sigma <- beta_forest$get_sigma()
    alpha_forest$set_sigma(sigma)
    
    R <- (y - Z * beta)
    alpha <- alpha_forest$do_gibbs(prognostic_X, R, prognostic_X, 1)
    
    
    # Save posterior draws on TEST set
    alpha_out[i, ] <- alpha
    beta_out[i, ]  <- beta
    sigma_out[i] <- sigma
  }
  mu_out <- alpha_out + t(Z * t(beta_out))
  
  # burn point
  burn <- floor(num_iter / 2)
  
  return(list(alpha = alpha_out[(burn+1):num_iter, ], beta = beta_out[(burn+1):num_iter,],
              sigma = sigma_out[(burn+1):num_iter], mu = mu_out[(burn+1):num_iter,]))
}


##############################################################################
# fits the propensity score and outcome regression functions
# BART = 1 uses bart package, otherwise BCF is used 
nonparam_nuisance <- function(X, Z, y, BART) {
  n <- length(Z)
  if (BART ==1) {
    probit_fit <- gbart(
      x.train = X,
      y.train = Z,
      type = 'pbart',
      ntree = 200,
      k = 2, 
      nskip = 2000,     # burn-in
      ndpost = 2000,    # posterior samples to keep
      keepevery = 1, # keep every draw
      numcut= 100
    )
    
    
    pi_hat_distribution <- pnorm(probit_fit$yhat.train)  # apply probit link
    pi_hat <- colMeans(pi_hat_distribution)   
    
    
    
    outcome_fit <- gbart(
      x.train = as.matrix(cbind(Z,X)),
      y.train = y,
      type = 'wbart',
      ntree = 200, 
      k = 2,
      nskip = 2000, 
      ndpost = 2000
    )
    
    # predict mu1(x)
    X1 <- as.matrix(cbind(Z = rep(1, n), X))
    mu1_hat <- predict(outcome_fit, newdata = X1)
    
    
    # predict mu0(x)
    X0 <- as.matrix(cbind(Z = rep(0, n), X))
    mu0_hat <- predict(outcome_fit, newdata = X0)
    
    
    
  }  else {
    regular_X_scaled <- apply(X, 2, normalize01)
    fitted_probit <- fit_probit(X=regular_X_scaled, Y=Z, num_tree=50, num_iter = 4000)
    
    pi_hat_distribution <- fitted_probit$pi_train
    pi_hat <- fitted_probit$pi_hat
    
    
    
    
    ##################### end of pi_hat
    #plot(pi_true, pi_hat, pch=20, col=rgb(0,0,1,0.3))
    #abline(0, 1, col="red", lwd=2)
    
    # adding pi_hat as a covariate to the prognostic forest
    prognostic_X <- cbind(X, pi_hat)
    
    
    
    #scaling the prognostic and treatment covariates
    prognostic_X_scaled <- apply(prognostic_X, 2, normalize01)
    treatment_X_scaled <- apply(X, 2, normalize01)
    
    
    
    y_scaled <- as.numeric(scale(y))
    
    y_mean <- mean(y)
    y_sd <- sd(y)
    
    
    # priors from Hahn 2020 section 5.2
    alpha_hypers <- Hypers(
      X = prognostic_X_scaled, 
      Y = y_scaled,
      beta = 2,
      gamma = 0.95,
      num_tree = 200
    )
    
    beta_hypers <- Hypers(
      X = treatment_X_scaled,
      Y = y_scaled,
      num_tree = 50,
      beta = 3,       
      gamma = 0.25,  
    )
    
    opts <- Opts(
      num_burn = 2500, 
      num_save = 2500
    )
    
    alpha_forest <- MakeForest(alpha_hypers, opts)
    beta_forest  <- MakeForest(beta_hypers, opts)
    
    Z_correct <- 0.5 - Z
    
    out <- fit_vc_bart(alpha_forest, beta_forest, y_scaled, prognostic_X_scaled, treatment_X_scaled, Z_correct, num_iter = 4000)
    mu0_hat <- (out$alpha + 0.5 * out$beta) * y_sd + y_mean
    mu1_hat <- (out$alpha - 0.5 * out$beta) * y_sd + y_mean
  }
  
  return(list(pi_hat = pi_hat, pi_hat_distribution = pi_hat_distribution, mu0_hat = mu0_hat, mu1_hat = mu1_hat))
}

# scale
# Quantile normalization to [0,1]
normalize01 <- function(x) {
  # scales by x and evaluates at x
  ecdf(x)(x)   # empirical CDF transform
}


###################################################################################
#
bayes_boot <- function(X, Z, y, B, nuisance_fit, delta_seq ) {
  
  n <- length(Z) 
  I <- length(delta_seq) 
  m_t_matrix <- matrix(NA, nrow = B, ncol = I)
  psi_delta_matrix <- matrix(NA, nrow = B, ncol = I)
  efficient_influence_matrix <- matrix(NA, nrow = B, ncol = I)
  
  
  
  pi_hat <- nuisance_fit$pi_hat_distribution
  mu0_hat <- nuisance_fit$mu0_hat
  mu1_hat <- nuisance_fit$mu1_hat

  
  # burn point
  num_iter <- 4000
  burn <- floor(num_iter / 2)
  size <- num_iter - burn
  
  for (b in 1:B) {
    # dirchlet distribution is the same as the standardized exponential distribution
    
    
    random_exp <- rexp(n)
    
    random_dirch <- random_exp / sum(random_exp)
    
    pi_hat_sample <- sample(1:size, 1)
    mu_sample <- sample(1:size, 1)
      
    # grabs a random pi_hat, mu1_hat and mu0_hat to account for variability in their models
    pi_hat_draw <- pi_hat[pi_hat_sample, ]
    mu0_hat_draw <- mu0_hat[mu_sample, ]
    mu1_hat_draw <- mu1_hat[mu_sample, ]
      
    for (i in 1:I) {
      delta <- delta_seq[i]
      # denominator
      denom <- delta * pi_hat_draw + (1 - pi_hat_draw)
      ratio <- (delta * pi_hat_draw * mu1_hat_draw + (1 - pi_hat_draw) * mu0_hat_draw) / denom
      
      psi_delta_matrix[b, i] <- sum(ratio * random_dirch)
      
      efficient_influence_value <- efficient_influence(Z=Z, y=y, pi_hat=pi_hat_draw,
                                                       delta=delta, m_t_1=mu1_hat_draw,m_t_0=mu0_hat_draw, ratio = ratio)
      
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




##################################################################################
# function: efficient influence
# takes in: Z, y, pi_hat_draw, delta, m_t, t
# returns efficient influence funciton psi_delta
efficient_influence <- function(Z, y, pi_hat, delta, m_t_1, m_t_0, ratio) {
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
  
  psi_hat <- score_term * shift_term * cum_weight_term + y_term - mean(ratio)
  return(psi_hat)
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


######## 
# returns if 95% credible interval contains true psi
coverage <- function(est_matrix, true_vec) {
  psi_hat <- colMeans(est_matrix)
  
  B <- nrow(est_matrix)
  I <- length(true_vec)
  max_vals <- numeric(B)
  
  
  pointwise_coverage <- numeric(I)
  
  for (i in 1:I) {
    eff_ll <- quantile(est_matrix[,i], 0.025)
    eff_ul <- quantile(est_matrix[,i], 0.975)
    if(true_vec[i] >= eff_ll & true_vec[i] <= eff_ul) {
      pointwise_coverage[i] <- 1
    } else {
      pointwise_coverage[i] <- 0
    }
  }
  
  # uniform coverage 
  for (b in 1:B) {
    diff_vec <- abs(est_matrix[b, ] - psi_hat)
    max_vals[b] <- max(diff_vec)
  }
  
  crit_val <- quantile(max_vals, 0.95)
  
  uniform_ll <- psi_hat - crit_val
  uniform_ul <- psi_hat + crit_val
  
  if(all(true_vec >= uniform_ll & true_vec <= uniform_ul)) {
    uniform_coverage <- 1
  } else {
    uniform_coverage <- 0
  }
  return(list( pointwise_coverage = pointwise_coverage, uniform_coverage = uniform_coverage))
}

compute_simulation <- function(n, I, J) {
  sim_seed = 70 + J
  set.seed(sim_seed)
  
  delta_seq <-delta_sequence((I=I))
  true_psi <- psi_true(Delta_seq = delta_seq) 
  
  # computes the true psi for each delta
  truth <- psi_true(Delta_seq = delta_seq)
  
  sim <- generate_data(n=n)
  Z <- sim[,1]
  y <- sim[,2]
  X <- sim[,3:6]
  
  nuisance_fit <- nonparam_nuisance(X = X, Z=Z, y = y, BART = 1)
  psi_matrix <- bayes_boot(X = X, Z = Z, y = y, B = 10000, nuisance_fit = nuisance_fit, delta_seq = delta_seq)
  psi_hat <- colMeans(psi_matrix$psi_delta_matrix)
  efficient_psi_hat <- colMeans(psi_matrix$efficient_influence_matrix)
  reg_coverage <- coverage(est_matrix = psi_matrix$psi_delta_matrix, true_vec = truth)
  reg_efficient_coverage <- coverage(est_matrix = psi_matrix$efficient_influence_matrix, true_vec = truth)
  ############# transformed data
  trans_sim <- data_transformation(dat = sim)
  trans_Z <- trans_sim[,1]
  trans_y <- trans_sim[,2]
  trans_X <- trans_sim[,3:6]
  
  trans_nuisance_fit <- nonparam_nuisance(X = trans_X, Z= trans_Z, y = trans_y, BART = 1)
  trans_psi_matrix <- bayes_boot(X = trans_X, Z = trans_Z, y = trans_y, B = 10000, nuisance_fit = trans_nuisance_fit, delta_seq = delta_seq)
  trans_psi_hat <- colMeans(trans_psi_matrix$psi_delta_matrix)
  trans_efficient_psi_hat <- colMeans(trans_psi_matrix$efficient_influence_matrix)
  trans_coverage <- coverage(est_matrix = trans_psi_matrix$psi_delta_matrix, true_vec = truth)
  trans_efficient_coverage <- coverage(est_matrix = trans_psi_matrix$efficient_influence_matrix, true_vec = truth)
  
  return(list(psi_hat = psi_hat, 
              efficient_psi_hat = efficient_psi_hat,
              reg_coverage = reg_coverage, 
              reg_efficient_coverage = reg_efficient_coverage,
              trans_psi_hat = trans_psi_hat,
              trans_efficient_psi_hat = trans_efficient_psi_hat,
              trans_coverage = trans_coverage,
              trans_efficient_coverage = trans_efficient_coverage))
}


sim_id <- as.numeric(commandArgs(TRUE))
data <- compute_simulation(n = 5000, I=100, J = sim_id)

filename <- paste0("sim", sim_id, ".rds")
saveRDS(data, file=filename)




