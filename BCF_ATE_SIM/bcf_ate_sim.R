library(BART)
library(SoftBart)



# used in generate_data
g <- function(x) {
  ifelse(x == 1, 2,
         ifelse(x == 2, -1, -4))
}



# Simulating data sets from BCF paper sim setup 
generate_data = function(n, tau, mu){
  
  sim_matrix <- matrix(NA, nrow = n, ncol = 7)
  x = matrix(rnorm(n * 5), nrow=n, ncol=5)
  x[,1] = rnorm(n)
  x[,2] = rnorm(n)
  x[,3] = rnorm(n)
  x[,4] = rbinom(n, 1, 0.5)
  x[,5] = sample(c(1,2,3), n, replace = TRUE)
  u = runif(n,0,1)
  sigma = 1
  # Tau -------- 
  if (tau == 'homogeneous'){tau_x = 3}
  if (tau == 'heterogeneous'){tau_x = 1 + 2*x[,2]*x[,4]} # x[,4] after talking to Eoghan. It makes sense!
  
  # Mu -----------
  if (mu == 'linear'){mu_x = 1 + g(x[,5]) + x[,1]*x[,3]}
  if (mu == 'nonlinear'){ mu_x = -6 + g(x[,5]) + 6*abs(x[,3] - 1)}
  
  
  # Pi and Pihat -------- 
  pi.x = 0.8*pnorm((3*mu_x/sd(mu_x)) - 0.5*x[,1]) + 0.05 + u/10
  z = rbinom(n, 1, pi.x)
  
  
  # Response variable
  y = mu_x + tau_x * z + rnorm(n, sd=sigma)
  # updating columns of sim_matrix
  sim_matrix[,1] <- z
  sim_matrix[,2] <- y
  sim_matrix[,3] <- x[,1]
  sim_matrix[,4] <- x[,2]
  sim_matrix[,5] <- x[,3]
  sim_matrix[,6] <- x[,4]
  sim_matrix[,7] <- x[,5]
  return(sim_matrix)
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
# computes the bayesain bootstrap estimate
bayes_boot <- function(X, Z, y, B, nuisance_fit) {
  
  n <- length(Z) 
  
  ate_est <- numeric(B)
  efficient_ate_est <- numeric(B)
  
  
  
 
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
    
    ate_est[b] <- sum(random_dirch * (mu1_hat_draw - mu0_hat_draw))
    
    w1 <- Z / pi_hat_draw
    w0 <- (1 - Z) / (1 - pi_hat_draw)
    
    w1_hajek <- w1 / sum(w1)
    w0_hajek <- w0 / sum(w0)
    
    # Centered EIF structure with Hájek weights
    eif_value <- w1_hajek * (y - mu1_hat_draw) - 
      w0_hajek * (y - mu0_hat_draw) +
      (mu1_hat_draw - mu0_hat_draw) - mean(mu1_hat_draw - mu0_hat_draw)
    
   
    
    
    efficient_ate_est[b] <- mean(mu1_hat_draw - mu0_hat_draw) + sum(random_dirch * eif_value)
                                                    
    
    
    ## progress update every 100 iterations
    if (b %% 100 == 0) {
      cat("Finished bootstrap", b, "of", B, "\n")
      flush.console() 
    }
    
  }
  
  
  return(list(ate_est = ate_est, efficient_ate_est = efficient_ate_est))
}






compute_simulation <- function(n, J, tau, mu, BART) {
  sim_seed = 70 + J
  set.seed(sim_seed)
  
  sim <- generate_data(n=n, tau = tau, mu = mu)
  Z <- sim[,1]
  y <- sim[,2]
  X <- sim[,3:7]
  
  nuisance_fit <- nonparam_nuisance(X = X, Z=Z, y = y, BART = BART)
  psi_vec <- bayes_boot(X = X, Z = Z, y = y, B = 10000, nuisance_fit = nuisance_fit)
  
  ate_est <- mean(psi_vec$ate_est)
  ate_est_lower <- quantile(psi_vec$ate_est, 0.025)
  ate_est_upper <- quantile(psi_vec$ate_est, 0.975)
  
  efficient_ate_est <- mean(psi_vec$efficient_ate_est)
  efficient_ate_est_lower <- quantile(psi_vec$efficient_ate_est, 0.025)
  efficient_ate_est_upper <- quantile(psi_vec$efficient_ate_est, 0.975)
  
  return(list(ate_est = ate_est,
              ate_est_lower = ate_est_lower,
              ate_est_upper = ate_est_upper,
              efficient_ate_est = efficient_ate_est, 
              efficient_ate_est_lower = efficient_ate_est_lower,
              efficient_ate_est_upper = efficient_ate_est_upper))
}
n = 250
sim_id <- as.numeric(commandArgs(TRUE))


# BART = 1
bart_hom_lin     <- compute_simulation(n = n, J = sim_id, tau = 'homogeneous', mu = 'linear', BART = 1)
bart_hom_nonlin  <- compute_simulation(n = n, J = sim_id, tau = 'homogeneous', mu = 'nonlinear', BART = 1)
bart_het_lin     <- compute_simulation(n = n, J = sim_id, tau = 'heterogeneous', mu = 'linear', BART = 1)
bart_het_nonlin  <- compute_simulation(n = n, J = sim_id, tau = 'heterogeneous', mu = 'nonlinear', BART = 1)

# BART = 0 indicates Horseshoe / BCP
bcp_hom_lin      <- compute_simulation(n = n, J = sim_id, tau = 'homogeneous', mu = 'linear', BART = 0)
bcp_hom_nonlin   <- compute_simulation(n = n, J = sim_id, tau = 'homogeneous', mu = 'nonlinear', BART = 0)
bcp_het_lin      <- compute_simulation(n = n, J = sim_id, tau = 'heterogeneous', mu = 'linear', BART = 0)
bcp_het_nonlin   <- compute_simulation(n = n, J = sim_id, tau = 'heterogeneous', mu = 'nonlinear', BART = 0)

# putting in a list
data <- list(
  bart_hom_lin    = bart_hom_lin,
  bart_hom_nonlin = bart_hom_nonlin,
  bart_het_lin    = bart_het_lin,
  bart_het_nonlin = bart_het_nonlin,
  bcp_hom_lin     = bcp_hom_lin,
  bcp_hom_nonlin  = bcp_hom_nonlin,
  bcp_het_lin     = bcp_het_lin,
  bcp_het_nonlin  = bcp_het_nonlin
)

filename <- paste0("sim", sim_id, ".rds")
saveRDS(data, file=filename)


