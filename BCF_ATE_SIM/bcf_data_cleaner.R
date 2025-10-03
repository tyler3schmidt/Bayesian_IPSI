
data_cleaner <- function(n, J) {
  hom_true_psi <- 3
  het_true_psi <- 1
  
  
  
  all_sims <- vector("list", J)  # preallocate a list
  
  for (i in 1:J) {
    filename <- paste0("sim", i, ".rds")
    all_sims[[i]] <- readRDS(filename)
  }
  
  #################################################### ate
  # data storage setup
  bart_hom_lin_ate_est     <- numeric(J)
  bart_hom_nonlin_ate_est  <- numeric(J)
  bart_het_lin_ate_est     <- numeric(J)
  bart_het_nonlin_ate_est  <- numeric(J)
  
  bcp_hom_lin_ate_est      <- numeric(J)
  bcp_hom_nonlin_ate_est   <- numeric(J)
  bcp_het_lin_ate_est      <- numeric(J)
  bcp_het_nonlin_ate_est   <- numeric(J)
  
  ####################################################### ate lower
  # data storage setup
  bart_hom_lin_ate_est_lower     <- numeric(J)
  bart_hom_nonlin_ate_est_lower  <- numeric(J)
  bart_het_lin_ate_est_lower     <- numeric(J)
  bart_het_nonlin_ate_est_lower  <- numeric(J)
  
  bcp_hom_lin_ate_est_lower      <- numeric(J)
  bcp_hom_nonlin_ate_est_lower   <- numeric(J)
  bcp_het_lin_ate_est_lower      <- numeric(J)
  bcp_het_nonlin_ate_est_lower   <- numeric(J)
  
  ####################################################### ate upper
  # data storage setup
  bart_hom_lin_ate_est_upper     <- numeric(J)
  bart_hom_nonlin_ate_est_upper  <- numeric(J)
  bart_het_lin_ate_est_upper     <- numeric(J)
  bart_het_nonlin_ate_est_upper  <- numeric(J)
  
  bcp_hom_lin_ate_est_upper      <- numeric(J)
  bcp_hom_nonlin_ate_est_upper   <- numeric(J)
  bcp_het_lin_ate_est_upper      <- numeric(J)
  bcp_het_nonlin_ate_est_upper   <- numeric(J)
  
  #**************************************************************************#
  for (j in 1:J) {
    
    # ATE
    bart_hom_lin_ate_est[j]     <- all_sims[[j]]$bart_hom_lin$ate_est
    bart_hom_nonlin_ate_est[j]  <- all_sims[[j]]$bart_hom_nonlin$ate_est
    bart_het_lin_ate_est[j]     <- all_sims[[j]]$bart_het_lin$ate_est
    bart_het_nonlin_ate_est[j]  <- all_sims[[j]]$bart_het_nonlin$ate_est
    
    bcp_hom_lin_ate_est[j]      <- all_sims[[j]]$bcp_hom_lin$ate_est
    bcp_hom_nonlin_ate_est[j]   <- all_sims[[j]]$bcp_hom_nonlin$ate_est
    bcp_het_lin_ate_est[j]      <- all_sims[[j]]$bcp_het_lin$ate_est
    bcp_het_nonlin_ate_est[j]   <- all_sims[[j]]$bcp_het_nonlin$ate_est
    
    ########################################################## ate lower
    bart_hom_lin_ate_est_lower[j]     <- all_sims[[j]]$bart_hom_lin$ate_est_lower
    bart_hom_nonlin_ate_est_lower[j]  <- all_sims[[j]]$bart_hom_nonlin$ate_est_lower
    bart_het_lin_ate_est_lower[j]     <- all_sims[[j]]$bart_het_lin$ate_est_lower
    bart_het_nonlin_ate_est_lower[j]  <- all_sims[[j]]$bart_het_nonlin$ate_est_lower
    
    bcp_hom_lin_ate_est_lower[j]      <- all_sims[[j]]$bcp_hom_lin$ate_est_lower
    bcp_hom_nonlin_ate_est_lower[j]   <- all_sims[[j]]$bcp_hom_nonlin$ate_est_lower
    bcp_het_lin_ate_est_lower[j]      <- all_sims[[j]]$bcp_het_lin$ate_est_lower
    bcp_het_nonlin_ate_est_lower[j]   <- all_sims[[j]]$bcp_het_nonlin$ate_est_lower
    
    ########################################################## ate upper
    bart_hom_lin_ate_est_upper[j]     <- all_sims[[j]]$bart_hom_lin$ate_est_upper
    bart_hom_nonlin_ate_est_upper[j]  <- all_sims[[j]]$bart_hom_nonlin$ate_est_upper
    bart_het_lin_ate_est_upper[j]     <- all_sims[[j]]$bart_het_lin$ate_est_upper
    bart_het_nonlin_ate_est_upper[j]  <- all_sims[[j]]$bart_het_nonlin$ate_est_upper
    
    bcp_hom_lin_ate_est_upper[j]      <- all_sims[[j]]$bcp_hom_lin$ate_est_upper
    bcp_hom_nonlin_ate_est_upper[j]   <- all_sims[[j]]$bcp_hom_nonlin$ate_est_upper
    bcp_het_lin_ate_est_upper[j]      <- all_sims[[j]]$bcp_het_lin$ate_est_upper
    bcp_het_nonlin_ate_est_upper[j]   <- all_sims[[j]]$bcp_het_nonlin$ate_est_upper
  }
  
  
  
  
  
  
  
  
  
  #######################################################################################
  #**************************************************************************************
  #######################################################################################
  
  
  #################################################### efficient ate
  # data storage setup
  bart_hom_lin_efficient_ate     <- numeric(J)
  bart_hom_nonlin_efficient_ate  <- numeric(J)
  bart_het_lin_efficient_ate     <- numeric(J)
  bart_het_nonlin_efficient_ate  <- numeric(J)
  
  bcp_hom_lin_efficient_ate      <- numeric(J)
  bcp_hom_nonlin_efficient_ate   <- numeric(J)
  bcp_het_lin_efficient_ate      <- numeric(J)
  bcp_het_nonlin_efficient_ate   <- numeric(J)
  
  ####################################################### efficient ate lower
  # data storage setup
  bart_hom_lin_efficient_ate_lower     <- numeric(J)
  bart_hom_nonlin_efficient_ate_lower  <- numeric(J)
  bart_het_lin_efficient_ate_lower     <- numeric(J)
  bart_het_nonlin_efficient_ate_lower  <- numeric(J)
  
  bcp_hom_lin_efficient_ate_lower      <- numeric(J)
  bcp_hom_nonlin_efficient_ate_lower   <- numeric(J)
  bcp_het_lin_efficient_ate_lower      <- numeric(J)
  bcp_het_nonlin_efficient_ate_lower   <- numeric(J)
  
  ####################################################### efficient ate upper
  # data storage setup
  bart_hom_lin_efficient_ate_upper     <- numeric(J)
  bart_hom_nonlin_efficient_ate_upper  <- numeric(J)
  bart_het_lin_efficient_ate_upper     <- numeric(J)
  bart_het_nonlin_efficient_ate_upper  <- numeric(J)
  
  bcp_hom_lin_efficient_ate_upper      <- numeric(J)
  bcp_hom_nonlin_efficient_ate_upper   <- numeric(J)
  bcp_het_lin_efficient_ate_upper      <- numeric(J)
  bcp_het_nonlin_efficient_ate_upper   <- numeric(J)
  
  #**************************************************************************#
  for (j in 1:J) {
    
    # Efficient ATE
    bart_hom_lin_efficient_ate[j]     <- all_sims[[j]]$bart_hom_lin$efficient_ate_est
    bart_hom_nonlin_efficient_ate[j]  <- all_sims[[j]]$bart_hom_nonlin$efficient_ate_est
    bart_het_lin_efficient_ate[j]     <- all_sims[[j]]$bart_het_lin$efficient_ate_est
    bart_het_nonlin_efficient_ate[j]  <- all_sims[[j]]$bart_het_nonlin$efficient_ate_est
    
    bcp_hom_lin_efficient_ate[j]      <- all_sims[[j]]$bcp_hom_lin$efficient_ate_est
    bcp_hom_nonlin_efficient_ate[j]   <- all_sims[[j]]$bcp_hom_nonlin$efficient_ate_est
    bcp_het_lin_efficient_ate[j]      <- all_sims[[j]]$bcp_het_lin$efficient_ate_est
    bcp_het_nonlin_efficient_ate[j]   <- all_sims[[j]]$bcp_het_nonlin$efficient_ate_est
    
    ########################################################## efficient ate lower
    bart_hom_lin_efficient_ate_lower[j]     <- all_sims[[j]]$bart_hom_lin$efficient_ate_est_lower
    bart_hom_nonlin_efficient_ate_lower[j]  <- all_sims[[j]]$bart_hom_nonlin$efficient_ate_est_lower
    bart_het_lin_efficient_ate_lower[j]     <- all_sims[[j]]$bart_het_lin$efficient_ate_est_lower
    bart_het_nonlin_efficient_ate_lower[j]  <- all_sims[[j]]$bart_het_nonlin$efficient_ate_est_lower
    
    bcp_hom_lin_efficient_ate_lower[j]      <- all_sims[[j]]$bcp_hom_lin$efficient_ate_est_lower
    bcp_hom_nonlin_efficient_ate_lower[j]   <- all_sims[[j]]$bcp_hom_nonlin$efficient_ate_est_lower
    bcp_het_lin_efficient_ate_lower[j]      <- all_sims[[j]]$bcp_het_lin$efficient_ate_est_lower
    bcp_het_nonlin_efficient_ate_lower[j]   <- all_sims[[j]]$bcp_het_nonlin$efficient_ate_est_lower
    
    ########################################################## efficient ate upper
    bart_hom_lin_efficient_ate_upper[j]     <- all_sims[[j]]$bart_hom_lin$efficient_ate_est_upper
    bart_hom_nonlin_efficient_ate_upper[j]  <- all_sims[[j]]$bart_hom_nonlin$efficient_ate_est_upper
    bart_het_lin_efficient_ate_upper[j]     <- all_sims[[j]]$bart_het_lin$efficient_ate_est_upper
    bart_het_nonlin_efficient_ate_upper[j]  <- all_sims[[j]]$bart_het_nonlin$efficient_ate_est_upper
    
    bcp_hom_lin_efficient_ate_upper[j]      <- all_sims[[j]]$bcp_hom_lin$efficient_ate_est_upper
    bcp_hom_nonlin_efficient_ate_upper[j]   <- all_sims[[j]]$bcp_hom_nonlin$efficient_ate_est_upper
    bcp_het_lin_efficient_ate_upper[j]      <- all_sims[[j]]$bcp_het_lin$efficient_ate_est_upper
    bcp_het_nonlin_efficient_ate_upper[j]   <- all_sims[[j]]$bcp_het_nonlin$efficient_ate_est_upper
  }
  
  
  
  
  ###################################################################################################### end of data cleaning 
  ###################################################################################################### 
  ###################################################################################################### 
  ###################################################################################################### 
  
  bias_results <- list(
    # Regular ATE
    bart_hom_lin        = mean(bart_hom_lin_ate_est - hom_true_psi),
    bart_hom_nonlin     = mean(bart_hom_nonlin_ate_est - hom_true_psi),
    bart_het_lin        = mean(bart_het_lin_ate_est - het_true_psi),
    bart_het_nonlin     = mean(bart_het_nonlin_ate_est - het_true_psi),
    
    bcp_hom_lin         = mean(bcp_hom_lin_ate_est - hom_true_psi),
    bcp_hom_nonlin      = mean(bcp_hom_nonlin_ate_est - hom_true_psi),
    bcp_het_lin         = mean(bcp_het_lin_ate_est - het_true_psi),
    bcp_het_nonlin      = mean(bcp_het_nonlin_ate_est - het_true_psi),
    
    # Efficient ATE
    bart_hom_lin_eff    = mean(bart_hom_lin_efficient_ate - hom_true_psi),
    bart_hom_nonlin_eff = mean(bart_hom_nonlin_efficient_ate - hom_true_psi),
    bart_het_lin_eff    = mean(bart_het_lin_efficient_ate - het_true_psi),
    bart_het_nonlin_eff = mean(bart_het_nonlin_efficient_ate - het_true_psi),
    
    bcp_hom_lin_eff     = mean(bcp_hom_lin_efficient_ate - hom_true_psi),
    bcp_hom_nonlin_eff  = mean(bcp_hom_nonlin_efficient_ate - hom_true_psi),
    bcp_het_lin_eff     = mean(bcp_het_lin_efficient_ate - het_true_psi),
    bcp_het_nonlin_eff  = mean(bcp_het_nonlin_efficient_ate - het_true_psi)
  )
  
  
  rmse_results <- list(
    # Regular ATE
    bart_hom_lin        = sqrt(mean((bart_hom_lin_ate_est - hom_true_psi)^2)),
    bart_hom_nonlin     = sqrt(mean((bart_hom_nonlin_ate_est - hom_true_psi)^2)),
    bart_het_lin        = sqrt(mean((bart_het_lin_ate_est - het_true_psi)^2)),
    bart_het_nonlin     = sqrt(mean((bart_het_nonlin_ate_est - het_true_psi)^2)),
    
    bcp_hom_lin         = sqrt(mean((bcp_hom_lin_ate_est - hom_true_psi)^2)),
    bcp_hom_nonlin      = sqrt(mean((bcp_hom_nonlin_ate_est - hom_true_psi)^2)),
    bcp_het_lin         = sqrt(mean((bcp_het_lin_ate_est - het_true_psi)^2)),
    bcp_het_nonlin      = sqrt(mean((bcp_het_nonlin_ate_est - het_true_psi)^2)),
    
    # Efficient ATE
    bart_hom_lin_eff    = sqrt(mean((bart_hom_lin_efficient_ate - hom_true_psi)^2)),
    bart_hom_nonlin_eff = sqrt(mean((bart_hom_nonlin_efficient_ate - hom_true_psi)^2)),
    bart_het_lin_eff    = sqrt(mean((bart_het_lin_efficient_ate - het_true_psi)^2)),
    bart_het_nonlin_eff = sqrt(mean((bart_het_nonlin_efficient_ate - het_true_psi)^2)),
    
    bcp_hom_lin_eff     = sqrt(mean((bcp_hom_lin_efficient_ate - hom_true_psi)^2)),
    bcp_hom_nonlin_eff  = sqrt(mean((bcp_hom_nonlin_efficient_ate - hom_true_psi)^2)),
    bcp_het_lin_eff     = sqrt(mean((bcp_het_lin_efficient_ate - het_true_psi)^2)),
    bcp_het_nonlin_eff  = sqrt(mean((bcp_het_nonlin_efficient_ate - het_true_psi)^2))
  )
  
  
  mae_results <- list(
    # Regular ATE
    bart_hom_lin        = mean(abs(bart_hom_lin_ate_est - hom_true_psi)),
    bart_hom_nonlin     = mean(abs(bart_hom_nonlin_ate_est - hom_true_psi)),
    bart_het_lin        = mean(abs(bart_het_lin_ate_est - het_true_psi)),
    bart_het_nonlin     = mean(abs(bart_het_nonlin_ate_est - het_true_psi)),
    
    bcp_hom_lin         = mean(abs(bcp_hom_lin_ate_est - hom_true_psi)),
    bcp_hom_nonlin      = mean(abs(bcp_hom_nonlin_ate_est - hom_true_psi)),
    bcp_het_lin         = mean(abs(bcp_het_lin_ate_est - het_true_psi)),
    bcp_het_nonlin      = mean(abs(bcp_het_nonlin_ate_est - het_true_psi)),
    
    # Efficient ATE
    bart_hom_lin_eff    = mean(abs(bart_hom_lin_efficient_ate - hom_true_psi)),
    bart_hom_nonlin_eff = mean(abs(bart_hom_nonlin_efficient_ate - hom_true_psi)),
    bart_het_lin_eff    = mean(abs(bart_het_lin_efficient_ate - het_true_psi)),
    bart_het_nonlin_eff = mean(abs(bart_het_nonlin_efficient_ate - het_true_psi)),
    
    bcp_hom_lin_eff     = mean(abs(bcp_hom_lin_efficient_ate - hom_true_psi)),
    bcp_hom_nonlin_eff  = mean(abs(bcp_hom_nonlin_efficient_ate - hom_true_psi)),
    bcp_het_lin_eff     = mean(abs(bcp_het_lin_efficient_ate - het_true_psi)),
    bcp_het_nonlin_eff  = mean(abs(bcp_het_nonlin_efficient_ate - het_true_psi))
  )
  
  
  int_length_results <- list(
    # Regular ATE
    bart_hom_lin        = mean(bart_hom_lin_ate_est_upper - bart_hom_lin_ate_est_lower),
    bart_hom_nonlin     = mean(bart_hom_nonlin_ate_est_upper - bart_hom_nonlin_ate_est_lower),
    bart_het_lin        = mean(bart_het_lin_ate_est_upper - bart_het_lin_ate_est_lower),
    bart_het_nonlin     = mean(bart_het_nonlin_ate_est_upper - bart_het_nonlin_ate_est_lower),
    
    bcp_hom_lin         = mean(bcp_hom_lin_ate_est_upper - bcp_hom_lin_ate_est_lower),
    bcp_hom_nonlin      = mean(bcp_hom_nonlin_ate_est_upper - bcp_hom_nonlin_ate_est_lower),
    bcp_het_lin         = mean(bcp_het_lin_ate_est_upper - bcp_het_lin_ate_est_lower),
    bcp_het_nonlin      = mean(bcp_het_nonlin_ate_est_upper - bcp_het_nonlin_ate_est_lower),
    
    # Efficient ATE
    bart_hom_lin_eff    = mean(bart_hom_lin_efficient_ate_upper - bart_hom_lin_efficient_ate_lower),
    bart_hom_nonlin_eff = mean(bart_hom_nonlin_efficient_ate_upper - bart_hom_nonlin_efficient_ate_lower),
    bart_het_lin_eff    = mean(bart_het_lin_efficient_ate_upper - bart_het_lin_efficient_ate_lower),
    bart_het_nonlin_eff = mean(bart_het_nonlin_efficient_ate_upper - bart_het_nonlin_efficient_ate_lower),
    
    bcp_hom_lin_eff     = mean(bcp_hom_lin_efficient_ate_upper - bcp_hom_lin_efficient_ate_lower),
    bcp_hom_nonlin_eff  = mean(bcp_hom_nonlin_efficient_ate_upper - bcp_hom_nonlin_efficient_ate_lower),
    bcp_het_lin_eff     = mean(bcp_het_lin_efficient_ate_upper - bcp_het_lin_efficient_ate_lower),
    bcp_het_nonlin_eff  = mean(bcp_het_nonlin_efficient_ate_upper - bcp_het_nonlin_efficient_ate_lower)
  )

  coverage_results <- list(
  # Regular ATE
  bart_hom_lin        = mean(bart_hom_lin_ate_est_lower <= hom_true_psi & bart_hom_lin_ate_est_upper >= hom_true_psi),
  bart_hom_nonlin     = mean(bart_hom_nonlin_ate_est_lower <= hom_true_psi & bart_hom_nonlin_ate_est_upper >= hom_true_psi),
  bart_het_lin        = mean(bart_het_lin_ate_est_lower <= het_true_psi & bart_het_lin_ate_est_upper >= het_true_psi),
  bart_het_nonlin     = mean(bart_het_nonlin_ate_est_lower <= het_true_psi & bart_het_nonlin_ate_est_upper >= het_true_psi),
  
  bcp_hom_lin         = mean(bcp_hom_lin_ate_est_lower <= hom_true_psi & bcp_hom_lin_ate_est_upper >= hom_true_psi),
  bcp_hom_nonlin      = mean(bcp_hom_nonlin_ate_est_lower <= hom_true_psi & bcp_hom_nonlin_ate_est_upper >= hom_true_psi),
  bcp_het_lin         = mean(bcp_het_lin_ate_est_lower <= het_true_psi & bcp_het_lin_ate_est_upper >= het_true_psi),
  bcp_het_nonlin      = mean(bcp_het_nonlin_ate_est_lower <= het_true_psi & bcp_het_nonlin_ate_est_upper >= het_true_psi),
  
  # Efficient ATE
  bart_hom_lin_eff    = mean(bart_hom_lin_efficient_ate_lower <= hom_true_psi & bart_hom_lin_efficient_ate_upper >= hom_true_psi),
  bart_hom_nonlin_eff = mean(bart_hom_nonlin_efficient_ate_lower <= hom_true_psi & bart_hom_nonlin_efficient_ate_upper >= hom_true_psi),
  bart_het_lin_eff    = mean(bart_het_lin_efficient_ate_lower <= het_true_psi & bart_het_lin_efficient_ate_upper >= het_true_psi),
  bart_het_nonlin_eff = mean(bart_het_nonlin_efficient_ate_lower <= het_true_psi & bart_het_nonlin_efficient_ate_upper >= het_true_psi),
  
  bcp_hom_lin_eff     = mean(bcp_hom_lin_efficient_ate_lower <= hom_true_psi & bcp_hom_lin_efficient_ate_upper >= hom_true_psi),
  bcp_hom_nonlin_eff  = mean(bcp_hom_nonlin_efficient_ate_lower <= hom_true_psi & bcp_hom_nonlin_efficient_ate_upper >= hom_true_psi),
  bcp_het_lin_eff     = mean(bcp_het_lin_efficient_ate_lower <= het_true_psi & bcp_het_lin_efficient_ate_upper >= het_true_psi),
  bcp_het_nonlin_eff  = mean(bcp_het_nonlin_efficient_ate_lower <= het_true_psi & bcp_het_nonlin_efficient_ate_upper >= het_true_psi)
)

  
  
  simulation_results <- list(
    bias       = bias_results,
    rmse       = rmse_results,
    mae        = mae_results,
    coverage   = coverage_results,
    int_length = int_length_results
  )
  
  return(simulation_results)
  
}


n3_data <- data_cleaner(n = 1000, J = 200)

filename <-"bcf_n3_data.rds"
saveRDS(n3_data, file=filename)
