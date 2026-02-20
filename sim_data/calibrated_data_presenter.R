################################################################################
# SET DIRECTORY AND LOAD DATA
################################################################################

setwd("/Users/tylerschmidt/Work/Research/calibrated_stochastic_int/sim_data")

frequentist_n1_data <- readRDS("frequentist_n1_completed_data.rds")
frequentist_n2_data <- readRDS("frequentist_n2_completed_data.rds")
frequentist_n3_data <- readRDS("frequentist_n3_completed_data.rds")

bart_n1_data <- readRDS("bart_n1_completed_data.rds")
bart_n2_data <- readRDS("bart_n2_completed_data.rds")
bart_n3_data <- readRDS("bart_n3_completed_data.rds")

softbart_n1_data <- readRDS("softbart_n1_completed_data.rds")
softbart_n2_data <- readRDS("softbart_n2_completed_data.rds")
softbart_n3_data <- readRDS("softbart_n3_completed_data.rds")

softbcf_n1_data <- readRDS("softbcf_n1_completed_data.rds")
softbcf_n2_data <- readRDS("softbcf_n2_completed_data.rds")
softbcf_n3_data <- readRDS("softbcf_n3_completed_data.rds")

################################################################################
# LOAD LIBRARIES
################################################################################

library(dplyr)
library(kableExtra)

################################################################################
# ORGANIZE DATA
################################################################################

data_list <- list(
  
  "500" = list(
    Frequentist = frequentist_n1_data,
    BART = bart_n1_data,
    SoftBART = softbart_n1_data,
    SoftBCF = softbcf_n1_data
  ),
  
  "1000" = list(
    Frequentist = frequentist_n2_data,
    BART = bart_n2_data,
    SoftBART = softbart_n2_data,
    SoftBCF = softbcf_n2_data
  ),
  
  "5000" = list(
    Frequentist = frequentist_n3_data,
    BART = bart_n3_data,
    SoftBART = softbart_n3_data,
    SoftBCF = softbcf_n3_data
  )
  
)

################################################################################
# FUNCTION TO COMPUTE METRICS
################################################################################

compute_metrics <- function(df, prefix){
  
  plug_bias <- mean(df[[paste0(prefix, "_bias")]])
  
  plug_rmse <- mean(df[[paste0(prefix, "_rmse")]])
  
  plug_cov <- mean(df[[paste0(prefix, "_uniform_coverage")]])
  
  plug_len <- mean(colMeans(df[[paste0(prefix, "_uniform_length")]]))
  
  
  eif_bias <- mean(df[[paste0(prefix, "_efficient_bias")]])
  
  eif_rmse <- mean(df[[paste0(prefix, "_efficient_rmse")]])
  
  eif_cov <- mean(df[[paste0(prefix, "_efficient_uniform_coverage")]])
  
  eif_len <- mean(colMeans(df[[paste0(prefix, "_efficient_uniform_length")]]))
  
  
  return(c(
    plug_bias,
    plug_rmse,
    plug_cov,
    plug_len,
    eif_bias,
    eif_rmse,
    eif_cov,
    eif_len
  ))
  
}

################################################################################
# BUILD TABLE DATAFRAME
################################################################################

build_table_df <- function(prefix){
  
  results <- data.frame()
  
  for(n in names(data_list)){
    
    for(method in names(data_list[[n]])){
      
      df <- data_list[[n]][[method]]
      
      vals <- compute_metrics(df, prefix)
      
      temp <- data.frame(
        
        n = n,
        
        Method = method,
        
        Bias_plugin = vals[1],
        RMSE_plugin = vals[2],
        Cov_plugin = vals[3],
        Len_plugin = vals[4],
        
        Bias_eif = vals[5],
        RMSE_eif = vals[6],
        Cov_eif = vals[7],
        Len_eif = vals[8]
        
      )
      
      results <- rbind(results, temp)
      
    }
    
  }
  
  return(results)
  
}

################################################################################
# FUNCTION TO CONVERT TO LATEX TABLE
################################################################################

make_latex_table <- function(df, caption){
  
  df %>%
    mutate(across(where(is.numeric), round, 3)) %>%
    
    kbl(
      format = "latex",
      booktabs = TRUE,
      align = "llcccccccc",
      escape = FALSE,
      caption = caption,
      col.names = c(
        "$n$", "Method",
        "Bias", "RMSE", "Cov", "Int Len",
        "Bias", "RMSE", "Cov", "Int Len"
      )
    ) %>%
    
    add_header_above(
      c(" " = 2, "Plug-in" = 4, "EIF" = 4)
    ) %>%
    
    collapse_rows(
      columns = 1,
      latex_hline = "major"
    )
  
}


################################################################################
# GENERATE TABLE DATA
################################################################################

ipsi_reg_df <- build_table_df("ipsi_reg")

ipsi_trans_df <- build_table_df("ipsi_trans")

static_reg_df <- build_table_df("static_reg")

static_trans_df <- build_table_df("static_trans")

################################################################################
# GENERATE LATEX TABLES
################################################################################

tab_ipsi_reg <- make_latex_table(
  ipsi_reg_df,
  "IPSI Regular Data Comparison of plug-in and EIF estimators across methods and sample sizes."
)

tab_ipsi_trans <- make_latex_table(
  ipsi_trans_df,
  "IPSI Transformed Data Comparison of plug-in and EIF estimators across methods and sample sizes."
)

tab_static_reg <- make_latex_table(
  static_reg_df,
  "Static Regular Data Comparison of plug-in and EIF estimators across methods and sample sizes."
)

tab_static_trans <- make_latex_table(
  static_trans_df,
  "Static Transformed Data Comparison of plug-in and EIF estimators across methods and sample sizes."
)

################################################################################
# SAVE TO FILE
################################################################################

cat(
  tab_ipsi_reg,
  "\n\n",
  tab_ipsi_trans,
  "\n\n",
  tab_static_reg,
  "\n\n",
  tab_static_trans,
  file = "simulation_tables.tex"
)

################################################################################
# PRINT TABLES IN CONSOLE
################################################################################

tab_ipsi_reg
tab_ipsi_trans
tab_static_reg
tab_static_trans
