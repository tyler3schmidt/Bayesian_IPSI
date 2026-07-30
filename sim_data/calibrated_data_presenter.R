################################################################################
# SET DIRECTORY AND LOAD DATA
################################################################################

setwd("/Users/tylerschmidt/Work/Git_Repos/Stochastic_Interventions/sim_data")

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
  
  plug_cov <- mean(df[[paste0(prefix, "_pointwise_coverage")]])
  
  plug_len <- mean(df[[paste0(prefix, "_pointwise_length")]])
  
  plug_un_cov <- mean(df[[paste0(prefix, "_uniform_coverage")]])
  
  plug_un_len <- mean(df[[paste0(prefix, "_uniform_length")]])
  
  eif_bias <- mean(df[[paste0(prefix, "_efficient_bias")]])
  
  eif_rmse <- mean(df[[paste0(prefix, "_efficient_rmse")]])
  
  eif_cov <- mean(df[[paste0(prefix, "_efficient_pointwise_coverage")]])
  
  eif_len <- mean(df[[paste0(prefix, "_efficient_pointwise_length")]])
  
  eif_un_cov <- mean(df[[paste0(prefix, "_efficient_uniform_coverage")]])
  
  eif_un_len <- mean(df[[paste0(prefix, "_efficient_uniform_length")]])
  
  return(c(
    plug_bias,
    plug_rmse,
    plug_cov,
    plug_len,
    plug_un_cov,
    plug_un_len,
    eif_bias,
    eif_rmse,
    eif_cov,
    eif_len,
    eif_un_cov,
    eif_un_len
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
        Cov_un_plugin = vals[5],
        Len_un_plugin = vals[6],
        
        Bias_eif = vals[7],
        RMSE_eif = vals[8],
        Cov_eif = vals[9],
        Len_eif = vals[10],
        Cov_un_eif = vals[11],
        Len_un_eif = vals[12]
        
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
      align = "llcccccccccccc",
      escape = FALSE,
      caption = caption,
      col.names = c(
        "$n$", "Method",
        "Bias", "RMSE", "Cov", "Int Len", "Un Cov", "Un Len",
        "Bias", "RMSE", "Cov", "Int Len", "Un Cov", "Un Len"
      )
    ) %>%
    
    add_header_above(
      c(" " = 2, "Plug-in" = 6, "EIF" = 6)
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




# pointwise only tables
################################################################################
# POINTWISE-ONLY TABLE
################################################################################

make_latex_table_pw <- function(df, caption){
  df %>%
    select(n, Method,
           Bias_plugin, RMSE_plugin, Cov_plugin, Len_plugin,
           Bias_eif, RMSE_eif, Cov_eif, Len_eif) %>%
    mutate(across(where(is.numeric), round, 3)) %>%
    kbl(
      format = "latex", booktabs = TRUE,
      align = "llcccccccc", escape = FALSE, caption = caption,
      col.names = c("$n$", "Method",
                    "Bias", "RMSE", "Cov", "Int Len",
                    "Bias", "RMSE", "Cov", "Int Len")
    ) %>%
    add_header_above(c(" " = 2, "Plug-in" = 4, "EIF" = 4)) %>%
    collapse_rows(columns = 1, latex_hline = "major")
}

tab_ipsi_reg_pw    <- make_latex_table_pw(ipsi_reg_df,    "IPSI Regular (pointwise only).")
tab_ipsi_trans_pw  <- make_latex_table_pw(ipsi_trans_df,  "IPSI Transformed (pointwise only).")
tab_static_reg_pw  <- make_latex_table_pw(static_reg_df,  "Static Regular (pointwise only).")
tab_static_trans_pw<- make_latex_table_pw(static_trans_df,"Static Transformed (pointwise only).")

cat(
  tab_ipsi_reg_pw, "\n\n", tab_ipsi_trans_pw, "\n\n",
  tab_static_reg_pw, "\n\n", tab_static_trans_pw,
  file = "simulation_tables_pointwise.tex"
)

tab_ipsi_reg_pw
tab_ipsi_trans_pw
tab_static_reg_pw
tab_static_trans_pw





################################################################################
# UNIFORM-ONLY TABLE (coverage and length only)
################################################################################

make_latex_table_un <- function(df, caption){
  df %>%
    select(n, Method,
           Cov_un_plugin, Len_un_plugin,
           Cov_un_eif, Len_un_eif) %>%
    mutate(across(where(is.numeric), round, 3)) %>%
    kbl(
      format = "latex", booktabs = TRUE,
      align = "llcccc", escape = FALSE, caption = caption,
      col.names = c("$n$", "Method",
                    "Un Cov", "Un Len",
                    "Un Cov", "Un Len")
    ) %>%
    add_header_above(c(" " = 2, "Plug-in" = 2, "EIF" = 2)) %>%
    collapse_rows(columns = 1, latex_hline = "major")
}

tab_ipsi_reg_un     <- make_latex_table_un(ipsi_reg_df,     "IPSI Regular (uniform only).")
tab_ipsi_trans_un   <- make_latex_table_un(ipsi_trans_df,   "IPSI Transformed (uniform only).")
tab_static_reg_un   <- make_latex_table_un(static_reg_df,   "Static Regular (uniform only).")
tab_static_trans_un <- make_latex_table_un(static_trans_df, "Static Transformed (uniform only).")

cat(
  tab_ipsi_reg_un, "\n\n", tab_ipsi_trans_un, "\n\n",
  tab_static_reg_un, "\n\n", tab_static_trans_un,
  file = "simulation_tables_uniform.tex"
)

tab_ipsi_reg_un
tab_ipsi_trans_un
tab_static_reg_un
tab_static_trans_un
