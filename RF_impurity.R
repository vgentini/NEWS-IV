################# TO DO LIST

set.seed(123)
library(highfrequency)
library(readxl)
library(xts)
library(data.table)
library(glmnet)
library(zoo)
library(HDeconometrics)
library(dplyr)
library(foreach)
library(doParallel)
library(MCS)
library(forecast)
library(ranger)
library(tidyr)
library(ggplot2)

# Vector of tickers to process
tickers <- c(
  "BBDC4", "CSNA3","GGBR4"
  ,"IBOV_BS"
  ,"ITUB4",
  "BOVA11", "PETR4")

forecast_horizon <-c(
  # "daily",
  # "lead_5"
  # , 
  "lead_10",
  "lead_21"
  )

# Base path up to each ticker’s folder
base_dir <- "~/FEA_USP/Mestrado/Novos/Tese/R/Correto_sqrt_ln/HAR_and_HIV_BASE_NOVA_HARQ/Impurity"

############# Function for regularização  #########################################



#Optionally, wrap the inner loop in a function and compile it
runRollingWindow <- function(model_data,
                             train_size,
                             time_index,
                             n_obs) {
  predictions_ols <- xts(rep(NA, n_obs - train_size), order.by = time_index)
  predictions_wls     <- predictions_ols
  # predictions_lasso   <- predictions_o ols
  # predictions_ridge   <- predictions_ols
  # predictions_adalasso<- predictions_ols
  
  
  # nonzero_coefs_lasso    <- vector("list", length = n_obs - train_size)
  # nonzero_coefs_adalasso <- vector("list", length = n_obs - train_size)
  
  for (i in train_size:(n_obs-1)) {
    
    # Extract and scale training data
    X_train_raw <- as.matrix(model_data[(i - train_size + 1):i, ])
    y_train <- as.matrix(target_variable[(i - train_size + 1):i])
    
    # Calculate scaling parameters
    # meanain <- colMeans(X_train_raw)
    # sdain <- apply(X_train_raw, 2, sd)
    
    # Scale features
    # Xain <- scale(X_train_raw, center = meanain, scale = sdain)
    train_df <- data.frame(y = y_train, X_train_raw)
    
    # Prepare test data
    X_test_raw <- model_data[i + 1, , drop = FALSE]
    # X_test <- scale(X_test_raw, center = meanain, scale = sdain)
    
    # OLS
    fit_ols <- lm(rv_var ~ ., data = train_df)
    predictions_ols[i - train_size + 1] <- predict(fit_ols, newdata = X_test_raw)
    
    # WLS
    weights <- 1 / (abs(fitted(fit_ols)))
    fit_wls <- lm(rv_var ~ ., data = train_df, weights = weights)
    predictions_wls[i - train_size + 1] <- predict(fit_wls, newdata = X_test_raw)
    
    X_test <- as.matrix(X_test_raw)
    
    # Ridge
    # fit_ridge <- ic.glmnet(Xain,
    #                        y_train,
    #                        alpha = 0,
    #                        crit = "bic",
    #                        standardize = FALSE)
    # predictions_ridge[i - train_size + 1] <- predict(fit_ridge, newdata = X_test)
    
    # Lasso
    # fit_lasso <- ic.glmnet(Xain,
    #                        y_train,
    #                        alpha = 1,
    #                        crit = "bic",
    #                        standardize = FALSE)
    # coef_lasso <- as.matrix(coef(fit_lasso))
    # nonzero_indices_lasso <- which(coef_lasso != 0)
    # nonzero_coefs_lasso[[i - train_size + 1]] <- rownames(coef_lasso)[nonzero_indices_lasso]
    # predictions_lasso[i - train_size + 1] <- predict(fit_lasso, newdata = X_test)
    
    
    # Adaptive Lasso: use Ridge coefficients to build penalty factor
    # first_step_coef <- coef(fit_ols)[-1]
    # penalty_factor <- abs(first_step_coef)^(-1)
    # adalasso <- ic.glmnet(Xain,
    #                       y_train,
    #                       crit = "bic",
    #                       penalty.factor = penalty_factor,
    #                       standardize = FALSE)
    # coef_adalasso <- as.matrix(coef(adalasso))
    # nonzero_indices_adalasso <- which(coef_adalasso != 0)
    # nonzero_coefs_adalasso[[i - train_size + 1]] <- rownames(coef_adalasso)[nonzero_indices_adalasso]
    # predictions_adalasso[i - train_size + 1] <- predict(adalasso, newdata = X_test)
    
    
    # Cleanup to prevent memory bloat
    rm(X_train,
       y_train,
       # fit_ridge,
       # fit_lasso,
       # adalasso,
       fit_ols,
       fit_wls
    )
    if (i %% 10 == 0) gc(full = TRUE, verbose = FALSE)  # Lightweight garbage collectionarbage collection
  }
  
  
  list(
    predictions_ols      = predictions_ols,
    # predictions_ridge    = predictions_ridge,
    # predictions_lasso    = predictions_lasso,
    # predictions_adalasso = predictions_adalasso,
    # nonzero_coefs_lasso      = nonzero_coefs_lasso,
    # nonzero_coefs_adalasso   = nonzero_coefs_adalasso,
    predictions_wls      = predictions_wls
    
  )
}







############# Function for Random Forest #########################
# Define the rolling window function
runRollingWindow_rf <- function(model_data,
                                train_size,
                                time_index,
                                n_obs,
                                always_split_vars = NULL) {  # Added parameter
  
  # Initialize predictions and importance storage
  n_predictions <- n_obs - train_size
  predictions_rf <- xts(rep(NA, n_predictions), order.by = time_index)
  importance_list <- list()
  
  # Combine data into a matrix
  full_data <- data.frame(target_variable, model_data)
  full_matrix <- as.matrix(full_data)
  
  # Loop through rolling windows
  for (i in train_size:(n_obs - 1)) {
    idx <- i - train_size + 1
    
    # Training data
    train_idx <- idx:i
    X_train_raw <- full_matrix[train_idx, -1, drop = FALSE]
    y_train <- full_matrix[train_idx, 1, drop = FALSE]
    X_test_raw <- full_matrix[i + 1, -1, drop = FALSE]
    
    # Create Bagging model (mtry = number of predictors)
    fit_rf <- ranger(
      x = X_train_raw, 
      y = y_train,
      num.trees = 1000,
      # mtry = ncol(X_train_raw),  # Key change: use all features for each split
      num.threads = 20,
      seed = 123,
      importance = "impurity",
      always.split.variables = always_split_vars
    )
    
    # Store predictions and importance
    predictions_rf[idx] <- predict(fit_rf, data = X_test_raw)$predictions
    importance_list[[idx]] <- fit_rf$variable.importance
    
    # Cleanup
    rm(fit_rf)
    if (i %% 10 == 0) gc(verbose = FALSE, full = TRUE)
  }
  
  # Calculate importance metrics
  importance_df <- do.call(rbind, importance_list)
  mean_importance <- colMeans(importance_df, na.rm = TRUE)
  total_importance <- sum(mean_importance, na.rm = TRUE)
  mean_importance_perc <- round((mean_importance / total_importance) * 100, 2)
  
  # Return results
  list(
    predictions_rf = predictions_rf,
    raw_importance = mean_importance,
    percentage_importance = mean_importance_perc
  )
}




################### LOOP throught all the tickers daily ####################

for (horizon in forecast_horizon) {
  for (tkr in tickers) { 
    # Create a new environment for this iteration
    setwd(paste0("~/FEA_USP/Mestrado/Novos/Tese/R/Normal/", tkr, "/Bases_RDS_IV_loop"))
    
  
    start_time <-Sys.time()
    
    assign(paste0("RV_", horizon), readRDS(paste0("RV_", horizon, ".rds")))
    assign(paste0("ln_RV_", horizon), readRDS(paste0("ln_RV_", horizon, ".rds")))
    assign(paste0("sqrt_RV_", horizon), readRDS(paste0("sqrt_RV_", horizon, ".rds")))
    
    HAR <- readRDS("HAR.rds")
    regressores_HAR <- readRDS("regressores_HAR.rds")
    regressores_HAR <- regressores_HAR[index(get(paste0("RV_", horizon)))]
    
    topics <- subset(regressores_HAR, select = -c(
      RV_lag1, RV_lag5, RV_lag21, IVAR_1, IVAR_5, IVAR_21
    ))
    
    topics_IV <- subset(regressores_HAR, select = -c(
      RV_lag1, RV_lag5, RV_lag21
    ))
    
    IV <- subset(regressores_HAR, select = c(
      IVAR_1, IVAR_5, IVAR_21
    ))
    
    regressores_HAR <- subset(regressores_HAR, select = -c(
      IVAR_1, IVAR_5, IVAR_21
    ))
    
    HARCJ <- readRDS("HARCJ.rds")
    LHARCJ <- readRDS("LHARCJ.rds")
    HARQ_1 <- readRDS("HARQ_1.rds")
    SHAR <- readRDS("SHAR.rds")
    
    ln_HAR <- readRDS("ln_HAR.rds")
    ln_regressores_HAR <- readRDS("ln_regressores_HAR.rds")
    ln_regressores_HAR <- ln_regressores_HAR[index(get(paste0("ln_RV_", horizon)))]
    
    ln_topics <- subset(ln_regressores_HAR, select = -c(
      RV_lag1, RV_lag5, RV_lag21, IVAR_1, IVAR_5, IVAR_21
    ))
    
    ln_topics_IV <- subset(ln_regressores_HAR, select = -c(
      RV_lag1, RV_lag5, RV_lag21
    ))
    
    ln_IV <- subset(ln_regressores_HAR, select = c(
      IVAR_1, IVAR_5, IVAR_21
    ))
    
    ln_regressores_HAR <- subset(ln_regressores_HAR, select = -c(
      IVAR_1, IVAR_5, IVAR_21
    ))
    
    ln_HARCJ <- readRDS("ln_HARCJ.rds")
    ln_LHARCJ <- readRDS("ln_LHARCJ.rds")
    ln_HARQ_1 <- readRDS("ln_HARQ_1.rds")
    ln_SHAR <- readRDS("ln_SHAR.rds")
    
    sqrt_HAR <- readRDS("sqrt_HAR.rds")
    sqrt_regressores_HAR <- readRDS("sqrt_regressores_HAR.rds")
    sqrt_regressores_HAR <- sqrt_regressores_HAR[index(get(paste0("sqrt_RV_", horizon)))]
    
    sqrt_topics <- subset(sqrt_regressores_HAR, select = -c(
      RV_lag1, RV_lag5, RV_lag21, IVAR_1, IVAR_5, IVAR_21
    ))
    
    sqrt_topics_IV <- subset(sqrt_regressores_HAR, select = -c(
      RV_lag1, RV_lag5, RV_lag21
    ))
    
    sqrt_IV <- subset(sqrt_regressores_HAR, select = c(
      IVAR_1, IVAR_5, IVAR_21
    ))
    
    sqrt_regressores_HAR <- subset(sqrt_regressores_HAR, select = -c(
      IVAR_1, IVAR_5, IVAR_21
    ))
    
    sqrt_HARCJ <- readRDS("sqrt_HARCJ.rds")
    sqrt_LHARCJ <- readRDS("sqrt_LHARCJ.rds")
    sqrt_HARQ_1 <- readRDS("sqrt_HARQ_1.rds")
    sqrt_SHAR <- readRDS("sqrt_SHAR.rds")
    
    
    
    
    # Access objects through loop_env
    rv_var <- get(paste0("RV_", horizon))
    ln_rv_var <- get(paste0("ln_RV_", horizon))
    sqrt_rv_var <- get(paste0("sqrt_RV_", horizon))
    
    
      ############# Carregar bases RV #############
      
      
      level <- list(
        # Regressores first
        # regressores_HAR = regressores_HAR,
        # topics = topics,
        # topics_IV = topics_IV,
        IV = IV,
        
        # Remaining variables
        HAR = HAR,
        HARCJ = HARCJ,
        LHARCJ = LHARCJ,
        HARQ_1 = HARQ_1,
        SHAR = SHAR
      )
      
      
      
      
      train_size <- sum(index(rv_var) <= as.POSIXct("2020-12-31 23:59", format="%Y-%m-%d %H:%M"))
      
      train_rv_var = rv_var[1:train_size, ]
      test_rv_var <- rv_var[(train_size + 1):nrow(rv_var), ] 
      
      gc()
      
      target_variable <- rv_var
      # 2a. If you know which column you want to rename (e.g. the first column):
      colnames(target_variable)[1] <- "rv_var"
      realized_variable <- test_rv_var
      
      # Precompute constant values outside the parallel loop
      n_obs <- nrow(target_variable)
      time_index <- index(target_variable[(train_size + 1):n_obs])
      
      ############# Regularização em level ############################
      results <- list()
      for(model_name in names(level)) {
        model_data <- level[[model_name]]
        out <- runRollingWindow(
          model_data = model_data,
          train_size = train_size,
          time_index = time_index,
          n_obs = n_obs
        )
        results[[model_name]] <- out
      }
      
      nonzero_coefs_lasso_results    <- do.call(cbind, lapply(results, `[[`, "nonzero_coefs_lasso"))
      nonzero_coefs_adalasso_results <- do.call(cbind, lapply(results, `[[`, "nonzero_coefs_adalasso"))
      gc()
      
      
      
      ############# Random Forest em level ##################
      level <- list(
        # Regressores first
        regressores_HAR = regressores_HAR,
        topics = topics,
        topics_IV = topics_IV,
        IV = IV,
        
        # Remaining variables
        HAR = HAR,
        HARCJ = HARCJ,
        LHARCJ = LHARCJ,
        HARQ_1 = HARQ_1,
        SHAR = SHAR
      )
      
      # Initialize results list
      results_rf <- list()
      
      # Main sequential execution
      for (model_name in names(level)) {
        model_data <- level[[model_name]]
        
        # Determine always.split.variables
        if (grepl("_IV$", model_name)) {
          always_vars <- c("IVAR_1", "IVAR_5", "IVAR_21")
        } else {
          model_suffixes <- c("HAR", "HARCJ", "LHARCJ", "HARQ_1", "SHAR")
          pattern <- paste0("_(", paste(model_suffixes, collapse = "|"), ")$")
          
          if (grepl(pattern, model_name)) {
            base_model <- regmatches(model_name, regexpr(pattern, model_name))
            base_model <- gsub("_", "", base_model)
            always_vars <- if (base_model %in% names(level)) colnames(level[[base_model]]) else NULL
          } else {
            always_vars <- NULL
          }
        }
        
        # Execute rolling window
        out <- runRollingWindow_rf(
          model_data = model_data,
          train_size = train_size,
          time_index = time_index,
          n_obs = n_obs,
          always_split_vars = always_vars
        )
        
        # Store results directly with model name as key
        results_rf[[model_name]] <- out
      }
      
      ############# Carregar bases ln_RV #############
      
      
      
      ln <- list(
        # Regressores first
        # ln_regressores_HAR    = ln_regressores_HAR,
        # ln_topics  = ln_topics,
        # ln_topics_IV = ln_topics_IV,
        # ln_IV      = ln_IV,
        
        # Remaining variables
        ln_HAR    = ln_HAR,
        ln_HARCJ  = ln_HARCJ,
        ln_LHARCJ = ln_LHARCJ,
        ln_HARQ_1 = ln_HARQ_1,
        ln_SHAR   = ln_SHAR
      )
      
      train_size <- sum(index(ln_rv_var) <= as.POSIXct("2020-12-31 23:59", format = "%Y-%m-%d %H:%M"))
      
      ln_train_rv_var <- ln_rv_var[1:train_size, ]
      ln_test_rv_var  <- ln_rv_var[(train_size + 1):nrow(ln_rv_var), ]
      
      gc()
      
      target_variable    <- ln_rv_var
      # 2a. If you know which column you want to rename (e.g. the first column):
      colnames(target_variable)[1] <- "rv_var"
      realized_variable <- ln_test_rv_var
      
      # Precompute constant values outside the parallel loop
      n_obs     <- nrow(target_variable)
      time_index <- index(target_variable[(train_size + 1):n_obs])
      
      ############# Regularização em ln ############################
      
      ln_results <- list()
      for(model_name in names(ln)) {
        model_data <- ln[[model_name]]
        out <- runRollingWindow(
          model_data = model_data,
          train_size = train_size,
          time_index = time_index,
          n_obs = n_obs
        )
        ln_results[[model_name]] <- out
      }
      
      ln_nonzero_coefs_lasso_results    <- do.call(cbind, lapply(ln_results, `[[`, "nonzero_coefs_lasso"))
      ln_nonzero_coefs_adalasso_results <- do.call(cbind, lapply(ln_results, `[[`, "nonzero_coefs_adalasso"))
      gc()
      
      
      
      
      ############# Random Forest em ln ##################
      ln <- list(
        # Regressores first
        ln_regressores_HAR    = ln_regressores_HAR,
        ln_topics  = ln_topics,
        ln_topics_IV = ln_topics_IV,
        ln_IV      = ln_IV,
        
        # Remaining variables
        ln_HAR    = ln_HAR,
        ln_HARCJ  = ln_HARCJ,
        ln_LHARCJ = ln_LHARCJ,
        ln_HARQ_1 = ln_HARQ_1,
        ln_SHAR   = ln_SHAR
      )
      
      # Initialize list to store results
      ln_results_rf <- list()
      
      # Loop through each model sequentially
      for (model_name in names(ln)) {
        model_data <- ln[[model_name]]
        
        # First check if model_name ends with _IV
        if (grepl("s_IV$", model_name)) {
          always_vars <- c("IVAR_1","IVAR_5","IVAR_21")
        } else {
          # Determine base model for always.split.variables
          model_suffixes <- c("HAR", "HARCJ", "LHARCJ", "HARQ_1", "SHAR")
          pattern <- paste0("(_ln|s)_(", paste(model_suffixes, collapse = "|"), ")$")
          
          if (grepl(pattern, model_name)) {
            base_model <- regmatches(model_name, regexpr(pattern, model_name))
            base_model <- gsub(pattern, "ln_\\2", base_model)
            
            if (base_model %in% names(ln)) {
              always_vars <- c(colnames(ln[[base_model]]))
            } else {
              always_vars <- NULL
            }
          } else {
            always_vars <- NULL
          }
        }
        
        # Execute rolling window
        out <- runRollingWindow_rf(
          model_data = model_data,
          train_size = train_size,
          time_index = time_index,
          n_obs = n_obs,
          always_split_vars = always_vars
        )
        
        # Store results directly in list with model name as key
        ln_results_rf[[model_name]] <- out
      }
      ############# Carregar bases sqrt_RV #############
      
      
      
      # Agrupar conjuntos de dados em lista única
      sqrt <- list(
        # sqrt_regressores_HAR      = sqrt_regressores_HAR,
        # sqrt_topics    = sqrt_topics,
        # sqrt_topics_IV = sqrt_topics_IV,
        sqrt_IV        = sqrt_IV,
        sqrt_HAR       = sqrt_HAR,
        sqrt_HARCJ     = sqrt_HARCJ,
        sqrt_LHARCJ    = sqrt_LHARCJ,
        sqrt_HARQ_1    = sqrt_HARQ_1,
        sqrt_SHAR      = sqrt_SHAR
      )
      
      # Definir tamanho de treinamento
      train_size <- sum(index(sqrt_rv_var) <= as.POSIXct("2020-12-31 23:59", format = "%Y-%m-%d %H:%M"))
      
      # Dividir em treino e teste
      sqrt_train_rv_var <- sqrt_rv_var[1:train_size, ]
      sqrt_test_rv_var  <- sqrt_rv_var[(train_size + 1):nrow(sqrt_rv_var), ]
      
      gc()
      
      target_variable    <- sqrt_rv_var
      # 2a. If you know which column you want to rename (e.g. the first column):
      colnames(target_variable)[1] <- "rv_var"
      realized_variable  <- sqrt_test_rv_var
      
      # Pré-computar constantes
      n_obs      <- nrow(target_variable)
      time_index <- index(target_variable[(train_size + 1):n_obs])
      
      ############# Regularização em sqrt #############
      
      sqrt_results <- list()
      for(model_name in names(sqrt)) {
        model_data <- sqrt[[model_name]]
        out <- runRollingWindow(
          model_data = model_data,
          train_size = train_size,
          time_index = time_index,
          n_obs = n_obs
        )
        sqrt_results[[model_name]] <- out
      }
      
      sqrt_nonzero_coefs_lasso_results    <- do.call(cbind, lapply(sqrt_results, `[[`, "nonzero_coefs_lasso"))
      sqrt_nonzero_coefs_adalasso_results <- do.call(cbind, lapply(sqrt_results, `[[`, "nonzero_coefs_adalasso"))
      gc()
      
      
      ############# Random Forest em sqrt #############
      sqrt <- list(
        sqrt_regressores_HAR      = sqrt_regressores_HAR,
        sqrt_topics    = sqrt_topics,
        sqrt_topics_IV = sqrt_topics_IV,
        sqrt_IV        = sqrt_IV,
        sqrt_HAR       = sqrt_HAR,
        sqrt_HARCJ     = sqrt_HARCJ,
        sqrt_LHARCJ    = sqrt_LHARCJ,
        sqrt_HARQ_1    = sqrt_HARQ_1,
        sqrt_SHAR      = sqrt_SHAR
      )
      
      sqrt_results_rf <- list()
      
      for (model_name in names(sqrt)) {
        model_data <- sqrt[[model_name]]
        
        # First check if model_name ends with _IV
        if (grepl("s_IV$", model_name)) {
          always_vars <- c("IVAR_1","IVAR_5","IVAR_21")
        } else {
          # Determine base model for always.split.variables
          model_suffixes <- c("HAR", "HARCJ", "LHARCJ", "HARQ_1", "SHAR")
          pattern <- paste0("(_sqrt|s)_(", paste(model_suffixes, collapse = "|"), ")$")
          
          if (grepl(pattern, model_name)) {
            base_model <- regmatches(model_name, regexpr(pattern, model_name))
            base_model <- gsub(pattern, "sqrt_\\1", base_model)
            
            if (base_model %in% names(sqrt)) {
              always_vars <- colnames(sqrt[[base_model]])
            } else {
              always_vars <- NULL
            }
          } else {
            always_vars <- NULL
          }
        }
        
        out <- runRollingWindow_rf(
          model_data = model_data,
          train_size = train_size,
          time_index = time_index,
          n_obs = n_obs,
          always_split_vars = always_vars
        )
        
        sqrt_results_rf[[model_name]] <- out
      }
      ######################  consolidar previsões #############
      # Map each transform name to its two lists (OLS/WLS) and RF list
      transforms <- list(
        level = list(base = results,      rf = results_rf),
        ln    = list(base = ln_results,   rf = ln_results_rf),
        sqrt  = list(base = sqrt_results, rf = sqrt_results_rf)
      )
      
      metrics <- c("predictions")
      all_results <- list()
      
      for (trans in names(transforms)) {
        base_list <- transforms[[trans]]$base
        rf_list   <- transforms[[trans]]$rf
        
        # for each metric, build and merge
        m_list <- list()
        for (m in metrics) {
          # collect all columns for OLS/WLS and RF
          ols <- do.call(cbind, lapply(base_list,   `[[`, paste0(m, "_ols")))
          wls <- do.call(cbind, lapply(base_list,   `[[`, paste0(m, "_wls")))
          rf  <- do.call(cbind, lapply(rf_list,     `[[`, paste0(m, "_rf")))
          
          # name the columns with the source + method
          colnames(ols) <- paste0(names(base_list), "_ols")
          colnames(wls) <- paste0(names(base_list), "_wls")
          colnames(rf)  <- paste0(names(rf_list),   "_rf")
          
          # merge into one xts
          m_list[[m]] <- merge.xts(ols, wls, rf)
        }
        
        all_results[[trans]] <- m_list
        gc()
      }
      
      # You now have:
      # all_results$level$predictions   — the `xts` for level-scale predictions
      # all_results$ln$predictions      — the `xts` for log-scale predictions
      # all_results$sqrt$predictions    — the `xts` for sqrt-scale predictions
      ############### salvar novos arquivos\############
      
      
      end_time <-Sys.time()
      print(end_time - start_time)
      
      out_file <- file.path(base_dir, tkr, "RData_TR", paste0("MCS_", tkr, "_", horizon, ".RData"))
      objs <- setdiff(ls(envir = globalenv()), c("base_dir", "forecast_horizon", "tickers",
                                                 "tkr", "horizon",
                                                 "runRollingWindow_rf",
                                                 "runRollingWindow"))
      
      
      # Save objects
      save(
        list = objs,
        file = out_file
      )
      message("Saved objects to: ", basename(out_file))
      
      # Remove the actual objects from loop_env (if needed)
      rm(list = objs)
      gc(full = TRUE, verbose = FALSE)
  }
}




