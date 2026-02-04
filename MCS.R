################# TO DO LIST

set.seed(123)
library(readxl)
library(xts)
library(data.table)
library(dplyr)
library(MCS)
library(tidyr)
library(parallel)



# Vector of tickers to process
tickers <- c(
             # "BBDC4"
             # ,
             # "CSNA3"
             # ,
             "GGBR4"
             ,
             "IBOV_BS"
             ,
             "ITUB4"
             ,
             "BOVA11"
             ,
             "PETR4"
             )

forecast_horizon <-c(
                     # "daily"
                     # ,
                     # "lead_5"
                     # ,
                     # "lead_10"
                     # ,
                     "lead_21"
                     )

# Base path up to each ticker’s folder
base_dir <- "~/FEA_USP/Mestrado/Novos/Tese/R/Correto_sqrt_ln/HAR_and_HIV_BASE_NOVA_HARQ/Impurity"


for (h in forecast_horizon) {
  
  # Determine k based on h
  k <- switch(h,
              "daily" = 1,
              "lead_5" = 5,
              "lead_10" = 10,
              "lead_21" = 21)
  
  for (tkr in tickers) {  # Move ticker loop inside h loop
    in_file <- file.path(base_dir, tkr, "RData_TR", paste0("MCS_", tkr, "_", h, ".RData"))
    # Rest of the code per ticker and h
    
    if (!file.exists(in_file)) {
      warning("File not found, skipping: ", in_file); next
    }
    loaded_names <- load(in_file)
    message("Loaded: ", basename(in_file))
    
    # Remove objects containing "MCS" from the global environment
    rm(list = grep("MCS", ls(), value = TRUE), envir = .GlobalEnv)
    
    
    # test_rv_var <- get(paste0("test_RV_", h))
    # ln_test_rv_var <- get(paste0("ln_test_RV_", h))
    # sqrt_test_rv_var <- get(paste0("sqrt_test_RV_", h))
    start_time <-Sys.time()
        ############# Importance em level #####################################################
    # Extract percentage_importance for all models
    percentage_importance_list <- lapply(results_rf, function(x) x$percentage_importance)
    
    # Ensure the list retains model names (inherited from results_rf)
    names(percentage_importance_list) <- names(results_rf)
    
    # Define the variables to sum for each model
    target_vars <- c(
      "RV_lag1", "RV_lag5", "RV_lag21",
      "C_lag1", "C_lag5", "C_lag21",
      "J_lag1", "J_lag5", "J_lag21", 
      "L_lag1", "L_lag5", "L_lag21",
      "RVQ_lag1", "JN_lag1", "JP_lag1", "BV_lag1"
    )
    
    # Calculate total percentage for target variables in each model
    model_importance <- lapply(percentage_importance_list, function(x) {
      sum(x[names(x) %in% target_vars], na.rm = TRUE)
    })
    
    # Keep the same model names structure
    names(model_importance) <- names(percentage_importance_list)
    
    # Sum IVAR for each model
    IV_importance <- sapply(percentage_importance_list, function(x) {
      sum(x[c("IVAR_1","IVAR_5","IVAR_21")], na.rm = TRUE)
    })
    
    # Convert model_importance from list to numeric vector (same structure as IV_importance)
    model_importance <- unlist(model_importance)
    
    # Create a data frame with model names, model_importance, IV_importance, and their difference
    importance_comparison <- data.frame(
      Model_Importance = round(model_importance, digits = 1),
      IV_Importance = round(IV_importance, digits = 1),
      Topics = round(100 - model_importance - IV_importance, digits=1)
    )
    
    # Print the result
    print(importance_comparison)
############ MCS  ################

cl <- makeCluster(5)
realized_variable <- test_rv_var

predictions <- all_results$level$predictions

na_columns <- colnames(predictions)[colSums(is.na(log(predictions))) > 0]
clean_predictions <- predictions[, !(colnames(log(predictions)) %in% na_columns)]


clean_predictions <- clean_predictions[, !grepl("HARQ|lasso|^IV|regressores_HAR_ols|topics_ols|topics_IV_ols|regressores_HAR_wls|topics_wls|topics_IV_wls", colnames(clean_predictions))]

QLIKE <-  log(clean_predictions) +
  (as.numeric(realized_variable) / clean_predictions)



AE <- abs(sweep(clean_predictions, 1, as.numeric(realized_variable), "-"))
SE <- (sweep(clean_predictions, 1, as.numeric(realized_variable), "-"))^2
PAE <- 100 * abs(sweep(clean_predictions, 1, as.numeric(realized_variable), "-"))/ as.numeric(realized_variable)
HSE <- 100 * (sweep(clean_predictions, 1, as.numeric(realized_variable), "-")/ as.numeric(realized_variable))^2





stopCluster(cl)
rm(cl)
gc()
cl <- makeCluster(5)
MCS_QLIKE <- MCSprocedure(cbind(QLIKE),
                          alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)
stopCluster(cl)
rm(cl)
gc()
cl <- makeCluster(5)
MCS_SE <- MCSprocedure(cbind(SE),
                       alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)

stopCluster(cl)
rm(cl)
gc()
cl <- makeCluster(5)
MCS_AE <- MCSprocedure(cbind(AE),
                       alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)


stopCluster(cl)
rm(cl)
gc()
cl <- makeCluster(5)
MCS_PAE <- MCSprocedure(cbind(PAE),
                        alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)

stopCluster(cl)
rm(cl)
gc()
cl <- makeCluster(5)
MCS_HSE <- MCSprocedure(cbind(HSE),
                        alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)



# Cleanup parallel cluster
stopCluster(cl)
rm(cl)
gc()




########### MCS regressores ################


cl <- makeCluster(5)
realized_variable <- test_rv_var

clean_predictions <- predictions[, !(colnames(log(predictions)) %in% na_columns)]
clean_predictions_regressores_HAR <- clean_predictions[, !grepl("HARQ|topics|IV|regressores_HAR_ols|regressores_HAR_wls", colnames(clean_predictions))]


QLIKE_regressores_HAR <-  log(clean_predictions_regressores_HAR) +
  (as.numeric(realized_variable) / clean_predictions_regressores_HAR)


AE_regressores_HAR <- abs(sweep(clean_predictions_regressores_HAR, 1, as.numeric(realized_variable), "-"))
SE_regressores_HAR <- (sweep(clean_predictions_regressores_HAR, 1, as.numeric(realized_variable), "-"))^2
PAE_regressores_HAR <- 100 * abs(sweep(clean_predictions_regressores_HAR, 1, as.numeric(realized_variable), "-"))/ as.numeric(realized_variable)
HSE_regressores_HAR <- 100 * (sweep(clean_predictions_regressores_HAR, 1, as.numeric(realized_variable), "-")/ as.numeric(realized_variable))^2

stopCluster(cl)
rm(cl)
gc()
cl <- makeCluster(5)
MCS_QLIKE_regressores_HAR <- MCSprocedure(cbind(QLIKE_regressores_HAR),
                               alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)
stopCluster(cl)
rm(cl)
gc()
cl <- makeCluster(5)
MCS_SE_regressores_HAR <- MCSprocedure(cbind(SE_regressores_HAR),
                            alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)

stopCluster(cl)
rm(cl)
gc()
cl <- makeCluster(5)
MCS_AE_regressores_HAR <- MCSprocedure(cbind(AE_regressores_HAR),
                            alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)


stopCluster(cl)
rm(cl)
gc()
cl <- makeCluster(5)
MCS_PAE_regressores_HAR <- MCSprocedure(cbind(PAE_regressores_HAR),
                             alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)

stopCluster(cl)
rm(cl)
gc()
cl <- makeCluster(5)


MCS_HSE_regressores_HAR <- MCSprocedure(cbind(HSE_regressores_HAR),
                             alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)



# Cleanup parallel cluster
stopCluster(cl)
rm(cl)
gc()
cl <- makeCluster(5)



############ MCS topics_IV ################


cl <- makeCluster(5)
realized_variable <- test_rv_var

clean_predictions <- predictions[, !(colnames(log(predictions)) %in% na_columns)]
clean_predictions_TopicsIV <- clean_predictions[
  ,
  !grepl("^HARQ|^IV|regressores_HAR|topics_IV_ols|topics_IV_wls", colnames(clean_predictions)) &
    # If it's a topics column, it must start with topics_IV
    (!grepl("topics", colnames(clean_predictions)) | grepl("^topics_IV", colnames(clean_predictions)))
]


QLIKE_TopicsIV <-  log(clean_predictions_TopicsIV) +
  (as.numeric(realized_variable) / clean_predictions_TopicsIV)


AE_TopicsIV <- abs(sweep(clean_predictions_TopicsIV, 1, as.numeric(realized_variable), "-"))
SE_TopicsIV <- (sweep(clean_predictions_TopicsIV, 1, as.numeric(realized_variable), "-"))^2
PAE_TopicsIV <- 100 * abs(sweep(clean_predictions_TopicsIV, 1, as.numeric(realized_variable), "-"))/ as.numeric(realized_variable)
HSE_TopicsIV <- 100 * (sweep(clean_predictions_TopicsIV, 1, as.numeric(realized_variable), "-")/ as.numeric(realized_variable))^2

gc()
MCS_QLIKE_TopicsIV <- MCSprocedure(cbind(QLIKE_TopicsIV),
                                   alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k)
stopCluster(cl)
rm(cl)
gc()
cl <- makeCluster(5)
MCS_SE_TopicsIV <- MCSprocedure(cbind(SE_TopicsIV),
                                alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k)

stopCluster(cl)
rm(cl)
gc()
cl <- makeCluster(5)
MCS_AE_TopicsIV <- MCSprocedure(cbind(AE_TopicsIV),
                                alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)


stopCluster(cl)
rm(cl)
gc()
cl <- makeCluster(5)
MCS_PAE_TopicsIV <- MCSprocedure(cbind(PAE_TopicsIV),
                                 alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)
stopCluster(cl)
rm(cl)
gc()
cl <- makeCluster(5)
MCS_HSE_TopicsIV <- MCSprocedure(cbind(HSE_TopicsIV),
                                 alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)



# Cleanup parallel cluster
stopCluster(cl)
rm(cl)
gc()




#
############ MCS topics ################

cl <- makeCluster(5)

realized_variable <- test_rv_var

clean_predictions <- predictions[, !(colnames(log(predictions)) %in% na_columns)]
clean_predictions_Topics <- clean_predictions[, !grepl("HARQ|regressores_HAR|IV|topics_wls|topics_ols", colnames(clean_predictions))]


QLIKE_Topics <- log(clean_predictions_Topics) +
  (as.numeric(realized_variable) / clean_predictions_Topics)


AE_Topics <- abs(sweep(clean_predictions_Topics, 1, as.numeric(realized_variable), "-"))
SE_Topics <- (sweep(clean_predictions_Topics, 1, as.numeric(realized_variable), "-"))^2
PAE_Topics <- 100 * abs(sweep(clean_predictions_Topics, 1, as.numeric(realized_variable), "-"))/ as.numeric(realized_variable)
HSE_Topics <- 100 * (sweep(clean_predictions_Topics, 1, as.numeric(realized_variable), "-")/ as.numeric(realized_variable))^2


gc()
MCS_QLIKE_Topics <- MCSprocedure(cbind(QLIKE_Topics),
                                 alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)
stopCluster(cl)
rm(cl)
gc()
cl <- makeCluster(5)
MCS_SE_Topics <- MCSprocedure(cbind(SE_Topics),
                              alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)

stopCluster(cl)
rm(cl)
gc()
cl <- makeCluster(5)
MCS_AE_Topics <- MCSprocedure(cbind(AE_Topics),
                              alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)


stopCluster(cl)
rm(cl)
gc()
cl <- makeCluster(5)
MCS_PAE_Topics <- MCSprocedure(cbind(PAE_Topics),
                               alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)
stopCluster(cl)
rm(cl)
gc()
cl <- makeCluster(5)
MCS_HSE_Topics <- MCSprocedure(cbind(HSE_Topics),
                               alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)


# Cleanup parallel cluster
stopCluster(cl)
rm(cl)
gc()





############ MCS IV ################

cl <- makeCluster(5)

realized_variable <- test_rv_var

clean_predictions <- predictions[, !(colnames(log(predictions)) %in% na_columns)]
clean_predictions_IV <- clean_predictions[, !grepl("HARQ|regressores_HAR|topics", colnames(clean_predictions))]


QLIKE_IV <- log(clean_predictions_IV) +
  (as.numeric(realized_variable) / clean_predictions_IV)


AE_IV <- abs(sweep(clean_predictions_IV, 1, as.numeric(realized_variable), "-"))
SE_IV <- (sweep(clean_predictions_IV, 1, as.numeric(realized_variable), "-"))^2
PAE_IV <- 100 * abs(sweep(clean_predictions_IV, 1, as.numeric(realized_variable), "-"))/ as.numeric(realized_variable)
HSE_IV <- 100 * (sweep(clean_predictions_IV, 1, as.numeric(realized_variable), "-")/ as.numeric(realized_variable))^2


gc()
MCS_QLIKE_IV <- MCSprocedure(cbind(QLIKE_IV),
                             alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)
stopCluster(cl)
rm(cl)
gc()
cl <- makeCluster(5)
MCS_SE_IV <- MCSprocedure(cbind(SE_IV),
                          alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)

stopCluster(cl)
rm(cl)
gc()
cl <- makeCluster(5)
MCS_AE_IV <- MCSprocedure(cbind(AE_IV),
                          alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)


stopCluster(cl)
rm(cl)
gc()
cl <- makeCluster(5)
MCS_PAE_IV <- MCSprocedure(cbind(PAE_IV),
                           alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)
stopCluster(cl)
rm(cl)
gc()
cl <- makeCluster(5)
MCS_HSE_IV <- MCSprocedure(cbind(HSE_IV),
                           alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)



# Cleanup parallel cluster
stopCluster(cl)
rm(cl)
gc()





        ############# Importance em ln #####################################################
# Extract percentage_importance for all models
ln_percentage_importance_list <- lapply(ln_results_rf, function(x) x$percentage_importance)

# Ensure the list retains model names (inherited from results_rf)
names(ln_percentage_importance_list) <- names(ln_results_rf)

# Define the variables to sum for each model
target_vars <- c(
  "RV_lag1", "RV_lag5", "RV_lag21",
  "C_lag1", "C_lag5", "C_lag21",
  "J_lag1", "J_lag5", "J_lag21", 
  "L_lag1", "L_lag5", "L_lag21",
  "RVQ_lag1", "JN_lag1", "JP_lag1", "BV_lag1"
)

# Calculate total percentage for target variables in each model
ln_model_importance <- lapply(ln_percentage_importance_list, function(x) {
  sum(x[names(x) %in% target_vars], na.rm = TRUE)
})

# Keep the same model names structure
names(ln_model_importance) <- names(ln_percentage_importance_list)

# Sum IVAR for each model
ln_IV_importance <- sapply(ln_percentage_importance_list, function(x) {
  sum(x[c("IVAR_1","IVAR_5","IVAR_21")], na.rm = TRUE)
})

# Convert model_importance from list to numeric vector (same structure as IV_importance)
ln_model_importance <- unlist(ln_model_importance)

# Create a data frame with model names, model_importance, IV_importance, and their difference
ln_importance_comparison <- data.frame(
  Model_Importance = round(ln_model_importance, digits = 1),
  IV_Importance = round(ln_IV_importance, digits = 1),
  Topics = round(100 - model_importance - IV_importance, digits=1)
)

# Print the result
print(ln_importance_comparison)
############ ln_MCS ################

cl <- makeCluster(5)

realized_variable <- ln_test_rv_var

ln_predictions <- all_results$ln$predictions
ln_predictions <- ln_predictions[, !grepl("HARQ|regressores_HAR_ols|regressores_HAR_wls|ln_IV|topics_wls|topics_ols", colnames(ln_predictions))]


ln_AE <- abs(sweep(ln_predictions, 1, as.numeric(realized_variable), "-"))
ln_SE <- (sweep(ln_predictions, 1, as.numeric(realized_variable), "-"))^2


gc()
ln_MCS_SE <- MCSprocedure(cbind(ln_SE),
                                 alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)

stopCluster(cl)
rm(cl)
gc()
cl <- makeCluster(5)
ln_MCS_AE <- MCSprocedure(cbind(ln_AE),
                                 alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)

# Cleanup parallel cluster
stopCluster(cl)
rm(cl)
gc()




############ ln_MCS regressores ################

cl <- makeCluster(5)

realized_variable <- ln_test_rv_var

ln_predictions <- all_results$ln$predictions
ln_predictions_regressores_HAR <- ln_predictions[, !grepl("HARQ|topics|IV|regressores_HAR_wls|regressores_HAR_ols", colnames(ln_predictions))]


ln_AE_regressores_HAR <- abs(sweep(ln_predictions_regressores_HAR, 1, as.numeric(realized_variable), "-"))
ln_SE_regressores_HAR <- (sweep(ln_predictions_regressores_HAR, 1, as.numeric(realized_variable), "-"))^2


gc()
ln_MCS_SE_regressores_HAR <- MCSprocedure(cbind(ln_SE_regressores_HAR),
                               alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)

stopCluster(cl)
rm(cl)
gc()
cl <- makeCluster(5)
ln_MCS_AE_regressores_HAR <- MCSprocedure(cbind(ln_AE_regressores_HAR),
                               alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)

# Cleanup parallel cluster
stopCluster(cl)
rm(cl)
gc()




############ ln_MCS topics_IV ################

cl <- makeCluster(5)

realized_variable <- ln_test_rv_var

ln_predictions <- all_results$ln$predictions
ln_predictions_TopicsIV <- ln_predictions[
  ,
  !grepl("^HARQ|ln_IV|regressores_HAR|topics_IV_wls|topics_IV_ols", colnames(ln_predictions)) &
    # If it's a topics column, keep only those starting with ln_topics_IV
    (!grepl("topics", colnames(ln_predictions)) | grepl("^ln_topics_IV", colnames(ln_predictions)))
]



ln_AE_TopicsIV <- abs(sweep(ln_predictions_TopicsIV, 1, as.numeric(realized_variable), "-"))
ln_SE_TopicsIV <- (sweep(ln_predictions_TopicsIV, 1, as.numeric(realized_variable), "-"))^2


gc()
ln_MCS_SE_TopicsIV <- MCSprocedure(cbind(ln_SE_TopicsIV), 
                                   alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)

stopCluster(cl)
rm(cl)
gc()
cl <- makeCluster(5)
ln_MCS_AE_TopicsIV <- MCSprocedure(cbind(ln_AE_TopicsIV), 
                                   alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)

# Cleanup parallel cluster
stopCluster(cl)
rm(cl)
gc()




############ ln_MCS topics ################

cl <- makeCluster(5)

realized_variable <- ln_test_rv_var

ln_predictions <- all_results$ln$predictions
ln_predictions_Topics <- ln_predictions[, !grepl("HARQ|regressores_HAR|IV|topics_wls|topics_ols", colnames(ln_predictions))]


ln_AE_Topics <- abs(sweep(ln_predictions_Topics, 1, as.numeric(realized_variable), "-"))
ln_SE_Topics <- (sweep(ln_predictions_Topics, 1, as.numeric(realized_variable), "-"))^2


gc()
ln_MCS_SE_Topics <- MCSprocedure(cbind(ln_SE_Topics),
                                 alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)

stopCluster(cl)
rm(cl)
gc()
cl <- makeCluster(5)
ln_MCS_AE_Topics <- MCSprocedure(cbind(ln_AE_Topics),
                                 alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)

# Cleanup parallel cluster
stopCluster(cl)
rm(cl)
gc()




############ ln_MCS IV ################

cl <- makeCluster(5)

realized_variable <- ln_test_rv_var

ln_predictions <- all_results$ln$predictions
ln_predictions_IV <- ln_predictions[, !grepl("HARQ|regressores_HAR|topics", colnames(ln_predictions))]


ln_AE_IV <- abs(sweep(ln_predictions_IV, 1, as.numeric(realized_variable), "-"))
ln_SE_IV <- (sweep(ln_predictions_IV, 1, as.numeric(realized_variable), "-"))^2


gc()
ln_MCS_SE_IV <- MCSprocedure(cbind(ln_SE_IV),
                             alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)

stopCluster(cl)
rm(cl)
gc()
cl <- makeCluster(5)
ln_MCS_AE_IV <- MCSprocedure(cbind(ln_AE_IV),
                             alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)

# Cleanup parallel cluster
stopCluster(cl)
rm(cl)
gc()


  
        ############# Importance em sqrt #####################################################
# Extract percentage_importance for all models
sqrt_percentage_importance_list <- lapply(sqrt_results_rf, function(x) x$percentage_importance)

# Ensure the list retains model names (inherited from results_rf)
names(sqrt_percentage_importance_list) <- names(sqrt_results_rf)

# Define the variables to sum for each model
target_vars <- c(
  "RV_lag1", "RV_lag5", "RV_lag21",
  "C_lag1", "C_lag5", "C_lag21",
  "J_lag1", "J_lag5", "J_lag21", 
  "L_lag1", "L_lag5", "L_lag21",
  "RVQ_lag1", "JN_lag1", "JP_lag1", "BV_lag1"
)

# Calculate total percentage for target variables in each model
sqrt_model_importance <- lapply(sqrt_percentage_importance_list, function(x) {
  sum(x[names(x) %in% target_vars], na.rm = TRUE)
})

# Keep the same model names structure
names(sqrt_model_importance) <- names(sqrt_percentage_importance_list)

# Sum IVAR for each model
sqrt_IV_importance <- sapply(sqrt_percentage_importance_list, function(x) {
  sum(x[c("IVAR_1", "IVAR_5", "IVAR_21")], na.rm = TRUE)
})

# Convert model_importance from list to numeric vector (same structure as IV_importance)
sqrt_model_importance <- unlist(sqrt_model_importance)

# Create a data frame with model names, model_importance, IV_importance, and their difference
sqrt_importance_comparison <- data.frame(
  Model_Importance = round(sqrt_model_importance, digits = 1),
  IV_Importance    = round(sqrt_IV_importance, digits = 1),
  Topics           = round(100 - sqrt_model_importance - sqrt_IV_importance, digits = 1)
)

# Print the result
print(sqrt_importance_comparison)

############ sqrt_MCS #############

cl <- makeCluster(5)

sqrt_predictions <- all_results$sqrt$predictions
sqrt_predictions <- sqrt_predictions[, !grepl("HARQ|^sqrt_IV|topics_ols|topics_IV_ols|topics_wls|topics_IV_wls|regressores_HAR_wls|regressores_HAR_ols", colnames(sqrt_predictions))]


realized_variable <- sqrt_test_rv_var

sqrt_SE <- (sweep(sqrt_predictions, 1, as.numeric(realized_variable), "-"))^2
sqrt_AE <- abs(sweep(sqrt_predictions, 1, as.numeric(realized_variable), "-"))





gc()
sqrt_MCS_SE <- MCSprocedure(cbind(sqrt_SE), alpha=0.1, B=10000, cl=cl, statistic="TR", k=k, verbose=TRUE)
stopCluster(cl)
rm(cl)
gc()
cl <- makeCluster(5)
sqrt_MCS_AE <- MCSprocedure(cbind(sqrt_AE), alpha=0.1, B=10000, cl=cl, statistic="TR", k=k, verbose=TRUE)

# Cleanup parallel cluster
stopCluster(cl)
rm(cl)
gc()

########### sqrt_MCS regressores_HAR ################

cl <- makeCluster(5)

realized_variable <- sqrt_test_rv_var

sqrt_predictions <- all_results$sqrt$predictions
sqrt_predictions_regressores_HAR <- sqrt_predictions[, !grepl("HARQ|topics|IV|regressores_HAR_wls|regressores_HAR_ols", colnames(sqrt_predictions))]


sqrt_AE_regressores_HAR <- abs(sweep(sqrt_predictions_regressores_HAR, 1, as.numeric(realized_variable), "-"))
sqrt_SE_regressores_HAR <- (sweep(sqrt_predictions_regressores_HAR, 1, as.numeric(realized_variable), "-"))^2


gc()
sqrt_MCS_SE_regressores_HAR <- MCSprocedure(cbind(sqrt_SE_regressores_HAR),
                                 alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)


stopCluster(cl)
rm(cl)
gc()
cl <- makeCluster(5)
sqrt_MCS_AE_regressores_HAR <- MCSprocedure(cbind(sqrt_AE_regressores_HAR),
                                 alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)

# Cleanup parallel cluster
stopCluster(cl)
rm(cl)
gc()




############ sqrt_MCS topics_IV ################

cl <- makeCluster(5)

realized_variable <- sqrt_test_rv_var

sqrt_predictions <- all_results$sqrt$predictions
sqrt_predictions_TopicsIV <- sqrt_predictions[
  ,
  !grepl("HARQ|sqrt_IV|regressores_HAR|topics_IV_wls|topics_IV_ols", colnames(sqrt_predictions)) &
    # Keep all non-topic cols, but for topics only those starting with sqrt_topics_IV
    (!grepl("topics", colnames(sqrt_predictions)) | grepl("^sqrt_topics_IV", colnames(sqrt_predictions)))
]



sqrt_AE_TopicsIV <- abs(sweep(sqrt_predictions_TopicsIV, 1, as.numeric(realized_variable), "-"))
sqrt_SE_TopicsIV <- (sweep(sqrt_predictions_TopicsIV, 1, as.numeric(realized_variable), "-"))^2


gc()
sqrt_MCS_SE_TopicsIV <- MCSprocedure(cbind(sqrt_SE_TopicsIV), 
                                     alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)


stopCluster(cl)
rm(cl)
gc()
cl <- makeCluster(5)
sqrt_MCS_AE_TopicsIV <- MCSprocedure(cbind(sqrt_AE_TopicsIV), 
                                     alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)

# Cleanup parallel cluster
stopCluster(cl)
rm(cl)
gc()



########### sqrt_MCS topics ################

cl <- makeCluster(5)

realized_variable <- sqrt_test_rv_var

sqrt_predictions <- all_results$sqrt$predictions
sqrt_predictions_Topics <- sqrt_predictions[, !grepl("HARQ|regressores_HAR|IV|topics_wls|topics_ols", colnames(sqrt_predictions))]


sqrt_AE_Topics <- abs(sweep(sqrt_predictions_Topics, 1, as.numeric(realized_variable), "-"))
sqrt_SE_Topics <- (sweep(sqrt_predictions_Topics, 1, as.numeric(realized_variable), "-"))^2


gc()
sqrt_MCS_SE_Topics <- MCSprocedure(cbind(sqrt_SE_Topics),
                                   alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)


stopCluster(cl)
rm(cl)
gc()
cl <- makeCluster(5)
sqrt_MCS_AE_Topics <- MCSprocedure(cbind(sqrt_AE_Topics),
                                   alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)

# Cleanup parallel cluster
stopCluster(cl)
rm(cl)
gc()



########### sqrt_MCS IV ################

cl <- makeCluster(5)

realized_variable <- sqrt_test_rv_var

sqrt_predictions <- all_results$sqrt$predictions
sqrt_predictions_IV <- sqrt_predictions[, !grepl("HARQ|regressores_HAR|topics", colnames(sqrt_predictions))]


sqrt_AE_IV <- abs(sweep(sqrt_predictions_IV, 1, as.numeric(realized_variable), "-"))
sqrt_SE_IV <- (sweep(sqrt_predictions_IV, 1, as.numeric(realized_variable), "-"))^2


gc()
sqrt_MCS_SE_IV <- MCSprocedure(cbind(sqrt_SE_IV),
                               alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)

stopCluster(cl)
rm(cl)
gc()
cl <- makeCluster(5)
sqrt_MCS_AE_IV <- MCSprocedure(cbind(sqrt_AE_IV),
                               alpha = 0.1, B = 10000, cl = cl, statistic = "TR", k=k, verbose=TRUE)

# Cleanup parallel cluster
stopCluster(cl)
rm(cl)
gc()


############### salvar novos arquivos\############

end_time <-Sys.time()
print(end_time - start_time)

objs <- setdiff(ls(envir = globalenv()), c("base_dir", "forecast_horizon", "tickers",
                                           "tkr", "h","k",
                                           # "test_rv_var","ln_test_rv_var","sqrt_test_rv_var",
                                           "start_time","end_time","cl"))


out_file <- file.path(base_dir, tkr, "RData_TR", paste0("MCS_", tkr, "_", h, ".RData"))
# Save objects from loop_env
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
  