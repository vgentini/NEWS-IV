library(xtable)
library(readxl)
library(xts)
library(data.table)
library(glmnet)
library(zoo)
library(dplyr)
library(tidyr)
library(ggplot2)
library(knitr)
library(kableExtra)
library(stringr)

# Vector of tickers_2_2 to process
tickers_2 <- c("BBDC4", "CSNA3", "GGBR4",
               # IBOV_BS",
               "BOVA11",
               "ITUB4", "PETR4")

forecast_horizon_2 <-c("daily", "lead_5", "lead_10", "lead_21")

# Base path up to each ticker’s folder
base_dir_2 <- "~/FEA_USP/Mestrado/Novos/Tese/R/Correto_sqrt_ln/HAR_and_HIV_BASE_NOVA_HARQ/Impurity"


for (h_2 in forecast_horizon_2) {
  
  # Determine k based on h
  k_2 <- switch(h_2,
              "daily" = 1,
              "lead_5" = 5,
              "lead_10" = 10,
              "lead_21" = 21)


for (tkr in tickers_2) {
  # 1. Locate and load
  in_file <- file.path(base_dir_2, tkr, "RData_TR", paste0("MCS_", tkr, "_", h_2, ".RData"))
  if (!file.exists(in_file)) {
    warning("File not found, skipping: ", in_file); next
  }
  loaded_names <- load(in_file)
  message("Loaded from ", tkr, ": ", paste(loaded_names, collapse = ", "))
  
  # 2. Copy all “MCS…” objects unchanged, suffixing with _<ticker>
  mcs_objs <- grep("MCS", loaded_names, value = TRUE)
  for (nm in mcs_objs) {
    new_nm <- paste0(nm, "_", tkr)
    assign(new_nm, get(nm), envir = .GlobalEnv)
    message("  • Copied MCS object: ", new_nm)
  }
  
  # 3. For each “results…” object with specific suffixes, compute column‐means and save only that
  # pattern <- "^(QLIKE|SE|AE|PAE|HSE)(?=.*results)"
  # res_objs <- grep(pattern, loaded_names, value = TRUE, ignore.case = TRUE)
  
  
  hits <- grep("^(QLIKE|SE|AE|PAE|HSE)",
               loaded_names, value = TRUE, ignore.case = TRUE, perl = TRUE)
  
  # 2) Remove anything that has 'topics' in it (case-insensitive)
  res_objs <- hits[!grepl("results", hits, ignore.case = TRUE)]
  
  
  for (nm in res_objs) {
    obj <- get(nm, envir = .GlobalEnv)
    if (is.data.frame(obj) || is.matrix(obj)) {
      cm <- colMeans(obj, na.rm = TRUE)
      
      # Identify HAR_ols reference column (excluding SHAR)
      tmp_idx <- grep("HAR_ols", names(cm), ignore.case = TRUE)
      har_idx <- tmp_idx[!grepl("SHAR|regressores", names(cm)[tmp_idx], ignore.case = TRUE)]
      
      if (length(har_idx) != 1) {
        warning("Expected exactly one HAR_ols (no SHAR) in ", nm,
                "; found ", length(har_idx), ". Skipping.")
        next
      }
      
      # Normalize metrics relative to HAR_ols
      cm <- cm / cm[har_idx]
      
      # Create new normalized object
      cm_nm <- paste0(nm, "_", tkr)
      assign(cm_nm, cm, envir = .GlobalEnv)
      message("  • Created normalized col‐mean object: ", cm_nm)
    } else {
      warning("Not a data.frame/matrix, skipping colMeans for: ", nm)
    }
  }
  
  # Clean up names by removing "results_" prefix
  objs_to_rename <- grep(paste0("results.*_", tkr, "$"), ls(.GlobalEnv), value = TRUE)
  
  
  # for (old_name in objs_to_rename) {
  #   new_name <- sub("results_", "", old_name)
  #   new_name <- sub("BOVA11", "BOVA11", new_name)
  #   if (exists(new_name, envir = .GlobalEnv)) {
  #     warning("Object ", new_name, " already exists. Skipping rename.")
  #     next
  #   }
  #   assign(new_name, get(old_name, envir = .GlobalEnv))
  #   rm(list = old_name, envir = .GlobalEnv)
  #   message("Renamed: ", old_name, " → ", new_name)
  # }
  
  # List all objects in the global environment
  all_objects <- ls(envir = .GlobalEnv)
  
  # Identify objects containing both "MCS" and "BOVA11"
  # target_objects <- all_objects[
  #   grepl("MCS", all_objects) & 
  #     grepl("BOVA11", all_objects)
  # ]
  
  # Rename each target object by replacing "BOVA11" with "BOVA11"
  # for (obj_name in target_objects) {
  #   new_name <- gsub("BOVA11", "BOVA11", obj_name)
  #   assign(new_name, get(obj_name, envir = .GlobalEnv), envir = .GlobalEnv)
  #   rm(list = obj_name, envir = .GlobalEnv)
  # }
  # 
  # 4. Clean up originals
  rm(list = loaded_names, envir = .GlobalEnv)
}


# --- PART 2: Create metric dataframes (FIXED) ---
tickers<- c("BBDC4", "CSNA3", "GGBR4","BOVA11","ITUB4", "PETR4")
all_objs <- ls()
metrics <- c("QLIKE_regressores_HAR","QLIKE_TopicsIV","QLIKE_Topics","QLIKE_IV","QLIKE",
             "SE_regressores_HAR","SE_TopicsIV","SE_Topics","SE_IV","SE",
             "AE_regressores_HAR","AE_TopicsIV","AE_Topics","AE_IV","AE",
             "PAE_regressores_HAR","PAE_TopicsIV","PAE_Topics","PAE_IV","PAE",
             "HSE_regressores_HAR","HSE_TopicsIV","HSE_Topics","HSE_IV","HSE",
             "sqrt_SE_regressores_HAR", "sqrt_SE_TopicsIV", "sqrt_SE_Topics", "sqrt_SE_IV", "sqrt_SE",
             "sqrt_AE_regressores_HAR", "sqrt_AE_TopicsIV", "sqrt_AE_Topics", "sqrt_AE_IV", "sqrt_AE",
             "ln_SE_regressores_HAR",   "ln_SE_TopicsIV",   "ln_SE_Topics",   "ln_SE_IV",   "ln_SE",
             "ln_AE_regressores_HAR",   "ln_AE_TopicsIV",   "ln_AE_Topics",   "ln_AE_IV",   "ln_AE"
)

df_list <- list()

for (obj_name in all_objs) {
  # Check each metric to find a valid prefix match
  found_metric <- NULL
  found_ticker <- NULL
  for (m in metrics) {
    # Check if the object name starts with the metric followed by an underscore
    prefix <- paste0(m, "_")
    if (startsWith(obj_name, prefix)) {
      # Extract the part after the metric prefix
      potential_rest <- sub(prefix, "", obj_name)
      # Split into parts by underscores
      parts <- unlist(strsplit(potential_rest, "_"))
      if (length(parts) >= 1) {
        # The last part is the potential ticker
        possible_ticker <- parts[length(parts)]
        if (possible_ticker %in% tickers) {
          found_metric <- m
          found_ticker <- possible_ticker
          break
        }
      }
    }
  }
  if (is.null(found_metric)) next  # Skip if no valid metric-ticker pair found
  
  metric_part <- found_metric
  tkr_part <- found_ticker
  
  metric_obj <- get(obj_name)
  
  if (isS4(metric_obj)) {
    vals <- tryCatch(as.numeric(metric_obj), error = function(e) numeric(0))
    nms <- tryCatch(names(metric_obj), error = function(e) character(0))
  } else {
    vals <- unname(metric_obj)
    nms <- names(metric_obj)
  }
  
  if (length(vals) > 0 && length(nms) > 0) {
    df_list[[obj_name]] <- data.frame(
      Model = nms,
      Value = vals,
      Metric = metric_part,
      Ticker = tkr_part,
      stringsAsFactors = FALSE
    )
  }
}

# Combine all dataframes into long format
all_data <- bind_rows(df_list)

# Split into individual dataframes by metric
metric_dfs <- split(all_data, all_data$Metric)

# --- PART 3: Create final dataframes with coloring and column ordering ---

# Define desired column order (excluding "Model")
column_order <- c("BOVA11", "BBDC4", "CSNA3", "GGBR4", "ITUB4", "PETR4")

# Loop through each metric dataframe and process
for (metric_name in names(metric_dfs)) {
  # Ensure the metric dataframe exists and has data
  df_metric <- metric_dfs[[metric_name]]
  if (is.null(df_metric) || nrow(df_metric) == 0) next
  
  # Pivot to wide format, handle duplicates by taking the first occurrence
  df_wide <- df_metric %>%
    select(-Metric) %>%
    pivot_wider(
      names_from = Ticker,
      values_from = Value,
      values_fn = list(Value = first)  # Handle duplicates
    )
  
  # Check if all required columns exist; if not, create them with NA
  missing_cols <- setdiff(column_order, names(df_wide))
  for (col in missing_cols) {
    df_wide[[col]] <- NA_real_
  }
  
  # Select and order columns, then format values
  df <- df_wide %>%
    select(Model, all_of(column_order)) %>%
    mutate(across(-Model, ~ {
      ifelse(is.na(.x), 
             "-", 
             sprintf("%.3f", .x))  # Format non-NA to 3 decimals
    }))
  
  # Apply cell coloring using ORIGINAL model names
  for (ticker_col in names(df)[-1]) {
    
    # Build mcs_obj_name based on metric_name prefix
    if (startsWith(metric_name, "sqrt_")) {
      # remove the leading "sqrt_"
      metric_base <- sub("^sqrt_", "", metric_name)
      
      # now build the object name with a single "sqrt_"
      mcs_obj_name <- paste0(
        "sqrt_",
        "MCS_",
        metric_base,
        "_",
        ticker_col
      )
    }
    else if (startsWith(metric_name, "ln_")) {
      mcs_obj_name <- paste0(
        "ln_",
        "MCS_",
        sub("^ln_", "", metric_name),
        "_",
        ticker_col
      )
    } else {
      mcs_obj_name <- paste0(
        "MCS_",
        metric_name,
        "_",
        ticker_col
      )
    }
    if (exists(mcs_obj_name)) {
      mcs_obj <- get(mcs_obj_name)
      mcs_models <- mcs_obj@Info$model.names
      bold_rows <- df$Model %in% mcs_models
      
      df[[ticker_col]] <- ifelse(
        bold_rows,
        paste0("\\cellcolor{gray!25}{", df[[ticker_col]], "}"),
        df[[ticker_col]]
      )
    }
  }
  
  # Rename Model column with corrected substitutions
  df <- df %>%
    mutate(Model = Model %>%
             str_remove_all('^(ln_|sqrt_)') %>%
             str_replace_all('regressores_HAR', 'HARX') %>%
             str_replace_all('topics_IV', 'NEWS IV') %>%
             str_replace_all('(^|_)topics($|_)', '\\1NEWS\\2') %>%
             # str_replace_all('HARQ_1', 'HARQ') %>%
             str_replace_all('_ols$', ' (OLS)') %>%
             str_replace_all('_wls$', ' (WLS)') %>%
             str_replace_all('_rf$', ' (RF)') %>%
             str_replace_all('_BG$',  ' (BG)')
    )
  
  # Assign processed df back to environment
  assign(paste0(metric_name, '_df'), df)
}

generate_custom_table <- function(df, metric) {
  # 1) Identify row indices for each category
  ols_indices <- grep("\\(OLS\\)", df[[1]])
  wls_indices <- grep("\\(WLS\\)", df[[1]])
  rf_indices  <- grep("^(HARX \\(RF\\)|NEWS \\(RF\\)|NEWS IV \\(RF\\))$", df[[1]])
  bg_indices  <- grep("\\(BG\\)", df[[1]])
  
  # Create a vector of all categorized indices to find uncategorized ones
  all_categorized_indices <- c(ols_indices, wls_indices, rf_indices, bg_indices)
  all_row_indices <- 1:nrow(df)
  other_indices <- setdiff(all_row_indices, all_categorized_indices) # Rows not in OLS, WLS, RF, or BG
  
  # 2) Reorder the data frame
  # Start with OLS, then WLS, then RF, then BG, then any other rows
  ordered_df <- df[c(ols_indices, wls_indices, rf_indices, bg_indices, other_indices), ]
  
  # 3) Recalculate the 'last' row indices based on the *newly ordered* data frame
  # These are now relative to the 'ordered_df'
  current_row_count <- 0
  
  last_ols <- if (length(ols_indices) > 0) {
    current_row_count <- current_row_count + length(ols_indices)
    current_row_count
  } else {
    NULL
  }
  
  last_wls <- if (length(wls_indices) > 0) {
    current_row_count <- current_row_count + length(wls_indices)
    current_row_count
  } else {
    NULL
  }
  
  last_rf <- if (length(rf_indices) > 0) {
    current_row_count <- current_row_count + length(rf_indices)
    current_row_count
  } else {
    NULL
  }
  
  last_bg <- if (length(bg_indices) > 0) {
    # Only add linespace if BG rows are not the absolute last rows in the table
    # (i.e., if there are 'other_indices' rows following them)
    if (length(other_indices) > 0) {
      current_row_count <- current_row_count + length(bg_indices)
      current_row_count
    } else {
      # If BG rows are the last, we don't need linespace after them,
      # but we still need to know their position for other logic if any.
      # However, for addlinespace, we set it to NULL.
      NULL # Or handle as per specific needs if BG is truly last.
      # For addlinespace, setting to NULL means no line will be added.
    }
  } else {
    NULL
  }
  
  
  # 4) Build the kable and styling, using the 'ordered_df'
  #    Add a LaTeX label based on metric, and enable threeparttable for footnotes:
  tbl <- kable(
    ordered_df, # Use the reordered data frame
    format    = "latex",
    booktabs  = TRUE,
    linesep   = "",              # disable automatic \addlinespace everywhere
    align     = "lcccccc",
    col.names = c(
      "",
      "\\textbf{BOVA11}",
      "\\textbf{BBDC4}",
      "\\textbf{CSNA3}",
      "\\textbf{GGBR4}",
      "\\textbf{ITUB4}",
      "\\textbf{PETR4}"
    ),
    caption   = paste0(metric, " results for h=",k_2),
    label     = paste0(metric, " ",h_2),
    escape    = FALSE
  ) %>%
    kable_styling(
      latex_options = c("hold_position", "threeparttable"),
      full_width    = FALSE
    ) %>%
    row_spec(0, bold = TRUE)      # header row bold
  
  # 5) Add \addlinespace after Last OLS, WLS, RF, and BG rows if found
  #    These are now based on the new counts in the ordered_df
  if (!is.null(last_ols) && last_ols < nrow(ordered_df)) { # Check if not the last row
    tbl <- tbl %>% row_spec(last_ols, extra_latex_after = "\\addlinespace[2em]")
  }
  if (!is.null(last_wls) && last_wls < nrow(ordered_df)) { # Check if not the last row
    tbl <- tbl %>% row_spec(last_wls, extra_latex_after = "\\addlinespace[2em]")
  }
  if (!is.null(last_rf) && last_rf < nrow(ordered_df)) { # Check if not the last row
    tbl <- tbl %>% row_spec(last_rf, extra_latex_after = "\\addlinespace[2em]")
  }
  # For last_bg, the logic in step 3 already considers if it's followed by 'other' rows.
  # If last_bg was set to a row number, it means there are rows after it (other_indices).
  if (!is.null(last_bg) && last_bg < nrow(ordered_df)) { # Ensure it's not the absolute last row
    tbl <- tbl %>% row_spec(last_bg, extra_latex_after = "\\addlinespace[2em]")
  }
  
  # 6) Add the tablenote via a threeparttable footnote
  tbl <- tbl %>%
    footnote(
      general = "\\textit{Note:} This table shows the loss ratios of various models compared to the HAR estimated through OLS. Cells highlighted indicate models that are part of the MCS at a 90\\% confidence level. Cells with missing values pertain to models that yielded forecasts with negative predictions",
      general_title = "",          # suppress “Note:” prefix if you like
      escape = FALSE,
      threeparttable = TRUE
    )
  
  return(tbl)
}



## R Script: Generate LaTeX Tables for All Metrics


setwd(file.path(base_dir_2, paste0("Tabelas_MCS_k_", k_2)))

for (metric in metrics) {
  # 1) build the data‐frame name
  df_name <- paste0(metric, "_df")
  if (!exists(df_name)) {
    warning(sprintf("Data frame '%s' not found. Skipping metric: %s", df_name, metric))
    next
  }
  df <- get(df_name)
  
  # 2) apply all your renaming rules in sequence
  metric_label <- metric %>%
    # HAR and topics first
    str_replace_all("regressores_HAR", "HARX") %>%
    str_replace_all("TopicsIV",      "NEWS IV") %>%
    str_replace_all("(^|_)Topics($|_)", "\\1NEWS\\2") %>%
    # MSE group: longest first
    str_replace_all("^ln_SE",    "MSE-LOG ") %>%
    str_replace_all("^sqrt_SE",  "MSE-SD ") %>%
    str_replace_all("^HSE",      "HMSE ") %>%
    str_replace_all("^SE",       "MSE ") %>%
    # MAE group: longest first
    str_replace_all("^ln_AE",     "MAE-LOG ") %>%
    str_replace_all("^sqrt_AE",   "MAE-SD ") %>%
    str_replace_all("^PAE",       "MAPE ") %>%
    str_replace_all("^AE",        "MAE ") %>%
    str_replace_all("^QLIKE",        "QLIKE ") %>%
    # finally, remove all hyphens
    str_replace_all("_", "")
  
  # 3) generate LaTeX using the cleaned‐up label
  latex_code <- generate_custom_table(df, metric_label)
  assign(paste0(metric, "_latex"), latex_code, envir = .GlobalEnv)
  
  # 4) write out the file, using the transformed name if you like
  out_file <- paste0(metric_label, ".tex")
  cat(latex_code, file = out_file)
  message(sprintf(
    "Written LaTeX table for '%s' (as '%s') to %s",
    metric, metric_label, out_file
  ))
}






# --- PART 2: Create metric dataframes (FIXED) ---
tickers<- c("BBDC4", "CSNA3", "GGBR4","BOVA11","ITUB4", "PETR4")
all_objs <- ls()
metrics <- c("sqrt_SE", "sqrt_AE")

df_list <- list()

for (obj_name in all_objs) {
  # Check each metric to find a valid prefix match
  found_metric <- NULL
  found_ticker <- NULL
  for (m in metrics) {
    # Check if the object name starts with the metric followed by an underscore
    prefix <- paste0(m, "_")
    if (startsWith(obj_name, prefix)) {
      # Extract the part after the metric prefix
      potential_rest <- sub(prefix, "", obj_name)
      # Split into parts by underscores
      parts <- unlist(strsplit(potential_rest, "_"))
      if (length(parts) >= 1) {
        # The last part is the potential ticker
        possible_ticker <- parts[length(parts)]
        if (possible_ticker %in% tickers) {
          found_metric <- m
          found_ticker <- possible_ticker
          break
        }
      }
    }
  }
  if (is.null(found_metric)) next  # Skip if no valid metric-ticker pair found
  
  metric_part <- found_metric
  tkr_part <- found_ticker
  
  metric_obj <- get(obj_name)
  
  if (isS4(metric_obj)) {
    vals <- tryCatch(as.numeric(metric_obj), error = function(e) numeric(0))
    nms <- tryCatch(names(metric_obj), error = function(e) character(0))
  } else {
    vals <- unname(metric_obj)
    nms <- names(metric_obj)
  }
  
  if (length(vals) > 0 && length(nms) > 0) {
    df_list[[obj_name]] <- data.frame(
      Model = nms,
      Value = vals,
      Metric = metric_part,
      Ticker = tkr_part,
      stringsAsFactors = FALSE
    )
  }
}

# Combine all dataframes into long format
all_data <- bind_rows(df_list)

# Split into individual dataframes by metric
metric_dfs <- split(all_data, all_data$Metric)

# --- PART 3: Create final dataframes with coloring and column ordering ---

# Define desired column order (excluding "Model")
column_order <- c("BOVA11", "BBDC4", "CSNA3", "GGBR4", "ITUB4", "PETR4")

# Loop through each metric dataframe and process
for (metric_name in names(metric_dfs)) {
  # Ensure the metric dataframe exists and has data
  df_metric <- metric_dfs[[metric_name]]
  if (is.null(df_metric) || nrow(df_metric) == 0) next
  
  # Pivot to wide format, handle duplicates by taking the first occurrence
  df_wide <- df_metric %>%
    select(-Metric) %>%
    pivot_wider(
      names_from = Ticker,
      values_from = Value,
      values_fn = list(Value = first)  # Handle duplicates
    )
  
  # Check if all required columns exist; if not, create them with NA
  missing_cols <- setdiff(column_order, names(df_wide))
  for (col in missing_cols) {
    df_wide[[col]] <- NA_real_
  }
  
  # Select and order columns, then format values
  df <- df_wide %>%
    select(Model, all_of(column_order)) %>%
    mutate(across(-Model, ~ {
      ifelse(is.na(.x), 
             "-", 
             sprintf("%.3f", .x))  # Format non-NA to 3 decimals
    }))
  
  # Apply cell coloring using ORIGINAL model names
  for (ticker_col in names(df)[-1]) {
    
    # Build mcs_obj_name based on metric_name prefix
    if (startsWith(metric_name, "sqrt_")) {
      # remove the leading "sqrt_"
      metric_base <- sub("^sqrt_", "", metric_name)
      
      # now build the object name with a single "sqrt_"
      mcs_obj_name <- paste0(
        "sqrt_",
        "MCS_",
        metric_base,
        "_",
        ticker_col
      )
    }
    else if (startsWith(metric_name, "ln_")) {
      mcs_obj_name <- paste0(
        "ln_",
        "MCS_",
        sub("^ln_", "", metric_name),
        "_",
        ticker_col
      )
    } else {
      mcs_obj_name <- paste0(
        "MCS_",
        metric_name,
        "_",
        ticker_col
      )
    }
    if (exists(mcs_obj_name)) {
      mcs_obj <- get(mcs_obj_name)
      mcs_models <- mcs_obj@Info$model.names
      bold_rows <- df$Model %in% mcs_models
      
      df[[ticker_col]] <- ifelse(
        bold_rows,
        paste0("\\cellcolor{gray!25}{", df[[ticker_col]], "}"),
        df[[ticker_col]]
      )
    }
  }
  
  # Rename Model column with corrected substitutions
  df <- df %>%
    mutate(Model = Model %>%
             str_remove_all('^(ln_|sqrt_)') %>%
             str_replace_all('regressores_HAR', 'HARX') %>%
             str_replace_all('topics_IV', 'NEWS IV') %>%
             str_replace_all('(^|_)topics($|_)', '\\1NEWS\\2') %>%
             # str_replace_all('HARQ_1', 'HARQ') %>%
             str_replace_all('_ols$', ' (OLS)') %>%
             str_replace_all('_wls$', ' (WLS)') %>%
             str_replace_all('_BG$',  ' (RF)')
    )
  
  # Assign processed df back to environment
  assign(paste0(metric_name, '_df'), df)
}



## R Script: Generate LaTeX Tables for All Metrics

setwd(file.path(base_dir_2, paste0("Tabelas_MCS_k_", k_2)))

for (metric in metrics) {
  # 1) build the data‐frame name
  df_name <- paste0(metric, "_df")
  if (!exists(df_name)) {
    warning(sprintf("Data frame '%s' not found. Skipping metric: %s", df_name, metric))
    next
  }
  df <- get(df_name)
  
  # 2) apply all your renaming rules in sequence
  metric_label <- metric %>%
    # HAR and topics first
    str_replace_all("regressores_HAR", "HARX") %>%
    str_replace_all("Topics_IV",      "NEWS IV") %>%
    str_replace_all("(^|_)Topics($|_)", "\\1NEWS\\2") %>%
    # MSE group: longest first
    str_replace_all("^ln_SE",    "MSE-LOG ") %>%
    str_replace_all("^sqrt_SE",  "MSE-SD ") %>%
    str_replace_all("^HSE",      "HMSE ") %>%
    str_replace_all("^SE",       "MSE ") %>%
    # MAE group: longest first
    str_replace_all("^ln_AE",     "MAE-LOG ") %>%
    str_replace_all("^sqrt_AE",   "MAE-SD ") %>%
    str_replace_all("^PAE",       "MAPE ") %>%
    str_replace_all("^AE",        "MAE ") %>%
    # finally, remove all hyphens
    str_replace_all("_", "")
  
  # 3) generate LaTeX using the cleaned‐up label
  latex_code <- generate_custom_table(df, metric_label)
  assign(paste0(metric, "_latex"), latex_code, envir = .GlobalEnv)
  
  # 4) write out the file, using the transformed name if you like
  out_file <- paste0(metric_label, ".tex")
  cat(latex_code, file = out_file)
  message(sprintf(
    "Written LaTeX table for '%s' (as '%s') to %s",
    metric, metric_label, out_file
  ))
}

rm(list=setdiff(ls(envir = .GlobalEnv, all.names = TRUE), c("tickers_2","tkr",
                                                    "forecast_horizon_2","k_2","h_2",
                                                    "base_dir_2")))
                                                    
}





