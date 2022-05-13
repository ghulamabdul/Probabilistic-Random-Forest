## Objective Identification: Station data
# Calculate objects for partial dependence plots

#### Settings ####
# Load package
library(lubridate)
library(dplyr)
library(ranger)
library(rpart)
library(ggplot2)
library(gridExtra)
library(pdp)

# Path of RF objects
data_rf_path <- "/.../"

# Path of R-Data
data_path <- "/.../"

# Path for PDP objects
data_pdp_path <- "/.../"

#### Load data ####
# Name of file
file_name <- paste0(data_path, "station_preds.RData")

# Load data frame
load(file = file_name)

#### Initiation ####
# Number of evaluations for PDP
n_eval_pdp <- 51

# Vector of labels
lab_vec <- sort(unique(df_storms[["lab"]]))

# Names of the winter storms
name_vec <- c("Niklas", "Susanna", 
              "Egon", "Thomas", "Xavier", "Herwart",
              "Burglind", "Friederike", "Fabienne",
              "Bennet", "Eberhard",
              "Sabine") 

# Dates of storms
date_vec <- unique(df_storms$date)

# Read out predictors
pred_vars <- names(df_storms)[!is.element(names(df_storms),
                                          c("date", "time", "lab", "lon", "lat", "wind", "alt"))]

# Check which columns include missing values
na_cols <- vector(length = ncol(df_storms))
for(i in 1:ncol(df_storms)){ na_cols[i] <- any(is.na(df_storms[,i])) }

# Limits for variables
xlim_ls <- list(
  "dpdt" = c(-5, 5),
  "dthdt" = c(-0.01, 0.01),
  "prec" = c(0, 8),
  "dwddt" = c(-75, 150),
  "mslp" = c(970, 1025),
  "vgl" = c(0.8, 3),
  "wd" = c(0, 359),
  "th" = c(0.94, 1.08)
)

# List for pdp-objects
pdp_ls <- grid_ls <- list()

# For-Loop over predictors
for(temp_var in pred_vars){ 
  # List for PDP plots
  pdp_ls[[temp_var]] <- list() 
  
  # Equidistant binning within given limits
  grid_ls[[temp_var]] <- as.data.frame(seq(from = min(xlim_ls[[temp_var]],
                                                      na.rm = TRUE),
                                           to = max(xlim_ls[[temp_var]],
                                                    na.rm = TRUE),
                                           length.out = n_eval_pdp))
  
  # Name binning accordingly
  colnames(grid_ls[[temp_var]]) <- temp_var
}

#### Partial dependence plots ####
# For-Loop over storms
for(temp_date in date_vec){
  # Name of storm
  temp_name <- name_vec[which(date_vec == temp_date)]
  
  # Output
  print(paste0("--- Storm: ", temp_name))
  
  # List for storm-specific pdp-data
  pdp_storm_ls <- list()
  
  # For-Loop over predictors
  for(temp_var in pred_vars){ pdp_storm_ls[[temp_var]] <- list() }
  
  # Test set
  i_test <- which(df_storms$date == temp_date)
  
  # Get training and test data
  df_train <- slice(df_storms, -i_test)
  df_test <- slice(df_storms, i_test)
  
  # Replace missing values with mean
  for(i in 1:ncol(df_storms)){ if(na_cols[i]){
    # Calculate mean on training set
    mu <- mean(df_train[,i], na.rm = TRUE)
    
    # Replace missing values
    df_train[is.na(df_train[,i]), i] <- mu 
    df_test[is.na(df_test[,i]), i] <- mu 
  }}
  
  # Name of file
  file_name <- paste0(data_rf_path, temp_name, ".RData")

  # Load data
  load(file = file_name)
  
  # For-Loop over labels
  for(i_lab in 1:length(lab_vec)){
    # Get label
    temp_lab <- lab_vec[i_lab]
    
    # Output
    print(paste0("--- Label: ", temp_lab))
    
    # For-Loop over predictor variables
    for(temp_var in pred_vars){
      # Generate PDP
      pdp_storm_ls[[temp_var]][[paste0(temp_lab)]] <- 
        pdp_ls[[temp_var]][[paste0(temp_name, temp_lab)]] <- 
        partial(object = est, 
                pred.data = df_test, 
                pred.var = temp_var, 
                which.class = i_lab,
                pred.grid = grid_ls[[temp_var]],
                type = "classification",
                prob = TRUE,
                parallel = FALSE
                )
    }
  }
  
  # Name of file
  file_save <- paste0(pdp_path, temp_name, ".RData")

  # Save storm-specific data
  save(file = file_save,
       list = c("pdp_storm_ls"))
}

# For-Loop over labels
for(i_lab in 1:length(lab_vec)){
  # Get label
  temp_lab <- lab_vec[i_lab]
  
  # For-Loop over predictor variables
  for(temp_var in pred_vars){ 
    # Make matrix of yhats
    temp_mtx <- matrix(nrow = length(date_vec),
                       ncol = n_eval_pdp)
    
    # For-Loop over storms
    for(temp_date in date_vec){
      # Get name of storm
      temp_name <- name_vec[which(date_vec == temp_date)]
      
      # Read out predictions
      temp_mtx[which(temp_date == date_vec),] <- 
        pdp_ls[[temp_var]][[paste0(temp_name, temp_lab)]][["yhat"]]*sum(df_storms[["date"]] == temp_date)
    }
    
    # Copy pdp object
    pdp_ls[[temp_var]][[paste0(temp_lab)]] <- pdp_ls[[temp_var]][[paste0(name_vec[1], temp_lab)]]
    
    # Overwrite predictions
    pdp_ls[[temp_var]][[paste0(temp_lab)]][["yhat"]] <- colSums(temp_mtx)/nrow(df_storms)
  }
}

#### Save data ####
# Name of file
file_save <- paste0(data_path, "station_pdp.RData")

# Save data
save(file = file_save,
     list = c("pdp_ls"))