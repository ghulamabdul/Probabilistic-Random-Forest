## Objective Identification: Station data
# Training of a random forest on all storm cases

#### Settings ####
# Load package
library(lubridate)
library(dplyr)
library(ranger)
library(rpart)

# Path of RF objects
data_rf_path <- "/.../"

# Path of R data
data_path <- "/.../"

#### Load data ####
# Name of file
# file_name <- paste0(data_path, "station_data.RData") # Not provided
file_name <- paste0(data_path, "station_fake_data.RData")

# Load data frame
load(file = file_name)

#### Initiation ####
# Vector of labels
lab_vec <- sort(unique(df_storms$lab))

# Read out predictors
pred_vars <- names(df_storms)[!is.element(names(df_storms),
                                   c("date", "time", "lab", "lon", "lat", "wind", "alt"))]

# Formula
pred_fm <- as.formula(paste0("lab ~ ", paste0(pred_vars,
                                                collapse = " + ")))

#### Estimation ####
# Check which columns include missing values
na_cols <- vector(length = ncol(df_storms))
for(i in 1:ncol(df_storms)){ na_cols[i] <- any(is.na(df_storms[,i])) }

# Entire set for training
df_train <- df_storms

# Replace missing values with mean
for(i in 1:ncol(df_storms)){ if(na_cols[i]){
  # Calculate mean on training set
  mu <- mean(df_train[,i], na.rm = TRUE)
  
  # Replace missing values
  df_train[is.na(df_train[,i]), i] <- mu 
}}
  
# Estimate
est <- ranger(formula = pred_fm,
              data = df_train[order(df_train$lab),],
              importance = "permutation",
              probability = TRUE,
              num.trees = 1000)

# Name of file
file_save <- paste0(data_rf_path, "total.RData")

# Save data
save(file = file_save,
     list = c("est"))