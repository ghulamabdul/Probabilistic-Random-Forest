## Objective Identification: COSMO-REA6
# Prediction via station-based random forests

#### Settings ####
# Load package
library(lubridate)
library(dplyr)
library(ranger)

# Path of preprocessed COSMO-data
data_cosmo_path <- "/.../"

# Path of R-data
data_path <- "/.../"

# Path of random forest files
data_rf_path <- "/.../"

#### Initiation ####
# Dates corresponding to storms
date_vec <- c("Niklas" = "20150331", 
              "Susanna"= "20160209",
              "Egon" = "20170112", 
              "Thomas" = "20170223", 
              "Xavier" = "20171005", 
              "Herwart" = "20171029", 
              "Burglind" = "20180103", 
              "Friederike" = "20180118", 
              "Fabienne"= "20180923", 
              "Bennet"= "20190304",
              "Eberhard" = "20190310")

# Labels considered
lab_vec <- c(0, 1, 2, 3, 5)

# Data frame for overall evaluation
df_total <- data.frame(lab = vector(length = 0),
                       storm = vector(length = 0))

# Probability for each label
for(i in lab_vec){ df_total[[paste0("p", i)]] <- vector(length = 0) }

#### Read out ####
# For-Loop over storms
for(temp_storm in names(date_vec)){
  #### Get data ####
  # Load gridded data
  file_name <- paste0(data_cosmo_path, "cosmo_data_", temp_storm, ".RData")
  
  # Load data
  load(file = file_name)
  
  # Load RF data
  file_name <- paste0(data_rf_path, temp_storm, ".RData") 
  
  # Load data
  load(file = file_name)
  rm(pred)
  
  #### Initiation ####
  # Rows for data frame
  i_rows <- nrow(df_total) + 1:nrow(df_storms)
  
  # Write in data frame
  df_total[i_rows, "storm"] <- temp_storm
  
  #### RF Prediction ####
  # Predict
  pred <- predict(object = est,
                  data = df_storms)
  
  #### Save data ####
  # Name of file
  file_save <- paste0(data_cosmo_path, "cosmo_preds_", temp_storm, ".RData")
  
  # Save data
  save(file = file_save,
       list = c("pred"))
  
  #### Read out ####
  # Probabilities and labels
  df_total[i_rows, paste0("p", lab_vec)] <- pred$predictions
  df_total[i_rows, "lab"] <- df_storms[["lab"]]
}

#### Save data ####
# Name of file
file_save <- paste0(data_path, "cosmo_preds.RData")

# Save data
save(file = file_save,
     list = c("df_total"))

#### Climatology: Calculation ####
# Initiate climatological forecasts
df_clim <- data.frame(storm = character(),
                      stringsAsFactors = FALSE)

# Probability for each label
for(i in lab_vec){ df_clim[[paste0("p", i)]] <- numeric() }

# For-Loop over storms
for(temp_storm in names(date_vec)){
  # Get index
  i <- nrow(df_clim) + 1
  
  # Read out
  df_clim[i, "storm"] <- temp_storm
  
  # Calculate frequencies
  df_clim[i, paste0("p", lab_vec)] <- sapply(lab_vec, function(i){
    mean(subset(df_total, storm != temp_storm)[["lab"]] == i) })
}

# Mean multivariate Brier Score
brier_score <- function(df_p, y_vec){
  ### Input
  #...df_p....Data Frame where each sample is a probability vector (nxk df)
  #...y_vec...Label of feature (n vector)
  
  # Calculate mean over samples
  res <- mean(sapply(1:nrow(df_p), function(i_row){
    # Brier Score of i_row-th sample
    sum((df_p[i_row,] - as.numeric(y_vec[i_row] == lab_vec))^2)
  }))
  
  # Return
  return(res)
}

# Variables for data frame
col_vec <- c("storm", "n", "method", "value")

# Data frame for scores
df_scores <- data.frame(matrix(nrow = 0,
                               ncol = length(col_vec)))
colnames(df_scores) <- col_vec

# For-Loop over scores, methods and labels
for(temp_storm in names(date_vec)){ for(temp_meth in c("p", "p_clim")){
  # Get indices of storm
  i_test <- which(df_total[["storm"]] == temp_storm)
  
  # Get current index
  i_row <- nrow(df_scores) + 1
  
  # Write in data frame
  df_scores[i_row, "storm"] <- temp_storm
  df_scores[i_row, "n"] <- length(i_test)
  df_scores[i_row, "method"] <- temp_meth
  
  # Get forecasts
  if(temp_meth == "p"){ df_p <- df_total[i_test, paste0("p", lab_vec)] }
  else{ 
    df_p <- df_clim[which(df_clim[["storm"]] == temp_storm), paste0("p", lab_vec)]
    df_p <- rbind(df_p[rep(1, length(i_test)),])
  }
  
  # Calculate scores
  df_scores[i_row, "value"] <-  brier_score(df_p = df_p,
                                            y_vec = df_total[i_test, "lab"])
}}

# File name
file_name <- paste0(data_path, "cosmo_clim.RData") 

# Save scores
save(file = file_name,
     list = c("df_scores"))

#### Climatology: Evaluation ####
# Overall skill
df_sub_p <- subset(df_scores, method == "p")
df_sub_clim <- subset(df_scores, method == "p_clim")
bss <- 1 - sum(df_sub_p[["n"]]*df_sub_p[["value"]])/sum(df_sub_clim[["n"]]*df_sub_clim[["value"]])

# Function for apply
fn_apply <- function(x){
  1 - df_scores[which((df_scores$storm == x) & (df_scores$method == "p")), "value"]/
    df_scores[which((df_scores$storm == x) & (df_scores$method == "p_clim")), "value"] }

# Calculate skill
bss_storms <- sapply(names(date_vec), fn_apply)

# Output
print(round(100*bss, 2))
print(round(100*sort(bss_storms, decreasing = TRUE), 2))