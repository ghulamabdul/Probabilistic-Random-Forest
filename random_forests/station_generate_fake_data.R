## Objective Identification: Station data
# Generation of fake data based on original data

#### Settings ####
# Path of data
data_path <- "/.../"

#### Load data ####
# Load original data
load(paste0(data_path, "station_data.RData")) # Not provided

#### Generation ####
# For-Loop over variables to change
for(temp_var in c("vgl", "dpdt", "dthdt", "dwddt", "prec", 
                  "wind", "wd", "mslp", "th", "lab")){
  # Get rows with missing values
  i_na <- is.na(df_storms[[temp_var]])
  
  # Round to d decimals (depending on variable)
  d <- 1 - is.element(temp_var, c("wd", "dwddt", "mslp", "lab")) + 
    3*is.element(temp_var, c("th", "dthdt"))
  
  # Values of variable (rounded)
  v_vec <- round(df_storms[[temp_var]], d)
  
  # Frequency of variable values
  n_vec <- table(v_vec)/length(v_vec)
  
  # Sort variable values
  v_vec <- sort(unique(v_vec))
  
  # Sample values for fake data
  df_storms[[temp_var]] <- sample(x = v_vec,
                                  size = nrow(df_storms),
                                  replace = TRUE,
                                  prob = n_vec)
  
  # Integer valued
  if(is.element(temp_var, c("wd", "dwddt", "mslp", "lab"))){
    df_storms[[temp_var]] <- as.integer(df_storms[[temp_var]])
  }
  
  # Recreate missing values
  df_storms[[temp_var]][i_na] <- NaN
}

#### Save ####
# Save fake data
save(file = paste0(data_path, "station_fake_data.RData"), 
     list = c("df_storms"))