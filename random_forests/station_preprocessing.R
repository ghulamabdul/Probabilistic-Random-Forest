## Objective Identification: Station data
# Data preprocessing

#### Settings ####
# Load package
library(lubridate)
library(dplyr)
library(ncdf4)

# Path of ncdf files (not provided)
data_in_path <- "/.../"

# Path for R-Data
data_out_path <- "/.../"

#### Initiation ####
# Dates corresponding to storms
date_vec <- c("20150331", "20160209",
              "20170112", "20170223", "20171005", "20171029",
              "20180103", "20180118", "20180923", 
              "20190304", "20190310",
              "20200202")

# Get all file names
file_names <- list.files(data_in_path)

#### Function ####
# Function to read out one file
fn_nc <- function(nc_name, nc_path = data_in_path){
  ### Input
  #...nc_name...Name of ncdf file (string)
  #...nc_path...Path to ncdf file (string)
  ### Output
  #...res...Data frame with forecast information
  ###

  # Open ncdf file
  fc_nc <- nc_open(paste0(nc_path, nc_name))
  
  # Get variable names
  fc_names <- names(fc_nc$var)
  fc_dim <- names(fc_nc$dim)
  
  # Columns of resulting data frame
  col_vec <- unique(c("date", "time", fc_names, "alt"))
  
  # Get label vector
  label_vec <- as.vector(ncvar_get(nc = fc_nc, 
                                   varid = "lab"))
  
  # Get altitude vector directly
  alt_vec <- as.vector(ncvar_get(nc = fc_nc,
                                 varid = "alt"))

  # Get v_tilde vector
  vgl_vec <- as.vector(ncvar_get(nc = fc_nc,
                                 varid = "vgl"))
  
  # Cut all v_tilde below 0.8 and NAs
  i_cut <- (vgl_vec < 0.8) | (alt_vec > 800) | is.na(label_vec) 
  i_cut[is.na(vgl_vec) | is.na(alt_vec)] <- TRUE
  
  # If everything is cut, return NULL
  if(sum(i_cut) == length(label_vec)){ return(NULL) }
  
  # Initiate resulting data frame
  res <- as.data.frame(matrix(nrow = length(label_vec) - sum(i_cut),
                              ncol = length(col_vec)))
  
  # Set names of variables
  names(res) <- col_vec
  
  # Read out date, time and member
  res[,"date"] <- substr(x = nc_name, 
                         start = 1, 
                         stop = 8)
  res[,"time"] <- as.numeric(substr(x = nc_name, 
                                    start = 10, 
                                    stop = 11))
  res[,"alt"] <- alt_vec[!i_cut]
  
  # For-Loop over variables
  for(temp_var in fc_names[fc_names != "alt"]){
    # Make matrix to vector (by columns)
    res[,temp_var] <- as.vector(ncvar_get(nc = fc_nc, 
                                          varid = temp_var))[!i_cut] }
  
  # Close nc
  nc_close(fc_nc)
  
  # Output
  return(res)
}

#### Read out and save ####
# Load files and write in data frame
df_storms <- do.call(rbind, lapply(file_names, fn_nc, nc_path = data_in_path))

# Name of file
file_name <- paste0(data_out_path, "station_data.RData") # Not provided

# Save
save(file = file_name,
     list = c("df_storms"))