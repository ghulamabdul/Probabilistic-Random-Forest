## Objective Identification: COSMO-REA6
# Data preprocessing

#### Settings ####
# Load package
library(lubridate)
library(dplyr)
library(ncdf4)

# Path of ncdf files
data_in_path <- "/.../"

# Path for R-data
data_out_path <- "/.../"

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
  
  # No index-variable
  fc_names <- fc_names[fc_names != "index"]
  
  # Columns of resulting data frame
  col_vec <- unique(c("date", "time", "lon", "lat", fc_names))
  
  # Get longitude vector
  lon_vec <- as.vector(ncvar_get(nc = fc_nc, 
                                 varid = "lon"))
  
  # Get latitude vector
  lat_vec <- as.vector(ncvar_get(nc = fc_nc, 
                                 varid = "lat"))
  
  # Get label vector
  label_vec <- as.vector(ncvar_get(nc = fc_nc, 
                                   varid = "lab"))
  
  # Get vgl vector
  vgl_vec <- as.vector(ncvar_get(nc = fc_nc,
                                 varid = "vgl"))
  
  # Cut all v_tilde below 0.8 and NAs
  i_cut <- (vgl_vec < 0.8) | is.na(label_vec) 
  i_cut[is.na(vgl_vec)] <- TRUE
  
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
  
  # Latitude and longitude
  res[,"lon"] <- lon_vec
  res[,"lat"] <- lat_vec

  # For-Loop over variables
  for(temp_var in fc_names){
    # Make matrix to vector (by columns)
    res[,temp_var] <- as.vector(ncvar_get(nc = fc_nc, 
                                          varid = temp_var))[!i_cut] }
  
  # Close nc
  nc_close(fc_nc)
  
  # Output
  return(res)
}

#### For-Loop over storms ####
for(temp_storm in names(date_vec)){
  # Get files of storm
  storm_files <- file_names[grepl(date_vec[temp_storm], file_names, fixed = TRUE)]
  
  # Load files and write in data frame
  df_storms <- do.call(rbind, lapply(storm_files, fn_nc, nc_path = data_in_path))
  
  # Name of file
  file_name <- paste0(data_out_path, "cosmo_data_", temp_storm, ".RData")
  
  # Save
  save(file = file_name,
       list = c("df_storms"))
}