## Objective Identification: COSMO-REA6
# Transform data in csv-files

#### Settings ####
# Load package
library(lubridate)
library(dplyr)
library(sp)
library(maps)
library(maptools)

# Path of COSMO-data
data_path <- "/.../"

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

# Vector of countries to consider
regions_vec <- c("austria", "belgium", "czech republic", "denmark", "france", 
                 "germany", "luxembourg",  "netherlands", "poland", "switzerland", "uk")

#### Functions ####
# Making a function that assigns a state, country to a given coordinate
lonlat_to_state_sp <- function(pointsDF) {
  # Prepare SpatialPolygons object with one SpatialPolygon
  # per state (plus DC, minus HI & AK)
  states <- map('world',
                regions = regions_vec, 
                fill = TRUE, col = "transparent", 
                plot = FALSE)
  IDs <- sapply(strsplit(states$names, ":"), function(x) x[1])
  states_sp <- map2SpatialPolygons(states, 
                                   IDs = IDs,
                                   proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  # Convert pointsDF to a SpatialPoints object 
  pointsSP <- SpatialPoints(pointsDF, 
                            proj4string = CRS("+proj=longlat +datum=WGS84"))
  
  # Use 'over' to get _indices_ of the Polygons object containing each point 
  indices <- over(pointsSP, states_sp)
  
  # Return the state names of the Polygons object containing each point
  stateNames <- sapply(states_sp@polygons, function(x) x@ID)
  stateNames[indices]
}

# Cut forecasts to relevant regions
cut_data <- function(my.data){
  ### Input
  #...my.data...Data-frame including longitude and latitude variables
  ### Output
  #...res...Data including only relevant lon-lat pairs
  ###
  
  # Cut outside region
  my.data <- my.data[(my.data[["lon"]] <= 20) & (my.data[["lon"]] >= -5) & 
                       (my.data[["lat"]] <= 60) & (my.data[["lat"]] >= 43) &
                       (my.data[["alt"]] <= 800),]
  
  # Only relevant countries
  check.points <- lonlat_to_state_sp(my.data[,c("lon", "lat")])
  res <- my.data[!is.na(check.points),]
  
  # Output
  return(res)
}

#### Read out ####
# For-Loop over files
for(temp_storm in names(date_vec)){
  # Load gridded data
  file_name <- paste0(data_path, "cosmo_data_", temp_storm, ".RData")
  
  # Save data
  load(file = file_name)
  
  # Cut irrelevant data
  df_storms <- df_storms[,c("date", "time", "lon", "lat", "alt", "lab")]
  
  # Load prediction data
  file_name <- paste0(data_path, "cosmo_preds_", temp_storm, ".RData")
  
  # Load data
  load(file = file_name)
  
  # Add predictions to data
  df_storms[,paste0("p", c(0:3, 5))] <- pred$predictions
  
  # Cut storms to relevant regions
  df_storms <- cut_data(my.data = df_storms)
  
  # Name of csv-file
  file_name <- paste0(data_path, "cosmo_preds_", temp_storm, ".csv")
  
  # Make csv
  write.csv(x = df_storms[,colnames(df_storms) != "alt"],
            file = file_name,
            row.names = FALSE)
}