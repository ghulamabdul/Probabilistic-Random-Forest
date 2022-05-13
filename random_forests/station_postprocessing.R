## Objective Identification: Station data
# Transform data in csv-files

#### Settings ####
# Path of R data
data_path <- "/.../"

#### Initiation ####
# Name of file
file_name <- paste0(data_path, "station_preds.RData")

# Load data frame
load(file = file_name)

# Get label vector
lab_vec <- sort(unique(df_total$lab))

# Name of csv-file
file_name <- paste0(data_path, "station_preds.csv")

# Make csv
write.csv(x = cbind(df_storms[,c("date", "time", "lon", "lat")], 
                    df_total[,c("storm", "lab", paste0("p", lab_vec))]),
          file = file_name,
          row.names = FALSE)