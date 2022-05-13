## Objective Identification: Station data
# Training and prediction via random forests
# Includes generation of an evaluation file and reliability diagrams

#### Settings ####
# Load package
library(lubridate)
library(dplyr)
library(ranger)
library(rpart)
library(ggplot2)
library(reliabilitydiag)
library(gridExtra)

# Path of RF objects
data_rf_path <- "/.../"

# Path of R data
data_path <- "/.../"

# Path for evaluation results
eval_path <- "/.../"

#### Load data ####
# Name of file
# file_name <- paste0(data_path, "station_data.RData") # Not provided
file_name <- paste0(data_path, "station_fake_data.RData")

# Load data frame
load(file = file_name)

#### Initiation ####
# Vector of labels
lab_vec <- sort(unique(df_storms$lab))

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

# Formula
pred_fm <- as.formula(paste0("lab ~ ", paste0(pred_vars,
                                                collapse = " + ")))

# Data frame for overall evaluation
df_total <- data.frame(lab = vector(length = nrow(df_storms)),
                       storm = vector(length = nrow(df_storms)))

# Probability for each label
for(i in lab_vec){ df_total[[paste0("p", i)]] <- df_total[[paste0("p_clim", i)]] <- 
  vector(length = nrow(df_storms)) }

# Make a list for permuted forecasts
fi_ls <- list()

# Generate matrix for each predictor
for(temp_var in pred_vars){
  # Data frame for overall evaluation
  fi_ls[[temp_var]] <- data.frame(lab = vector(length = nrow(df_storms)),
                                  storm = vector(length = nrow(df_storms)))
  
  # Probability for each label
  for(i in lab_vec){ fi_ls[[temp_var]][[paste0("p", i)]] <- 
    vector(length = nrow(df_storms)) }
}

# Check which columns include missing values
na_cols <- vector(length = ncol(df_storms))
for(i in 1:ncol(df_storms)){ na_cols[i] <- any(is.na(df_storms[,i])) }

# Vector of scores to consider
sr_vec <- c("mbs", "bs", "mcb", "dsc", "unc")

# Variables for data frame
col_vec <- c("storm", "n", "hits", "method", "metric", "lab", "value")

# Data frame for scores
df_scores <- data.frame(matrix(nrow = 0,
                               ncol = length(col_vec)))
colnames(df_scores) <- col_vec

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

#### Start summary file ####
# Name of .txt file for output
file_sink <- paste0(eval_path, "rf_eval.txt")

# Initialize .txt-files
sink(file = file_sink, 
     type = "output")

# New line
cat(" ------ Random Forest: Cross-Validation over storms ")
cat("\n \n ------ Handling of missing values: Replace with mean")

# End sink
sink()

#### Estimation and prediction ####
# For-Loop over storms
for(temp_date in date_vec){
  # Name of storm
  temp_name <- name_vec[which(date_vec == temp_date)]
  
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
  
  # Estimate
  est <- ranger(formula = pred_fm,
                data = df_train[order(df_train$lab),],
                importance = "permutation",
                probability = TRUE,
                num.trees = 1000)

  # Predict
  pred <- predict(object = est,
                  data = df_test)

  # Name of file
  file_save <- paste0(data_rf_path, temp_name, ".RData")

  # Save data
  save(file = file_save,
       list = c("est", "pred"))
  
  # Climatology
  p_clim <- sapply(lab_vec, function(x) mean(df_train$lab == x))
  
  # Read out (Every label present in at least two storms)
  df_total[,paste0("p_clim", lab_vec)][i_test,] <- 
    matrix(data = p_clim,
           ncol = length(lab_vec),
           nrow = length(i_test),
           byrow = TRUE)
  df_total[,paste0("p", lab_vec)][i_test,] <- pred$predictions
  df_total[["lab"]][i_test] <- df_test$lab
  df_total[["storm"]][i_test] <- temp_name
  
  # Make predictions for each variable
  for(temp_var in pred_vars){
    # Copy test set for permutation
    df_perm <- df_test

    # Permute variable in test set
    df_perm[[temp_var]] <- sample(x = df_perm[[temp_var]],
                                  size = nrow(df_perm),
                                  replace = FALSE)

    # Predict
    pred <- predict(object = est,
                    data = df_perm)

    # Read out
    fi_ls[[temp_var]][,paste0("p", lab_vec)][i_test,] <- pred$predictions
    fi_ls[[temp_var]][["lab"]][i_test] <- df_perm$lab
    fi_ls[[temp_var]][["storm"]][i_test] <- temp_name
  }
  
  # Make lists
  ls_rd <- ls_plot <- list()
  
  # Loop over labels
  for(i in lab_vec){
    # Make reliability diagram of predictions
    ls_rd[[paste0("p", i)]] <- reliabilitydiag(p = df_total[[paste0("p", i)]][i_test],
                                               y = as.numeric(df_total[["lab"]][i_test] == i))
    
    # Make reliability diagram of climatology
    ls_rd[[paste0("p_clim", i)]] <- reliabilitydiag(p = df_total[[paste0("p_clim", i)]][i_test],
                                                    y = as.numeric(df_total[["lab"]][i_test] == i), 
                                                    region.level = NA)
    
    # Plot
    ls_plot[[paste0("p", i)]] <- ggplot2::autoplot(
      ls_rd[[paste0("p", i)]],
      params_histogram = list(breaks = (0:50)/50,
                              fill = NA,
                              colour = "black")) + 
      ggplot2::ggtitle(lab = paste0("Label: ", i, 
                                      "; n: ", format(length(i_test), big.mark = ","), 
                                      "; hits: ", format(sum(df_total[["lab"]][i_test] == i), big.mark = ","))) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) 
  }
  
  # For-Loop over scores, methods and labels
  for(temp_meth in c("p", "p_clim")){ for(temp_sr in sr_vec){ for(i_lab in lab_vec){
    # Skip for multivariate BS
    if((temp_sr == "mbs") & (i_lab > 0)){ next }
    
    # Get current index
    i_row <- nrow(df_scores) + 1
    
    # Write in data frame
    df_scores[i_row, "storm"] <- temp_name
    df_scores[i_row, "n"] <- nrow(df_test)
    df_scores[i_row, "method"] <- temp_meth
    df_scores[i_row, "metric"] <- temp_sr
    
    # Special case: Multivariate BS
    if(temp_sr == "mbs"){
      df_scores[i_row, "metric"] <- temp_sr
      df_scores[i_row, "value"] <- 
        brier_score(df_p = df_total[i_test, paste0(temp_meth, lab_vec)],
                    y_vec = df_total[i_test, "lab"])
    }else{
      df_scores[i_row, "lab"] <- i_lab
      df_scores[i_row, "hits"] <- sum(df_test[["lab"]] == i_lab)
      
      # Get index of score
      if(temp_sr == "bs"){ i_rd <- 2 }
      else if(temp_sr == "mcb"){ i_rd <- 3 }
      else if(temp_sr == "dsc"){ i_rd <- 4 }
      else if(temp_sr == "unc"){ i_rd <- 5 }
      
      # Read out of reliability diagram
      df_scores[i_row, "value"] <- as.numeric(summary(ls_rd[[paste0(temp_meth, i_lab)]])[i_rd])
    }
  }}}
  
  # Open .txt-files
  sink(file = file_sink, 
       type = "output",
       append = TRUE)
  
  # Which storm?
  cat(paste0("\n \n --- Storm: ", temp_name))
  
  # Occurrences
  cat(paste0("\n Table - Train: \n"))
  print(table(df_train$lab))
  cat(paste0("\n Table - Test: \n"))
  print(table(df_test$lab))
  
  # Evaluation
  cat(paste0("\n Mean Brier score of climatology: ", 
             round(subset(df_scores, (storm == temp_name) &
                            (method == "p_clim") & (metric == "mbs"))[["value"]], 4)))
  cat(paste0("\n Mean Brier score: ", 
             round(subset(df_scores, (storm == temp_name) &
                            (method == "p") & (metric == "mbs"))[["value"]], 4)))
  
  # Importance
  cat(paste0("\n \n Importance mode: ", est$importance.mode))
  cat("\n \n Variable importance: \n ")
  print(round(sort(est$variable.importance, decreasing = TRUE), 4))
  
  # New line
  cat("\n\n Scores of reliability diagrams")
  
  # For-Loop over labels
  for(i in lab_vec){
    # Scores of reliability diagrams
    cat(paste0("\n i = ", i, ". Climatology: \n"))
    print(round(as.matrix(summary(ls_rd[[paste0("p_clim", i)]])[,-1]), 4))
    cat(paste0("\n i = ", i, ". Random Forest: \n"))
    print(round(as.matrix(summary(ls_rd[[paste0("p", i)]])[,-1]), 4))
  }
  
  # New line
  cat("\n")
  
  # End sink
  sink()
  
  # Name of PDF
  file_pdf <- paste0(eval_path, "rf_rd_", temp_name, ".pdf")

  # Plot together
  pdf_plot <- marrangeGrob(grobs = ls_plot[is.element(names(ls_plot), paste0("p", lab_vec))],
                           layout_matrix = matrix(1:(4*ceiling(length(lab_vec)/4)),
                                                  ncol = 4,
                                                  byrow = TRUE),
                           top = NULL)

  # Save PDF
  ggplot2::ggsave(filename = file_pdf,
                  plot = pdf_plot,
                  width = 20,
                  height = 5*ceiling(length(lab_vec)/4))
}

#### Overall diagrams and scores ####
# Make lists
ls_rd <- ls_plot <- list()

# Loop over labels
for(i in lab_vec){ 
  # Make reliability diagram of predictions
  ls_rd[[paste0("p", i)]] <- reliabilitydiag(p = df_total[[paste0("p", i)]],
                                             y = as.numeric(df_total[["lab"]] == i))
  
  # Make reliability diagram of climatologies
  ls_rd[[paste0("p_clim", i)]] <- reliabilitydiag(p = df_total[[paste0("p_clim", i)]],
                                                  y = as.numeric(df_total[["lab"]] == i))
  
  # Plot
  ls_plot[[paste0("p", i)]] <- ggplot2::autoplot(
    ls_rd[[paste0("p", i)]],
    params_histogram = list(breaks = (0:50)/50,
                            fill = NA,
                            colour = "black")) + 
    ggplot2::ggtitle(lab = paste0("Label: ", i, 
                                    "; n: ", format(nrow(df_total), big.mark = ","), 
                                    "; hits: ", format(sum(df_total[["lab"]] == i), big.mark = ","))) +
    ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))  
}

# For-Loop over scores, methods and labels
for(temp_meth in c("p", "p_clim")){ for(temp_sr in sr_vec){ for(i_lab in lab_vec){
  # Skip for multivariate BS
  if((temp_sr == "mbs") & (i_lab > 0)){ next }
  
  # Get current index
  i_row <- nrow(df_scores) + 1
  
  # Write in data frame
  df_scores[i_row, "storm"] <- "Total"
  df_scores[i_row, "n"] <- nrow(df_total)
  df_scores[i_row, "method"] <- temp_meth
  df_scores[i_row, "metric"] <- temp_sr
  
  # Special case: Multivariate BS
  if(temp_sr == "mbs"){
    df_scores[i_row, "metric"] <- temp_sr
    df_scores[i_row, "value"] <- 
      brier_score(df_p = df_total[,paste0(temp_meth, lab_vec)],
                  y_vec = df_total[,"lab"])
  }else{
    df_scores[i_row, "lab"] <- i_lab
    df_scores[i_row, "hits"] <- sum(df_total[["lab"]] == i_lab)
    
    # Get index of score
    if(temp_sr == "bs"){ i_rd <- 2 }
    else if(temp_sr == "mcb"){ i_rd <- 3 }
    else if(temp_sr == "dsc"){ i_rd <- 4 }
    else if(temp_sr == "unc"){ i_rd <- 5 }
    
    # Read out of reliability diagram
    df_scores[i_row, "value"] <- as.numeric(summary(ls_rd[[paste0(temp_meth, i_lab)]])[i_rd])
  }
}}}

# Open .txt-files
sink(file = file_sink, 
     type = "output",
     append = TRUE)

# Which storm?
cat("\n \n ------ Total")

# Occurrences
cat(paste0("\n Table: \n"))
print(table(df_total$lab))

# Evaluation
cat(paste0("\n Mean Brier score of climatology: ", 
           round(subset(df_scores, (storm == "Total") &
                          (method == "p_clim") & (metric == "mbs"))[["value"]], 4)))
cat(paste0("\n Mean Brier score: ", 
           round(subset(df_scores, (storm == "Total") &
                          (method == "p") & (metric == "mbs"))[["value"]], 4)))

# New line
cat("\n\n Scores of reliability diagrams")

# For-Loop over labels
for(i in lab_vec){
  # Scores of reliability diagrams
  cat(paste0("\n i = ", i, ". Climatology: \n"))
  print(round(as.matrix(summary(ls_rd[[paste0("p_clim", i)]])[,-1]), 4))
  cat(paste0("\n i = ", i, ". Random Forest: \n"))
  print(round(as.matrix(summary(ls_rd[[paste0("p", i)]])[,-1]), 4))
}

# New line
cat("\n")

# End sink
sink()

# Name of PDF
file_pdf <- paste0(eval_path, "rf_rd_storms.pdf")

# Plot together
pdf_plot <- marrangeGrob(grobs = ls_plot[is.element(names(ls_plot), paste0("p", lab_vec))],
                         layout_matrix = matrix(1:(4*ceiling(length(lab_vec)/4)),
                                                ncol = 4,
                                                byrow = TRUE),
                         top = NULL)

# Save PDF
ggplot2::ggsave(filename = file_pdf,
                plot = pdf_plot,
                width = 20,
                height = 5*ceiling(length(lab_vec)/4))

#### Save data ####
# Name of file
file_save <- paste0(data_path, "station_preds.RData")

# Save data
save(file = file_save,
     list = c("df_storms", "df_total", "fi_ls", "df_scores"))

#### Calculate Feature Importance ####
# Make list for BS
bs_ls <- list()

# Function for binary BS
bs <- function(p, y){
  ###-----------------------------------------------------------------------------
  ###Input
  #p...........n probability forecasts (n vector)
  #y...........n binary observations (n vector)
  ###-----------------------------------------------------------------------------
  ###Output
  #res...Mean brier score (scalar)
  ###-----------------------------------------------------------------------------
  
  #### Calculation ####
  return(mean((p-y)^2))
}

# For-Loop over variables (and default)
for(temp_var in c("def", names(fi_ls))){
  # Get forecasts
  if(temp_var == "def"){ fi_preds <- df_total }
  else{ fi_preds <- fi_ls[[temp_var]] }
  
  # Create one-vs-all names
  temp <- paste0("p", lab_vec)
  
  # Pairwise column names
  for(i0 in lab_vec){ for(i1 in lab_vec[lab_vec > i0]){
    temp <- c(temp, paste0("p", i0, i1))
  }}
  
  # Create BS matrix
  bs_ls[[temp_var]] <- matrix(nrow = 1 + length(date_vec),
                              ncol = length(temp))
  
  # Name columns
  colnames(bs_ls[[temp_var]]) <- temp
  row.names(bs_ls[[temp_var]]) <- c("all", name_vec)
  
  # Calculate Brier scores
  for(temp_row in row.names(bs_ls[[temp_var]])){
    # Get subset of forecasts
    if(temp_row == "all"){ temp_df <- fi_preds }
    else{ temp_df <- subset(fi_preds, storm == temp_row) }
    
    # For-Loop over labels
    for(i_lab in lab_vec){
      # Calculate one-vs-all BS
      bs_ls[[temp_var]][temp_row, paste0("p", i_lab)] <- 
        bs(p = temp_df[,paste0("p", i_lab)],
           y = (temp_df[,"lab"] == i_lab))
      
      # Calculate pairwise BS
      for(j_lab in lab_vec[lab_vec > i_lab]){
        # Pair
        temp_pair <- c(i_lab, j_lab)
        
        # Get relevant subset
        df_cp <- subset(temp_df, is.element(lab, temp_pair))[,c("lab", paste0("p", temp_pair))]
        
        # If no pairs included -> BS of 0
        if(nrow(df_cp) == 0){
          # Set BS to zero
          bs_ls[[temp_var]][temp_row, paste0("p", i_lab, j_lab)] <- 0
          
          # Go to next pair
          next
        }
        
        # If both probabilities are zero, both are equally likely -> 0.5
        df_cp[which(rowSums(df_cp[,paste0("p", temp_pair)]) == 0), paste0("p", temp_pair)] <- rep(1/2, 2)
        
        # Calculate conditional probabilities
        df_cp[,paste0("p", temp_pair)] <- 
          df_cp[,paste0("p", temp_pair)]/rowSums(df_cp[,paste0("p", temp_pair)])
        
        # Calculate BS
        bs_ls[[temp_var]][temp_row, paste0("p", i_lab, j_lab)] <-
          bs(p = df_cp[,paste0("p", i_lab)],
             y = (df_cp[,"lab"] == i_lab))
      }
    }
  }
}

# Calculate (relative) BS importance
fi_bs_ls <- fi_bs_rel_ls <- list()

# For-Loop over variables
for(temp_var in pred_vars){
  # Calculate BS permutation importance
  fi_bs_ls[[temp_var]] <- bs_ls[[temp_var]] - bs_ls[["def"]] 
  
  # Calculate relative BS permutation importance
  fi_bs_rel_ls[[temp_var]] <- fi_bs_ls[[temp_var]]/bs_ls[["def"]] 
  
  # Set to zero if default and BS is equal to zero
  fi_bs_rel_ls[[temp_var]][(bs_ls[["def"]] == 0) & (bs_ls[[temp_var]] == 0)] <- 0
}

# Save importance
save(file = paste0(data_path, "station_fi.RData"),
     list = c("fi_bs_ls", "fi_bs_rel_ls", "bs_ls"))
