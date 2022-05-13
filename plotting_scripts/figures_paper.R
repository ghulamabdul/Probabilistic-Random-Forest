## Objective Identification: Station data
# Training and prediction via random forests

#### Settings ####
# Load package
library(lubridate)
library(dplyr)
library(ranger)
library(rpart)
library(ggplot2)
library(reliabilitydiag)
library(gridExtra)

# Path of R-data
data_path <- "/.../"

# Path of COSMO data
data_cosmo_path <- "/.../"

# Path for evaluation results
pdf_path <- "/.../"

#### Initiation ####
# Names of the winter storms
name_vec <- c("Niklas", "Susanna", 
              "Egon", "Thomas", "Xavier", "Herwart",
              "Burglind", "Friederike", "Fabienne",
              "Bennet", "Eberhard",
              "Sabine")

# Labels
lab_vec <- c(0, 1, 2, 3, 5)

# Feature of label
lab_names <- c("0" = "NF",
               "1" = "WJ",
               "2" = "CFC",
               "3" = "CJ",
               "5" = "CS")

# Colors of features
lab_cols <- c("NF" = "#808080", # grey
              "WJ" = "#FF0000", # red
              "CFC" = "#008000", # green
              "CJ" = "#0000FF", # blue
              "CS" = "#FFD700") # gold

# Name of predictor variables
var_names <- c("mslp" = "Mean Sea Level Pressure",
               "dpdt" = "Mean Sea Level Pressure Tendency",
               "th" = "Normalised Potential Temperature",
               "dthdt" = "Normalised Potential Temperature Tendency",
               "prec" = "Precipitation",
               "vgl" = "Normalised Wind Speed",
               "wd" = "Wind Direction",
               "dwddt" = "Wind Direction Tendency")

# Abbreviations of predictor variables
var_abbr <- c( "mslp" = "\"p\"",
               "dpdt" = "Delta * \"p\"",
               "th" = "tilde(theta)",
               "dthdt" = "Delta * tilde(theta)",
               "prec" = "\"RR\"",
               "vgl" = "tilde(v)",
               "wd" = "\"d\"",
               "dwddt" = "Delta * \"d\"")

# In italic
var_abbr[] <- paste0("italic(", var_abbr, ")")

# Order predictor variables
pred_vars <- names(var_names)

# Colors of variables
var_cols <- RColorBrewer::brewer.pal(n = 8,
                                     name = "Dark2")
var_cols <- var_cols[c(7, 6, 4, 2,
                       3, 8, 1, 5)]
names(var_cols) <- pred_vars

#### Stations: Load data ####
# Name of file
file_name <- paste0(data_path, "station_preds.RData")

# Load data frame
load(file = file_name)

# Storms
date_vec <- unique(df_storms[["date"]])

#### Stations: Reliability diagrams: One vs. all ####
# Make lists
ls_rd <- ls_plot <- list()

# Loop over labels
for(i in lab_vec){ 
  # Make reliability diagram of predictions
  ls_rd[[paste0(i)]] <- reliabilitydiag(p = df_total[[paste0("p", i)]],
                                        y = as.numeric(df_total[["lab"]] == i))
}

# Loop over labels
for(i in lab_vec){ 
  # Make plots
  ls_plot[[paste0(i)]] <- ggplot2::autoplot(
    ls_rd[[paste0(i)]],
    params_histogram = list(breaks = (0:50)/50,
                            fill = NA,
                            colour = "black")) +
    ggtitle(paste0("(", letters[which(i == lab_vec)], ") 1 = ", lab_names[paste0(i)])) +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          plot.title = element_text(face = "bold",
                                    size = 25)) +
    annotate("label", 
             label = paste0("Samples: ", 
                            format(sum(ls_rd[[paste0(i)]][["p"]][["bins"]][["n"]]), big.mark = ","), 
                            "\n", "Hits: ", 
                            format(sum(ls_rd[[paste0(i)]][["p"]][["cases"]][["y"]]), big.mark = ",")), 
             x = 0, 
             y = 0.92,
             hjust = 0,
             size = 6.5,
             label.padding = unit(0.55, "lines"))
}

# Plot together
pdf_plot <- marrangeGrob(grobs = ls_plot,
                         layout_matrix = matrix(1:5,
                                                ncol = 5,
                                                byrow = TRUE),
                         top = NULL)

# Name of PDF
file_pdf <- paste0(pdf_path, "reliability_one_vs_all.pdf")

# Save PDF
ggplot2::ggsave(filename = file_pdf,
                plot = pdf_plot,
                width = 25,
                height = 5)

#### Stations: Reliability diagrams: all-pairs ####
# Make matrix of all pairs
mtx_pairs <- as.matrix(expand.grid(lab_vec, lab_vec))
mtx_pairs <- mtx_pairs <- mtx_pairs[mtx_pairs[,1] > mtx_pairs[,2],]

# List for plots
ls_rd <- ls_plot <- list()

# For-Loop over all pairs
for(i_pair in 1:nrow(mtx_pairs)){
  # Get corresponding pair
  temp_pair <- mtx_pairs[i_pair, c(2, 1)]
  
  # Get relevant subset
  df_cp <- subset(df_total, is.element(lab, temp_pair))[,c("lab", paste0("p", temp_pair))]
  
  # If both probabilities are zero, both are equally likely -> 0.5
  df_cp[which(rowSums(df_cp[,paste0("p", temp_pair)]) == 0), paste0("p", temp_pair)] <- rep(1/2, 2)
  
  # Calculate conditional probabilities
  df_cp[,paste0("p", temp_pair)] <- 
    df_cp[,paste0("p", temp_pair)]/rowSums(df_cp[,paste0("p", temp_pair)])
  
  # Make reliability diagram of predictions
  ls_rd[[paste0(i_pair)]] <- reliabilitydiag(p = df_cp[[paste0("p", temp_pair[2])]],
                                             y = as.numeric(df_cp[["lab"]] == temp_pair[2]))
}

# For-Loop over all pairs
for(i_pair in 1:nrow(mtx_pairs)){
  # Get corresponding pair
  temp_pair <- mtx_pairs[i_pair, c(2, 1)]
  
  # Get names
  temp_name0 <- lab_names[paste0(temp_pair[1])]
  temp_name1 <- lab_names[paste0(temp_pair[2])]
  
  # Plot
  ls_plot[[paste0(i_pair)]] <- ggplot2::autoplot(
    ls_rd[[paste0(i_pair)]],
    params_histogram = list(breaks = (0:50)/50,
                            fill = NA,
                            colour = "black")) +
    ggtitle(paste0("(", letters[i_pair], ") 0 = ", temp_name0, ", 1 = ", temp_name1)) +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          plot.title = element_text(face = "bold",
                                    size = 25)) +
    annotate("label", 
             label = paste0("Samples: ", 
                            format(sum(ls_rd[[paste0(i_pair)]][["p"]][["bins"]][["n"]]), big.mark = ","), 
                            "\n", "Hits: ", 
                            format(sum(ls_rd[[paste0(i_pair)]][["p"]][["cases"]][["y"]]), big.mark = ",")), 
             x = 0, 
             y = 0.92,
             hjust = 0,
             size = 6.5,
             label.padding = unit(0.55, "lines"))
}

# Plot together
pdf_plot <- marrangeGrob(grobs = ls_plot,
                         layout_matrix = matrix(1:10,
                                                ncol = 5,
                                                byrow = TRUE),
                         top = NULL)

# Name of PDF
file_pdf <- paste0(pdf_path, "reliability_all_pairs.pdf")

# Save PDF
ggplot2::ggsave(filename = file_pdf,
                plot = pdf_plot,
                width = 25,
                height = 10)


#### Stations: Predictor importance: BS One vs. All ####
# Load BS importance
load(file = paste0(data_path, "station_fi.RData"))

# Make data frame
df_pi <- data.frame(storm = character(),
                    var = character(),
                    var0 = character(),
                    feature = character(),
                    pi = numeric(),
                    stringsAsFactors = FALSE)


# Start with first row
i <- 1

# For-Loop over variables, storms and features
for(temp_storm in name_vec){ for(temp_pred in pred_vars){ for(temp_lab in lab_vec){
  # # Read out importance
  df_pi[i, "storm"] <- temp_storm
  df_pi[i, "var"] <- var_abbr[temp_pred]
  df_pi[i, "var0"] <- temp_pred
  df_pi[i, "feature"] <- lab_names[paste0(temp_lab)]
  df_pi[i, "pi"] <- fi_bs_ls[[temp_pred]][temp_storm, paste0("p", temp_lab)]
  
  # Next row
  i <- i + 1
}}}

# Order of variables
var_order <- c(1, 4, 2, 3, 5:8)
var_order <- 1:8

# Features and variables as factor
df_pi[["feature"]] <- factor(x = df_pi[["feature"]],
                             levels = lab_names)
df_pi[["var"]] <- factor(x = df_pi[["var"]],
                         levels = var_abbr)
df_pi[["var0"]] <- factor(x = df_pi[["var0"]],
                         levels = pred_vars[var_order])

# Panel per predictor: For-Loop over scale of y-axis
for(temp_scale in c("free", "fixed")){
  # Make plot
  pdf_plot <- 
    ggplot(df_pi, 
           aes(x = feature, y = pi,
               color = feature, fill = feature)) + 
    facet_wrap("var",
               ncol = 4,
               scales = temp_scale,
               labeller = label_parsed) +
    geom_hline(yintercept = 0,
               color = "lightgrey",
               linetype = "dashed") +
    geom_boxplot(alpha = 0.5,
                 outlier.size = 1.7) + 
    scale_color_manual(values = lab_cols) +
    scale_fill_manual(values = lab_cols) +
    theme_bw() +
    theme(legend.position = "top", 
          legend.title = element_blank(),
          legend.text = element_text(size = 20),
          legend.spacing.x = unit(1.0, 'cm'),
          legend.key.width = unit(1.5, 'cm'),
          axis.text.x = element_text(size = 18),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          strip.text.x = element_text(size = 20),
          strip.text.y = element_text(size = 20)) +
    guides(colour = guide_legend(nrow = 1,
                                 label.position = "bottom"),
           fill = guide_legend(nrow = 1,
                               label.position = "bottom")) + 
    xlab("Wind Feature") + 
    ylab("Brier Score Permutation Importance")
  
  # Name of file
  file_pdf <- paste0(pdf_path, "pred_importance_one_vs_all_", temp_scale, ".pdf")
  
  # Save plot
  ggplot2::ggsave(filename = file_pdf,
                  plot = pdf_plot,
                  height = 10,
                  width = 20,
                  scale = 0.8)
}

# Panel per wind feature: For-Loop over scale of y-axis
for(temp_scale in c("free", "fixed")){
  # Make plot
  pdf_plot <- 
    ggplot(df_pi, 
           aes(x = var, y = pi,
               color = var0, fill = var0)) + 
    facet_wrap("feature",
             ncol = 5,
             scales = temp_scale,
             labeller = label_parsed) +
    geom_hline(yintercept = 0,
             color = "lightgrey",
             linetype = "dashed") +
    geom_boxplot(alpha = 0.5,
                 outlier.size = 1.7) + 
    scale_color_manual(values = var_cols[var_order],
                       labels = as.expression(parse(text = var_abbr))) +
    scale_fill_manual(values = var_cols[var_order],
                      labels = as.expression(parse(text = var_abbr))) +
    scale_x_discrete(labels = as.expression(parse(text = var_abbr))) +
    theme_bw() +
    theme(legend.position = "top", 
          legend.title = element_blank(),
          legend.text = element_text(size = 20),
          legend.spacing.x = unit(1.0, 'cm'),
          legend.key.width = unit(1.5, 'cm'),
          axis.text.x = element_text(size = 16,
                                     vjust = 0),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 18,
                                      vjust = 0),
          axis.title.y = element_text(size = 18),
          strip.text.x = element_text(size = 20),
          strip.text.y = element_text(size = 20)) +
    # guides(colour = guide_legend(nrow = 1)) +
    guides(colour = guide_legend(nrow = 1,
                                 label.position = "bottom"),
           fill = guide_legend(nrow = 1,
                               label.position = "bottom")) + 
    xlab("Predictor Variable") + 
    ylab("Brier Score Permutation Importance")
  
  # Name of file
  file_pdf <- paste0(pdf_path, "pred_importance_one_vs_all_feature_", temp_scale, ".pdf")
  
  # Save plot
  ggplot2::ggsave(filename = file_pdf,
                  plot = pdf_plot,
                  height = 6,
                  width = 20,
                  scale = 0.8)
}


#### Stations: Predictor importance: BS all-pairs ####
# Load BS importance
load(file = paste0(data_path, "station_fi.RData"))

# Get possible pairs
mtx_pairs <- as.matrix(expand.grid(lab_vec, lab_vec))
mtx_pairs <- mtx_pairs[mtx_pairs[,1] > mtx_pairs[,2],]
pair_vec <- apply(mtx_pairs, 1, function(x) paste0(rev(x), collapse = ""))

# Names of pairs
pair_names <- sapply(pair_vec, function(x){
  paste0(lab_names[substr(x, 1, 1)], " vs. ",lab_names[substr(x, 2, 2)]) })
  
# Make data frame
df_pi <- data.frame(storm = character(),
                    var = character(),
                    var0 = character(),
                    pair = character(),
                    pi = numeric(),
                    stringsAsFactors = FALSE)


# Start with first row
i <- 1

# For-Loop over variables, storms and features
for(temp_storm in name_vec){ for(temp_pred in pred_vars){ for(temp_pair in pair_vec){
  # # Read out importance
  df_pi[i, "storm"] <- temp_storm
  df_pi[i, "var"] <- var_abbr[temp_pred]
  df_pi[i, "var0"] <- temp_pred
  df_pi[i, "pair"] <- pair_names[temp_pair]
  df_pi[i, "pi"] <- fi_bs_ls[[temp_pred]][temp_storm, paste0("p", temp_pair)]
  
  # Next row
  i <- i + 1
}}}

# Order of variables
var_order <- c(1, 4, 2, 3, 5:8)
var_order <- 1:8

# Variables and pairs as factor
df_pi[["var"]] <- factor(x = df_pi[["var"]],
                         levels = var_abbr)
df_pi[["var0"]] <- factor(x = df_pi[["var0"]],
                          levels = pred_vars[var_order])
df_pi[["pair"]] <- factor(x = df_pi[["pair"]],
                         levels = pair_names)

# For-Loop over scale of y-axis
for(temp_scale in c("free", "fixed")){
  # Make plot
  pdf_plot <- 
    ggplot(df_pi, 
           aes(x = var, y = pi,
               color = var0, fill = var0)) + 
    facet_wrap("pair",
               ncol = 5,
               scales = temp_scale) +
    geom_hline(yintercept = 0,
               color = "lightgrey",
               linetype = "dashed") +
    geom_boxplot(alpha = 0.5,
                 outlier.size = 1.7) + 
    scale_color_manual(values = var_cols[var_order],
                       labels = as.expression(parse(text = var_abbr))) +
    scale_fill_manual(values = var_cols[var_order],
                      labels = as.expression(parse(text = var_abbr))) +
    scale_x_discrete(labels = as.expression(parse(text = var_abbr))) +
    theme_bw() +
    theme(legend.position = "top", 
          legend.title = element_blank(),
          legend.text = element_text(size = 20),
          legend.spacing.x = unit(1.0, 'cm'),
          legend.key.width = unit(1.5, 'cm'),
          axis.text.x = element_text(size = 16,
                                     vjust = 0),
          axis.text.y = element_text(size = 18),
          axis.title.x = element_text(size = 18,
                                      vjust = 0),
          axis.title.y = element_text(size = 18),
          strip.text.x = element_text(size = 20),
          strip.text.y = element_text(size = 20)) +
    # guides(colour = guide_legend(nrow = 2)) + 
    guides(colour = guide_legend(nrow = 1,
                                 label.position = "bottom"),
           fill = guide_legend(nrow = 1,
                               label.position = "bottom")) + 
    xlab("Predictor Variable") + 
    ylab("Brier Score Permutation Importance")
  
  # Name of file
  file_pdf <- paste0(pdf_path, "pred_importance_all_pairs_", temp_scale, ".pdf")
  
  # Save plot
  ggplot2::ggsave(filename = file_pdf,
                  plot = pdf_plot,
                  height = 10,
                  width = 20,
                  scale = 0.8)
}

#### Stations: PDP: Average curves ####
# x-Limits for variables
x_lim <- list(
  "vgl" = c(0.8, 2.5),
  "wd" = c(90, 360),
  "dwddt" = c(-50, 100),
  "prec" = c(0, 6),
  "dpdt" = c(-2.5, 3.5),
  "th" = c(0.96, 1.05)
)

# Load PDP data
load(file = paste0(data_path, "station_pdp.RData"))

# For-Loop over variables and labels
for(temp_var in names(x_lim)){ for(temp_lab in lab_vec){
  # Read out data frame
  temp_df <- pdp_ls[[temp_var]][[paste0(temp_lab)]]
  
  # Indices of relevant rows
  i_rows <- (x_lim[[temp_var]][1] <= temp_df[[temp_var]]) & 
    (temp_df[[temp_var]] <= x_lim[[temp_var]][2])
  
  # Cut rows
  temp_df <- temp_df[i_rows,]
  
  # Save data frame
  pdp_ls[[temp_var]][[paste0(temp_lab)]] <- temp_df
}}

# Make data frame
df_pi <- data.frame(var = character(),
                    feature = character(),
                    x = numeric(),
                    p = numeric(),
                    stringsAsFactors = FALSE)


# Start with first row
i <- 1

# For-Loop over variables, storms and features
for(temp_pred in pred_vars){ for(temp_lab in lab_vec){
  # Number of evaluations
  n_grid <- nrow(pdp_ls[[temp_pred]][[paste0(temp_lab)]])
  
  # For-Loop over evaluations
  for(j in 1:n_grid){
    # # Read out importance
    df_pi[i, "feature"] <- lab_names[paste0(temp_lab)]
    df_pi[i, "var"] <- var_abbr[temp_pred]
    df_pi[i, "x"] <- pdp_ls[[temp_pred]][[paste0(temp_lab)]][[temp_pred]][j]
    df_pi[i, "p"] <- pdp_ls[[temp_pred]][[paste0(temp_lab)]][["yhat"]][j]
    
    # Next row
    i <- i + 1
  }
}}

# Features and variables as factor
df_pi[["feature"]] <- factor(x = df_pi[["feature"]],
                             levels = lab_names)
df_pi[["var"]] <- factor(x = df_pi[["var"]],
                         levels = var_abbr)

# Make plot
pdf_plot <- 
  ggplot(df_pi, 
         aes(x = x, y = p, color = feature)) + 
  facet_wrap("var",
             ncol = 4,
             scales = "free_x",
             labeller = label_parsed) +
  geom_line(size = 1) +
  scale_color_manual(values = lab_cols) +
  theme_bw() +
  theme(legend.position = "top", 
        legend.title = element_blank(),
        legend.text = element_text(size = 18),
        legend.spacing.x = unit(1.0, 'cm'),
        legend.key.width = unit(1.5, 'cm'),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 12),
        axis.title.x = element_text(size = 15),
        axis.title.y = element_text(size = 15),
        strip.text.x = element_text(size = 18),
        strip.text.y = element_text(size = 18)) +
  guides(colour = guide_legend(nrow = 1,
                               label.position = "bottom",
                               override.aes = list(size = 5))) + 
  ylab("Probability") + 
  xlab("Value of Predictor Variable")

# Save plot
ggplot2::ggsave(filename = paste0(pdf_path, "pdp.pdf"),
                plot = pdf_plot,
                height = 10,
                width = 20,
                scale = 0.7)

#### COSMO-REA: Load data ####
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

# Load data
df_total <- do.call(rbind, lapply(names(date_vec), function(x) 
  read.csv(paste0(data_cosmo_path, "cosmo_preds_", x, ".csv"))))

# Transform dates to names
df_total[["storm"]] <- names(date_vec)[sapply(df_total[["date"]], function(x) which(date_vec == x))]

# Cut columns
df_total <- df_total[,!is.element(colnames(df_total), c("date", "time", "lon", "lat"))]

# Storms
date_vec <- unique(df_storms[["date"]])

#### COSMO-REA: Reliability diagrams: One vs. all #####
# Make lists
ls_rd <- ls_plot <- list()

# Loop over labels
for(i in lab_vec){ 
  # Make reliability diagram of predictions
  ls_rd[[paste0(i)]] <- reliabilitydiag(p = df_total[[paste0("p", i)]],
                                        y = as.numeric(df_total[["lab"]] == i))
}

# Loop over labels
for(i in lab_vec){ 
  # Make plots
  ls_plot[[paste0(i)]] <- ggplot2::autoplot(
    ls_rd[[paste0(i)]],
    params_histogram = list(breaks = (0:50)/50,
                            fill = NA,
                            colour = "black")) +
    ggtitle(paste0("(", letters[which(i == lab_vec)], ") 1 = ", lab_names[paste0(i)])) +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          plot.title = element_text(face = "bold",
                                    size = 25)) +
    annotate("label", 
             label = paste0("Samples: ", 
                            format(sum(ls_rd[[paste0(i)]][["p"]][["bins"]][["n"]]), big.mark = ","), 
                            "\n", "Hits: ", 
                            format(sum(ls_rd[[paste0(i)]][["p"]][["cases"]][["y"]]), big.mark = ",")), 
             x = 0, 
             y = 0.92,
             hjust = 0,
             size = 6.5,
             label.padding = unit(0.55, "lines"))
}

# Plot together
pdf_plot <- marrangeGrob(grobs = ls_plot,
                         layout_matrix = matrix(1:5,
                                                ncol = 5,
                                                byrow = TRUE),
                         top = NULL)

# Name of PDF
file_pdf <- paste0(pdf_path, "cosmo_rea_reliability_one_vs_all.pdf")

# Save PDF
ggplot2::ggsave(filename = file_pdf,
                plot = pdf_plot,
                width = 25,
                height = 5)

#### COSMO-REA: Reliability diagrams: all-pairs #####
# Make matrix of all pairs
mtx_pairs <- as.matrix(expand.grid(lab_vec, lab_vec))
mtx_pairs <- mtx_pairs <- mtx_pairs[mtx_pairs[,1] > mtx_pairs[,2],]

# List for plots
ls_rd <- ls_plot <- list()

# For-Loop over all pairs
for(i_pair in 1:nrow(mtx_pairs)){
  # Get corresponding pair
  temp_pair <- mtx_pairs[i_pair, c(2, 1)]
  
  # Get relevant subset
  df_cp <- subset(df_total, is.element(lab, temp_pair))[,c("lab", paste0("p", temp_pair))]
  
  # If both probabilities are zero, both are equally likely -> 0.5
  df_cp[which(rowSums(df_cp[,paste0("p", temp_pair)]) == 0), paste0("p", temp_pair)] <- rep(1/2, 2)
  
  # Calculate conditional probabilities
  df_cp[,paste0("p", temp_pair)] <- 
    df_cp[,paste0("p", temp_pair)]/rowSums(df_cp[,paste0("p", temp_pair)])
  
  # Make reliability diagram of predictions
  ls_rd[[paste0(i_pair)]] <- reliabilitydiag(p = df_cp[[paste0("p", temp_pair[2])]],
                                             y = as.numeric(df_cp[["lab"]] == temp_pair[2]))
}

# For-Loop over all pairs
for(i_pair in 1:nrow(mtx_pairs)){
  # Get corresponding pair
  temp_pair <- mtx_pairs[i_pair, c(2, 1)]
  
  # Get names
  temp_name0 <- lab_names[paste0(temp_pair[1])]
  temp_name1 <- lab_names[paste0(temp_pair[2])]
  
  # Plot
  ls_plot[[paste0(i_pair)]] <- ggplot2::autoplot(
    ls_rd[[paste0(i_pair)]],
    params_histogram = list(breaks = (0:50)/50,
                            fill = NA,
                            colour = "black")) +
    ggtitle(paste0("(", letters[i_pair], ") 0 = ", temp_name0, ", 1 = ", temp_name1)) +
    theme(axis.text.x = element_text(size = 16),
          axis.text.y = element_text(size = 16),
          axis.title.x = element_text(size = 18),
          axis.title.y = element_text(size = 18),
          plot.title = element_text(face = "bold",
                                    size = 25)) +
    annotate("label", 
             label = paste0("Samples: ", 
                            format(sum(ls_rd[[paste0(i_pair)]][["p"]][["bins"]][["n"]]), big.mark = ","), 
                            "\n", "Hits: ", 
                            format(sum(ls_rd[[paste0(i_pair)]][["p"]][["cases"]][["y"]]), big.mark = ",")), 
             x = 0, 
             y = 0.92,
             hjust = 0,
             size = 6.5,
             label.padding = unit(0.55, "lines"))
}

# Plot together
pdf_plot <- marrangeGrob(grobs = ls_plot,
                         layout_matrix = matrix(1:10,
                                                ncol = 5,
                                                byrow = TRUE),
                         top = NULL)

# Name of PDF
file_pdf <- paste0(pdf_path, "cosmo_rea_reliability_all_pairs.pdf")

# Save PDF
ggplot2::ggsave(filename = file_pdf,
                plot = pdf_plot,
                width = 25,
                height = 10)