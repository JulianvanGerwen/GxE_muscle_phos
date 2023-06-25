###Background
#Here I store code for processing blood glucose data during my insulin injections

###Initialise
library(tidyverse)
library(purrr)

###Function: Load blood glucose data
#directory: Where data is stored
load_bloodgluc_data <- function(directory){
  bloodgluc_data <- read_csv(directory) %>% as.data.frame()
  rownames(bloodgluc_data) <- bloodgluc_data$Mouse
  bloodgluc_data <- bloodgluc_data[, grep("min", colnames(bloodgluc_data))]
  colnames(bloodgluc_data) <- paste("bloodgluc", colnames(bloodgluc_data), sep = "_")
  #Add delta from baseline
  frombaseline_cols <- colnames(bloodgluc_data)
  bloodgluc_data[, gsub("bloodgluc", "bloodgluc_frombase", frombaseline_cols)] <- 
    sweep(bloodgluc_data[, frombaseline_cols],
          1,
          "-",
          STATS = bloodgluc_data[, grep("-1min", colnames(bloodgluc_data))])
  return(bloodgluc_data)
}



###Function: Extract slope from linear fit of blood glucose data
bloodgluc_slope <- function(bloodgluc_vec,
                            times){
  #Make dataframe
  model_df <- data.frame("bloodgluc" = bloodgluc_vec,
                         "time" = times)
  #Fit model
  model <- lm(bloodgluc ~ time, model_df)
  return(model$coefficients[2])
}

###Function: Add AOC and slope to bloodgluc data
process_bloodgluc_data <- function(bloodgluc_data,
                                   bloodgluc_times = c(1, 5, 7.5, 10)){
  #Set up data
  bloodgluc_cols <- paste("bloodgluc_", bloodgluc_times, "min", sep = "")
  bloodgluc_data_num <- bloodgluc_data[, bloodgluc_cols]
  
  #Calculate AOC
  bloodgluc_data$bloodgluc_AOC <- apply(bloodgluc_data_num,
                                        1,
                                        FUN = function(x) GTT_AUC_calculator(x,
                                                                             times = bloodgluc_times,
                                                                             AOC = TRUE))
  #Calculate slope
  bloodgluc_data$bloodgluc_slope <- apply(bloodgluc_data_num,
                                          1,
                                          FUN = function(x) bloodgluc_slope(x,
                                                                            times = bloodgluc_times))
  return(bloodgluc_data)
}



















