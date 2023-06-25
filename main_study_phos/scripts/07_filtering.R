###Background
#Code to filter phos data

###Initialise
library(tidyverse)
library(purrr)

#####Extremely low intensity values#####
###Function: target values to filter if they are below a threshold (threshold_low), and their difference from ppeptide mean is above a threshold (threshold_diff))
#threshold_diff must be positive
filter_low_and_belowbulk <- function(data,
                                     num_cols,
                                     threshold_low,
                                     threshold_diff){
  #low threshold
  threshold_low_m <- data[, num_cols] < threshold_low
  #diff threshold
  medians <- apply(data[, num_cols], 1, FUN = "median", na.rm = TRUE)
  threshold_diff_m <- sweep(data[, num_cols],
                            1,
                            STATS = medians,
                            FUN = function(x, y){y - x > threshold_diff})
  #Combine
  threshold_m <- threshold_low_m & threshold_diff_m
  threshold_m[is.na(threshold_m)] <- FALSE
  threshold_df <- data.frame(threshold_m)
  colnames(threshold_df) <- paste(colnames(threshold_df), "_filt", sep = "")
  threshold_df$filt <- rowSums(threshold_m, na.rm = TRUE) > 0
  return(threshold_df)
}


















































