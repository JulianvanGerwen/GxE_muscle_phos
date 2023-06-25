###background
#Functions for processing blood count data

###Initialise
library(tidyverse)


###Function: Turn long count data into wide
#data_long: count data in long form. Columns are Mouse, Time, DPM
#volcounted: The volume of blood counted in uL. 5 by default
blood_data_towide <- function(data_long, volcounted = 5*0.7){
  data_long <- data_long %>%
    #Format Time
    mutate(Time = paste(Time, "min", sep = "")) %>%
    #Normalise DPMs
    mutate(DPM = DPM/volcounted)
  #Pivot separately for DPM, Timepoint, and Include
  pivot_cols <- c("DPM", "Timepoint", "Include")
  data_wide_list <- map(pivot_cols,
                        function(x){
                          data_long[, c("Mouse", "Time", x)] %>%
                            pivot_wider(names_from = Time, values_from = x)
                        })
  names(data_wide_list) <- pivot_cols
  return(data_wide_list)
}


###Function: single linear fit of log count to time
#time: Numeric vector of time points
#logcount: numeric log count data
single_logcount_fit <- function(time, logcount){
  fit_data <- data.frame("time" = as.numeric(time),
                         "log_count" = as.numeric(logcount))
  #Fit linear model and perform correlation
  fit <- lm(log_count ~ time, data = fit_data)
  corr <- cor.test(fit_data$time, fit_data$log_count)
  output_vec <- c("intercept" = coefficients(fit)[1],
                  "slope" = coefficients(fit)[2],
                  "Rsq" = corr$estimate^2)
  names(output_vec) <- c("intercept", "slope", "Rsq")
  return(output_vec)
}

###Function: linear fits of multiple log count data to time
#data_wide_list: List of data frames as output by blood_data_towide, include a dataframe DPM_log of log counts
mult_logcount_fit <- function(data_wide_list){
  #Prep data
  data_wide_list_forfit <- data_wide_list[c("DPM_log", "Timepoint", "Include")] %>%
    map(~ t(select(., -Mouse))) 
  #Loop over rows
  fit_data <- list()
  for (i in 1:ncol(data_wide_list_forfit$DPM_log)){
    to_use <- data_wide_list_forfit$Include[, i]
    temp_fit <- single_logcount_fit(data_wide_list_forfit$Timepoint[to_use, i],
                                    data_wide_list_forfit$DPM_log[to_use, i])
    fit_data[[i]] <- temp_fit
  }
  fit_data <- fit_data %>%
    purrr::reduce(rbind) %>%
    as.data.frame() %>%
    mutate(Mouse = data_wide_list$DPM_log$Mouse)
  #Add rate and initial
  fit_data$rate <- -fit_data$slope
  fit_data$initial <- exp(fit_data$intercept)
  return(fit_data)
}


###Function: log wide blood data and fit
#data_wide: wide data as outputed by blood_data_towide
blood_logandfit <- function(data_wide_list){
  #Log data
  time_cols <- colnames(data_wide_list$DPM)[grep("min$", colnames(data_wide_list$DPM))]
  data_wide_list$DPM_log <- data_wide_list$DPM 
  data_wide_list$DPM_log[, time_cols] <- log(data_wide_list$DPM_log[, time_cols])
  data_wide_list$DPM_log <- rename_with(data_wide_list$DPM_log, ~ gsub("min", "min_log", .))
  #Fit
  data_wide_list$fit_data <- mult_logcount_fit(data_wide_list)
  #Combine into one df
  #Rename Time and Include
  data_wide_list$Timepoint <- rename_with(data_wide_list$Timepoint, ~gsub("min", "min_timepoint", .))
  data_wide_list$Include <- rename_with(data_wide_list$Include, ~gsub("min", "min_include", .))
  #Condense
  mouse_IDs <- data_wide_list$DPM$Mouse[order(data_wide_list$DPM$Mouse)]
  data_wide_condensed <- map(data_wide_list,
                             function(data){
                               data <- data[order(data$Mouse), ] %>%
                                 select(-Mouse)
                             }) %>%
    purrr::reduce(cbind) %>%
    mutate(Mouse = mouse_IDs)
  return(data_wide_condensed)
}

###Function: Take long blood data, put it to wide, log it, and fit
#data_long: count data in long form. Columns are Mouse, Time, DPM
#volcounted: The volume of blood counted in uL. 5 by default
blood_raw_to_logandfit <- function(data_long,
                                   volcounted = 5*0.7,
                                   rename = TRUE,
                                   ...){
  #Wide data
  data_wide_list <- blood_data_towide(data_long = data_long, volcounted = volcounted, ...)
  #Fit and get rate
  data_wide_fit <- blood_logandfit(data_wide_list)
  #Add rownames and make colnames unique if desired
  if (rename){
    data_wide_fit <- as.data.frame(data_wide_fit)
    rownames(data_wide_fit) <- data_wide_fit$Mouse
    data_wide_fit <- select(data_wide_fit, -Mouse)
    colnames(data_wide_fit) <- paste("bloodcounts", colnames(data_wide_fit), sep = "_")
  }
  return(data_wide_fit)
}



####Below is code for filtering blood counts 

###Function: Modift the Include column on raw bloodcount data based on a filtering function that uses DPM column
#data: the data to use
#filt_func: A function that takes in a vector of counts, ordered by 1min, 5min, 7.5min, 10min, and returns a boolean vector of timepoints to include
filter_bloodcount_data <- function(data, filt_func){
  #Order data
  data <- data %>% .[order(.$Mouse, .$Time), ]
  for (mouse in unique(data$Mouse)){
    #Run filter
    include_temp <- filt_func(data$DPM[data$Mouse == mouse])
    #Match back into data
    data$Include[data$Mouse == mouse] <- include_temp
  }
  return(data)
}

#Filt B: Remove first timepoint if it is smaller than second by a great enough margin
filtB_func <- function(counts){
  if (counts[2] > counts[1] & (counts[2] - counts[1]) > 0.5*(counts[2] - counts[3])){
    return(c(F, T, T, T))
  } else {
    return(rep(T, 4))
  }
}













