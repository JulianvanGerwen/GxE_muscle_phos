###Background
#Here I store global functions that are useful across R projects e.g. in pheno data and phospho data

###Initialise
library(tidyverse)
library(purrr)
library(ggplot2)


#####Data processing#####
###Function: Take a vector of mouse names (e.g. NOD_CHOW_ins_1) or conditions (e.g. NOD_CHOW_ins) and turn them into a data frame with Strain, Diet, Ins, and all combinations as factors
extract_factors_from_vector <- function(vector){
  #Get list of strain, diet, and ins
  StrainDietIns_list <- strsplit(vector, "_") %>%
    map(function(x){
      x <- x[1:3]
      names(x) <- c("Strain", "Diet", "Ins")
      return(x)
    })
  factors <- c("Strain", "Diet", "Ins")
  #Produce combinations
  combinations <- map(1:3, ~combn(factors, ., simplify = FALSE)) %>%
    purrr::reduce(c)
  #Loop over combinations to fuse
  level_combinations <- map(combinations,
                            function(x){
                              map(StrainDietIns_list, ~paste(.[x], collapse = "_")) %>% unlist
                            })
  names(level_combinations) <- map(combinations, ~paste(., collapse = ""))
  #Turn into data frame
  level_combinations_df <- purrr::reduce(level_combinations, cbind) %>%
    as.data.frame()
  colnames(level_combinations_df) <- names(level_combinations)
  #Make each column a factor
  for (col in colnames(level_combinations_df)){
    level_combinations_df[, col] <- factor(level_combinations_df[, col],
                                           levels = all_levels[[col]])
  }
  return(level_combinations_df)
}

#Function: Bestow data where rownames are mice with Strain, Diet, and Ins cols
extract_factors_from_rownames <- function(data){
  factor_data <- extract_factors_from_vector(rownames(data))
  #Remove existing colnames
  factor_data <- factor_data[, which(colnames(factor_data) %in%
                                       colnames(data) == FALSE)]
  data <- cbind(data, factor_data)
  return(data)
}


#####Visualisation#####

###Function: Standard error of the mean
sem <- function(x){return(sd(x, na.rm = TRUE)/sqrt(length(x)))}


###Function: Summarise data over StrainDietIns levels. Generally used to prepare a corr_plot
#xaxis and yaxis: Columns in the data to summarise
#data: Data containign xaxis and yaxis, as well as Strain, Diet, Ins, and all of their combinations
summarise_for_corrplot <- function(data, xaxis, yaxis,
                                   summary_func = function(x){mean(x, na.rm = TRUE)},
                                   error_func = sem){
  #Rename
  colnames(data)[colnames(data) == yaxis] <- "yaxis"
  colnames(data)[colnames(data) == xaxis] <- "xaxis"
  #Summarise
  data <- group_by(data, StrainDietIns) %>%
    summarise(Strain = Strain[1], Diet = Diet[1], Ins = Ins[1],
              StrainDiet = StrainDiet[1], StrainIns = StrainIns[1], DietIns = DietIns[1],
              yaxis_sum = summary_func(yaxis), xaxis_sum = summary_func(xaxis),
              yaxis_error = error_func(yaxis), xaxis_error = error_func(xaxis))
  #Rename back
  colnames(data)[colnames(data) == "yaxis_sum"] <- yaxis
  colnames(data)[colnames(data) == "xaxis_sum"] <- xaxis
  return(data)
}

#Function: Point plot for correlation of two columns across Strain, Diet, and Ins
#Colours by Strain, shape for Diet (CHOW is circle, HFD is square), bas points have transparent inside
corrplot_general <- function(data, include_legend = T,
                          yaxis, xaxis,
                          ylab = NULL, xlab = NULL,
                          point_size = 1, alpha = 0.8,
                          summarise = FALSE, errorbars = FALSE,
                          trendline = FALSE, trendline_data = "all",
                          include_basins = TRUE,
                          ...){
  #Make axis labels if specific
  if (is.null(xlab)){
    xlab <- xaxis
  }
  if (is.null(ylab)){
    ylab <- yaxis
  }
  #relabel
  colnames(data)[colnames(data) == xaxis] <- "xaxis"
  colnames(data)[colnames(data) == yaxis] <- "yaxis"
  #Summarise if desired
  errorbar_plot1 <- NULL
  errorbar_plot2 <- NULL
  unsumm_data <- data
  if (summarise){
    data <- summarise_for_corrplot(data, xaxis = "xaxis", yaxis = "yaxis", ...)
    if (errorbars){
      errorbar_plot1 <- geom_errorbar(data = data,
                                      aes(y = yaxis, x = xaxis, 
                                          ymin = yaxis - yaxis_error, ymax = yaxis + yaxis_error,
                                          colour = Strain),
                                      size = 0.15, width = 0, alpha = alpha)
      errorbar_plot2 <- geom_errorbarh(data = data,
                                       aes(y = yaxis, 
                                           xmin = xaxis - xaxis_error, xmax = xaxis + xaxis_error,
                                           colour = Strain),
                                       size = 0.15, height = 0, alpha = alpha)
    }
  }
  #Trendline if desired
  if (trendline){
    #Set up data
    if (trendline_data == "all"){
      trend_data <- unsumm_data
    } else {
      trend_data <- data
    }
    #Make plot
    trendline_plot <- geom_smooth(data = trend_data, aes(x = xaxis, y = yaxis),
                                  method = lm, colour = "black", se = TRUE,
                                  size = 0.235,
                                  alpha = 0.15)
  } else {
    trendline_plot <- NULL
  }
  
  #Set up geom_points
  #If bas and ins data are included, make bas alpha'd (hollow)
  if (include_basins){
    #Ins geom
    point_plot_1 <- geom_point(data = data[data$Ins == "ins", ],
                               aes(x = xaxis, y= yaxis, color = Strain, fill = Strain, shape = DietIns),
                               stroke = 0.2,
                               size = point_size,
                               alpha = alpha)
    #Bas geom
    point_plot_2 <- geom_point(data = data[data$Ins == "bas", ],
                               aes(x = xaxis, y= yaxis, color = Strain, fill = Strain, shape = DietIns),
                               stroke = 0.2,
                               size = point_size)
    #Scale shape
    scale_shape <- scale_shape_manual(values = c("CHOW_bas" = 21, 
                                                 "CHOW_ins" = 16, 
                                                 "HFD_bas" = 22, 
                                                 "HFD_ins" = 15))
  } else {
    #HFD geom
    point_plot_1 <- geom_point(data = data[data$Diet == "HFD", ],
                               aes(x = xaxis, y= yaxis, color = Strain, fill = Strain, shape = Diet),
                               stroke = 0.2,
                               size = point_size,
                               alpha = alpha)
    #CHOW geom
    point_plot_2 <- geom_point(data = data[data$Diet == "CHOW", ],
                               aes(x = xaxis, y= yaxis, color = Strain, fill = Strain, shape = Diet),
                               stroke = 0.2,
                               size = point_size)
    #Scale shape
    scale_shape <- scale_shape_manual(values = c("CHOW" = 21, 
                                                 "HFD" = 16))
  }
  
  #Plot
  output_plot <- ggplot() +
    #Trend line 
    trendline_plot +
    #errorbar 
    errorbar_plot1 + errorbar_plot2 +
    #INS geom
    point_plot_1 +
    #BAS geom
    point_plot_2 +
    #Scale aesthetics
    scale_shape +
    scale_colour_manual(values = fivestraincols_main) +
    scale_fill_manual(values = alpha(fivestraincols_main, 0.2)) +
    comfy_theme(include_legend = include_legend)  +
    labs(x = xlab, y = ylab)
  return(output_plot)
}

#####Stats#####

###Function: Convert p values to stars
ps_to_stars <- function(p, star = "*",
                        sig_thresholds = c(0.05, 0.01, 0.001)){
  p_num <- sum(sig_thresholds > p, na.rm = TRUE)
  p_stars <- rep(star, p_num) %>% paste(collapse = "")
  return(p_stars)
}

ps_to_stars_vec <- function(p, ...){
  return(sapply(p, function(x){ ps_to_stars(p = x, ...) }))
}

ps_to_stars_df <- function(p, ...){
  stars_df <- p
  for (col in colnames(p)){
    stars_df[, col] <- ps_to_stars_vec(p[, col], ...)
  }
  return(stars_df)
}




