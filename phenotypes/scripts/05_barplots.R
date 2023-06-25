###Background
#Here I store barplot visualisations

###Initialise
library(ggplot2)
library(tidyverse)

###SEM
sem <- function(vec){
  return(sd(vec)/sqrt(length(vec)))
}


###Function: Make a barplot of strains dodged by diet. CHOW is hollow and HFD is full
#value_col: The column to plot
#yaxis_lab: Label for the yaxis
#padding: padding parameter for dodging
#errorbar_include: Adds SEM errorbars
#points_include: Adds data points
barplot_diets <- function(data,
                          value_col,
                          yaxis_lab = "Value",
                          padding = 0.2,
                          errorbar_include = FALSE,
                          points_include = TRUE){
  #Prep data
  #Add factors and rename columns
  data <- data %>%
    mutate(Diet_Strain = paste(Diet, Strain, sep = "_"))
  Diet_Strain_levels <- cross(list(levels(data$Strain), levels(data$Diet))) %>%
    map(~ paste(.[[2]], .[[1]], sep = "_")) %>%
    unlist()
  data <- mutate(data, Diet_Strain = factor(Diet_Strain, levels = Diet_Strain_levels))
  data$Value <- data[, value_col]
  #summarise data
  data_summarised <- data %>%
    group_by(Strain, Diet) %>%
    summarise(se = sem(Value), Value = mean(Value))
  
  #Prep colours
  Diet_Strain_cols <- c(rep("#ffffff", 5),
                        fivestraincols_main)
  names(Diet_Strain_cols) <- Diet_Strain_levels
  
  #Make error bar plot
  error_bar_plot <- geom_errorbar(data = data_summarised,
                                  aes(ymin = Value, ymax = Value + se,
                                      fill = Strain,
                                      group = Diet), 
                                  position = position_dodge2(preserve = "single", padding = padding), 
                                  size = 0.2)
  if (!errorbar_include){ error_bar_plot <- NULL }
  
  #Make point plot
  point_plot <- geom_point(size = 0.5,
                           stroke = 0.25,
                           position = position_dodge(
                             #Width controls how much to dodge. Empirically I found that 0.9 matches the default parameters for geom_bar with position_dodge2
                             width = 0.9))
  if (!points_include){ point_plot <- NULL}
  
  #Plot data
  output_plot <- ggplot(data,
                        aes(y = Value,
                            x = Strain,
                            fill = Diet_Strain,
                            colour = Strain,
                            group = Diet)) +
    geom_bar(stat = "summary", fun = "mean", 
             size = 0.5,
             alpha = 0.65,
             #use position_dodge2 to dodge columsn based on the aesthetic group
             #padding governs how close dodged columns are to each other
             position = position_dodge2(preserve = "single", padding = padding)) +
    error_bar_plot +
    point_plot +
    #Fill all chows as white
    scale_fill_manual(values = Diet_Strain_cols) +
    scale_colour_manual(values = fivestraincols_main) +
    comfy_theme(include_xaxis = FALSE,
                rotate_x_text = TRUE) +
    theme(axis.title.x = element_blank()) +
    labs(y = yaxis_lab)
  return(output_plot)
}

###Function: Save barplot diet plots with SEM or points
barplot_diets_SEMpoints <- function(data,
                                    file,
                                    width, height,
                                    ...){
  #SEM
  barplot_diets(data = data,
                errorbar_include = T, points_include = F,
                ...)
  ggsave_pdfpng(paste(file, "_sum", sep = ""),
                width = width, height = height)
  #points
  barplot_diets(data = data,
                errorbar_include = F, points_include = T,
                ...)
  ggsave_pdfpng(paste(file, "_alldat", sep = ""),
                width = width, height = height)
}



###Function: Prepare data for barplot_ins
#data: data that contains columns Strain, Diet, Ins, and a column specified by value_col
barplot_ins_prep <- function(data, value_col){
  #Prep data
  #Add factors and rename columns
  data <- data %>%
    mutate(Strain_Diet_Ins = paste(Strain, Diet, Ins, sep = "_"),
           Strain_Diet = paste(Strain, Diet, sep = " "))
  Strain_Diet_Ins_levels <- cross(list(levels(data$Ins), 
                                       levels(data$Diet),
                                       levels(data$Strain))) %>%
    map(~ paste(rev(unlist(.)), collapse = "_")) %>%
    unlist()
  Strain_Diet_levels <- cross(list(levels(data$Diet),
                                   levels(data$Strain))) %>%
    map(~ paste(rev(unlist(.)), collapse = " ")) %>%
    unlist()
  data <- mutate(data, Strain_Diet_Ins = factor(Strain_Diet_Ins, levels = Strain_Diet_Ins_levels),
                 Strain_Diet = factor(Strain_Diet, levels = Strain_Diet_levels))
  data$Value <- data[, value_col]
  
  #summarise data
  data_summarised <- data %>%
    group_by(Strain, Diet, Ins) %>%
    summarise(se = sem(Value), Value = mean(Value),
              #Retain categorical columns so I can add to ggplot object
              Strain_Diet = Strain_Diet[1])
  
  #Prep colours
  #White for all basal, coloured by strain otherwise
  Strain_Diet_Ins_cols <- sapply(Strain_Diet_Ins_levels,
                                 function(x){
                                   x_split <- strsplit(x, "_")[[1]]
                                   if (x_split[3] == "bas"){
                                     return("#FFFFFF")
                                   } else {
                                     return(fivestraincols_main[x_split[1]])
                                   }
                                 })
  names(Strain_Diet_Ins_cols) <- Strain_Diet_Ins_levels
  
  #Output as list
  output_list <- list("data" = data,
                      "data_summarised" = data_summarised,
                      "Strain_Diet_Ins_cols" = Strain_Diet_Ins_cols)
  return(output_list)
}

###Function: Make a barplot of strains and diets dodged by insulin. Bas is hollow and Ins is full
#data: data that contains columns Strain, Diet, Ins, and a column specified by value_col
#value_col: The column to plot
#yaxis_lab: Label for the yaxis
#padding: padding parameter for dodging
#errorbar_include: Adds SEM errorbars
#points_include: Adds data points
#facet_diet: Specifies whether diets are grouped within each strain (FALSE) or faceted horizontally into two plots (TRUE)
barplot_ins <- function(data,
                        value_col,
                        yaxis_lab = "Value",
                        padding = 0.2,
                        errorbar_include = FALSE,
                        points_include = TRUE,
                        facet_diet = FALSE){
  #Prep data
  data_prepped <- barplot_ins_prep(data = data, value_col = value_col)
  data <- data_prepped$data
  data_summarised <- data_prepped$data_summarised
  Strain_Diet_Ins_cols <- data_prepped$Strain_Diet_Ins_cols
  
  #Make error bar plot
  error_bar_plot <- geom_errorbar(data = data_summarised,
                                  aes(ymin = Value, ymax = Value + se,
                                      fill = Strain,
                                      group = Ins), 
                                  position = position_dodge2(preserve = "single", padding = padding), 
                                  size = 0.2)
  if (!errorbar_include){ error_bar_plot <- NULL }
  
  #Make point plot
  point_plot <- geom_point(size = 0.5,
                           stroke = 0.25,
                           position = position_dodge(
                             #Width controls how much to dodge. Empirically I found that 0.9 matches the default parameters for geom_bar with position_dodge2
                             width = 0.9))
  if (!points_include){ point_plot <- NULL}
  
  #Make changes based on ordering of plot
  if (facet_diet){
    facet <- facet_grid(cols = vars(Diet))
    x_aes <- "Strain"
  } else {
    facet <- NULL
    x_aes <- "Strain_Diet"
  }
  
  #Plot data
  output_plot <- ggplot(data,
                        aes_string(y = "Value",
                                   x = x_aes,
                                   fill = "Strain_Diet_Ins",
                                   colour = "Strain",
                                   group = "Ins")) +
    geom_bar(stat = "summary", fun = "mean", 
             size = 0.5,
             alpha = 0.65,
             #use position_dodge2 to dodge columsn based on the aesthetic group
             #padding governs how close dodged columns are to each other
             position = position_dodge2(preserve = "single", padding = padding)) +
    error_bar_plot +
    point_plot +
    #Fill all chows as white
    scale_fill_manual(values = Strain_Diet_Ins_cols) +
    scale_colour_manual(values = fivestraincols_main) +
    comfy_theme(include_xaxis = FALSE,
                rotate_x_text = TRUE) +
    theme(axis.title.x = element_blank()) +
    labs(y = yaxis_lab) +
    facet
  return(output_plot)
}

###Function: Save barplot diet plots with SEM or points
barplot_ins_SEMpoints <- function(data,
                                    file,
                                    width, height,
                                    ...){
  #SEM
  barplot_ins(data = data,
                errorbar_include = T, points_include = F,
                ...)
  ggsave_pdfpng(paste(file, "_sum", sep = ""),
                width = width, height = height)
  #points
  barplot_ins(data = data,
                errorbar_include = F, points_include = T,
                ...)
  ggsave_pdfpng(paste(file, "_alldat", sep = ""),
                width = width, height = height)
}



###Function: Prepare stat_summary to make significance labels to add to barplot_diets
#data: data to go into barplot_diets
#value_col: COlumn to be plotted in barplot
#stat_summary: Summary of t-test statistics as output by ttests_mult
#pcol: column of stat_summary to use for pvalues
#show_stars: If TRUE, show significance stars, if FALSE, show p values
#star_types: significance star types for diet and strain comparisons
#buffer: what proportion of the yaxis is added to stars so that they do not touch data points
barplot_diets_siglabs_prep <- function(data,
                                       value_col,
                                       stat_summary,
                                       pcol = "adj_p",
                                       show_stars = TRUE,
                                       star_types = c("diet" = "*", "strain" = "#"),
                                       buffer = 0.075){
  ##Getting test context
  #Extract strains and diets, which is used for plotting
  stat_summary[, c("Strain_1", "Diet_1", "Strain_2", "Diet_2")] <- strsplit(stat_summary$comparisons, "vs") %>%
    map(function(x){
      x_split <- strsplit(x, "_")
      x_trimmed <- map(x_split, ~.[1:2])
      return(unlist(x_trimmed))
    }) %>%
    purrr::reduce(rbind)
  #Assing Strain and Diet to the second, which is used for plotting
  stat_summary <- mutate(stat_summary, 
                         Strain = Strain_2, Diet = Diet_2,
                         #Make merged columns so that aesthetics inherited from the barplot are present here
                         Diet_Strain = paste(Diet_2, Strain_2, sep = "_"))
  #Figure out which are strain and diet comparisons
  stat_summary$comp_type <- "strain"
  stat_summary$comp_type[stat_summary$Strain_1 == stat_summary$Strain_2] <- "diet"
  stat_summary$dodging <- TRUE
  stat_summary$dodging[stat_summary$comp_type == "diet"] <- FALSE
  
  #Make significance labels
  #Get significance levels
  stat_summary$sig_levels <- (stat_summary[, pcol] < 0.05) + (stat_summary[, pcol] < 0.01) + (stat_summary[, pcol] < 0.001)
  #Loop over tests to make labels
  stat_summary$sig_label <- NA
  for (i in 1:nrow(stat_summary)){
    #sig labels
    if (stat_summary[i, pcol] < 0.05 & !is.na(stat_summary[i, pcol])){
      comp_type <- stat_summary$comp_type[i]
      stat_summary$sig_label[i] <- paste(rep(star_types[comp_type], stat_summary$sig_levels[i]),
                                         collapse = "")
    }
  }
  #Set sig labels to pvalues if desired
  if (!show_stars){stat_summary$sig_label <- stat_summary$p}
  
  #Get yvalues
  #Get max of abs values, and mean in each strain and diet
  stat_summary[, c("absmax_value_1", "mean_value_1", "absmax_value_2", "mean_value_2")] <-
    list(stat_summary[, c("Strain_1", "Diet_1")],
         stat_summary[, c("Strain_2", "Diet_2")]) %>%
    #Loop over rows
    map(~apply(.,
               1,
               function(x){
                 dat_tosummarise <- data[data$Strain == x[1] &
                                           data$Diet == x[2], value_col]
                 abs_max <- max(abs(dat_tosummarise), na.rm = TRUE)
                 mean <- mean(dat_tosummarise, na.rm = TRUE)
                 return(c(abs_max, mean))
               })) %>%
    purrr::reduce(rbind) %>%
    t
  #Set sign 
  stat_summary$y_sign <- sign(stat_summary$mean_value_1 + stat_summary$mean_value_2)
  #Pick y value depending on comparison type
  stat_summary$y <- apply(stat_summary[, c("absmax_value_1", "absmax_value_2")],
                          1, max)
  stat_summary$y[stat_summary$comp_type == "strain"] <- stat_summary$absmax_value_2[stat_summary$comp_type == "strain"]
  #Expand y values to avoid collisions. Here I assume the plot contains 0
  y_max <- max(data[, value_col], na.rm = TRUE)
  y_min <- min(c(0, min(data[, value_col], na.rm = TRUE)))
  y_height <- y_max - y_min
  #Add buffer to y heights, and add twice if diet comparison
  stat_summary$y <- stat_summary$y + buffer*y_height
  stat_summary$y[stat_summary$comp_type == "diet"] <- stat_summary$y[stat_summary$comp_type == "diet"] + buffer*y_height*0.75
  #Multiply by sign in case they need to be flipped
  stat_summary$y <- stat_summary$y*stat_summary$y_sign
  return(stat_summary)
}

###Function: Make significance labels from a pre-prepared summary of stats for use on dodged plots
#stat_summary: A processed summary of stats as output by barplot_diets_siglabs_prep, for example
#dodged_size: size of labels for labels that will be dodged
#undodged_size: size of labels for labels that will not be dodged
#x_aes, group_aes: names of aesthetics for x-axis and group (must be contained in stat_summary)
siglabels_dodged_make <- function(stat_summary,
                                  dodged_size = 4.5/3,
                                  undodged_size = 6/3,
                                  x_aes, group_aes){
  output_plots <- list(
    "dodged" = geom_text(data = stat_summary[stat_summary$dodging == TRUE, ],
                         aes_string(x = x_aes,
                                    group = group_aes,
                                    label = "sig_label",
                                    y = "y"),
                         colour = "#000000",
                         size = dodged_size,
                         position = position_dodge(width = 0.9)),
    "undodged" = geom_text(data = stat_summary[stat_summary$dodging == FALSE, ],
                           aes_string(x = x_aes,
                                      group = group_aes,
                                      label = "sig_label",
                                      y = "y"),
                           colour = "#000000",
                           size = undodged_size,
                           position = position_dodge(width = 0.9)),
    "expand_lims" = scale_y_continuous(expand = expansion(mult = c(0, 0.15)))
  )
  return(output_plots)
}

###Function: Wrapper for barplot_diets that adds siglabs on top
barplot_diets_siglabs <- function(data,
                                  value_col,
                                  stat_summary,
                                  ...){
  #Prep stat_summary
  stat_summary <- barplot_diets_siglabs_prep(data = data,
                                             value_col = value_col,
                                             stat_summary = stat_summary)
  #Plot
  barplot <- barplot_diets(data = data,
                           value_col = value_col,
                           ...)
  siglabs <- siglabels_dodged_make(stat_summary,
                                   x_aes = "Strain", group_aes = "Diet")
  output_plot <- barplot + 
    siglabs[[1]] + siglabs[[2]]
  return(output_plot)
}

###Function: Save barplot diet plots with SEM or points, and with or without sig labels
barplot_diets_SEMpoints_siglabs <- function(data,
                                            stat_summary,
                                    file,
                                    width, height,
                                    ...){
  #Make directory to deposit plots
  new_dir <- file
  dir.create(new_dir)
  filename <- rev(strsplit(file, "/")[[1]])[1]
  file <- paste(new_dir, "/", filename, sep = "")
  
  #SEM
  barplot_diets(data = data,
                errorbar_include = T, points_include = F,
                ...)
  ggsave_pdfpng(paste(file, "_sum", sep = ""),
                width = width, height = height)
  barplot_diets_siglabs(data = data,
                        stat_summary = stat_summary,
                errorbar_include = T, points_include = F,
                ...)
  ggsave_pdfpng(paste(file, "_sum_siglabs", sep = ""),
                width = width, height = height)
  #points
  barplot_diets(data = data,
                errorbar_include = F, points_include = T,
                ...)
  ggsave_pdfpng(paste(file, "_alldat", sep = ""),
                width = width, height = height)
  barplot_diets_siglabs(data = data,
                        stat_summary = stat_summary,
                        errorbar_include = F, points_include = T,
                        ...)
  ggsave_pdfpng(paste(file, "_alldat_siglabs", sep = ""),
                width = width, height = height)
}


###Function: Prepare stat_summary to make significance labels to add to barplot_ins
#data: data to go into barplot_diets
#value_col: COlumn to be plotted in barplot
#stat_summary: Summary of t-test statistics as output by ttests_mult
#pcol: column of stat_summary to use for pvalues
#show_stars: If TRUE, show significance stars, if FALSE, show p values
#star_types: significance star types for diet and strain comparisons
#buffer: what proportion of the yaxis is added to stars so that they do not touch data points
#fill_basorins_stats: If only bas or only ins stats are used, this fills up empty tests with the other condition to ensure dodging works
barplot_ins_siglabs_prep <- function(data,
                                     value_col,
                                     stat_summary,
                                     pcol = "adj_p",
                                     show_stars = TRUE,
                                     star_types = c("diet" = "*", "strain" = "#"),
                                     buffer = 0.075,
                                     fill_basorins_stats = TRUE){
  #Fill in empty conditions from bas tests or ins tests so that dodging works
  if (fill_basorins_stats){
    Ins_1 <- sapply(stat_summary$comparisons, function(x) strsplit(x, "vs|_")[[1]][3])
    rownames(stat_summary) <- stat_summary$comparisons
    #Find if the tests are bas or ins
    if ("ins" %in% Ins_1){
      stat_summary[gsub("ins", "bas", rownames(stat_summary)), ] <- NA
    } else if ("bas" %in% Ins_1){
      stat_summary[gsub("bas", "ins", rownames(stat_summary)), ] <- NA
    }
    stat_summary$comparisons <- rownames(stat_summary)
  }
  ##Getting test context
  #Extract strains and diets, which is used for plotting
  stat_summary[, c("Strain_1", "Diet_1", "Ins_1", "Strain_2", "Diet_2", "Ins_2")] <- strsplit(stat_summary$comparisons,
                                                                                              "vs|_") %>%
    purrr::reduce(rbind)
  #Assing Strain, Diet, and Ins to the second, which is used for plotting
  stat_summary <- mutate(stat_summary, 
                         Strain = Strain_2, Diet = Diet_2, Ins = Ins_2,
                         #Make merged columns so that aesthetics inherited from the barplot are present here
                         Strain_Diet = paste(Strain_2, Diet_2, sep = " "),
                         Strain_Diet_Ins = paste(Strain_2, Diet_2, Ins_2, sep = "_"))
  #Figure out which are strain, diet, and ins comparisons
  stat_summary$comp_type <- "strain"
  stat_summary$comp_type[stat_summary$Strain_1 == stat_summary$Strain_2 &
                           stat_summary$Ins_1 == stat_summary$Ins_2] <- "diet"
  stat_summary$comp_type[stat_summary$Strain_1 == stat_summary$Strain_2 &
                           stat_summary$Diet_1 == stat_summary$Diet_2] <- "ins"
  #Only ins comparisons are not dodged
  stat_summary$dodging <- TRUE
  stat_summary$dodging[stat_summary$comp_type == "ins"] <- FALSE
  
  #Make significance labels
  #Get significance levels
  stat_summary$sig_levels <- (stat_summary[, pcol] < 0.05) + (stat_summary[, pcol] < 0.01) + (stat_summary[, pcol] < 0.001)
  #Loop over tests to make labels
  stat_summary$sig_label <- NA
  for (i in 1:nrow(stat_summary)){
    #sig labels
    if (stat_summary[i, pcol] < 0.05 & !is.na(stat_summary[i, pcol])){
      comp_type <- stat_summary$comp_type[i]
      stat_summary$sig_label[i] <- paste(rep(star_types[comp_type], stat_summary$sig_levels[i]),
                                         collapse = "")
    }
  }
  #Set sig labels to pvalues if desired
  if (!show_stars){stat_summary$sig_label <- stat_summary$p}
  
  
  #Get yvalues
  #Get max in each strain and diet
  #Get yvalues
  #Get max of abs values, and mean in each strain and diet
  stat_summary[, c("absmax_value_1", "mean_value_1", "absmax_value_2", "mean_value_2")] <-
    list(stat_summary[, c("Strain_1", "Diet_1", "Ins_1")],
         stat_summary[, c("Strain_2", "Diet_2", "Ins_2")]) %>%
    #Loop over rows
    map(~apply(.,
               1,
               function(x){
                 dat_tosummarise <- data[data$Strain == x[1] &
                                           data$Diet == x[2] &
                                           data$Ins == x[3], value_col]
                 abs_max <- max(abs(dat_tosummarise), na.rm = TRUE)
                 mean <- mean(dat_tosummarise, na.rm = TRUE)
                 return(c(abs_max, mean))
               })) %>%
    purrr::reduce(rbind) %>%
    t
  #Set sign 
  stat_summary$y_sign <- sign(stat_summary$mean_value_1 + stat_summary$mean_value_2)
  #Pick y value depending on comparison type
  #If ins comparsons (not dodged), take the max of the two conditions
  stat_summary$y <- apply(stat_summary[, c("absmax_value_1", "absmax_value_2")],
                          1, max)
  stat_summary$y[stat_summary$comp_type %in% c("strain", "diet")] <- 
    stat_summary$absmax_value_2[stat_summary$comp_type %in% c("strain", "diet")]
  #Expand y values to avoid collisions. Here I assume the plot contains 0 on the yaxis
  y_max <- max(data[, value_col], na.rm = TRUE)
  y_min <- min(c(0, min(data[, value_col], na.rm = TRUE)))
  y_height <- y_max - y_min
  #Add buffer to y heights, and add twice if diet comparison
  stat_summary$y <- stat_summary$y + buffer*y_height
  stat_summary$y[stat_summary$comp_type == "diet"] <- stat_summary$y[stat_summary$comp_type == "diet"] + buffer*y_height*0.75
  #Multiply by sign
  stat_summary$y <- stat_summary$y*stat_summary$y_sign
  return(stat_summary)
}

###Function: Wrapper for barplot_diets that adds siglabs on top
barplot_ins_siglabs <- function(data,
                                value_col,
                                stat_summary,
                                ...){
  #Prep stat_summary
  stat_summary <- barplot_ins_siglabs_prep(data = data,
                                           value_col = value_col,
                                           stat_summary = stat_summary)
  #Plot
  barplot <- barplot_ins(data = data,
                         value_col = value_col,
                         ...)
  siglabs <- siglabels_dodged_make(stat_summary,
                                   x_aes = "Strain_Diet", group_aes = "Ins")
  output_plot <- barplot + 
    siglabs[[1]] + siglabs[[2]]
  return(output_plot)
}







