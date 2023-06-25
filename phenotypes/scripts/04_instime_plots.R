###Background
#Functions for plotting insulin time courses e.g. blood glucose, blood tracer

###Data format
#Long data needs to have the columns: Time, Value (become yaxis), Mouse, Ins

###Initialise
library(ggplot2)
library(tidyverse)

###Function: Turn wide data into long data
#data: Phenotypic data in wide format
#numcols: Colnames or indices for numerical data to be sent to long. These colnames must contain time followed by min e.g. 10min
instime_widetolong <- function(data,
                               numcols,
                               time_unit = "min"){
  data_long <- pivot_longer(data,
                            cols = colnames(data[, numcols]),
                            names_to = "Time",
                            values_to = "Value") %>%
    mutate(Time = sapply(Time, function(x){
      x_split <- strsplit(x, "_")[[1]]
      x_time <- x_split[grep(time_unit, x_split)]
      x_time <- gsub(time_unit, "", x_time)
    })) %>%
    mutate(Time = as.numeric(Time))
  return(data_long)
}



###Funtion: Line plot
#data: Data in long format, as described above
#yaxis_name: Name for yaxis in the plot
#colour_Ins: Boolean, indicates whether to colour bas (grey) and ins(black)
#facet_Ins: Boolean, indicates whether to use faceting to plot bas and ins values in two plots side-by-side
#yaxis_includes_0: Boolean, indicates whether yaxis must contain 0
#line_size: Width of the line
#xaxis_breaks: Where to place ticks on the x axis
instime_lineplot <- function(data,
                             yaxis_name = "DPM",
                             alpha_Ins = TRUE,
                             diet_shape = FALSE,
                             diet_shape_vals = c(16, 15),
                             yaxis_includes_0 = TRUE,
                             line_size = 0.4,
                             xaxis_breaks = c(0, 5, 10),
                             include_geom_point = TRUE,
                             include_legend = TRUE,
                             point_size = medium_pointsize){
  
  #Prep colour
  used_levels <- levels(data$Strain)[which(levels(data$Strain) %in% data$Strain)]
  used_colours <- fivestraincols_main[used_levels]
  scale_colour <- scale_colour_manual(values = used_colours)
  
  #Prep alpha
  if (alpha_Ins){
    scale_alpha <- scale_alpha_manual(values = c(0.5, 1))
    alpha_col <- "Ins"
  } else {
    scale_alpha <- scale_alpha_manual(values = c(1, 1))
    alpha_col <- "Diet"
  }
  
  #Prep shape
  if (diet_shape){
    scale_shape <- scale_shape_manual(values = diet_shape_vals)
  } else {
    scale_shape <- scale_shape_manual(values = c(16, 16))
  }
  
  #Prep axis limit
  if (yaxis_includes_0){
    scale_y <- expand_limits(y = 0)
  } else {
    scale_y <- NULL
  }
  
  #Prep xaxis_breaks
  if (is.null(xaxis_breaks)){
    scale_x <- NULL
  } else {
    scale_x <- scale_x_continuous(breaks = xaxis_breaks,
                                  labels = xaxis_breaks)
  }
  
  #Include geom_point
  if (include_geom_point){
    point_plot <- geom_point(fill = "#ffffff",
                             stroke = line_size*1.488,
                             size = point_size)
  } else {
    point_plot <- NULL
  }
  
  #Set up legend
  if (include_legend){
    legend_theme <- theme(legend.text = element_text(colour = "black",
                                                     size = 7),
                          legend.title = element_text(colour = "black",
                                                      size = 7),
                          legend.key = element_blank(),
                          legend.key.size = unit(0.25, 'cm'))
  } else {
    legend_theme <- theme(legend.position = "none")
  }
  
  output_plot <- ggplot(data = data,
                        aes_string(x = "Time", y = "Value", group = "Mouse", colour = "Strain", 
                                   alpha = alpha_col, shape = "Diet")) +
    geom_line(size = line_size) +
    point_plot +
    theme(panel.background = element_blank(),
          panel.grid = element_blank(),
          axis.line.x = element_line(colour = "black", size = 0.235),
          axis.line.y = element_line(colour = "black", size = 0.235),
          axis.ticks.x = element_line(size = 0.235, colour = "black"),
          axis.ticks.y = element_line(size = 0.235, colour = "black"),
          axis.text.x = element_text(colour = "black", size = 7),
          axis.text.y = element_text(colour = "black", size = 7),
          axis.title.x = element_text(colour = "black", size = 7),
          axis.title.y = element_text(colour = "black", size = 7),
          plot.title = element_text(colour = "black", size = 7, hjust = 0.5),
          strip.text = element_text(colour = "black", size = 7),
          strip.background = element_blank()) +
    legend_theme +
    scale_alpha +
    scale_colour +
    scale_shape +
    scale_y +
    scale_x +
    labs(y = yaxis_name, x = "Time (min)")
  return(output_plot)
}



###Function: Line plot of either blood glucose or blood tracer counts
#pheno_data: Phenotype data in standard form (produced by combine_phenotype_data)
#measure: The phenotype to display. 
##"bloodgluc": Blood glucose
##"bloodcounts": Blood tracer
#summarise: Whether or not to plot mean +- SEM
#xaxis_breaks_manual: overrides xaxis_breaks generated by default for each measure
instime_lineplot_bloodmeasure <- function(pheno_data,
                                          measure = "bloodgluc",
                                          summarise = FALSE,
                                          xaxis_breaks_manual = NULL,
                                          line_size = 0.4,
                                          point_size = medium_pointsize,
                                          errorbar_widener = 0.1,
                                          ...){
  #Prep data
  numcols <- intersect(grep(paste(measure, "_-*\\d+", sep = ""), colnames(pheno_data)),
                       grep("min$", colnames(pheno_data)))
  pheno_data <- pheno_data[, c("Mouse", "Strain", "Diet", "Ins",
                               colnames(pheno_data)[numcols])]
  pheno_data_long <- instime_widetolong(pheno_data,
                                        grep(measure, colnames(pheno_data)))
  #Configure based on measure (mainly accounts for GTT)
  if (measure == "postdiet_GTTgluc" |
      measure == "baseline_GTTgluc"){
    xaxis_breaks <- c(0, 30, 60, 90)
    include_ins <- FALSE
    alpha_ins <- FALSE
  } else if (measure == "postdiet_GTTins" |
             measure == "baseline_GTTins"){
    xaxis_breaks <- c(0, 15)
    include_ins <- FALSE
    alpha_ins <- FALSE
  } else {
    xaxis_breaks <- c(0, 5, 10)
    include_ins <- TRUE
    alpha_ins <- TRUE
    errorbar_widener <- 0.15
  }
  #Override xaxis_breaks if desires
  if (!is.null(xaxis_breaks_manual)){
    xaxis_breaks <- xaxis_breaks_manual
    scale_xaxis <- expand_limits(x = c(min(xaxis_breaks_manual), max(xaxis_breaks_manual)))
  } else{
    scale_xaxis <- NULL
  }
  
  #Summarise if necessary
  if (summarise){
    pheno_data_long <- summarise_long_data(pheno_data_long,
                                           include_ins = include_ins)
    #Calculate width based on xaxis_breaks
    errorbar_width <- errorbar_widener*(max(xaxis_breaks) - min(xaxis_breaks))
    errorbar_plot <- geom_errorbar(aes(ymin = Value - se, ymax = Value + se), 
                                   #width is width of the horizontal error bar
                                   width = errorbar_width, size = 0.15)
  } else{
    errorbar_plot <- NULL
  }
  
  #Prep naming
  if (measure == "bloodgluc" |
      measure == "postdiet_GTTgluc" |
      measure == "baseline_GTTgluc"){
    yaxis_name <- "Blood glucose (mM)"
  } else if (measure == "bloodgluc_frombase"){
    yaxis_name <- "Blood glucose change from baseline (mM)"
  } else if (measure == "bloodcounts"){
    yaxis_name <- "Blood 2DG (DPM/uL)"
  } else if (measure == "postdiet_GTTins" |
             measure == "baseline_GTTins"){
    yaxis_name <- "Blood insulin (ng/mL)"
  } else {
    print("Incorrect value of measure")
    return(NULL)
  }
  #Plot
  output_plot <- instime_lineplot(pheno_data_long,
                                  yaxis_name = yaxis_name,
                                  xaxis_breaks = xaxis_breaks,
                                  alpha_Ins = alpha_ins,
                                  include_geom_point = FALSE,
                                  line_size = line_size,
                                  point_size = point_size,
                                  ...) +
    errorbar_plot +
    geom_point(fill = "#ffffff",
               stroke = line_size*1.488,
               size = point_size) +
    scale_xaxis
  return(output_plot)
}


###Function: Line plot of body weight
#pheno_data: Phenotype data in standard form (produced by combine_phenotype_data)
#summarise: Whether or not to plot mean +- SEM
#xaxis_breaks_manual: overrides xaxis_breaks generated by default for each measure
lineplot_bodyweight <- function(pheno_data,
                                summarise = FALSE,
                                xaxis_breaks_manual = NULL,
                                line_size = 0.4,
                                point_size = medium_pointsize,
                                errorbar_widener = 0.1,
                                ...){
  #Prep data
  numcols <- grep("bodyweight.+week$", colnames(pheno_data))
  pheno_data <- pheno_data[, c("Mouse", "Strain", "Diet", "Ins",
                               colnames(pheno_data)[numcols])]
  pheno_data_long <- instime_widetolong(pheno_data,
                                        numcols = grep("bodyweight.+week$", colnames(pheno_data)),
                                        time_unit = "week")
  #Set xaxis_breaks
  xaxis_breaks <- c(0, 2, 4, 6)
  #Override xaxis_breaks if desires
  if (!is.null(xaxis_breaks_manual)){
    xaxis_breaks <- xaxis_breaks_manual
    scale_xaxis <- expand_limits(x = c(min(xaxis_breaks_manual), max(xaxis_breaks_manual)))
  } else{
    scale_xaxis <- NULL
  }
  
  #Summarise if necessary
  if (summarise){
    pheno_data_long <- summarise_long_data(pheno_data_long,
                                           include_ins = FALSE)
    #Calculate width based on xaxis_breaks
    errorbar_width <- errorbar_widener*(max(xaxis_breaks) - min(xaxis_breaks))
    errorbar_plot <- geom_errorbar(aes(ymin = Value - se, ymax = Value + se), 
                                   #width is width of the horizontal error bar
                                   width = errorbar_width, size = 0.25)
  } else{
    errorbar_plot <- NULL
  }
  #Plot
  output_plot <- instime_lineplot(pheno_data_long,
                                  yaxis_name = "Bodyweight (g)",
                                  xaxis_breaks = xaxis_breaks,
                                  alpha_Ins = FALSE,
                                  include_geom_point = FALSE,
                                  line_size = line_size,
                                  yaxis_includes_0 = FALSE,
                                  point_size = point_size,
                                  ...) +
    errorbar_plot +
    geom_point(fill = "#ffffff",
               stroke = line_size*1.488,
               size = point_size) +
    scale_xaxis +
    labs(x = "Weeks on diet")
  return(output_plot)
}



###Function: Standard error of the mean
sem <- function(vec){
  return(sd(vec)/sqrt(length(vec)))
}


###Function: Extract mean and SEM from long data
#data: long data as output by instime_widetolong
summarise_long_data <- function(data,
                                include_ins = TRUE){
  #Summarise based on whether or not to include Ins
  if (include_ins){
    data_summarised <- data %>%
      group_by(Strain, Diet, Ins, Time)
  } else {
    data_summarised <- data %>%
      group_by(Strain, Diet, Time)
  }
  data_summarised <- data_summarised %>%
    summarise(se = sem(Value), Value = mean(Value), Mouse = Mouse[1])
  return(data_summarised)
}











