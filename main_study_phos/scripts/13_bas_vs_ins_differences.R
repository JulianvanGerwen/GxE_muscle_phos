###Background
#Here I store code to help pinning down differences to bas and ins changes

###Initialise
library(tidyverse)
library(purrr)
library(ggplot2)

#####Filtering for differences#####
###Function: Filters for basins changes
#data: phos data to filter
#enhsupp: subset of data that is enhsupp columns
#bas_FC_data: subset of data that is bas FC columns
#strains: The strains tested. Order must match that of enhsupp and FC_data
#FC_threshold: Threshold of FCs to filter
basins_change_filterer <- function(data,
                                   enhsupp,
                                   bas_FC_data, ins_FC_data,
                                   strains,
                                   FC_threshold = 0.58){
  ###Build mask to transpose data so all comparisons are the same
  ##Get insreg data
  insreg_dir <- data$insreg_dir
  ##Make mask
  #This inverts FC data so that everything is like up-regulated supp
  mask <- map(1:length(insreg_dir), function(i){
    enhsupp_row <- as.character(enhsupp[i, ])
    if (insreg_dir[i] == "up"){
      enhsupp_row[which(enhsupp_row == "suppressed")] <- 1
      enhsupp_row[which(enhsupp_row == "enhanced")] <- -1
    } else if (insreg_dir[i] == "down"){
      enhsupp_row[which(enhsupp_row == "suppressed")] <- -1
      enhsupp_row[which(enhsupp_row == "enhanced")] <- 1
    }
    return(as.numeric(enhsupp_row))
  }) %>% purrr::reduce(rbind)
  
  ###Get FC data and apply mask
  bas_FC_data_mask <- bas_FC_data*mask
  ins_FC_data_mask <- ins_FC_data*mask
  
  ###Filter FCs and declare bas and ins changes
  bas_changes <- bas_FC_data_mask > FC_threshold
  ins_changes <- ins_FC_data_mask < -FC_threshold
  ##Combine
  bas_changes[which(bas_changes == TRUE)] <- "bas"
  bas_changes[which(bas_changes != "bas")] <- NA
  ins_changes[which(ins_changes == TRUE)] <- "ins"
  ins_changes[which(ins_changes != "ins")] <- NA
  basins_changes <- matrix(NA, nrow = nrow(bas_changes), ncol = ncol(bas_changes))
  for (i in 1:nrow(basins_changes)){
    for (j in 1:ncol(basins_changes)){
      vec <- c(bas_changes[i, j], ins_changes[i, j])
      vec <- vec[which(!is.na(vec))]
      if (length(vec) > 0){
        basins_changes[i, j] <- paste(vec, collapse = ";")
      }
    }
  }
  ##Add "undetermined"
  basins_changes[is.na(basins_changes) & !is.na(enhsupp)] <- "undetermined"
  ##Turn to df
  colnames(basins_changes) <- strains
  rownames(basins_changes) <- rownames(data)
  basins_changes <- as.data.frame(basins_changes)
  ##REMOVED (Add total col)
  basins_total <- apply(basins_changes, 1, function(vec){
    vec <- vec[which(!is.na(vec))]
    vec <- unique(vec)
    if (length(vec) == 0){
      return(NA)
    } else if (length(vec) == 1 & vec[1] == "undetermined"){
      return("undetermined")
    } else {
      vec <- vec[which(vec != "undetermined")]
      vec <- strsplit(vec, ";") %>% unlist %>% unique %>% sort
      return(paste(vec, collapse = ";"))
    }
  })
  #basins_changes$total <- basins_total
  return(basins_changes)
}

#####Summary#####
###Function: Summarise basins changes. 
#data: basins_changes output
basins_change_summariser <- function(data){
  levels <- c("undetermined", "bas;ins", "ins", "bas")
  basins_summ <- map(colnames(data), function(col){
    df <- data.frame("level" = levels,
                     "num" = map(levels, ~length(which(data[, col] == .))) %>% unlist)
    df$Strain <- col
    return(df)
  }) %>% purrr::reduce(rbind) %>%
    mutate(Strain = factor(Strain, levels = colnames(data)),
           level = factor(level, levels = levels))
  return(basins_summ)
}

#####Visualisations#####
###Function: Barplot to summarise basins differences
basins_summ_bplot <- function(basins_summ, 
                              ylab = "Phosphopeptides",
                              bplot_cols){
  output_plot <- ggplot(basins_summ,
                        aes(x = Strain, y = num,
                            colour = Strain, fill = Strain,
                            alpha = level, pattern = level,
                            pattern_fill = Strain, pattern_colour = Strain)) +
    #geoms
    geom_col_pattern(position = "stack", width = 0.8,
                     linewidth = 0.235/2) +
    #scales
    scale_colour_manual(values = bplot_cols) +
    scale_fill_manual(values = bplot_cols) +
    scale_pattern_colour_manual(values = bplot_cols) +
    scale_pattern_fill_manual(values = bplot_cols) +
    scale_alpha_manual(values = c(0, 1, 0.3, 0)) +
    scale_pattern_manual(values = c("stripe", rep("none", 3))) +
    #themes
    comfy_theme(include_xaxis = FALSE,
                rotate_x_text = TRUE) +
    labs(y = ylab) +
    theme(axis.title.x = element_blank())
  return(output_plot)
}


































