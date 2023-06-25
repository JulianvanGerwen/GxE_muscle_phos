###Background
#Here I store code for basic visualisation e.g. barplots, correlation point plots

###Initialise
library(tidyverse)
library(ggplot2)
library(purrr)


#####Correlation plots#####

#Function: Point plot for correlation of two columns across Strain, Diet, and Ins
#Colours by Strain, shape for Diet (CHOW is circle, HFD is square), bas points have transparent inside
#As of 20221013 this is a wrapper for a more general corrplot implemented in scripts/02_basic_visualisations.R
corrplot_phos <- corrplot_general

###Function: Correlation plot of two ppeptides with one another
#num_cols: The data to use
corrplot_phos_ppeptides <- function(data, xaxis, yaxis, num_cols = UIDs, ...){
  #Prepare data
  data <- data[c(xaxis, yaxis), num_cols]
  data_t <- t(data) %>% extract_factors_from_rownames
  #Plot
  output_plot <- corrplot_phos(data = data_t,
                               yaxis = yaxis,
                               xaxis = xaxis,
                               ...)
  return(output_plot)
}

###Function: Correlation plot of a ppeptide witha phenotype
corrplot_phos_pheno <- function(phos_data,
                                pheno_data,
                                phos_row,
                                pheno_col,
                                swap_axes = FALSE,
                                num_cols,
                                phos_UIDs,
                                ...){
  #Combine data
  phos_pheno_data <- combine_phos_pheno_data(phos_data = phos_data,
                                             pheno_data = pheno_data,
                                             pheno_cols = c(pheno_col),
                                             phos_rows = c(phos_row),
                                             phos_num_cols = num_cols,
                                             phos_UIDs = phos_UIDs)
  #Make plot
  if (swap_axes){
    xaxis <- pheno_col
    yaxis <- phos_row
  } else {
    xaxis <- phos_row
    yaxis <- pheno_col
  }
  output_plot <- corrplot_phos(phos_pheno_data,
                               xaxis = xaxis,
                               yaxis = yaxis,
                               ...)
  return(output_plot)
}

#####Boxplots#####

###Function: Boxplot of GxE muscle phos data
#insresp: If TRUE, plots insulin data normalied to basal in each strain-diet condition
boxplot_GxEphos <- function(data,
                            protein,
                            is_raw = TRUE,
                            CHOW_only = FALSE,
                            insresp = FALSE,
                            ...){
  #Set up conditions to use
  if (insresp){
    conditions <- paste(StrainDietIns_levels, "_basnorm", sep = "")
  } else {
    conditions <- StrainDietIns_levels
  }
  
  #Plot
  if (CHOW_only){
    output_plot <- general_boxplot(data = data,
                                   protein = protein,
                                   conditions = conditions[grep("CHOW", conditions)],
                                   treatments = Strain_levels,
                                   treatment_colours = fivestraincols_main,
                                   is_raw = is_raw,
                                   yaxis_lab = "log2 Intensity",
                                   ...)
  } else {
    output_plot <- general_boxplot(data = data,
                                   protein = protein,
                                   conditions = conditions,
                                   treatments = Strain_levels,
                                   treatment_colours = fivestraincols_main,
                                   is_raw = is_raw,
                                   yaxis_lab = "log2 Intensity",
                                   ...)
  }
  
  return(output_plot)
}



#####Heatmaps#####
###basins heatmap
#Objects used to make basins heatmap
basins_FC_cols <- map(c("bas", "ins"),
                      ~StrainDietIns_levels[grep(., StrainDietIns_levels)]) %>%
  transpose %>%
  map(~ paste(.[[1]], "_vs_", .[[2]], "_logFC", sep = "")) %>%
  unlist
basins_FC_xnames <- strsplit(basins_FC_cols, "_") %>%
  map(~paste(.[1:2], collapse = " ")) %>%
  unlist
###Function: Heatmap for ins vs bas logFCs
heatmap_basins_GxEphos <- function(data,
                                   order_column = "C57Bl6J_CHOW_bas_vs_C57Bl6J_CHOW_ins_logFC",
                                   ...){
  fc_sig_heatmap(data,
                 fc_cols_wo_pvals = basins_FC_cols,
                 order_column = order_column,
                 x_axis_names = basins_FC_xnames,
                 legend_title = "log2 ins/bas",
                 ...)
}

###Heatmap for any of Ins, Diet, and Strain comparisons
#Set up FC cols
Ins_FC_cols <- map(c("bas", "ins"),
                   ~StrainDietIns_levels[grep(., StrainDietIns_levels)]) %>%
  transpose %>%
  map(~ paste(.[[1]], "_vs_", .[[2]], "_logFC", sep = "")) %>%
  unlist
Diet_FC_cols <- map(c("CHOW", "HFD"),
                    function(diet){
                      diet_levels <- map(c("bas", "ins"),
                                         ~StrainDietIns_levels[grep(paste(diet, ., sep = "_"),
                                                                    StrainDietIns_levels)]) %>%
                        unlist
                      return(diet_levels)
                    }) %>%
  transpose %>%
  map(~ paste(.[[1]], "_vs_", .[[2]], "_logFC", sep = "")) %>%
  unlist
Strain_FC_cols <- map(StrainDietIns_levels[-grep("C57Bl6J", StrainDietIns_levels)],
                      function(x){
                        x_split <- strsplit(x, "_")[[1]]
                        x_B6 <- paste(c("C57Bl6J", x_split[-1]), collapse = "_")
                        return(c(x_B6, x))
                      }) %>%
  map(~ paste(.[[1]], "_vs_", .[[2]], "_logFC", sep = "")) %>%
  unlist
Strain_FC_cols_order <- strsplit(Strain_FC_cols, "_") %>% map(~.[2:3]) %>% 
  purrr::reduce(rbind) %>% as.data.frame
Strain_FC_cols_order$FC_col <- Strain_FC_cols
Strain_FC_cols_order <- Strain_FC_cols_order %>% .[order(.$V2, .$V1), ]
Strain_FC_cols <- Strain_FC_cols_order$FC_col
all_FC_cols <- list("Ins" = Ins_FC_cols,
                    "Diet" = Diet_FC_cols,
                    "Strain" = Strain_FC_cols)
#Set up xaxis names
all_xnames <- all_FC_cols
all_xnames$Ins <- strsplit(all_xnames$Ins, "_") %>%
  map(~ paste(.[1], .[2], "ins/bas", sep = " ")) %>%
  unlist
all_xnames$Diet <- strsplit(all_xnames$Diet, "_") %>%
  map(~ paste(.[1], .[3], "HFD/CHOW", sep = " ")) %>%
  unlist
all_xnames$Strain <- strsplit(all_xnames$Strain, "_") %>%
  map(~ paste(.[2],  " ", .[3], " ", .[5], "/C57Bl6J", sep = "")) %>%
  unlist

###Function: Heatmap for ins vs bas logFCs
heatmap_GxEphos <- function(data,
                            order_column = NULL,
                            comparisons = c("Ins"),
                            CHOW_only = FALSE,
                            insresp = FALSE,
                            ...){
  #Set up FC cols and xaxis names
  FC_cols <- unlist(all_FC_cols[comparisons])
  xnames <- unlist(all_xnames[comparisons])
  #Change cols and xnames if insresp
  if (insresp){
    FC_cols <- colnames(data) %>% .[grep("basnorm_\\d+$", .)]
    xnames <- gsub("basnorm_", "", FC_cols)
  }
  #Restrict cols and names if CHOW_only
  if (CHOW_only){
    FC_cols <- FC_cols[grep("CHOW", FC_cols)]
    xnames <- xnames[grep("CHOW", xnames)]
  }
  #set up order_column if not specified
  if (is.null(order_column)){
    order_column <- FC_cols[1]
  }
  fc_sig_heatmap(data,
                 fc_cols_wo_pvals = FC_cols,
                 order_column = order_column,
                 x_axis_names = xnames,
                 legend_title = "log2FC",
                 ...)
}

###Function: Clustered heatmaps for GxE muscle phos
heatmap_GxEphos_clustered <- function(data, clust_cols, ...){
  #Cluster for ordering
  dist_m <- dist(data[, clust_cols])
  hclust <- hclust(dist_m)
  data <- data[hclust$order, ]
  data$order_col <- 1:nrow(data)
  output_plot <- heatmap_GxEphos(data = data,
                                 order_column = "order_col",
                                 ...)
  return(output_plot)
}