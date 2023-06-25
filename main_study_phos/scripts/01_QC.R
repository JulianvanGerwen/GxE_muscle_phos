###Background
#Here I store code for QC (PCAs, quantification in each sample, etc)

###Initialise
library(ggplot2)
library(tidyverse)
library(purrr)




#####Quantification#####

###Function: Number of unique quantified items, as specified by id_col
quant_nums_onecol <- function(data,
                              id_col,
                              num_cols){
  #Get quant nums in each sample
  quant_nums <- map(num_cols,
                    ~length(unique(data[!is.na(data[, .]) > 0, id_col])))
  names(quant_nums) <- num_cols
  #Get total quant
  quant_nums$total <- length(unique(data[rowSums(!is.na(data[, num_cols])) > 0, id_col]))
  return(unlist(quant_nums))
}

###Function: Number of unique quantified items, as specified by id_cols
quant_nums_multcols <- function(data,
                                id_cols,
                                num_cols){
  #Get quants
  quants <- map(id_cols,
                ~quant_nums_onecol(data = data,
                                   id_col = .,
                                   num_cols = num_cols))
  names <- names(quants[[1]])
  quant_df <- purrr::reduce(quants, cbind) %>% as.data.frame
  rownames(quant_df) <- names
  colnames(quant_df) <- id_cols
  quant_df$sample <- rownames(quant_df)
  return(quant_df)
}

###Function: Make df showing quantification of ppeptides, psites, and pprtoeins in each sample and overall
quant_df_maker <- function(data,
                           num_cols){
  #Make quant df
  quant_df <- quant_nums_multcols(data,
                                  id_cols = c("standard_name", "uniprot_site", "uniprot"),
                                  num_cols = num_cols)
  #Add levels
  quant_df[, c("Strain", "Diet", "Ins")] <- strsplit(rownames(quant_df), "_") %>%
    purrr::reduce(rbind)
  quant_df[rownames(quant_df) == "total", c("Strain", "Diet", "Ins")] <- "total"
  quant_df <- mutate(quant_df,
                     Strain = factor(Strain, levels = c(Strain_levels, "total")),
                     Diet = factor(Diet, levels = c(Diet_levels, "total")))
  quant_df$sample <- factor(quant_df$sample, levels = c(UIDs, "total"))
  #Rename
  quant_df <- rename(quant_df,
                     "Phosphopeptides" = "standard_name",
                     "Phosphosites" = "uniprot_site",
                     "Phosphoproteins" = "uniprot")
  return(quant_df)
}

###Function: Barplot of items quantified in each sample
#quant_df: df as output by quant_df_maker
quant_barplot <- function(quant_df,
                          ycol = "Phosphopeptides",
                          file,
                          width = 4, height = 2){
  #Colours
  colours <- fivestraincols_main
  colours["total"] <- "grey"
  #Make plot
  ggplot(quant_df, 
         aes_string(x = "sample", y = ycol, fill = "Strain")) +
    geom_col(width = 0.7) +
    comfy_theme(include_xaxis = FALSE,
                rotate_x_text = TRUE,
                x_text_angle = 90) +
    theme(axis.title.x = element_blank()) +
    scale_fill_manual(values = colours) +
    labs(y = ycol)
  ggsave_pdfpng(file,
                width = width, height = height)
}


###Function: Histogram of number of ppeptides quantification each number of times
quant_hist <- function(data, num_cols){
  quant_vec <- rowSums(!is.na(data[, num_cols]))
  quant_df <- data.frame(quant = quant_vec)
  hist <- ggplot(quant_df,
                 aes(x = quant)) +
    geom_histogram(fill = "#000000") +
    comfy_theme() +
    labs(x = "# quantified samples",
         y = "# phosphopeptides")
  return(hist)
}

###Function: List of conditions each pppeptide is quantified in
quant_condition_list_maker <- function(data,
                                 num_cols,
                                 conditions){
  #Boolean list for quantification in each condition
  condition_quant_bool <-
    map(conditions,
        function(condition){
          condition_cols <- num_cols[grep(condition, num_cols)]
          quant_df <- data.frame(!is.na(data[, condition_cols]))
          quant_vec <- rowSums(quant_df) > 0
        }) %>%
    transpose
  condition_quants <- map(condition_quant_bool,
                          ~conditions[unlist(.)])
  return(condition_quants)
}

###Function: Number of Conditions each ppeptide is in
#factors: The factors to be used for Conditions (e.g. Strain, Diet)
quant_condition_df_maker <- function(data,
                               num_cols,
                               factors){
  #Make df with quantified conditions per ppeptide
  ppeptide_quant_df <- data.frame(standard_name = data$standard_name)
  ppeptide_quant_df$quant_samples_list <- quant_condition_list_maker(data,
                                                                     num_cols = num_cols,
                                                                     conditions = num_cols)
  #Quantified conditions
  #Loop over factors to get quant list
  factors <- c("Strain", "StrainDiet", "StrainDietIns")
  factor_quant_list <- map(factors,
                           ~quant_condition_list_maker(data,
                                                 num_cols,
                                                 all_levels[[.]]))
  #Add to df by flattening list
  factor_quant_list_flat <- map(factor_quant_list,
                                function(list){map(list, ~paste(., collapse = ";")) %>%
                                    unlist})
  #Add to df then turn back into list
  ppeptide_quant_df[, paste("quant_", factors, sep = "")] <- purrr::reduce(factor_quant_list_flat,
                                                                           cbind)
  for (factor in factors){
    colnames(ppeptide_quant_df) <- gsub(factor, "Value", colnames(ppeptide_quant_df))
    ppeptide_quant_df <- mutate(ppeptide_quant_df,
                                quant_Value_list = strsplit(quant_Value, ";"))
    colnames(ppeptide_quant_df) <- gsub("Value", factor, colnames(ppeptide_quant_df))
  }
  return(ppeptide_quant_df)
}


#####PCA, correlation, and clustering#####

###Function: Plot of PC1 and PC2 from GxE muscle phos data
pca_plot <- function(data,
                     num_cols,
                     PCs = c(1, 2)){
  ###PCA
  #Run PCA and make df
  pca <- prcomp(t(na.omit(data[, num_cols])), center = TRUE, scale. = FALSE)
  pca_df <- as.data.frame(pca$x, stringsAsFactors = FALSE)
  pca_df <- extract_factors_from_rownames(pca_df)
  pca_df$Diet_Ins <- factor(paste(pca_df$Diet, pca_df$Ins, sep = "_"),
                            levels = c("CHOW_bas", "CHOW_ins",
                                       "HFD_bas", "HFD_ins"))
  #Make PC labels
  PC_labels <- paste("PC", PCs, sep = "")
  #Plot
  pca_plot <- corrplot_phos(pca_df, yaxis = PC_labels[1], xaxis = PC_labels[2]) +
    labs(x = paste(c(PC_labels[1], " (",
                     100*summary(pca)$importance[2, PCs[1]],
                     "%)"),
                   collapse = ""),
         y = paste(c(PC_labels[2], " (",
                     100*summary(pca)$importance[2, PCs[2]],
                     "%)"),
                   collapse = "")) +
    comfy_theme()
  return(pca_plot)
}

###Function: Make correlation matrix of samples with each other, and give summary stats (mean, max, median, min)
corr_m_maker <- function(data, num_cols){
  corr_m <- cor(data[, num_cols], use = "pairwise.complete.obs")
  #Get summary stats
  corr_m_summary <- corr_m
  for (i in 1:nrow(corr_m_summary)){
    corr_m_summary[i, i] <- NA
  }
  sum_stats <- c("max" = max(corr_m_summary, na.rm = TRUE),
                 "min" = min(corr_m_summary, na.rm = TRUE),
                 "mean" = mean(corr_m_summary, na.rm = TRUE),
                 "median" = median(corr_m_summary, na.rm = TRUE))
  return(list("corr_m" = corr_m,
              "sum_stats" = sum_stats))
}

library(corrplot)
library(RColorBrewer)
###Function: Plot correlation matrix, displaying all values below a threhsold as black
corr_m_plotter <- function(corr_m,
                           black_fraction = 7,
                           total_fraction = 8,
                           text_size = 0.5*7/6){
  corrplot(corr_m,
           method = "color",
           type = "upper",
           order = "hclust",
           #Text parameters
           tl.cex = text_size,
           tl.col = "#000000",
           #legend paramters,
           cl.cex = 0.5*7/6,
           cl.align.text = "l",
           col = c(rep("black", black_fraction*100),
                   colorRampPalette(brewer.pal(9, "Blues"))((total_fraction - black_fraction)*100)),
           title = "Sample correlation")
}

library(dendextend)
###Function: Clustered dendogram by hierarchical clustering
#colours: Vector of colours to be used
#colour_ns: Number of times each colour appears
#When multiplied out by colour_ns
clustered_dendogram <- function(data,
                                num_cols,
                                colours,
                                label_size = 0.2*7/6){
  
  ###Clustered dendogram
  #Build dendogram
  dist_m <- dist(t(data[, num_cols]))
  h_clust <- hclust(dist_m)
  h_dend <- as.dendrogram(h_clust)
  #Set up colours
  colours <- colours[order.dendrogram(h_dend)]
  labels_colors(h_dend) <- colours
  labels_cex(h_dend) <- label_size
  return(h_dend)
}


