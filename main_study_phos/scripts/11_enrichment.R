###Background
#Here I store functions for enrichment analysis e.g. kinase enrichment

###Initialise
library(tidyverse)
library(purrr)

#####Kinase enrichment and further analysis#####

##Function: Heatmap of all kinase substrates from kin_list
#Maps using PSP_SITE_GRP_ID by default
heatmap_kinase <- function(data, kin_list, kinase, heatmap_func, 
                           promisc_threshold = NULL, ...){
  #Remove promiscuous psites if desired
  if (!is.null(promisc_threshold)){
    promisc_keep_bool <- promisc_psites_trim(stat_names = data$PSP_SITE_GRP_ID,
                                             kin_sub_list = kin_list,
                                             promisc_threshold = promisc_threshold)
    data <- data[promisc_keep_bool, ]
  }
  #Reduce to kinase substrates
  data <- subset(data, PSP_SITE_GRP_ID %in% kin_list[[kinase]])
  #Hmap
  output_plot <- heatmap_func(data = data,
                              ...)
  return(output_plot)
}

#Function:  -log10 pval, signed by FC direction
#inf_val: -log10 pvalue value to give when p = 0
signed_nlog10_pval <- function(FCs,
                               pvals,
                               inf_val = NULL){
  
  #Calcualte -log10pvals
  nlog10_pvals <- -log10(pvals)
  #if inf_val is not supplied, set to the max
  if (is.null(inf_val)){
    temp_vals <- nlog10_pvals[which(nlog10_pvals != Inf)]
    inf_val <- max(temp_vals, na.rm = TRUE)
  }
  nlog10_pvals[which(nlog10_pvals == Inf)] <- inf_val
  #Sign
  signed_nlog10_pvals <- sign(FCs)*nlog10_pvals
  return(signed_nlog10_pvals)
}


###Function: Pivot KSEA data to longer
#data: KSEA output data
#value_cols: columns to be plotted as yaxis values
#pval_cols: the columns of p-values that correspond to value_cols
KSEA_pivot_longer <- function(data,
                              value_cols = NULL,
                              pval_cols = NULL){
  #Set up value_cols and pval_cols
  #Set to _ES and _p if nothing specified
  if (is.null(value_cols)){
    value_cols <- colnames(data) %>% .[grep("_ES$", .)]
  }
  if (is.null(pval_cols)){
    pval_cols <- colnames(data) %>% .[grep("_p$", .)]
  }
  #Add kinase
  data$kinase <- rownames(data)
  #Melt
  data_long <- pivot_longer(data[, c("kinase", value_cols)], cols = value_cols,
                            names_to = "sample",
                            values_to = "val")
  data_long_p <- pivot_longer(data[, c("kinase", pval_cols)], cols = pval_cols,
                              names_to = "sample",
                              values_to = "p")
  #Tidy and add columns
  data_long <- data_long %>%
    mutate(kinase = factor(kinase, levels = unique(kinase)),
           sample = factor(sample, levels = unique(sample)),
           p = data_long_p$p)
  #Add Strain and Diet cols
  StrainDiet_list <- strsplit(as.character(data_long$sample), "_") %>% map(~.[c(1, 2)])
  data_long <- mutate(data_long, 
                      StrainDiet = factor(map(StrainDiet_list, ~paste(., collapse = "_")),
                                          levels = all_levels$StrainDiet),
                      Strain = factor(map(StrainDiet_list, ~.[1]),
                                      levels = all_levels$Strain),
                      Diet = factor(map(StrainDiet_list, ~.[2]),
                                    levels = all_levels$Diet))
  return(data_long)
}


###Function: Prepare KSEA data for barplot
#data: KSEA output data
#value_cols: columns to be plotted as yaxis values
#pval_cols: the columns of p-values that correspond to value_cols
#order_col: column to order kinases by
#sig_cutoff: p-value cutoff to classify as significantly regulated
#sig_kinases_only: if TRUE, only retain kinases that are significantly regulated in at least one condition
KSEA_bplot_data_prepper <- function(data,
                                    value_cols,
                                    pval_cols,
                                    order_col,
                                    order_decreasing = TRUE,
                                    sig_cutoff = 0.05,
                                    sig_kinases_only = FALSE){
  #Add kinase
  data$kinase <- rownames(data)
  #Order data
  data <- data[order(data[, order_col], decreasing = order_decreasing), ]
  #Melt
  data_long <- KSEA_pivot_longer(data = data,
                                 value_cols = value_cols,
                                 pval_cols = pval_cols)
  
  #Restrict to kinases significant in at least one sample
  data_long$sig <- FALSE
  data_long[which(data_long$p < sig_cutoff), "sig"] <- TRUE
  if (sig_kinases_only){
    kinases_to_keep_bool <- map(levels(data_long$kinase),
                                function(kin){
                                  sig <- subset(data_long, kinase == kin) %>% .$sig
                                  return(sum(sig) > 0)
                                }) %>% unlist
    kinases_to_keep <- levels(data_long$kinase)[kinases_to_keep_bool]
    data_long <- data_long[which(data_long$kinase %in% kinases_to_keep), ] %>%
      mutate(kinase = factor(kinase, levels = unique(kinase)))
  }
  return(data_long)
}

###Function: Simple dodged KSEA barplot
#data_long: Prepared KSEA data as output by KSEA_bplot_data_prepper
#colours: a vector of colours corresponding to each sample
KSEA_bplot_simple_plotter <- function(data_long,
                                      colours,
                                      padding = 0,
                                      width = 0.5,
                                      ylab = "val"){
  #Prep colours
  names(colours) <- levels(data_long$sample)
  
  #Plot
  output_plot <- ggplot(data_long,
                        aes(x = kinase, y = val, fill = sample)) +
    geom_bar(stat = "identity",
             size = 0.235,
             width = width,
             alpha = 1,
             #use position_dodge2 to dodge columsn based on the aesthetic group
             #padding governs how close dodged columns are to each other
             position = position_dodge2(preserve = "single", padding = padding)) +
    geom_hline(yintercept = 0,
               colour = "black", size = 0.235) +
    comfy_theme(include_xaxis = FALSE,
                rotate_x_text = TRUE, x_text_angle = 90,
                include_legend = FALSE) +
    scale_fill_manual(values = colours) +
    labs(y = ylab)
  return(output_plot)
}

###Function: KSEA barplot simple, Wrapper for functions above
KSEA_bplot_simple <- function(data,
                                   value_cols,
                                   pval_cols,
                                   order_col,
                              order_decreasing = TRUE,
                                   sig_cutoff = 0.05,
                                   sig_kinases_only = FALSE,
                                   colours,
                                   ...){
  #Prep data
  data_long <- KSEA_bplot_data_prepper(data = data,
                                       value_cols = value_cols,
                                       pval_cols = pval_cols,
                                       order_col = order_col,
                                       order_decreasing = order_decreasing,
                                       sig_cutoff = sig_cutoff,
                                       sig_kinases_only = sig_kinases_only)
  #Plot
  output_plot <- KSEA_bplot_simple_plotter(data_long = data_long,
                                                colours = colours,
                                                ...)
  return(output_plot)
}



###Function: Dodged KSEA barplot with alpha for insignificant enrichment
#data_long: Prepared KSEA data as output by KSEA_bplot_data_prepper
#colours: a vector of colours corresponding to each sample
KSEA_bplot_sigalpha_plotter <- function(data_long,
                                        colours,
                                        padding = 0,
                                        width = 0.5,
                                        insig_alpha = 0.4,
                                        ylab = "val"){
  #Prep colours
  names(colours) <- levels(data_long$sample)
  
  #Plot
  output_plot <- ggplot(data_long,
                        aes(x = kinase, y = val, fill = sample,
                            alpha = sig)) +
    geom_bar(stat = "identity",
             size = 0.235,
             width = width,
             #use position_dodge2 to dodge columsn based on the aesthetic group
             #padding governs how close dodged columns are to each other
             position = position_dodge2(preserve = "single", padding = padding)) +
    geom_hline(yintercept = 0,
               colour = "black", size = 0.235) +
    comfy_theme(include_xaxis = FALSE,
                rotate_x_text = TRUE, x_text_angle = 90,
                include_legend = FALSE) +
    scale_fill_manual(values = colours) +
    scale_alpha_manual(values = c(insig_alpha, 1)) +
    labs(y = ylab)
  return(output_plot)
}

###Function: KSEA barplot simple, Wrapper for functions above
KSEA_bplot_sigalpha <- function(data,
                              value_cols,
                              pval_cols,
                              order_col,
                              order_decreasing = TRUE,
                              sig_cutoff = 0.05,
                              sig_kinases_only = FALSE,
                              colours,
                              ...){
  #Prep data
  data_long <- KSEA_bplot_data_prepper(data = data,
                                       value_cols = value_cols,
                                       pval_cols = pval_cols,
                                       order_col = order_col,
                                       order_decreasing = order_decreasing,
                                       sig_cutoff = sig_cutoff,
                                       sig_kinases_only = sig_kinases_only)
  #Plot
  output_plot <- KSEA_bplot_sigalpha_plotter(data_long = data_long,
                                           colours = colours,
                                           ...)
  return(output_plot)
}

###Function: Dodged KSEA barplot that is hollow for insignificant enrichment
#data_long: Prepared KSEA data as output by KSEA_bplot_data_prepper
#colours: a vector of colours corresponding to each sample
KSEA_bplot_sighollow_plotter <- function(data_long,
                                         colours,
                                         width = 0.5,
                                         padding = 0.1,
                                         ylab = "val"){
  #Prep colours
  names(colours) <- levels(data_long$sample)
  
  #Prep sig data (set all else to NA)
  data_long_sig <- data_long
  data_long_sig[which(data_long_sig$sig == FALSE), "val"] <- NA
  
  #Plot
  output_plot <- ggplot() +
    #Outline for all data
    geom_bar(data = data_long,
             aes(x = kinase, y = val, colour = sample),
             fill = "white",
             stat = "identity",
             size = 0.235,
             width = width,
             alpha = 1,
             #use position_dodge2 to dodge columsn based on the aesthetic group
             #padding governs how close dodged columns are to each other
             position = position_dodge2(preserve = "single", padding = padding)) +
    
    #Fill for sig data
    geom_bar(data = data_long_sig,
             aes(x = kinase, y = val, fill = sample),
             stat = "identity",
             width = width,
             alpha = 1,
             position = position_dodge2(preserve = "single", padding = padding)) +
    geom_hline(yintercept = 0,
               colour = "black", size = 0.235) +
    comfy_theme(include_xaxis = FALSE,
                rotate_x_text = TRUE, x_text_angle = 90,
                include_legend = FALSE) +
    scale_colour_manual(values = colours) +
    scale_fill_manual(values = colours) +
    labs(y = ylab)
  return(output_plot)
}

###Function: KSEA barplot simple, Wrapper for functions above
KSEA_bplot_sighollow <- function(data,
                              value_cols,
                              pval_cols,
                              order_col,
                              order_decreasing = TRUE,
                              sig_cutoff = 0.05,
                              sig_kinases_only = FALSE,
                              colours,
                              ...){
  #Prep data
  data_long <- KSEA_bplot_data_prepper(data = data,
                                       value_cols = value_cols,
                                       pval_cols = pval_cols,
                                       order_col = order_col,
                                       order_decreasing = order_decreasing,
                                       sig_cutoff = sig_cutoff,
                                       sig_kinases_only = sig_kinases_only)
  #Plot
  output_plot <- KSEA_bplot_sighollow_plotter(data_long = data_long,
                                           colours = colours,
                                           ...)
  return(output_plot)
}

###Function: Dodged KSEA barplot that is plots individual points
#data_long: Prepared KSEA data as output by KSEA_bplot_data_prepper
#value_col_suffix: a suffix which, when removed from value_cols, will give the conditions to go on x-axsi
#colours: a vector of colours corresponding to each condition
KSEA_bplot_indivpoints_plotter <- function(data_long,
                                           value_col_suffix,
                                           colours,
                                           padding = 0,
                                           width = 0.5,
                                           point_size = 0.3,
                                           ylab = "val"){
  #Prep colours
  names(colours) <- levels(data_long$condition)
  
  ##Prep data further
  #Add condition column (the broader conditions that individual samples fit into)
  conditions <- gsub(value_col_suffix, "", data_long$sample) %>% unique
  data_long <- mutate(data_long,
                      condition = factor(gsub(value_col_suffix, "",
                                              sample),
                                         levels = conditions))
  #Summarise data
  data_long_sum <- group_by(data_long,
                            kinase, condition) %>%
    summarise(val = median(val, na.rm = TRUE))
  
  #Plot
  output_plot <- ggplot() +
    #Barplot
    geom_bar(data = data_long_sum,
             aes(x = kinase, y = val, fill = condition),
             stat = "identity", 
             size = 0.235,
             colour = NA,
             width = width,
             alpha = 0.65,
             #use position_dodge2 to dodge columsn based on the aesthetic group
             #padding governs how close dodged columns are to each other
             position = position_dodge2(preserve = "single", padding = padding)) +
    #x-axis
    geom_hline(yintercept = 0,
               colour = "black", size = 0.235) +
    #point plot
    geom_point(data = data_long,
               aes(x = kinase, y = val, colour = condition),
               size = point_size,
               stroke = 0.25,
               position = position_dodge(
                 #Width controls how much to dodge.
                 width = width)) +
    #Themes
    comfy_theme(include_xaxis = FALSE,
                rotate_x_text = TRUE, x_text_angle = 90,
                include_legend = FALSE) +
    scale_fill_manual(values = colours) +
    scale_colour_manual(values = colours) +
    labs(y = ylab)
  return(output_plot)
}

###Function: KSEA barplot with individual points. Wrapper for functions above
KSEA_bplot_indivpoints <- function(data,
                                   value_cols,
                                   pval_cols,
                                   order_col,
                                   order_decreasing = TRUE,
                                   sig_cutoff = 0.05,
                                   sig_kinases_only = FALSE,
                                   value_col_suffix,
                                   colours,
                                   ...){
  #Prep data
  data_long <- KSEA_bplot_data_prepper(data = data,
                                       value_cols = value_cols,
                                       pval_cols = pval_cols,
                                       order_col = order_col,
                                       order_decreasing = order_decreasing,
                                       sig_cutoff = sig_cutoff,
                                       sig_kinases_only = sig_kinases_only)
  #Plot
  output_plot <- KSEA_bplot_indivpoints_plotter(data_long = data_long,
                                                value_col_suffix = value_col_suffix,
                                                colours = colours,
                                                ...)
  return(output_plot)
}

##Function: Bplot KSEA barplots for insbas and insresp
#KSEA_list: List of KSEA outputs where first argument is insbas and second is insresp
#kinases: vector of kinases to display. They are plotted in this order
#dir: Directory and prefix for saving file
#CHOW_only: If TRUE, only plot CHOW values
KSEA_insbas_bplotter <- function(KSEA_list,
                                 kinases = NULL,
                                 dir,
                                 CHOW_only = FALSE,
                                 width = 4, height = 2.5,
                                 ...){
  ##Prep data and set up parameters
  KSEA_list_kinases <- map(KSEA_list, function(data){
    if (is.null(kinases)){
      kinases <- rownames(data)
    }
    data <- data[kinases, ]
    data$kinase_ordered <- factor(kinases, levels = kinases)
    return(data)
  })
  ylab <- "log2(ins/bas) enrichment"
  if (CHOW_only){
    KSEA_colours <- fivestraincols_main
  } else {
    KSEA_colours <- map(fivestraincols_main, ~rep(., 2)) %>% unlist
  }
  
  ##Set up colnames
  cols_touse <- map(KSEA_list_kinases,
                    function(data){
                      cols <- colnames(data)
                      if (CHOW_only){
                        cols <- cols[grep("CHOW", cols)]
                      }
                      return(cols)
                    })
  
  ##Plot
  #insresp
  KSEA_bplot_indivpoints(KSEA_list_kinases$ins_resp,
                         value_cols = cols_touse$ins_resp %>% .[grep("_ES$", .)],
                         pval_cols = cols_touse$ins_resp %>% .[grep("_p$", .)],
                         value_col_suffix = "_\\d+_ES$",
                         order_col = "kinase_ordered",
                         order_decreasing = FALSE,
                         sig_kinases_only = FALSE,
                         colours = KSEA_colours,
                         ylab = ylab,
                         ...)
  ggsave_pdfpng(file = paste(dir, "insresp_ES", sep = ""),
                width = width, height = height)
}

###Function: Prepare KSEA data with individual values for summary heatmap. Takes median in each condition
KSEA_summary_hmap_dataprep <- function(KSEA_data){
  ##Pivot data
  value_col_suffix <- "_\\d+_ES$"
  data_long <- KSEA_bplot_data_prepper(data = KSEA_data,
                                       value_cols = colnames(KSEA_data) %>% .[grep("_ES$", .)],
                                       pval_cols = colnames(KSEA_data) %>% .[grep("_p$", .)],
                                       order_col = colnames(KSEA_data) %>% .[grep("_ES$", .)] %>% .[1],
                                       order_decreasing = TRUE,
                                       sig_cutoff = 0.05,
                                       sig_kinases_only = FALSE)
  ##Prep data further
  #Add condition column (the broader conditions that individual samples fit into)
  conditions <- gsub(value_col_suffix, "", data_long$sample) %>% unique
  data_long <- mutate(data_long,
                      condition = factor(gsub(value_col_suffix, "",
                                              sample),
                                         levels = conditions))
  ##Summarise data
  data_long_sum <- group_by(data_long,
                            kinase, condition) %>%
    summarise(val = median(val, na.rm = TRUE))
  ##Pivot to wider
  data_wide_sum <- pivot_wider(data_long_sum,
                               names_from = "condition",
                               values_from = "val",
                               id_cols = c("kinase")) %>%
    as.data.frame
  rownames(data_wide_sum) <- data_wide_sum$kinase
  data_wide_sum <- dplyr::select(data_wide_sum, -kinase)
  ##Rename columns to prepare for heatmap_KSEA
  conditions <- colnames(data_wide_sum)
  colnames(data_wide_sum) <- paste(conditions, "_ES", sep = "")
  data_wide_sum[, paste(conditions, "_p", sep = "")] <- NA
  ##Outptu
  return(data_wide_sum)
}

###Function: Heatmap of summarised KSEA data from individual values
KSEA_summary_hmap <- function(KSEA_data,
                              conditions,
                              ...){
  #Summarise data
  data_sum <- KSEA_summary_hmap_dataprep(KSEA_data = KSEA_data)
  #Plot
  output_plot <- heatmap_KSEA(data_sum,
                              conditions = conditions,
                              ...)
  return(output_plot)
}

###Function: Heatmaps of KSEA data from insbas (includes pvals) and insresp (summarises as median)
KSEA_insbas_hmapper <- function(KSEA_list,
                                 kinases = NULL,
                                conditions,
                                 dir,
                                 width = 4, height = 2.5){
  #insresp heatmap
  insresp_hmap <- KSEA_summary_hmap(KSEA_data = KSEA_list$ins_resp,
                                    conditions = conditions)
  ggsave_pdfpng(paste(dir, "_insresp", sep = ""),
                width = width, height = height)
}





#####Visualisations#####
###Function: Heatmap of multiple overrepresentation analysis that plots log2 enrichment score and significance
#All log2 enrichment < 0 is shown as 0
overrep_hmap <- function(data,
                         sig_threshold = 0.05,
                         ...){
  ##Get conditions, cols, etc
  ES_log2_cols <- colnames(data) %>% .[grep("ES_log2$", .)]
  ES_cols <- colnames(data) %>% .[grep("ES$", .)]
  adj_p_cols <- colnames(data) %>% .[grep("adj_pval_comb$", .)]
  conditions <- gsub("_ES_log2", "", ES_log2_cols)
  ##Subset for only sig
  sig_rows <- rowSums(data[, adj_p_cols] < sig_threshold, na.rm = TRUE) > 0
  data_sig <- data[sig_rows, ]
  ##Set ES_log2 < 0 to 0
  data_sig[data_sig == -Inf] <- 0
  data_sig[data_sig < 0] <- 0
  ##Make heatmap
  
  print(data_sig)
  
  output_plot <- fc_sig_heatmap(data = data_sig,
                                fc_cols_w_pvals = ES_log2_cols,
                                pval_cols = adj_p_cols,
                                x_axis_names = conditions,
                                order_column = ES_log2_cols[1],
                                yaxis_naming = "keepall",
                                clust_cols = ES_log2_cols,
                                legend_title = "log2 Enrichment",
                                sig_threshold = sig_threshold,
                                ...)
  return(output_plot)
}














