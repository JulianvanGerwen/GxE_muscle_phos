###Background
#Here I keep all my main visualisation functions


###Packages
library(reshape2)
library(ggplot2)
library(RColorBrewer)
library(ggrepel)
library(scales)
library(rlist)
library(gridExtra)
library(gtable)
library(grid)
library(ggVennDiagram)
library(ggpubr)
library(tidyverse)

###comfy_theme
comfy_theme <- function(include_xaxis = TRUE, include_yaxis = TRUE,
                        include_legend = TRUE,
                        rotate_x_text = FALSE, x_text_angle = 45){
    #Base theme
    base_theme <- theme(panel.background = element_blank(),
            panel.grid = element_blank(),
            axis.line.x = element_line(colour = "black", size = 0.235),
            axis.line.y = element_line(colour = "black", size = 0.235),
            axis.ticks.x = element_line(size = 0.235,
                                        colour = "black"),
            axis.ticks.y = element_line(size = 0.235,
                                        colour = "black"),
            axis.text.x = element_text(colour = "black",
                                       size = 7),
            axis.text.y = element_text(colour = "black",
                                       size = 7),
            axis.title.x = element_text(colour = "black",
                                        size = 7),
            axis.title.y = element_text(colour = "black",
                                        size = 7),
            plot.title = element_text(colour = "black",
                                      size = 7,
                                      hjust = 0.5),
            legend.text = element_text(colour = "black",
                                       size = 7),
            legend.title = element_text(colour = "black",
                                        size = 7),
            legend.key = element_blank(),
            legend.key.size = unit(0.25, 'cm'),
            strip.text = element_text(size = 7, colour = "black"),
            strip.background = element_blank())
    #Build list for additional arguments
    if (include_legend){
      legend_position <- NULL
    }else{
      legend_position <- "none"
    }
    if (include_xaxis){
      axis_line_x <- NULL
      axis_ticks_x <- NULL
    }else{
      axis_line_x <- element_blank()
      axis_ticks_x <- element_blank()
    }
    if (include_yaxis){
      axis_line_y <- NULL
      axis_ticks_y <- NULL
    }else{
      axis_line_y <- element_blank()
      axis_ticks_y <- element_blank()
    }
    #Rotate x text if desired
    if (rotate_x_text){
      #if 90 degrees
      if (x_text_angle == 90){
        axis_text_x <- element_text(colour = "black",
                                    size = 7,
                                    hjust = 1, vjust = 0.5, angle = 90)
      } else if (x_text_angle == 45){
        axis_text_x <- element_text(colour = "black",
                                    size = 7,
                                    hjust = 1, vjust = 1, angle = 45)
      }
      
    } else{
      axis_text_x <- NULL
    }
    theme_args <- list("legend.position" = legend_position,
                       "axis.line.x" = axis_line_x,
                       "axis.ticks.x" = axis_ticks_x,
                       "axis.line.y" = axis_line_y,
                       "axis.ticks.y" = axis_ticks_y,
                       "axis.text.x" = axis_text_x)
    theme_args <- theme_args[vapply(theme_args, function(x) {!is.null(x)}, logical(1))]
    
    #Make secondary theme
    secondary_theme <- do.call(function(...) theme(...), theme_args)
    return(base_theme + secondary_theme)
}








###Function: Volcano plot given pval and FC columns (can be adj_pval, and/or rval)
#Up to date 20210926
#By default, colours by pval and FC cutoff
#User can specificy a categorical column to colour points by, through the argument "colour_col". Then, the levels of this categorical column and the colours must be specified (colour_col_levels, colour_col_colours). Note that points will be placed in the oredr of the colour_col_levels levels, so order wisely 
#Colours significant (pval and FC) proteins
#label_col specifies were to get labels from (if desired), and rows_to_label points at what rows to label (needs to be row names)
volcano_plot_coloured <- function(data,
                                  pval_col,
                                  FC_col,
                                  pval_cutoff = 0.05,
                                  FC_cutoff = 0.58,
                                  inner_alpha_val = 0.3,
                                  outer_alpha_val = 0.5,
                                  size = 1,
                                  sig_colour = "red",
                                  colour_col = FALSE,
                                  colour_col_levels,
                                  colour_col_colours,
                                  x_lim = "max",
                                  y_lim = "max",
                                  remove_Inf = TRUE,
                                  return_df = FALSE,
                                  show_legend = FALSE,
                                  label_col = FALSE,
                                  rows_to_label = c()){
  
  ##Set up data
  plot_data <- data
  colnames(plot_data)[colnames(plot_data) == FC_col] <- "FC"
  colnames(plot_data)[colnames(plot_data) == pval_col] <- "pval"
  plot_data$nlog_pval <- -log10(plot_data$pval)
  plot_data$sig <- FALSE
  plot_data$sig[which(abs(plot_data$FC) > FC_cutoff &
                        plot_data$pval < pval_cutoff)] <- TRUE
  
  ##If colour col
  if (colour_col != FALSE){
    
    plot_data[, "colour_col"] <- data[, colour_col]
    
    #Tag so I colour by the right name
    colour_by <- "colour_col"
    
    #Adjusting colours
    colour_adjusting <- scale_colour_manual(values = alpha(colour_col_colours,
                                                           outer_alpha_val)) 
    fill_adjusting <-  scale_fill_manual(values = alpha(colour_col_colours,
                                                        alpha = inner_alpha_val)) 
  } else {
    
    colour_by <- "sig"
    colour_adjusting <- scale_colour_manual(values = alpha(c("black",
                                                             sig_colour),
                                                           outer_alpha_val)) 
    fill_adjusting <- scale_fill_manual(values = alpha(c("black",
                                                         sig_colour),
                                                       alpha = inner_alpha_val))
  }
  
  #If label_col
  if (label_col != FALSE){
    
    plot_data[, "label_col"] <- data[, label_col]
    plot_data[which(rownames(plot_data) %in%
                      rows_to_label == FALSE), "label_col"] <- ""
    
    label_plot <- geom_text_repel(aes_string(x = "FC",
                                             y = "nlog_pval",
                                             colour = colour_by,
                                             label = "label_col"),
                                  size = 1.875,
                                  min.segment.length = 0,
                                  max.overlaps = Inf) 
  } else {
    
    label_plot <- NULL
  }

  
  
  
  ##Remove Inf if desired
  if (remove_Inf){
    
    plot_data <- plot_data[which(plot_data$nlog_pval !=
                                   Inf), ]
  }
  
  ##Get x_lim and y_lim if using max
  if (x_lim == "max"){
    
    x_lim <- max(abs(plot_data$FC[which(is.na(plot_data$pval) == FALSE)]),
                 na.rm = TRUE)
  }  
  
  if (y_lim == "max"){
    
    y_lim <- max(plot_data$nlog_pval[which(plot_data$nlog_pval != Inf)])
  }
  
  ##Order df
  
  #If colour by sig:
  if (colour_col == FALSE){
    
    plot_data <- plot_data[order(plot_data$sig), ] 
  } else {
    
    #Make colour_col a factor
    plot_data$colour_col <- factor(plot_data$colour_col,
                                   levels = colour_col_levels)
    plot_data <- plot_data[order(plot_data$colour_col), ]
  }
  
  
  
  ##Decide whether to keep legend
  if (show_legend == FALSE){
    
    legend_theme <- theme(legend.position = "none")
  } else {
    
    legend_theme = NULL
  }
  
  
  
  
  
  ##Output
  
  #return_df
  if (return_df){
    
    return(plot_data)
    
    #Colour by sig
  } else {
    output_plot <- ggplot(data = plot_data,
                          aes_string(x = "FC",
                                     y = "nlog_pval",
                                     colour = colour_by,
                                     fill = colour_by))+
      geom_point(size = size,
                 shape = 21,
                 stroke = size*0.3) +
      label_plot +
      colour_adjusting +
      fill_adjusting +
      comfy_theme() +
      legend_theme +
      xlim(-x_lim, x_lim) +
      ylim(0, y_lim)
    return(output_plot) 
  } 
}



###Funcion that converts unique psites formatted formally into nice format
#mult = TRUE: gene_uniprot_site_multiplicity to gene site Pmultplicity (or gene uniprot site Pmultiplicity if there are double ups)
#mult = FALSE: gene_uniprot_site to gene site (or gene uniprot site if there are double ups)
psite_to_nice_name <- function(psites, mult = TRUE){
  name_df <- as.data.frame(t(sapply(psites, function(x) strsplit(x, "_")[[1]])))
  #If multiplicity is present
  if (mult){
    names(name_df) <- c("gene", "uniprot", "site", "mult")
    name_df$mult_nice <- paste("P", name_df$mult, sep = "")
    name_df$site_mult_nice <- paste(name_df$site, name_df$mult_nice, sep = " ")
    #If no multiplicity
  } else {
    names(name_df) <- c("gene", "uniprot", "site")
    name_df$site_mult_nice <- name_df$site
  }
  name_df$nice_name_old <- apply(name_df[, c("gene", "site_mult_nice")],
                                 1, function(x) paste(x, collapse = " "))
  #Include uniprots for sites that have same gene, site, and multiplicity, or no gene
  name_df$nice_name <- name_df$nice_name_old
  for (i in 1:nrow(name_df)){
    if (length(which(name_df$nice_name_old == name_df$nice_name[i])) > 1 |
        name_df$gene[i] %in% c("NA", "")){
      name_df$nice_name[i] <- paste(c(name_df[i, c("gene", "uniprot", "site_mult_nice")]), collapse = " ")
    }
  }
  #Strip of NA
  name_df$nice_name <- gsub("NA ", "", name_df$nice_name)
  return(name_df$nice_name)
}

###Functions for naming yaxis items in fc_sig_heatmap
gene_uniprot_site_mult_naming <- function(data){
  data$standard_name <- psite_to_nice_name(rownames(data))
  return(data)
}
gene_uniprot_naming <- function(data){
  data$standard_name <- sapply(rownames(data),
                               function(x) strsplit(x, "_")[[1]][1])
  #Double up genenames:
  data$old_standard_name <- data$standard_name
  for (i in 1:nrow(data)){
    
    if (length(which(data$old_standard_name ==
                     data$standard_name[i])) > 1){
      
      data$standard_name[i] <- paste(strsplit(rownames(data)[i],
                                              "_")[[1]],
                                     collapse = " ")
    }
  }
  return(data)
}


###Heatmap function
#Up to date 20211020
##Plots FC on heatmap
##Plots significance as point in tile
#Arguments:
#data: your data
#fc_cols_w_pvals: FC columns that have corresponding pval columns
#fc_cols_wo_pvals: FC columns without corresponding pval columns
#pval_cols: p-value columns
#shape2_pval_cols: p-value columns where you want a different shape for the significance point
#order_column: The FC column to order by
#is_decreasing: How you order
#return_df: Returns the dataframe that is put into ggplot
#x_axis_names: names for x_axis, correspond to fc_cols
#gap_column_index: index of column where you want a gap
#scale_colours_manually: do you want to scale the gradient manually? If so, you can change where the middle colour falls with mid_colour_threshold, and where the max colour falls with high_colour_threhsold
#mid_colour_threshold: See scale_colours_manually Typically you would want to highlight the changes between 0 and mid_colour_threshold
#high_colour_threshold: See scale_colours_manually
#colour_scheme: What colour scheme to use. Accepts "original" (quite gross), or brewer palette codes e.g. "RdBu". To reverse pallet, enter "rev_RdBu"
#middle_brewer_white: If using brewer palette, sets middle colour to white (FFFFFF) if true
#yaxis_naming: How rownames will be processed to get y-axis names. "gene_uniprot_site_mult" makes names "gene site Pmult", "gene_uniprot" makes name "gene", but "gene uniprot" for double-up genes, "psite_and_prot" applies both "gene_uniprot_site_mult" and "gene_uniprot", and "keepall" leaves rownames as is
#legend_title: Title for legend
#tile_border_size: Size of white border for tiles. 0.1 is good. If FALSE, no border
#shapes: Shapes to use for dots
#show_shape_legend: Whether or not to show shape in legend
#aspect_ratio_fixed: if true, force tiles to be squares
#clust_cols: If not NULL, these are columns used to cluster the rows of the hmap
fc_sig_heatmap <- function(data,
                          fc_cols_w_pvals = c(),
                          fc_cols_wo_pvals = c(),
                          pval_cols = c(),
                          sig_threshold = 0.05,
                          shape2_pval_cols = c(),
                          order_column,
                          is_decreasing = TRUE,
                          return_df = FALSE,
                          x_axis_names = c(),
                          gap_column_index = FALSE,
                          scale_colours_manually = FALSE,
                          mid_colour_threshold = 1.5,
                          high_colour_threshold = "max",
                          colour_scheme = "rev_RdBu",
                          middle_brewer_white = TRUE,
                          yaxis_naming = "gene_uniprot_site_mult",
                          legend_title = "log2 INS/BASAL",
                          tile_border_size = 0.05,
                           shapes = c(19, 17),
                          show_shape_legend = TRUE,
                          aspect_ratio_fixed = TRUE,
                          clust_cols = NULL){
  
  
  ##Cluster if desired
  if (!is.null(clust_cols)){
    dist_m <- dist(data[, clust_cols])
    hclust <- hclust(dist_m)
    data <- data[hclust$order, ]
    data$order_col <- 1:nrow(data)
    order_column <- "order_col"
  }
  
  
  ##Prep for melts
  
  #Set up data
  
  #yaxis_naming
  
  #Apply naming functions
  if (yaxis_naming == "gene_uniprot_site_mult"){
    data <- gene_uniprot_site_mult_naming(data)
  } else if (yaxis_naming == "gene_uniprot"){
    data <- gene_uniprot_naming(data)
  } else if (yaxis_naming == "psite_and_prot"){
    underscore_counts <- sapply(strsplit(rownames(data), "_"),
                                function(x) length(x))
    psite_indices <- which(underscore_counts == 4)
    prot_indices <- which(underscore_counts == 2)
    data$standard_name <- ""
    if (length(psite_indices) > 0){
      data$standard_name[psite_indices] <- gene_uniprot_site_mult_naming(data[psite_indices, ])$standard_name
    }
    if (length(prot_indices) > 0){
      data$standard_name[prot_indices] <- gene_uniprot_naming(data[prot_indices, ])$standard_name
    }
    
  } else if (yaxis_naming == "keepall"){
    
    data$standard_name <- rownames(data)
  }
  
  
  data <- data[order(data[, order_column],
                     decreasing = is_decreasing == FALSE), ]
  data$standard_name <- factor(data$standard_name,
                               levels = unique(data$standard_name))
  
  ##Gap column
  #Make FC and pval columns with 0 FC (appear white) and non-significant pvalue (no dot)
  #data$gap_FC_column <- 0
  #data$gap_pval_column <- 1
  ##Insert gap column in fc_cols and pval_cols. use gap_column_w_pvals to decide where
  #if (gap_column_w_pvals){
  #  
  #  fc_cols
  #}
  
  ##w_pvals melt
  
  #Check that we have pval_cols
  
  if (length(pval_cols) > 0){
    
    #FC melt
    FC_melt_df_w_pvals <- melt(data,
                               id.vars = "standard_name",
                               measure.vars = fc_cols_w_pvals)
    #p_val melt
    pval_melt_df_w_pvals <- melt(data,
                                 id.vars = "standard_name",
                                 measure.vars = pval_cols)
    pval_melt_df_w_pvals$value <- pval_melt_df_w_pvals$value < sig_threshold
    
    #dot_shape column
    pval_melt_df_w_pvals$dot_shape <- "dot"
    pval_melt_df_w_pvals$dot_shape[which(pval_melt_df_w_pvals$variable %in%
                                           shape2_pval_cols)] <- "shape2"
    pval_melt_df_w_pvals$dot_shape <- factor(pval_melt_df_w_pvals$dot_shape,
                                             levels = c("dot", "shape2"))
    
    #Combine
    post_melt_df_w_pvals <- FC_melt_df_w_pvals
    colnames(post_melt_df_w_pvals)[which(colnames(post_melt_df_w_pvals) == "value")] <- "FC_value"
    post_melt_df_w_pvals$sig_value <- pval_melt_df_w_pvals$value
    post_melt_df_w_pvals$dot_shape <- pval_melt_df_w_pvals$dot_shape
    
    #NAs
    post_melt_df_w_pvals$sig_value[is.na(post_melt_df_w_pvals$sig_value)] <- FALSE
  } else {
    
    post_melt_df_w_pvals <- NULL
    
  }
  
  
  
  ##wo_pvals melt
  
  #Check that we have fc_cols_wo_pvals
  
  if (length(fc_cols_wo_pvals) > 0){
    
    #FC melt
    FC_melt_df_wo_pvals <- melt(data,
                                id.vars = "standard_name",
                                measure.vars = fc_cols_wo_pvals)
    colnames(FC_melt_df_wo_pvals)[which(colnames(FC_melt_df_wo_pvals) == "value")] <- "FC_value"
    #Add FALSE pvals
    FC_melt_df_wo_pvals$sig_value <- FALSE
    FC_melt_df_wo_pvals$dot_shape <- "shape"
    
    #post_melt_df_wo_pvals
    post_melt_df_wo_pvals <- FC_melt_df_wo_pvals
  } else {
    
    post_melt_df_wo_pvals <- NULL
    
  }
  
  
  ##Combine w_pvals and wo_pvals
  
  post_melt_df_combined <- rbind(post_melt_df_w_pvals,
                                 post_melt_df_wo_pvals)
  
  
  ##x-axis names
  
  #check that we have x_axis_names
  if (length(x_axis_names)> 0){
    levels(post_melt_df_combined$variable) <- x_axis_names
  }
  
  
  
  ##Add gap_column
  
  if (gap_column_index != FALSE){
    #Make gap_column_df
    gap_column_df <- post_melt_df_combined[1:length(unique(post_melt_df_combined$standard_name)), ]
    gap_column_df$FC_value <- NA
    gap_column_df$sig_value <- FALSE
    gap_column_df$variable <- ""
    #Calculate new levels
    new_levels <- c(levels(post_melt_df_combined$variable)[1:gap_column_index],
                    "",
                    levels(post_melt_df_combined$variable)[-c(1:gap_column_index)])
    #Add to post_melt_df_combined
    post_melt_df_combined <- rbind(post_melt_df_combined,
                                   gap_column_df)
    #Reorder levels
    post_melt_df_combined$variable <- factor(post_melt_df_combined$variable,
                                             levels = new_levels)
  }
  
  ##Set aspect ratio
  if (aspect_ratio_fixed){
    aspect_ratio_fixed <- coord_fixed(ratio = 1) 
  } else {
    aspect_ratio_fixed <- NULL
  }
  
  
  ##Make colour scaling
  
  #If using max for high
  if (high_colour_threshold == "max"){
    high_colour_threshold <- max(abs(post_melt_df_combined$FC_value),
                                 na.rm = TRUE)
  }
  
  #If using original colours
  if (colour_scheme == "original"){
    
    #If scaling manually:
    if (scale_colours_manually){
      
      gradient <- scale_fill_gradientn(colors = c("#1A09FF",
                                                  "#B28AFF",
                                                  "white",
                                                  "#FF9E82",
                                                  "#FF1C0B"),
                                       values = rescale(c(-high_colour_threshold, 
                                                          -mid_colour_threshold,
                                                          0,
                                                          mid_colour_threshold, 
                                                          high_colour_threshold)),
                                       limits = c(-high_colour_threshold, high_colour_threshold),
                                       name = legend_title)
      
    } else {
      
      gradient <- scale_fill_gradient2(low = "blue",
                                       mid = "white",
                                       high = "red",
                                       name = legend_title)
    }
    
    
    #If using brewer
  } else{
    
    ##Extract brewer colours
    colour_scheme <- strsplit(colour_scheme, "_")[[1]]
    
    if (length(colour_scheme) == 1){
      
      colours <- brewer.pal(9, colour_scheme[1])
    } else if (colour_scheme[1] == "rev"){
      
      colours <- rev(brewer.pal(9, colour_scheme[2]))
    }
    
    if (middle_brewer_white){
      
      colours[5] <- "#FFFFFF"
    }
    
    
    #If scaling manually:
    if (scale_colours_manually){
      
      
      gradient <- scale_fill_gradientn(colours = colours,
                                       values = rescale(c(-high_colour_threshold,
                                                          -0.5*high_colour_threshold -0.5*mid_colour_threshold,
                                                          -mid_colour_threshold,
                                                          -0.5*mid_colour_threshold,
                                                          0,
                                                          0.5*mid_colour_threshold,
                                                          mid_colour_threshold,
                                                          0.5*high_colour_threshold + 0.5*mid_colour_threshold,
                                                          high_colour_threshold)),
                                       limits = c(-high_colour_threshold, high_colour_threshold),
                                       name = legend_title)
      
    } else {
      
      gradient <- scale_fill_gradientn(colours = colours,
                                       values = rescale(c(-high_colour_threshold,
                                                          -0.75*high_colour_threshold,
                                                          -0.5*high_colour_threshold,
                                                          -0.25*high_colour_threshold,
                                                          0,
                                                          0.25*high_colour_threshold,
                                                          0.5*high_colour_threshold,
                                                          0.75*high_colour_threshold,
                                                          high_colour_threshold)),
                                       limits = c(-high_colour_threshold, high_colour_threshold),
                                       name = legend_title)
    }
  }
  
  
  #Set up tile
  if (tile_border_size == FALSE){
    
    tile = geom_tile()
  } else {
    
    tile = geom_tile(colour = "#a5a5a5",
                     size = tile_border_size)
  }
  
  
  #Scale for shapes
  if (show_shape_legend){
    shape_scale <- scale_shape_manual(values = shapes)
  } else {
    shape_scale <- scale_shape_manual(values = shapes,
                                      guide = "none")
  }
  
  ##Heatmap
  output_plot <- ggplot(post_melt_df_combined,
                        aes(x = variable,
                            y = standard_name,
                            fill = FC_value)) + 
    gradient + 
    tile +
    aspect_ratio_fixed +
    geom_point(data = post_melt_df_combined[post_melt_df_combined$sig_value == TRUE,],
               aes(x = variable,
                   y = standard_name,
                   shape = dot_shape),
               size = 0.5) +
    shape_scale +
    #labs(x = "Condition",
    #     y = "Phosphopeptide") +
    theme(panel.grid = element_blank(),
          panel.background = element_blank(),
          axis.text.x = element_text(angle = 90,
                                     hjust = 1,
                                     vjust = 0.5,
                                     colour = "black",
                                     size = 7),
          axis.text.y = element_text(colour = "black",
                                     size = 7),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks = element_blank(),
          plot.title = element_text(hjust = 0.5,
                                    colour = "black",
                                    size = 7),
          legend.text = element_text(colour = "black",
                                     size = 7),
          legend.title = element_text(colour = "black",
                                      size = 7),
          legend.key.size = unit(0.15,
                                 units = "inch"))
  
  if (return_df){
    
    return(post_melt_df_combined)
    
  } else {
    
    return(output_plot)
    
  }
  
}



###Function: Horizontal barplot that shows -log10 pvalues of pathways from ernichment analysis
#data: Data frame containing pvalue (pval_col), where rownames are pathways
enrichment_barplot_pval <- function(data,
                                    pval_col,
                                    colour_col = NULL,
                                    colours = rev(brewer.pal(9, "RdBu")),
                                    limits_from0 = FALSE,
                                    p_cutoff = 0.05,
                                    sig_only = TRUE,
                                    sig_line = FALSE,
                                    xlab = "-log10 adj. p-value"){
  ##Prepare data
  data <- data.frame(data)
  data$nlog10_pval <- -log10(data[, pval_col])
  data$pathway <- rownames(data)
  data <- data[order(data$nlog10_pval, decreasing = FALSE), ] %>%
    mutate(pathway = factor(pathway, levels = unique(pathway)))
  #Filter for signficiant only
  if (sig_only){
    data <- data[which(data[, pval_col] < p_cutoff), ]
  }
  ##Set up signifiance line
  if (sig_line){
    sig_line_plot <- geom_vline(xintercept = -log10(p_cutoff),
                                colour = "black", alpha = 0.5,
                                size = 0.235)
  } else {
    sig_line_plot <- NULL
  }
  ##Set up colouring
  if (!is.null(colour_col)){
    col_plot <- geom_col(width = 0.8)
    max <- max(abs(data[, colour_col]), na.rm = TRUE)
    if (limits_from0){
      scale_fill <- scale_fill_gradientn(colors = colours,
                                         limits = c(0, max))
    } else {
      scale_fill <- scale_fill_gradientn(colors = colours,
                                         limits = c(-max, max))
    }
  } else {
    colour_col <- "pathway"
    col_plot <- geom_col(fill = "black", width = 0.8)
    scale_fill <- NULL
  }
  ##Barplot
  output_plot <- ggplot(data, aes_string(x = "nlog10_pval", y = "pathway",
                                         fill = colour_col)) +
    col_plot +
    sig_line_plot +
    scale_fill +
    comfy_theme(include_yaxis = FALSE) +
    theme(axis.title.y = element_blank()) +
    labs(x = xlab)
  return(output_plot)
}


###Function to make a nice ggvenndiagram.
ggVennDiagram_nice <- function(list, 
                               colours, 
                               fill_limits = NULL,
                               label_alpha = 0, 
                               label = "both"){
  
  return(ggVennDiagram(list, 
                       label_alpha = label_alpha,
                       label = label) +
           scale_fill_gradientn(colours = colours,
                                limits = fill_limits) +
           scale_colour_manual(values = rep("#bebebe",
                                            length(list))) +
           theme(plot.title = element_text(hjust = 0.5)))
}


###FUnction: ggsave as a pdf and png
ggsave_pdfpng <- function(file, width, height){
  ggsave(paste(file, ".pdf", sep = ""),
         width = width, height = height)
  ggsave(paste(file, ".png", sep = ""),
         width = width, height = height)
}

###Function: Save as a pdf and png
save_pdfpng <- function(plot, file, width, height, res = 300){
  pdf(paste(file, ".pdf", sep = ""),
      width = width, height = height)
  plot(plot)
  dev.off()
  png(paste(file, ".png", sep = ""),
      width = width, height = height, units = "in", res = res)
  plot(plot)
  dev.off()
}




















