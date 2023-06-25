####Background
#Store functions for identifying insulin regulation and differences therein

###initiliase
library(tidyverse)
library(purrr)
library(multcomp)

#####Processing#####
#Made quant df
condition_quanter <- function(data,
                              num_cols = phos_num_cols){
  factors <- c("StrainDietIns", "StrainDiet", "DietIns", "Strain", "Diet", "Ins")
  quant_df <- map(unlist(all_levels[factors]),
                  function(level){
                    cols <- num_cols[grep(level, num_cols)]
                    quant_vec <- rowSums(!is.na(data[, cols]))
                    return(quant_vec)
                  }) %>%
    purrr::reduce(cbind) %>% as.data.frame
  colnames(quant_df) <- unlist(all_levels[factors])
  return(quant_df)
}


#####Filtering#####
####Functions
###Function: filter a StrainDiet condition so that bas and ins are quantified more than a threshold
#args: list of bas_quant and ins_quant (number of times ppeptides quantified in bas and ns)
basins_quant_StrainDietfilter <- function(args,
                                          bas_threshold = 2,
                                          ins_threshold = 2){
  filter_bool <- args$bas_quant >= bas_threshold & 
    args$ins_quant >= ins_threshold
  return(filter_bool)
}
###Apply StrainDiet filter
apply_StrainDietfilter <- function(filter_func,
                                   arguments,
                                   ppeptides,
                                   StrainDiet_threshold = 8,
                                   ...){
  #Apply filter
  filter_bool_list <- map(arguments, 
                          function(arg){filter_func(args = arg, ...)})
  filter_bool_df <- purrr::reduce(filter_bool_list, cbind) %>% as.data.frame
  colnames(filter_bool_df) <- names(arguments)
  rownames(filter_bool_df) <- ppeptides
  #Overall filter by combining StrainDiet results
  filter_bool_df$combined <- rowSums(filter_bool_df[, all_levels$StrainDiet]) >=
    StrainDiet_threshold
  return(filter_bool_df)
}

###Function: Take filter_df (columns are bools indicating when StrainDiets have passed a filter)
#and use to make numerical data (num_df) NA when filter not passed
StrainDiet_filter_makeNAs <- function(filter_df, num_df){
  #Align data by rownames
  common_rownames <- intersect(rownames(filter_df), rownames(num_df))
  filter_df <- filter_df[common_rownames, ]
  num_df <- num_df[common_rownames, ]
  #Expand filter_df to match structure of num_df
  filter_df_expanded <- num_df
  conditions <- colnames(filter_df)
  for (col in colnames(filter_df_expanded)){
    condition_match <- map(conditions, ~length(grep(., col)) > 0) %>% unlist
    filter_df_expanded[, col] <- filter_df[, condition_match]
  }
  #Set values to NA when filter not passes
  num_df[!filter_df_expanded] <- NA
  return(num_df)
}


#####Stats#####

###Function: Run tests on phospho data across ppeptides. Tests use individual sample values
#test_func: The test function to use. 
#Takes in data (transposed phos data with ppeptides as columns)  and value_col: the column to test
#results_format: vector of names of the results output
ppeptide_stat_runner <- function(data,
                                 num_cols = UIDs,
                                 test_func,
                                 results_format = c("p", "F")){
  #Get transposed data
  ppeptides <- rownames(data)
  data_t <- transpose_phos_data(data = data,
                                num_cols = num_cols)
  #Loop over ppeptides and apply tests
  #Build possible test function to handle errors
  test_func_pos <- possibly(test_func, otherwise = NA)
  
  test_results <- map(ppeptides,
                      ~test_func_pos(data = data_t,
                                     value_col = .))
  #Turn NULL results into NAs
  test_results[map(test_results, is.null) %>% unlist] <- NA
  #Add a dummy to top of test resutls to trick it indo reducing
  test_results <- c(list(rep(0, length(results_format))), test_results)
  test_results_df <- purrr::reduce(test_results, rbind) %>% as.data.frame %>%
    .[-1, ]
  rownames(test_results_df) <- ppeptides
  colnames(test_results_df) <- results_format
  return(test_results_df)
}


###Function: ppeptide ~ Strain*Diet*Ins, takes pvlaue and F vlaue for Ins main effect
#data: Transposed phos data where column are ppeptides. Must contain columns for Strain, Diet, and Ins
ins_3wayaov <- function(data, value_col){
  #Rename
  data <- dplyr::rename(data, "Value" = value_col)
  #Run model
  ins_aov <- aov(Value ~ Strain*Diet*Ins, data = data)
  summary_stats <- c("p" = summary(ins_aov)[[1]][["Pr(>F)"]][3],
                     "F" = summary(ins_aov)[[1]][["F value"]][3])
  return(summary_stats)
}


###Function: ppeptide ~ Strain*Diet
#data: Transposed phos data where column are ppeptides. Must contain columns for Strain and Diet
strainxdiet_aov <- function(data, value_col){
  #Rename
  data <- dplyr::rename(data, "Value" = value_col)
  #Run model
  ins_aov <- aov(Value ~ Strain*Diet, data = data)
  #Extract results
  ps <- summary(ins_aov)[[1]][[("Pr(>F)")]][1:3]
  names(ps) <- paste(c("Strain", "Diet", "StrainxDiet"), "_p", sep = "")
  Fs <- summary(ins_aov)[[1]][[("F value")]][1:3]
  names(Fs) <- paste(c("Strain", "Diet", "StrainxDiet"), "_F", sep = "")
  return(c(ps, Fs))
}

###Function: ppeptide ~ Strain in CHPOW
#data: Transposed phos data where column are ppeptides. Must contain columns for Strain and Diet
strain_CHOW_aov <- function(data, value_col){
  #Rename
  data <- dplyr::rename(data, "Value" = value_col)
  #Run model
  ins_aov <- aov(Value ~ Strain, data = data)
  #Extract results
  results <- c(summary(ins_aov)[[1]][[("Pr(>F)")]][1],
               summary(ins_aov)[[1]][[("F value")]][1])
  names(results) <- c("p", "F")
  return(results)
}


#Function: ttest that handles errors
ttest_possibly <- possibly(t.test, otherwise = NULL)


#Function: perform ttests comparing each Strain to C57Bl6J
##data: Transposed phos data where column are ppeptides. Must contain columns for Strain
Strain_C57Bl6J_ttests <- function(data, value_col){
  #Loop over Strains that are not C57Bl6J
  other_strains <- levels(data$Strain)[-which(levels(data$Strain) == "C57Bl6J")]
  test_results <- map(other_strains,
                      function(strain){
                        ttest <- ttest_possibly(data[data$Strain == "C57Bl6J", value_col],
                                                data[data$Strain == strain, value_col])
                        #Abort if error
                        if (is.null(ttest)){
                          results <- c("p" = NA, "t" = NA)
                        } else {
                          results <- c("p" = ttest$p.value,
                                       "t" = ttest$statistic)
                          names(results) <- c("p", "t")
                        }
                        return(results)
                      })
  #Rename
  test_results <- map2(test_results, other_strains, 
                       function(results, name){
                         names(results) <- paste(name, names(results), sep = "_")
                         return(results)
                       })
  output <- unlist(test_results)
  return(output)
}

#Function: perform ttests comparing CHOW to HFD in each strain
##data: Transposed phos data where column are ppeptides. Must contain columns for Strain
StrainxDiet_Diet_ttests <- function(data, value_col){
  ##Diet comparisons
  diet_test_results <- map(all_levels$Strain,
                           function(strain){
                             ttest <- ttest_possibly(data[data$Strain == strain &
                                                            data$Diet == "CHOW", value_col],
                                                     data[data$Strain == strain &
                                                            data$Diet == "HFD", value_col])
                             #Abort if error
                             if (is.null(ttest)){
                               results <- c("p" = NA, "t" = NA)
                             } else {
                               results <- c("p" = ttest$p.value,
                                            "t" = ttest$statistic)
                               names(results) <- c("p", "t")
                             }
                             return(results)
                           })
  #Rename and unlist
  diet_test_results <- map2(diet_test_results, all_levels$Strain,
                            function(results, name){
                              names(results) <- paste(name, "_HFD - ",
                                                      name, "_CHOW_",
                                                      names(results), sep = "")
                              return(results)
                            }) %>% unlist
  return(diet_test_results)
}

#####Filtering for insulin-regulated sites#####
###Function: filter for insulin regulation by FC and qvalue cutoff. Generates volcano plot
#suffix: suffix to add to the end of the volcano plot
#if q_col is NULL, then adjusted p-values for you using p_col
#If suffix is NULL, it will not make a volcano plot
insreg_FC_filterer <- function(data,
                               FC_cutoff = 0.58,
                               FC_col,
                               q_cutoff = 0.05,
                               q_col,
                               p_col = NULL,
                               directory = "output/images/analysis/insulin_regulation/insreg/FC_methods/",
                               suffix,
                               ...){
  #Adjust pvalues if not done
  if (is.null(q_col)){
    q_col <- gsub("_p$", "_q", p_col)
    qvalues <- qvalue(data[, p_col])
    data[, q_col] <- qvalues$qvalues
  }
  #Filter
  data$insreg_dir <- "unregulated"
  data$insreg_dir[which(data[, q_col] < q_cutoff &
                          data[, FC_col] > FC_cutoff)] <- "up"
  data$insreg_dir[which(data[, q_col] < q_cutoff &
                          data[, FC_col] < -FC_cutoff)] <- "down"
  data$insreg_bool <- data$insreg_dir != "unregulated"
  data_insreg <- data[data$insreg_bool, ]
  data$insreg_dir <- factor(data$insreg_dir, levels = c("up", "down", "unregulated"))
  #Summary and visualisations
  insreg_summary <- table(data$insreg_dir)
  #Volcano plot, if desired
  if (!is.null(suffix)){
    volcano_plot_coloured(data = data,
                          pval_col = q_col,
                          FC_col = FC_col,
                          colour_col = "insreg_bool",
                          colour_col_levels = c(FALSE, TRUE),
                          colour_col_colours = c("grey", "black"),
                          ...) +
      labs(x = "ins/bas log2FC", y = "-log10 q-value")
    ggsave_pdfpng(paste(directory, "vplot_insreg_", suffix,
                        sep = ""), width = 3, height = 3)
  }
  #Output
  output_list <- list("data" = data,
                      "data_insreg" = data_insreg,
                      "summary" = insreg_summary)
  return(output_list)
}

#Function: Value of a vector with greatest value in the same direction as median
directed_max_func <- function(vec){
  if (sum(!is.na(vec)) == 0){
    return(NA)
  } else {
    #Get sign of summary_val
    median <- median(vec, na.rm = TRUE)
    sign <- sign(median)
    #Swap so values we care about are positive
    vec <- sign*vec
    #Get maximum and resign
    max <- max(vec, na.rm = TRUE)*sign
    return(max)
  }
}

###Function: insreg_FC_filterer that uses directed max of StrainDiet ins/bas logFCs
insreg_FC_filterer_directedMax <- function(data,
                               FC_cutoff = 0.58,
                               FC_columns,
                               q_cutoff = 0.05,
                               q_col,
                               p_col = NULL,
                               directory = "output/images/analysis/insulin_regulation/insreg/FC_methods/",
                               suffix,
                               ...){
  #Make directed max FC
  data$ins_bas_logFC_StrainDietmax <- 
    apply(data[, FC_columns], 1, directed_max_func)
  #Proceed with insreg_FC_filterer
  output_list <- insreg_FC_filterer(data = data,
                                    FC_cutoff = FC_cutoff,
                                    FC_col = "ins_bas_logFC_StrainDietmax",
                                    q_cutoff = q_cutoff,
                                    q_col = q_col,
                                    p_col = p_col,
                                    directory = directory,
                                    suffix = suffix,
                                    ...)
  return(output_list)
}

###function: StrainxDiet anovas for insreg differences. Only extract Diet and SxD terms
#data needs: basnorm_cols, insreg_bool
insregdiff_Diet_strainxdiet_aov <- function(data){
  #Set up numerical columns. basnormalised ins data by default
  num_columns <- colnames(data)[grep("basnorm_\\d+$", colnames(data))]
  #Run stats
  data_regdiff <- ppeptide_stat_runner(subset(data, insreg_bool == TRUE),
                                       num_cols = num_columns,
                                       test_func = strainxdiet_aov,
                                       results_format = c("Strain_p", "Diet_p", "StrainxDiet_p",
                                                          "Strain_F", "Diet_F", "StrainxDiet_F")) %>%
    .[, -grep("Strain_", colnames(.))]
  #q-value adjustment across all columns
  conditions <- c("Diet", "StrainxDiet")
  pcols <- paste(conditions, "_p", sep = "")
  qcols <- paste(conditions, "_q", sep = "")
  sig_cols <- paste(conditions, "_sig", sep = "")
  #qvalues
  qvalues <- qvalue(data_regdiff[, pcols])
  data_regdiff[, qcols] <- qvalues$qvalues
  #Significance
  for (i in 1:length(pcols)){
    data_regdiff[, sig_cols[i]] <- FALSE
    data_regdiff[which(data_regdiff[, qcols[i]] < 0.05), sig_cols[i]] <- TRUE
  }
  #Make a column for which posthoc tests will be performed
  #This is:
  #StrainxDiet always if StrainxDiet is significant. Otherwise:
  #Diet if Diet is signifciant only
  
  #Classify each as S, D, SxD, S+D
  data_regdiff$posthoc_class <- "none"
  data_regdiff$posthoc_class[data_regdiff$Diet_sig] <- "Diet"
  data_regdiff$posthoc_class[data_regdiff$StrainxDiet_sig] <- "StrainxDiet"
  data_regdiff$posthoc_class <- factor(data_regdiff$posthoc_class,
                                       levels = c("none", "StrainxDiet", "Diet"))
  return(data_regdiff)
}


####FC filtering

###Function: return result of a piecewise function specific by points.
#final_slope is the slope of the line before and after final points
piecewise_func <- function(x, points, final_slope){
  #Deal with extremeties
  if (x > points[[1]][1]){
    y <- final_slope*(x - points[[1]][1]) + points[[1]][2]
  } else if (x <= points[[length(points)]][1]){
    y <- final_slope*(x - points[[length(points)]][1]) + points[[length(points)]][2]
  } else {
    #Otherwise loop over points
    for (i in 1:(length(points) - 1)){
      point1 <- points[[i]]
      point2 <- points[[i + 1]]
      #Check if in threshold
      if (x > point2[1] & x <= point1[1]){
        gradient <- (point2[2] - point1[2])/(point2[1] - point1[1])
        y <- gradient*(x - point1[1]) + point1[2]
      }
    }
  }
  return(y)
}

###Functions: Upper and lower boundaries for FC_filter_adaplin
get_upper_boundary <- function(x, buff = 0.58){
  #Divide through by buffer for simplicity
  x <- x/buff
  y <- piecewise_func(x,
                      points = list(
                        c(2, 3),
                        c(0, 1/2),
                        c(-1/2, 0),
                        c(-3, -2)),
                      final_slope = 1)
  return(y*buff)
}
get_lower_boundary <- function(x, buff = 0.58){
  #Divide through by buffer for simplicity
  x <- x/buff
  y <- piecewise_func(x,
                      points = list(
                        c(3, 2),
                        c(1/2, 0),
                        c(0, -1/2),
                        c(-2, -3)),
                      final_slope = 1)
  return(y*buff)
}

###Function: FC filter by the adaptive linear method
#Accept a test_FC if its absoute difference to CTRL_FC is > 0.58. When CTRL_FC is smaller, this is relaxed
FC_filter_adaplin <- function(CTRL_FC, test_FC, ...){
  #Get boundaries
  test_boundary_upper <- get_upper_boundary(CTRL_FC, ...)
  test_boundary_lower <- get_lower_boundary(CTRL_FC, ...)
  #Compare to boundaries
  if (test_FC > test_boundary_upper){
    dir <- "up"
  } else if (test_FC < test_boundary_lower){
    dir <- "down"
  } else {
    dir <- NA
  }
  return(dir)
}

###Function: Compare vectors of FCs from test condition to CTRL conditio
#totest_list: list of boolean vectors indicating whether each observatin should be tested
#testFC_list: list of numerical vectors which are the test FCs
#CTRLFC_list: list of numerical vectors which are the CTRL FCs
#FC_comp_func: Function to compare CTRL and test FCs. Must take arguments CTRL_FC and test_FC
#All lists must be named with the comparisons they correspond to
FC_comparer <- function(totest_list, testFC_list, CTRLFC_list, FC_comp_func){
  #Loop over comparisons
  FC_comp_list <- map(names(totest_list),
                      function(comp){
                        #Loop over ppeptides
                        FC_comp_vec <- NULL
                        for (i in 1:length(totest_list[[comp]])){
                          if (totest_list[[comp]][i]){
                            FC_comp <- FC_comp_func(CTRL_FC = CTRLFC_list[[comp]][i],
                                                    test_FC = testFC_list[[comp]][i])
                            FC_comp_vec[i] <- FC_comp
                          } else {
                            FC_comp_vec[i] <- NA
                          }
                        }
                        return(FC_comp_vec)
                      })
  FC_comp_df <- purrr::reduce(FC_comp_list, cbind) %>% as.data.frame
  #Name columns and add bool
  colnames(FC_comp_df) <- paste(names(totest_list), "_change_dir", sep = "")
  FC_comp_df[, gsub("_dir$", "_bool", colnames(FC_comp_df))] <- 
    !is.na(FC_comp_df[, grep("_dir$", colnames(FC_comp_df))])
  return(FC_comp_df)
}

###Call differences as enhanced or suppressed
#Function: Convert up and down to 1 and -1
updown_to_plusmin <- function(vec){
  plusmin_vec <- rep(0, length(vec))
  plusmin_vec[which(vec == "up")] <- 1
  plusmin_vec[which(vec == "down")] <- -1
  plusmin_vec <- as.numeric(plusmin_vec)
  return(plusmin_vec)
}

###Function: Call insulin regulation differences (columns specified by insregdiff_dir_cols) as 
#enhanced or suppressed based on the overall insulin regulation
#data must contain the column insreg_dir
enhsupp_insreg_caller <- function(data, insregdiff_dir_cols){
  #Set up insreg_dir_plusmin
  insreg_dir_plusmin <- updown_to_plusmin(data$insreg_dir)
  #Loop over columns
  for (col in insregdiff_dir_cols){
    #Turn insregdiff_dir into 1 and -1 
    insregdiff_dir_plusmin <- updown_to_plusmin(data[, col])
    #Multiply columns
    outcome_plusmin <- insregdiff_dir_plusmin*insreg_dir_plusmin
    #Get outcome
    outcome_vec <- rep(NA, length(outcome_plusmin))
    outcome_vec[which(outcome_plusmin == 1)] <- "enhanced"
    outcome_vec[which(outcome_plusmin == -1)] <- "suppressed"
    #ASsign to data
    data[, gsub("_dir", "_enhsupp", col)] <- outcome_vec
  }
  return(data)
}


##Summary columns for overall enhsupp
###Function: Combine multiple columns showing enhanced and suppressed info by uniqueness
enhsupp_combiner <- function(data, columns){
  #Loop over columns
  combined_vec <- map(1:nrow(data), function(index){
    enhsupp_vec <- as.character(data[index, columns])
    #Remove NAs
    enhsupp_vec <- enhsupp_vec[!is.na(enhsupp_vec)]
    #Reduce
    enhsupp_vec <- sort(unique(enhsupp_vec))
    enhsupp_vec_pasted <- paste(enhsupp_vec, collapse = ";")
    if (enhsupp_vec_pasted == ""){enhsupp_vec_pasted <- NA}
    return(enhsupp_vec_pasted)
  }) %>% unlist
  return(combined_vec)
}


###Summaries

###Function: summarise insregdiffs at each combination of conditions (e.g. strains)
#changed_conds_col: Column that contains combinations of changed conditions separated by ;
#enhsupp_template: template where the word BLANK will be substitued for a specific strain to get enhsupp col
insregdiff_combs_summariser <- function(data,
                                        changed_conds_col,
                                        enhsupp_template){
  #Get combinations ready
  comb_table <- table(data[, changed_conds_col]) %>% as.data.frame %>%
    .[order(.$Freq, decreasing = TRUE), ] %>%
    dplyr::rename(combination = Var1) %>% mutate(combination = as.character(combination))
  #Loop over combinations to summarise
  summary_list <- map(comb_table$combination,
                      function(combination){
                        #Summarise for a single combination
                        #Reduce data
                        data <- data[which(data[, changed_conds_col] == combination), ]
                        #Make enhsupp column
                        enhsupp_cols <- strsplit(combination, ";")[[1]] %>%
                          map(~gsub("BLANK", ., enhsupp_template)) %>% unlist
                        data$enhsupp <- enhsupp_combiner(data, enhsupp_cols)
                        #Make summary table
                        summary <- table(data$enhsupp, data$insreg_dir) %>% as.data.frame %>%
                          dplyr::rename(enhsupp = Var1, insreg_dir = Var2) %>% subset(insreg_dir != "unregulated")
                        summary$combination <- combination    
                        return(summary)
                      })
  #Combine
  summary_df <- map(summary_list, ~mutate(., enhsupp = as.character(enhsupp),
                                          insreg_dir = as.character(insreg_dir))) %>%
    purrr::reduce(rbind) %>%
    mutate(enhsupp = factor(enhsupp, levels = c("enhanced", "enhanced;suppressed", "suppressed")),
           insreg_dir = factor(insreg_dir, levels = c("up", "down")),
           combination = factor(combination, levels = comb_table$combination)) %>%
    #Add insregdir_enhsupp
    .[order(.$insreg_dir, .$enhsupp), ] %>%
    mutate(insregdir_enhsupp = paste(insreg_dir, enhsupp, sep = "_")) %>%
    mutate(insregdir_enhsupp = factor(insregdir_enhsupp, levels = unique(insregdir_enhsupp))) %>%
    #Add column for number of conditions changed
    mutate(num_conds = strsplit(as.character(combination), ";") %>% map(length) %>% unlist) %>%
    mutate(num_conds = factor(num_conds, levels = sort(unique(num_conds))))
  return(summary_df)
}

###Function: summarise insregdiffs for each strain
#enhsupp_template: template where the word BLANK will be substitued for a specific strain to get enhsupp col
insregdiff_strains_summariser <- function(data, strains, enhsupp_col_template){
  #Loop over Strains to summarise
  summary_list <- map(strains,
                      function(strain){
                        #Get enhsupp column for change
                        enhsupp_col <- gsub("BLANK", strain, enhsupp_col_template)
                        #Reduce data to specific strain
                        data <- data[!is.na(data[, enhsupp_col]), ]
                        data$enhsupp <- data[, enhsupp_col]
                        #Make summary table
                        summary <- table(data$enhsupp, data$insreg_dir) %>% as.data.frame %>%
                          dplyr::rename(enhsupp = Var1, insreg_dir = Var2) %>% subset(insreg_dir != "unregulated")
                        summary$Strain <- strain    
                        return(summary)
                      })
  #Combine
  summary_df <- map(summary_list, ~mutate(., enhsupp = as.character(enhsupp),
                                          insreg_dir = as.character(insreg_dir))) %>%
    purrr::reduce(rbind) %>%
    mutate(enhsupp = factor(enhsupp, levels = c("enhanced", "enhanced;suppressed", "suppressed")),
           insreg_dir = factor(insreg_dir, levels = c("up", "down")),
           Strain = factor(Strain, levels = strains)) %>%
    #Add insregdir_enhsupp
    .[order(.$insreg_dir, .$enhsupp), ] %>%
    mutate(insregdir_enhsupp = paste(insreg_dir, enhsupp, sep = "_")) %>%
    mutate(insregdir_enhsupp = factor(insregdir_enhsupp, levels = unique(insregdir_enhsupp))) %>%
    #Reorder
    .[order(.$Strain), ]
  return(summary_df)
}


#Function: Barplot of insultion regultaion summary
#labels: Labels for colours
insreg_barplot <- function(data, colour_aes, colours, labels,
                           x_aes = "type"){
  #Get levels that are present
  present_levels_bool <- map(levels(data[, colour_aes]), 
                             ~. %in% data[, colour_aes]) %>% unlist
  present_levels <- levels(data[, colour_aes])[present_levels_bool]
  
  #Plot
  output_plot <- ggplot(data, aes_string(x = x_aes, y = "Freq", fill = colour_aes)) +
    geom_col(position = "stack", width = 0.8) +
    scale_fill_discrete(name = "Insulin regulation",
                        labels = labels[present_levels],
                        type = colours[present_levels]) +
    comfy_theme(include_xaxis = FALSE, rotate_x_text = TRUE) +
    theme(axis.title.x = element_blank()) +
    labs(y = "Phosphopeptides")
  return(output_plot)
}

#Function: insreg_barplot with updown, enhsupp, and updownenhsupp colours
insreg_barplot_3cols <- function(data, x_aes,
                                 filename, 
                                 width, height){
  #updown
  insreg_barplot(data, colour_aes = "insreg_dir",
                 x_aes = x_aes,
                 colour = ins_updown_colours, labels = c("Upregulated", "Downregulated"))
  ggsave_pdfpng(paste(filename, "_updown", sep = ""),
                width = width, height = height)
  #enhsupp
  insreg_barplot(data, colour_aes = "enhsupp",
                 x_aes = x_aes,
                 colour = enhsupp_colours, labels = enhsupp_labels)
  ggsave_pdfpng(paste(filename, "_enhsupp", sep = ""),
                width = width, height = height)
  #updown_enhsupp
  insreg_barplot(data, colour_aes = "insregdir_enhsupp",
                 x_aes = x_aes,
                 colour = updown_enhsupp_colours, labels = updown_enhsupp_labels)
  ggsave_pdfpng(paste(filename, "_updownenhsupp", sep = ""),
                width = width, height = height)
}





