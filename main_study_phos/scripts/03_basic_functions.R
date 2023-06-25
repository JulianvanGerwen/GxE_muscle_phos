###Background
#Here I store simple, broad-use functions for phospho data

###Initialise
library(tidyverse)
library(purrr)


#####Data processing#####


###Function: Combine phospho and pheno data
#phos_num_cols: Numerical cols of phospho to use. Must be mouse IDs in the form "NOD_CHOW_ins_1"
combine_phos_pheno_data <- function(phos_data,
                                    pheno_data,
                                    pheno_cols,
                                    phos_rows,
                                    phos_num_cols = UIDs,
                                    phos_UIDs = NULL){
  ##Transpose and combine
  phos_data_prep <- t(phos_data[phos_rows, phos_num_cols]) %>% as.data.frame
  pheno_data_prep <- pheno_data
  rownames(pheno_data_prep) <- pheno_data_prep$UID
  pheno_data_prep <- dplyr::select(pheno_data_prep, pheno_cols) %>%
    as.data.frame()
  #Set phospho UIDs if not specified
  if (!is.null(phos_UIDs)){
    rownames(phos_data_prep) <- phos_UIDs
  }
  #Get common rownames
  common_UIDs <- intersect(rownames(phos_data_prep), rownames(pheno_data_prep))
  #Trim down and rename columns in case they have gotten lost (happens when only one col)
  phos_data_common <- as.data.frame(phos_data_prep[common_UIDs, ])
  colnames(phos_data_common) <- phos_rows
  pheno_data_common <- as.data.frame(pheno_data_prep[common_UIDs, ])
  colnames(pheno_data_common) <- pheno_cols
  phos_pheno_data <- cbind(phos_data_common,
                           pheno_data_common)
  rownames(phos_pheno_data) <- common_UIDs
  #Add factors
  phos_pheno_data <- extract_factors_from_rownames(phos_pheno_data)
  return(phos_pheno_data)
}

###Function: Transpose numerical phos data so columns are ppeptides and rows are mice
transpose_phos_data <- function(data, num_cols = UIDs){
  data_t <- as.data.frame(t(data[, num_cols]))
  data_t <- extract_factors_from_rownames(data_t)
  return(data_t)
}

#Function: Map columns from a smaller dataframe into a larger one of which its rows are a subset
#data_source: smaller data that you are mapping from
#data_target: Larger data that you are mapping into
#cols_to_map: columns from the smaller data that you want to map
#factors_to_preserve: List of length 2 vectors containing factors and base_level. 
#These indicate factors where you want to prescribe a specific level (base_level) instead of NA
map_data_into_larger <- function(data_source,
                                 data_target,
                                 cols_to_map,
                                 factors_to_preserve = NULL){
  ##Perform mapping
  #Trim data_source to not have more rows than data_target
  common_rows <- intersect(rownames(data_source), rownames(data_target))
  data_source <- data_source[common_rows, ]
  #Expand data_source to have same rows as data_target
  data_source[setdiff(rownames(data_target), rownames(data_source)), ] <- NA
  data_source <- data_source[rownames(data_target), ]
  #Map in
  data_target[, cols_to_map] <- NA
  data_target[, cols_to_map] <- data_source[, cols_to_map]
  #Preserve boolean columns and factors
  if (!is.null(factors_to_preserve)){
    factors_to_preserve_list <- transpose(factors_to_preserve) %>% map(unlist)
    names(factors_to_preserve_list) <- c("factors", "base_level")
  } else {
    factors_to_preserve_list <- NULL
  }
  for (col in cols_to_map){
    col_type <- typeof(data_target[, col])
    #Preserve factors if desired
    if (col %in% factors_to_preserve_list$factors){
      base_level <- factors_to_preserve_list$base_level[which(factors_to_preserve_list$factors == col)]
      data_target[is.na(data_target[, col]), col] <- base_level
    } else if (col_type == "logical"){
      data_target[is.na(data_target[, col]), col] <- FALSE
    }
  }
  return(data_target)
}

#####Summary functions#####
###Function: Cross summary of data with two factors, including totals where one factor is summed along levels of the other
summariser_twoway <- function(data, Var1, Var2){
  summary_df <- table(data[, Var1],
                      data[, Var2]) %>% as.data.frame
  #Add totals
  var1_total <- tapply(summary_df$Freq, summary_df$Var1, sum)
  var1_total_df <- data.frame("Var1" = names(var1_total), "Var2" = "total", "Freq" = var1_total)
  var2_total <- tapply(summary_df$Freq, summary_df$Var2, sum)
  var2_total_df <- data.frame("Var1" = "total", "Var2" = names(var2_total), "Freq" = var2_total)
  total_df <- rbind(var1_total_df, var2_total_df)
  summary_df_totals <-  rbind(summary_df, var1_total_df, var2_total_df)
  #Rename
  colnames(summary_df_totals) <- c(Var1, Var2, "Freq")
  return(summary_df_totals)
}

###Function: Cross summary of data with two factors. For each Var1 levels, includes proprotion of Var2 levels within
summariser_twoway_prop <- function(data, Var1, Var2){
  temp_table <- table(data[, Var1], data[, Var2]) %>% as.data.frame
  #Get proportions
  Var1_totals <- group_by(temp_table, Var1) %>% summarise(total = sum(Freq))
  temp_table$prop_of_Var1 <- 0
  for (level in levels(temp_table$Var1)){
    temp_table[temp_table$Var1 == level, "prop_of_Var1"] <- temp_table[temp_table$Var1 == level, "Freq"]/
      as.numeric(Var1_totals[Var1_totals$Var1 == level, "total"])
  }
  Var2_totals <- group_by(temp_table, Var2) %>% summarise(total = sum(Freq))
  temp_table$prop_of_Var2 <- 0
  for (level in levels(temp_table$Var2)){
    temp_table[temp_table$Var2 == level, "prop_of_Var2"] <- temp_table[temp_table$Var2 == level, "Freq"]/
      as.numeric(Var2_totals[Var2_totals$Var2 == level, "total"])
  }
  #Rename
  colnames(temp_table) <- gsub("Var1", Var1, colnames(temp_table))
  colnames(temp_table) <- gsub("Var2", Var2, colnames(temp_table))
  return(temp_table)
}














