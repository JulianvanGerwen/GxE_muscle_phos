###Background
#This script contains functions for loading and combining phenotype data

###Initialising
library(tidyverse)


####Id data

###Function: Load id data into R
#directory: Where data is stored
#columns: Columns to select
process_id_data <- function(directory,
                            columns = c("Cull code", "Strain", "Diet", "Bas or Ins", "Cage", "Tail mark")){
  to_select <- unique(c("Cull code", "Strain", "Diet", "Bas or Ins", "Cage", "Tail mark", columns))
  id_data <- read_csv(directory) %>%
    select(to_select) %>%
    rename("Mouse" = "Cull code", "Ins" = "Bas or Ins", "Tailmark" = "Tail mark") %>%
    mutate(Strain = factor(Strain, levels = c("C57Bl6J", "NOD", "BXH9", "CAST", "BXD34")), 
           Diet = as.factor(Diet), 
           Ins = as.factor(Ins),
           cage_tailmark = paste(Cage, Tailmark, sep = "_")) %>%
    #Order by Strain, Diet, and Ins
    .[order(.$Strain, .$Diet, .$Ins, .$Mouse), ] %>%
    #Make Condition column
    mutate(Condition = as.factor(paste(Strain, Diet, Ins, sep = "_"))) %>%
    mutate(Condition = factor(Condition, levels = unique(Condition))) %>%
    as.data.frame()
  #Make columns for StrainDiet etc
  for (comb in c("Strain_Diet", "Strain_Ins", "Diet_Ins")){
    factors <- strsplit(comb, "_")[[1]]
    col <- paste(factors, collapse = "")
    id_data[, col] <- factor(paste(id_data[, factors[1]], id_data[, factors[2]],
                                   sep = "_"),
                             levels = all_levels[[col]])
  }
  
  rownames(id_data) <- id_data$Mouse
  return(id_data)
}




####Combining data

###Function: Combine phenotype data
#id_data: Data containing identifiers for mice - strain, diet, ins/bas, code. Rownames are code
#phenotype_data_list: list of phenotypic data where Rownames are code
combine_phenotype_data <- function(id_data, phenotype_data_list){
  ##Prep id_data
  id_data <- id_data[order(id_data$Mouse), ]
  id_rownames <- rownames(id_data)
  ##Order phenotype data so same as id
  phenotype_data_list <- lapply(phenotype_data_list,
                                function(x) x[id_rownames, ])
  ##Combine
  combined_data <- id_data
  for (i in 1:length(phenotype_data_list)){
    temp_data <- phenotype_data_list[[i]] %>%
      .[, colnames(.) %in% colnames(combined_data) == FALSE]
    combined_data <- cbind(combined_data, temp_data)
  }
  return(combined_data)
}
