###Background
#Here I store code for mapping annotation database info (e.g. phosphosite plus) into my data

###Intialise
library(tidyverse)
library(purrr)

#####Processing#####
###Function: Condense data so one column becomes only unqiue values, and ohter columns are squished into lists
condense_data_byUID <- function(data,
                                UID_col){
  #Set up columns to collapse
  collapse_cols <- colnames(data)[colnames(data) != UID_col]
  #Prep data by ensuring not a tibble
  data <- as.data.frame(data)
  #List of data for each UID
  UIDs <- unique(data[, UID_col])
  data_list <- map(UIDs,
                   ~data[data[, UID_col] == ., collapse_cols])
  #Condense each df
  data_condensed <- map(data_list,
                        function(data){
                          map(colnames(data),
                              ~data[, .])
                        }) %>%
    transpose %>%
    purrr::reduce(cbind) %>% as.data.frame
  colnames(data_condensed) <- collapse_cols
  rownames(data_condensed) <- UIDs
  data_condensed[, UID_col] <- UIDs
  return(data_condensed)
}

#####Phosphosite plus#####
###Function: Search PSP annotation data with SITE_GRP_IDs
#each SITE_GRP_ID could match to multiple PSP annotation entries. These are concatenated as vectors
#Returns a dataframe where each column is a list of these concatenated entries
psp_annotation_search <- function(SITE_GRP_IDs,
                                  annotation_data,
                                  cols_to_keep){
  #Loop over SITE_GRP_IDs to get list of annotation data
  data_list <- map(SITE_GRP_IDs,
                   ~subset(annotation_data, SITE_GRP_ID == .))
  column_list <- map(cols_to_keep,
                     function(col){
                       map(data_list, ~.[, col])
                     })
  column_list_df <- purrr::reduce(column_list, cbind) %>%
    as.data.frame
  colnames(column_list_df) <- cols_to_keep
  return(column_list_df)
}

###Function: Search multiple PSP annotation databases with SITE_GRP_IDs
psp_annotation_search_mult <- function(SITE_GRP_IDs,
                                       annotation_data_list,
                                       cols_to_keep_list){
  #Loop over annotation databases
  mapped_annotation_list <- map(names(annotation_data_list),
                                function(name){
                                  mapped_annotation_df <- 
                                    psp_annotation_search(SITE_GRP_IDs = SITE_GRP_IDs,
                                                          annotation_data = annotation_data_list[[name]],
                                                          cols_to_keep = cols_to_keep_list[[name]])
                                  #Label colnames with data type
                                  colnames(mapped_annotation_df) <- 
                                    paste(name, colnames(mapped_annotation_df), sep = "_")
                                  return(mapped_annotation_df)
                                })
  mapped_annotation_comb <- purrr::reduce(mapped_annotation_list, cbind)
  mapped_annotation_comb$SITE_GRP_ID <- SITE_GRP_IDs
  return(mapped_annotation_comb)
}

###Function: kinase, regulatory site, and disease site annotation from SITE_GRP_IDs
#PSP databases need to be preloaded and named specifically
psp_kinregdisease_search <- function(SITE_GRP_IDs){
  #Set up
  annotation_data_list <- list("PSPkin" = psp_kin_sub_data,
                               "PSPreg" = psp_reg_site_data,
                               "PSPdisease" = psp_disease_site_data)
  cols_to_keep_list <- list("PSPkin" = c("KINASE", "KIN_ORGANISM", "SUB_ORGANISM", "IN_VIVO_RXN", "IN_VITRO_RXN"),
                            "PSPreg" = c("ORGANISM", "ON_FUNCTION", "ON_PROCESS", "ON_PROT_INTERACT", "ON_OTHER_INTERACT", "PMIDs"),
                            "PSPdisease" = c("DISEASE", "ALTERATION"))
  #Search
  annotation_data <- psp_annotation_search_mult(SITE_GRP_IDs = SITE_GRP_IDs,
                                                annotation_data_list = annotation_data_list,
                                                cols_to_keep_list = cols_to_keep_list)
  #Columns indicating whether annotation is present
  annotation_data$PSPkin_annotated <- map(annotation_data$PSPkin_KIN_ORGANISM, ~length(.) > 0) %>% unlist
  annotation_data$PSPreg_annotated <- map(annotation_data$PSPreg_ORGANISM, ~length(.) > 0) %>% unlist
  annotation_data$PSPdisease_annotated <- map(annotation_data$PSPdisease_DISEASE, ~length(.) > 0) %>% unlist
  return(annotation_data)
}

###Generalised SITE_GRP_ID mapping
#target_data: Data to map into. Must have column PSP_SITE_GRP_ID 
#annotation_data: Data to take annotaiton from. Must have column SITE_GRP_ID
#cols_to_map: columns to map from annotation data to target_data
#col_renames: Names to give columns once mapped
map_by_SITE_GRP_ID <- function(target_data,
                               annotation_data,
                               cols_to_map,
                               col_renames){
  target_data[, col_renames] <- NA
  for (i in 1:nrow(target_data)){
    SITE_GRP_ID <- target_data$PSP_SITE_GRP_ID[i]
    if (!is.na(SITE_GRP_ID) & SITE_GRP_ID %in% annotation_data$SITE_GRP_ID){
      target_data[i, col_renames] <- 
        annotation_data[which(annotation_data$SITE_GRP_ID == SITE_GRP_ID),
                        cols_to_map]
    }
  }
  return(target_data)
}

###Function: Get all regultory sites on kinases that are present in phos_data
#phos_data: Some derivation of phos_data_proc
phosdata_get_kinaseregsites <- function(phos_data,
                                        kinases,
                                        pspdata = psp_data,
                                        pspreg_data = psp_reg_data,
                                        pspkin_sub_data = psp_kin_sub_data_merge){
  #Get all regulatory sites on kinases
  kin_regsites <- psp_get_kinaseregsites(kinases = kinases,
                                         psp_data = pspdata,
                                         psp_reg_data = pspreg_data,
                                         psp_kin_sub_data = pspkin_sub_data)
  #Map regulatory site info to phos_data
  phos_data_reg <- annotation_search_mapback(data = phos_data,
                                             annotation_groups = c("PSPreg"))
  #Map to regulatory sites on kinases
  phos_data_kin_regsites <- map(kin_regsites$list,
                                ~subset(phos_data_reg, 
                                        PSP_SITE_GRP_ID %in% .))
  return(phos_data_kin_regsites)
}





#####Summaries#####
###Function: Count how many ppeptides, psites, and pproteins have annotation
#data must have the columns standard_name, uniprot_site, uniprot
#indicator_cols: Named vector of boolean columns indicating which rows have annotation
#names are incorporated into the summary table to assist aesthetic plotting
annotation_summary_table <- function(data,
                                     indicator_cols){
  #Make total site column
  data$total <- rep(TRUE, nrow(data))
  indicator_cols <- c("Total" = "total",
                      indicator_cols)
  #Loop over indicator cols
  summary_list <- map(indicator_cols,
                      function(col){
                        bool <- data[, col]
                        counts <- map(c("standard_name",
                                        "uniprot_site",
                                        "uniprot"),
                                      ~length(unique(data[, .][bool]))) %>%
                          unlist
                        return(counts)
                      })
  summary_df <- purrr::reduce(summary_list, rbind) %>% as.data.frame
  colnames(summary_df) <- c("Phosphopeptides",
                            "Phosphosites",
                            "Phosphoproteins")
  rownames(summary_df) <- indicator_cols
  summary_df$annotation <- factor(names(indicator_cols),
                                  levels = names(indicator_cols))
  return(summary_df)
}

#####Searching annotation#####

###Function: Search for annotation uisng ppetpides
#annotation_groups: Prefixes for annotation that you want to see. By default uses all. Removes those with no annotation
#ppeptides_restricted: If TRUE, only return ppeptides with annotaiton
annotation_search <- function(ppeptides,
                              annotation_data = phos_data_annotation,
                              annotation_groups = NULL,
                              ppeptides_restricted = TRUE,
                              columns_restricted = TRUE){
  #Get annotation group column names
  #Use all annotation unless specified
  if (is.null(annotation_groups)){
    annotation_groups <- colnames(phos_data_annotation)[grep("_annotated$",
                                                             colnames(phos_data_annotation))] %>%
      gsub("_annotated", "", .)
  }
  annotation_indicators <- paste(annotation_groups, "_annotated", sep = "")
  
  #Subset annotation data by selected ppeptides
  annotation_data <- annotation_data[ppeptides, ]
  #Retain annotation groups with annotation, if desired
  if (columns_restricted){
    annotation_present <- map(annotation_indicators, ~sum(annotation_data[, .]) > 0) %>% unlist
    annotation_groups_present <- annotation_groups[annotation_present]
  } else {
    annotation_groups_present <- annotation_groups
  }
  cols_present <- map(annotation_groups_present, ~grep(., colnames(annotation_data))) %>% unlist %>% unique
  
  #Retain ppeptides with annotation if desired
  if (ppeptides_restricted){
    if (length(annotation_indicators) == 1){
      ppeptides_present <- rownames(annotation_data)[annotation_data[, annotation_indicators] == TRUE]
    } else {
      ppeptides_present <- rownames(annotation_data)[rowSums(annotation_data[, annotation_indicators]) > 0]
    }
    
  } else {
    ppeptides_present <- ppeptides
  }
  #Retain annotation
  annotation_data <- annotation_data[ppeptides_present, cols_present]
  return(annotation_data)
}

###Function: Annotation search, but also maps back into original data
annotation_search_mapback <- function(data,
                              annotation_data = phos_data_annotation,
                              annotation_groups = NULL,
                              columns_restricted = TRUE){
  annotation_data <- annotation_search(ppeptides = rownames(data),
                                       annotation_data = annotation_data,
                                       annotation_groups = annotation_groups,
                                       columns_restricted = columns_restricted,
                                       ppeptides_restricted = FALSE)
  data <- cbind(data, annotation_data)
  return(data)
}

###Function: Search for PSP reg sites on a gene
#data: data to search
#gene: gene to search. Can be truncated. Is searched for in the rownames of data
gene_regsites_search <- function(data, gene){
  #Restrict to sites on the kinase
  data_kin <- data[grep(gene, rownames(data)), ]
  #Search for regulatoru sites
  data_kin <- annotation_search_mapback(data_kin, 
                                        annotation_groups = c("PSPreg")) %>%
    subset(PSPreg_annotated == TRUE)
  return(data_kin)
}














