###Background
#Functions for processing tissue including calculating rate from blood tracer curve

#Initialise
library(tidyverse)




###Function: Tissue data to wide
#tissue_long must have columsn: Mouse, Ins, Tissue, p2DG_DPM. The latter is p2DG_DPM normalised to mg counted
tissue_data_towide <- function(tissue_long){
  tissue_wide <-  tissue_long %>%
    #Rename tissues as tissue_DPM
    mutate(Tissue = paste(Tissue, "DPM", sep = "_")) %>%
    #To wider
    pivot_wider(names_from = Tissue, values_from = p2DG_DPM)
  return(tissue_wide)
}




###Function: Calculate rate of tissue uptake from blood data and tissue DPMs
#blood_logintercept: logged yintercept of tracer disapperance from blood
#blood_rate: Rate of tracer disappearance from blood
#Time: Time point at which tissue was taken out
#tissue_DPM: Tissue p2DG DPM/mg of counted tissue
tissue_rate <- function(blood_logintercept,
                        blood_rate,
                        time = 10,
                        tissue_DPM){
  tissue_rate <- (tissue_DPM * blood_rate) /
    (exp(blood_logintercept) * (1 - exp(-blood_rate*time)))
  return(tissue_rate)
}

###Function: Calculate rate of tissue uptake from blood data and tissue DPMs. Multiple tissues
#blood_data: Processed data as output by blood_logandfit
#tissue_data: Wide tissue data as output by tissue_data_towide
#time: time point at which tissue was taken out
mult_tissue_rate <- function(tissue_data,
                             blood_data_intercept,
                             blood_data_rate,
                             time = 10){
  tissue_cols <- colnames(tissue_data)
  tissue_rate_cols <- paste(tissue_cols, "_rate", sep = "")
  tissue_data[, tissue_rate_cols] <- NA
  for (i in 1:length(tissue_cols)){
    tissue_data[, tissue_rate_cols[i]] <- 
      tissue_rate(tissue_DPM = tissue_data[, tissue_cols[i]],
                  blood_logintercept = blood_data_intercept,
                  blood_rate = blood_data_rate,
                  time = time)
  }
  return(tissue_data[, tissue_rate_cols])
}



#Average replicates
###Function to average replicates while preserving other data
average_replicates <- function(data){
  data[, c("CPM_avg", "DPM_avg")] <- NA
  data$keep <- FALSE
  for (i in 1:nrow(data)){
    #Get averages
    Mouse_ID <- data$Mouse[i]
    temp <- apply(data[data$Mouse == Mouse_ID, c("CPM", "DPM")],
                  2, FUN = "mean")
    data[i, c("CPM_avg", "DPM_avg")] <- as.list(temp)
    #Keep row if it is the first with that mouse ID
    if (i == sort(which(data$Mouse == Mouse_ID))[1]){
      data$keep[i] <- TRUE
    }
  }
  data <- data %>%
    subset(keep == TRUE) %>%
    select(-c(keep, CPM, DPM)) %>%
    rename(CPM = CPM_avg, DPM = DPM_avg)
  return(data)
}

###Function: Turn raw tissue count data into p2DG data
tissuecounts_raw_to_pD2G <- function(data,
                                     tissue = "sol"){
  #Rename
  data <- data %>%
    rename(total_or_flow = "Total or flow",
           counted_vol = "Volume counted (uL)",
           lysate_vol = "Lysate volume (uL)",
           count_time = "Count Time",
           CPM = "CPMA",
           DPM = "DPM1") %>%
    mutate(total_or_flow = factor(total_or_flow, levels = c("total", "flow")))
  #Average replicates
  data <- split(data, data$total_or_flow) %>%
    map(average_replicates) %>%
    purrr::reduce(rbind)
  #Correct for total lysate volume
  data <- mutate(data,
                 DPM_totallysate = (DPM/counted_vol)*lysate_vol)
  #Get p2DG
  data_p2DG <- data %>%
    select(Mouse, total_or_flow, DPM_totallysate) %>%
    pivot_wider(names_from = total_or_flow,
                values_from = DPM_totallysate) %>%
    mutate(p2DG = total - flow) %>%
    as.data.frame
  #Remame
  rownames(data_p2DG) <- data_p2DG$Mouse
  torename_cols <- which(colnames(data_p2DG) != "Mouse")
  colnames(data_p2DG)[torename_cols] <- paste(tissue, colnames(data_p2DG)[torename_cols], sep = "_")
  #Return
  return(data_p2DG)
}

###Function: normalise tissue count data
#data_p2DG: tissue count data as output by tissuecounts_raw_to_pD2G
tissuecounts_normalise <- function(data_p2DG,
                                   bloodcount_data,
                                   tissueweight_data,
                                   prot_conc_data,
                                   tissue = "sol"){
  #Combine data
  data_list <- list("bloodcount_data" = bloodcount_data,
                    "tissueweight_data" = tissueweight_data,
                    "prot_conc_data" = prot_conc_data,
                    "data_p2DG" = data_p2DG)
  data_comb <- map(data_list,
                   function(data){
                     data <- data[rownames(data_p2DG), ] %>%
                       .[, colnames(.) != "Mouse"]
                   }) %>%
    purrr::reduce(cbind)
  #Normalise
  #To tissue or protein
  p2DG_col <- paste(tissue, "p2DG", sep = "_")
  p2DG_weightnorm_col <- paste(tissue, "p2DG_weightnorm", sep = "_")
  p2DG_protnorm_col <- paste(tissue, "p2DG_protnorm", sep = "_")
  data_comb[, p2DG_weightnorm_col] <- data_comb[, p2DG_col] / data_comb[, paste(tissue, "mass", sep = "_")]
  data_comb[, p2DG_protnorm_col] <- data_comb[, p2DG_col] / data_comb[, paste(tissue, "prot_amount", sep = "_")]
  #Using blood counts
  tissue_rates <- mult_tissue_rate(tissue_data = data_comb[, c(p2DG_weightnorm_col, p2DG_protnorm_col)],
                                   blood_data_intercept = data_comb$bloodcounts_intercept,
                                   blood_data_rate = data_comb$bloodcounts_rate,
                                   time = 10)
  data_comb <- cbind(data_comb, tissue_rates)
  #Return
  return(data_comb)
}


