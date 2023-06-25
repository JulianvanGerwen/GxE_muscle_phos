###Background
#Here I record scripts for processing GTT and echo data, as well as other weight data e.g. ground tissue weights

###Initialise
library(tidyverse)



###Function: Extract the last n weights from a vector that may include NAs at the end
extract_weights <- function(vec,
                            n = 7){
  vec <- as.numeric(vec)
  extracted_vec <- vec[!is.na(vec)] %>%
    rev() %>%
    .[1:n] %>%
    rev()
  return(extracted_vec)
}


###Function: Process weight data, as stored in a monitroing sheet format
#file: Filename
#id_data: ID dataset as output by process_id_data
#n: Number of weights weeks to take. Counts from the back
process_weight_data <- function(file,
                                id_data,
                                n = 7){
  #Read in
  weights_monitoring <- read_csv("data/2022_inbred_weights_JvG.csv") %>%
    rename(tailmark = "Tail mark") %>%
    #Make Cage_tailmark so I can map into id_data
    mutate(cage_tailmark = paste(Cage, tailmark, sep = "_"))
  
  print("Initial read")
  print(weights_monitoring)
  
  #Restrict to mice in id_data and give Mouse column (cull code)
  weights_monitoring <- weights_monitoring %>%
    .[which(.$cage_tailmark %in% id_data$cage_tailmark), ] %>%
    mutate(Mouse = sapply(cage_tailmark,
                          function(x) id_data$Mouse[which(id_data$cage_tailmark == x)])) %>%
    as.data.frame()
  rownames(weights_monitoring) <- weights_monitoring$Mouse
  
  print("Restricted")
  print(weights_monitoring)
  
  #Extract 6 weeks of diet
  weights_extracted <- weights_monitoring %>%
    .[, grep("2022", colnames(.))] %>%
    apply(1,
          function(x) extract_weights(x, n = n)) %>%
    t() %>%
    as.data.frame()
  colnames(weights_extracted) <- paste("bodyweight_", 0:(n - 1), "week", sep = "")
  weights_extracted$bodyweight_final <- weights_extracted[, paste("bodyweight_", (n-1), "week", sep = "")]
  
  print("Processed")
  print(weights_extracted)
  
  weights_extracted$Mouse <- rownames(weights_extracted)
  
  print("Mouse reassigned")
  print(weights_extracted)
  
  
  
  return(weights_extracted)
}




###Function: Calculate GTT AUC
#gluc_data: vector of glucose data
#times: timepoints 
GTT_AUC_calculator <- function(gluc_data,
                               times,
                               AOC = FALSE){
  #Calculate area of the curve of desired
  if (AOC){
    gluc_data <- gluc_data - gluc_data[1]
  }
  trapezoids <- NULL
  for (i in 1:(length(times) - 1)){
    h <- times[i + 1] - times[i]
    trapezoids[i] <- h*(gluc_data[i + 1] + gluc_data[i])/2
  }
  return(sum(trapezoids))
}

###Function: Calculate matsuda index
#GTT_gluc: GTT glucose data (mM) ordered by timepoint
#GTT_ins: GTT insulin data (ng/mL) ordered by timepoint (0min, 15min)
matsuda_calculator <- function(GTT_gluc, GTT_ins){
  #Convert to appropriate units (glucose: mg/dL, insulin: mU/L)
  GTT_gluc <- GTT_gluc*18
  GTT_ins <- GTT_ins*24.8
  #Calculate matsuda
  GTT_gluc_mean <- mean(GTT_gluc)
  GTT_ins_mean <- mean(GTT_ins)
  denom <- sqrt(GTT_gluc[1]*GTT_ins[1]*GTT_gluc_mean*GTT_ins_mean)
  matsuda <- 10000/denom
  return(matsuda)
}

###Function: Calculate HOMAIR
#GTT_gluc: GTT glucose data (mM) ordered by timepoint
#GTT_ins: GTT insulin data (ng/mL) ordered by timepoint (0min, 15min)
HOMAIR_calculator <- function(GTT_gluc, GTT_ins){
  #Conver to appropriate units (insulin: mU/L)
  GTT_ins <- GTT_ins*24.8
  #Caluclate homa-ir
  HOMAIR <- GTT_gluc[1]*GTT_ins[1]/22.5
  return(HOMAIR)
}


###Function: Process GTT and echo data
#bodyweight_data: Data of bodyweights with a column bodyweight_final which is used for adiposity. Rownames must be Mouse code
process_GTTecho_data <- function(file,
                                 bodyweight_data){
  #Load in and set up colnames
  GTTecho_data <- read_csv(file) 
  #GTT AUCs
  GTT_labs <- c("postdiet")
  for (lab in GTT_labs){
    #Get glucose and insulin data
    GTTgluc_data <- GTTecho_data %>%
      .[, grep(lab, colnames(.))] %>%
      .[, grep("GTTgluc", colnames(.))] %>%
      .[, grep("min$", colnames(.))]
    GTTins_data <- GTTecho_data %>%
      .[, grep(lab, colnames(.))] %>%
      .[, grep("GTTins", colnames(.))] %>%
      .[, grep("min$", colnames(.))]
    
    #Calculate AUC, AOC
    temp_AUC <- apply(GTTgluc_data, 1, function(x){GTT_AUC_calculator(as.numeric(x), c(0, 15, 30, 45, 60, 90))})
    GTTecho_data[, paste(lab, "_GTTgluc_AUC", sep = "")] <- temp_AUC
    temp_AOC <- apply(GTTgluc_data, 1, function(x){GTT_AUC_calculator(as.numeric(x), c(0, 15, 30, 45, 60, 90), AOC = TRUE)})
    GTTecho_data[, paste(lab, "_GTTgluc_AOC", sep = "")] <- temp_AOC
    #Calculate matsuda
    #Gluc and insulin data to list
    GTTgluc_list <- map(1:nrow(GTTgluc_data), ~as.numeric(GTTgluc_data[., ]))
    GTTins_list <- map(1:nrow(GTTins_data), ~as.numeric(GTTins_data[., ]))
    #Get matsuda and HOMAIR
    temp_matsuda <- map2(GTTgluc_list, GTTins_list, matsuda_calculator) %>% unlist()
    GTTecho_data[, paste(lab, "_matsuda", sep = "")] <- temp_matsuda
    temp_HOMAIR <- map2(GTTgluc_list, GTTins_list, HOMAIR_calculator) %>% unlist()
    GTTecho_data[, paste(lab, "_HOMAIR", sep = "")] <- temp_HOMAIR
    
  }
  #Adiposity and leanosity
  GTTecho_data$bodyweight_final <- bodyweight_data[GTTecho_data$Mouse, "bodyweight_final"]
  GTTecho_data$bodyweight_prediet <- bodyweight_data[GTTecho_data$Mouse, "bodyweight_0week"]
  GTTecho_data$postdiet_adiposity <- GTTecho_data$postdiet_fat / GTTecho_data$bodyweight_final * 100
  GTTecho_data$postdiet_leanosity <- GTTecho_data$postdiet_lean / GTTecho_data$bodyweight_final * 100
  #Remove bodyweight cols
  GTTecho_data <- GTTecho_data %>%
    .[, colnames(.) %in% c("bodyweight_final", "bodyweight_prediet") == FALSE]
  
  #Set rownames as Mouse
  GTTecho_data <- as.data.frame(GTTecho_data)
  rownames(GTTecho_data) <- GTTecho_data$Mouse
  return(GTTecho_data)
}

###Function: Process ground tissue weight data
process_tissueweight_data <- function(file){
  tissueweight_data <- read_csv(file) %>%
    rename(date_ground = "Date ground", mass = "Ground mass (mg)") %>%
    mutate(Tissue = as.factor(Tissue),
           date_ground = as.factor(date_ground)) %>%
    pivot_wider(names_from = Tissue,
                values_from = c(date_ground, mass),
                names_sep = "__") %>%
    rename_with(function(x){
      x_split <- strsplit(x, "__")
      x_pasted <- map(x_split, ~ paste(rev(.), collapse = "_")) %>%
        unlist()
      return(x_pasted)
    }) %>%
    as.data.frame()
  rownames(tissueweight_data) <- tissueweight_data$Mouse
  return(tissueweight_data)
}

###Function: Process prot concentration data. Mainly adds yield
#tissues: a vector of tissues
process_prot_conc_data <- function(prot_conc,
                                   tissueweight_data,
                                   tissues){
  #Calculate yield
  yield_data <- map(tissues,
                    function(tissue){
                      yield <- NULL
                      for (i in 1:nrow(prot_conc)){
                        yield[i] <- as.numeric(prot_conc[i, paste(tissue, "_prot_amount", sep = "")])/
                          as.numeric(tissueweight_data[which(tissueweight_data$Mouse == prot_conc$Mouse[i]), 
                                                       paste(tissue, "_mass", sep = "")])
                      }
                      return(yield)
                    }) %>%
    purrr::reduce(cbind) %>% as.data.frame()
  colnames(yield_data) <- paste(tissues, "_yield", sep = "")
  prot_conc <- cbind(prot_conc, yield_data)
  #Make rownames
  rownames(prot_conc) <- prot_conc$Mouse
  #Return
  return(prot_conc)
}



