###Background
#Here I store functions for correlating phospho data with GU data

###Initialise
library(tidyverse)
library(purrr)


#####Data processing#####
###Function: Combine phospho  data with glucose uptake data in insulin-stimulated mice
combine_phos_GU_data <- function(phos_data, pheno_data,
                                  phos_rows,
                                  phos_num_cols,
                                 phos_UIDs,
                                 ins_only = TRUE){
  if (ins_only){
    pheno_data <- subset(pheno_data, Ins == "ins")
  }
  data <- combine_phos_pheno_data(phos_data = phos_data,
                                  pheno_data = pheno_data,
                                  pheno_cols = c("sol_p2DG_weightnorm_rate"),
                                  phos_rows = phos_rows,
                                  phos_num_cols = phos_num_cols,
                                  phos_UIDs = phos_UIDs)
  return(data)
}

#####Performing correlations#####

#Function: correlate phosphorylation with GU. Returns correlation r and p, and linear model B and p
#phos_GU_data: combined phospho and GU data, as output by combine_phosbasnorm_GU_data
phos_GU_correlator <- function(phos_GU_data,
                               ppeptide,
                               cor_method = "pearson"){
  #Correlation
  cor_test <- cor.test(phos_GU_data[, ppeptide], phos_GU_data$sol_p2DG_weightnorm_rate,
                       method = cor_method)
  #Linear model
  colnames(phos_GU_data)[colnames(phos_GU_data) == ppeptide] <- "ppeptide"
  lin_mod <- lm(sol_p2DG_weightnorm_rate ~ ppeptide, data = phos_GU_data)
  colnames(phos_GU_data)[colnames(phos_GU_data) == "ppeptide"] <- ppeptide
  #output
  cor_output <- c(cor_test$estimate, cor_test$p.value,
                  summary(lin_mod)$coefficients[2, 1])
  names(cor_output) <- c("cor_r", "cor_p",
                         "lm_beta")
  return(cor_output)
}

#Function: Correlate phospho  data with GU data for all ppetpides in data
#phos_GU_data: combined phospho and GU data, as output by combine_phosbasnorm_GU_data
phos_GU_correlator_mult_fromcomb <- function(phos_GU_data, ppeptides,
                                             cor_method = "pearson"){
  #Run tests
  tests <- map(ppeptides, ~phos_GU_correlator(phos_GU_data = phos_GU_data,
                                              ppeptide = .,
                                              cor_method = cor_method)) %>%
    purrr::reduce(rbind) %>% as.data.frame
  rownames(tests) <- ppeptides
  #Adjustt pvalues
  qvalue_possibly <- possibly(qvalue, otherwise = list("qvalues" = NA))
  tests$cor_q <- qvalue_possibly(tests$cor_p)$qvalues
  #Order
  tests <- tests[order(abs(tests$cor_r), decreasing = TRUE), ]
  return(tests)
}


#Function: Correlate phospho basnorm data with GU data for all ppetpides in data
#phos_data: phos data containing basnorm cols
#pheno_data: Pheno_data containing sol_p2DG_weightnorm_rate as a column
phos_GU_correlator_mult <- function(phos_data, pheno_data,
                                   phos_num_cols, phos_UIDs,
                                   ins_only = TRUE,
                                   cor_method = "pearson"){
  #Set up
  ppeptides <- rownames(phos_data)
  #Combine data
  phos_GU <- combine_phos_GU_data(phos_data = phos_data, 
                                         pheno_data = pheno_data,
                                         phos_rows = ppeptides,
                                  phos_num_cols = phos_num_cols,
                                  phos_UIDs = phos_UIDs,
                                  ins_only = ins_only)
  #Run tests
  tests <- phos_GU_correlator_mult_fromcomb(phos_GU_data = phos_GU, ppeptides = ppeptides,
                                            cor_method = cor_method)
  return(tests)
}



#####Visualisation#####
###Function: Correlation plot of phospho basnorm data with soleus GU
corrplot_phosbasnorm_solGU <- function(phos_data,
                                       pheno_data,
                                       ppeptide,
                                       num_cols = phos_num_cols_basnormfilt,
                                       phos_UIDs = gsub("basnormfilt_", "", phos_num_cols_basnormfilt),
                                       ...){
  output_plot <- corrplot_phos_pheno(phos_data = phos_data,
                                     pheno_data = subset(pheno_data, Ins == "ins"),
                                     phos_row = ppeptide,
                                     pheno_col = "sol_p2DG_weightnorm_rate",
                                     swap_axes = FALSE,
                                     num_cols = num_cols,
                                     phos_UIDs = phos_UIDs,
                                     include_basins = FALSE,
                                     ...) +
    labs(y = "2DG uptake (Ki/mg tissue)")
  return(output_plot)
}

###Function: Correlation plot of phospho ins data with soleus GU
corrplot_phosins_solGU <- function(phos_data,
                                       pheno_data,
                                       ppeptide,
                                       ...){
  output_plot <- corrplot_phos_pheno(phos_data = phos_data,
                                     pheno_data = subset(pheno_data, Ins == "ins"),
                                     phos_row = ppeptide,
                                     pheno_col = "sol_p2DG_weightnorm_rate",
                                     swap_axes = FALSE,
                                     num_cols = phos_num_cols,
                                     phos_UIDs = phos_num_cols,
                                     ...) +
    labs(y = "2DG uptake (Ki/mg tissue)")
  return(output_plot)
}

###Function: Correlation plot of phospho ins data with soleus GU
corrplot_phos_solGU <- function(phos_data,
                                   pheno_data,
                                   ppeptide,
                                   ...){
  output_plot <- corrplot_phos_pheno(phos_data = phos_data,
                                     pheno_data = pheno_data,
                                     phos_row = ppeptide,
                                     pheno_col = "sol_p2DG_weightnorm_rate",
                                     swap_axes = FALSE,
                                     num_cols = phos_num_cols,
                                     phos_UIDs = phos_num_cols,
                                     ...) +
    labs(y = "2DG uptake (Ki/mg tissue)")
  return(output_plot)
}














