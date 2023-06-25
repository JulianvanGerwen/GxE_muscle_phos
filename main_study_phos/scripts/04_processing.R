###Background
#Here I store code to process phos data e.g. filter, calculate summaries, etc.

###Initialise
library(tidyverse)
library(purrr)

#####Normalisation#####
median_normaliser <- function(image_directory,
                              data,
                              num_cols){
  #Boxplot before normalising
  pdf(paste(image_directory,
            "/boxplot_before_norm.pdf",
            sep = ""))
  boxplot(data[, num_cols], las = 2, pch = 20, cex = 0.5, cex.axis = 0.5)
  dev.off()
  #normalise
  medians <- apply(data[, num_cols], 2, FUN = "median", na.rm = TRUE)
  overall_median <- median(medians)
  diffs <- medians - overall_median
  data[, num_cols] <- sweep(data[, num_cols], 2, FUN = "-", STATS = diffs)
  #Boxplot after normalising
  pdf(paste(image_directory,
            "/boxplot_after_norm.pdf",
            sep = ""))
  boxplot(data[, num_cols], las = 2, pch = 20, cex = 0.5, cex.axis = 0.5)
  dev.off()
  return(data)
}



######Summary statistics e.g. means and FCs#####

###Function: average GxE phos data in each possible conditon
#summary_fun: The summary function to use. "mean", by default
#factors: A list of levels for specific factors e.g. levels for Strain, StrainDiet, etc.
#transposed: indicates if data has been transposed or still needs to be
#ppeptides: If data has been transposed, this specifies the columns that correspond to ppeptides
condition_averager <- function(data,
                               summary_fun = "mean",
                               num_cols,
                               factors = all_levels,
                               transposed = FALSE,
                               ppeptides = NULL){
  #Transpose data if not already done
  if (transposed == FALSE){
    #Get ppeptides
    ppeptides <- rownames(data)
    #Tranpose
    data_t <- transpose_phos_data(data, num_cols = num_cols)
  } else {
    data_t <- data
  }
  
  #Loop over factors
  all_factors_summary_list <- map2(factors,
                                   names(factors),
                                   function(levels, factor){
                                     factor_summary <- map(levels,
                                                           ~apply(data_t[data_t[, factor] == ., ppeptides],
                                                                  2,
                                                                  FUN = summary_fun,
                                                                  na.rm = TRUE)) %>%
                                       purrr::reduce(rbind) %>% as.data.frame
                                     rownames(factor_summary) <- levels
                                     return(factor_summary)
                                   })
  #Turn into data frame
  ppeptide_summary_df <- purrr::reduce(all_factors_summary_list, rbind) %>%
    t %>%
    as.data.frame
  colnames(ppeptide_summary_df) <- paste(colnames(ppeptide_summary_df), summary_fun, sep = "_")
  return(ppeptide_summary_df)
}


###Function: Make logFCs from means. By default, subtracts the first mean from teh second
#comparisons: List of Conditions to compare
#summary_func: The summary function initially used
#data: Data containing columns that are the conditions specified in comparisons pasted to the summary_func, separated by _ e.g. NOD_CHOW_ins_mean
condition_FC_maker <- function(data,
                               summary_fun = "mean",
                               comparisons){
  #Get FCs
  #Loop over comparisons
  FCs_df <- map(comparisons,
                function(comparison){
                  comparison <- paste(comparison, summary_fun, sep = "_")
                  FCs <- data[, comparison[2]] - data[, comparison[1]]
                  return(FCs)
                }) %>%
    #Turn into df
    purrr::reduce(cbind) %>%
    as.data.frame
  rownames(FCs_df) <- rownames(data)
  colnames(FCs_df) <- map(comparisons,
                          ~paste(.[1], "_vs_", .[2], "_logFC", sep = "")) %>%
    unlist
  return(FCs_df)
}

###Function: get means across conditions and logFCs between conditions for GxEphos data
GxEphos_mean_FC_maker <- function(data, 
                                  summary_fun = "mean",
                                  num_cols){
  #Get means
  data_means <- condition_averager(data = data, 
                                   summary_fun = summary_fun,
                                   num_cols = num_cols)
  #Get FCs
  #Set up comparisons
  Ins_comparisons <- map(c("bas", "ins"),
                         ~StrainDietIns_levels[grep(., StrainDietIns_levels)]) %>%
    transpose %>%
    map(unlist)
  Diet_comparisons <- map(c("CHOW", "HFD"),
                          ~StrainDietIns_levels[grep(., StrainDietIns_levels)]) %>%
    transpose %>%
    map(unlist)
  Strain_comparisons <- map(StrainDietIns_levels[-grep("C57Bl6J", StrainDietIns_levels)],
                            function(condition){
                              condition_split <- strsplit(condition, "_")[[1]]
                              B6_condition <- paste(c("C57Bl6J", condition_split[-1]),
                                                    collapse = "_")
                              return(c(B6_condition, condition))
                            })
  all_comparisons <- c(Ins_comparisons, Diet_comparisons, Strain_comparisons)
  data_FCs <- condition_FC_maker(data_means,
                                 comparisons = all_comparisons,
                                 summary_fun = summary_fun)
  return(cbind(data, data_means, data_FCs))
}





#More efficient strain diet summariser
summarise_StrainDiet_efficient <- function(data,
                                           pheno_cols,
                                           summary_func = function(x){mean(x, na.rm = TRUE)},
                                           ...){
  #list of summaries 
  summary_list <- split(data, data$StrainDiet) %>%
    map(~apply(.[, pheno_cols], 2, FUN = summary_func))
  #Make a df
  summary_df <- purrr::reduce(summary_list, rbind) %>%
    as.data.frame
  colnames(summary_df) <- pheno_cols
  rownames(summary_df) <- names(summary_list)
  summary_df[, c("Strain", "Diet")] <- strsplit(rownames(summary_df), "_") %>%
    purrr::reduce(rbind)
  return(summary_df)
}

#Function: Normalise ins data by subtracting mean of bas data
#basins_cols: Columns to summarise
#basins_func: Function used to combine bas and ins data e.g. delta, FC
normalise_StrainDiet_basins <- function(data,
                                        basins_cols,
                                        basins_func = function(bas, ins){ins - bas},
                                        basins_func_name = "delta",
                                        ...){
  #Add function name to end of basins_cols
  basins_cols <- paste(basins_cols, basins_func_name, sep = "_")
  rename_indices <- which(colnames(data) %in% c("Strain", "Diet", "Ins",
                                                "StrainDiet", "StrainIns", "DietIns") == FALSE)
  colnames(data)[rename_indices] <- paste(colnames(data)[rename_indices], 
                                          basins_func_name,
                                          sep = "_")
  #Split data by ins and summarise basal
  basins_summary_list <- split(data, data$Ins) 
  basins_summary_list$bas <- summarise_StrainDiet_efficient(basins_summary_list$bas,
                                                  pheno_cols = basins_cols,
                                                  ...)
  #Make a summarised bas df that matches ins
  basins_summary_list$bas_summary <- map(rownames(basins_summary_list$bas),
                                         function(row){
                                           #get bas data
                                           vector <- as.numeric(basins_summary_list$bas[row, ])
                                           #expand bas data
                                           n <- length(grep(row, rownames(basins_summary_list$ins)))
                                           df <- map(1:n, ~vector) %>% purrr::reduce(rbind) %>% as.data.frame
                                         }) %>%
    purrr::reduce(rbind) %>% as.data.frame
  rownames(basins_summary_list$bas_summary) <- rownames(basins_summary_list$ins)
  colnames(basins_summary_list$bas_summary) <- basins_cols
  
  #Normalise ins using bas
  basins_summary_list$ins_norm <- 
    basins_summary_list$ins[, c("Strain", "Diet", basins_cols)]
  basins_summary_list$ins_norm[, basins_cols] <- 
    basins_func(basins_summary_list$bas_summary[, basins_cols],
                basins_summary_list$ins[, basins_cols])
  return(basins_summary_list$ins_norm)
}

###Function: Normalise individual insulin samples to bas mean, to get insulin response
GxEphos_insresponse_maker <- function(data,
                                      summary_func = function(x){mean(x, na.rm = TRUE)},
                                      ...){
  #Transpose data
  data_t <- transpose_phos_data(data, ...)
  ppeptides <- rownames(data)
  #Get insulin responses (deltas)
  data_delta <- normalise_StrainDiet_basins(data_t,
                                            basins_cols = ppeptides,
                                            summary_func = summary_func,
                                            ...)
  #Untranspose
  colnames(data_delta) <- gsub("_delta", "", colnames(data_delta))
  data_delta <- as.data.frame(t(data_delta[, ppeptides]))
  colnames(data_delta) <- gsub("ins", "ins_basnorm", colnames(data_delta))
  return(data_delta)
}


###Function: CV calculator. Takes log2 data and first exponentiates
CV_calculator <- function(vec){
  vec_exp <- 2**vec
  sd <- sd(vec_exp, na.rm = TRUE)
  mean <- mean(vec_exp, na.rm = TRUE)
  CV <- sd/mean
  return(CV)
}

###Function: Get CVs for all ppeptides in phos data
GxEphos_CVs <- function(data, num_cols){
  CVs <- apply(data[, num_cols],
               1,
               CV_calculator)
  data$CV <- CVs
  return(data)
}

###Function: Number of quantified replicates per ppetpide
GxEphos_quant_calculator <- function(data, num_cols){
  quants <- rowSums(!is.na(data[, num_cols]))
  data$num_quant_samples <- quants
  return(data)
}


#####Sequence window#####
library(seqinr)

###Function: Read a fasta file into a list with uniprots as names
read_fasta_to_list <- function(file){
  fasta <- read.fasta(file, seqtype = "AA", set.attributes = FALSE)
  names(fasta) <- unlist(lapply(names(fasta),
                                function(x) strsplit(x, "\\|")[[1]][2]))
  return(fasta)
}

###Function: Get sequence window of specified length for uniprot_site and fasta database
get_sequence_window <- function(uniprot_fastas,
                                uniprot_site,
                                window_length = 31){
  #Get uniprot and site
  uniprot <- strsplit(uniprot_site, "_")[[1]][1]
  position <- as.numeric(strsplit(uniprot_site, "_[STY]")[[1]][2])
  
  #Pull out fasta
  fasta <- uniprot_fastas[[which(names(uniprot_fastas) == uniprot)]]
  radius <- round((window_length - 1)/2)
  
  #Add _ to both sides equal to the radius, so we deal with sequences that are too short
  fasta <- c(rep("_", radius), fasta, rep("_", radius))
  position <- position + radius
  sequence_window <- paste(fasta[(position - radius):(position + radius)], collapse = "")
  return(sequence_window)
}



