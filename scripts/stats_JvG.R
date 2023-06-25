###Background
#Here I keep functions for stats


###Packages
library(multcomp)
library(qvalue)
library(rlist)


#####ANOVAs and posthoc tests#####
###Function: Extract pvals and t statistics from a posthoc test
#posthoc_results: Object output by glht
#comparisons: Comparisons of the form "level1 - level2". Must contain all comparisons present in posthoc_results
posthoc_extract_stats <- function(posthoc_results,
                                  comparisons){
  posthoc_summary <- summary(posthoc_results)
  #Get tstats and comparisons used (sometimes shorter than what you expect)
  ts <- posthoc_summary$test$tstat
  results_comparisons <- names(ts)
  #Get pvalues
  pvals <- posthoc_summary$test$pvalues[1:length(results_comparisons)]
  names(pvals) <- results_comparisons
  #Turn into vector
  summary_m <- matrix(NA, ncol = length(comparisons), nrow = 2)
  colnames(summary_m) <- comparisons
  rownames(summary_m) <- c("p", "t")
  for (col in results_comparisons){
    summary_m[, col] <- c(pvals[col], ts[col])
  }
  summary_vec <- map(c("p", "t"),
                     function(stat){
                       vec <- summary_m[stat, ]
                       names(vec) <- paste(colnames(summary_m), stat, sep = "_")
                       return(vec)
                     }) %>% unlist
  return(summary_vec)
}

#####pvalue adjusstment#####
###Function: Adjust multiple columns of pvalues by qvalue method. Adjusts each column separately
qval_adjustment_mult <- function(data, pcols, return_sig = TRUE){
  #Get qvalues
  qvalues_df <- map(pcols, function(pcol){
    #Proceed if not all NA
    if (sum(!is.na(data[, pcol])) > 0){
      #Histogram
      hist(data[, pcol], breaks = 20)
      #Adjust
      qvalues <- qvalue(data[, pcol])
      return(qvalues$qvalues)
    } else {
      return(data[, pcol])
    }
   
  }) %>% purrr::reduce(cbind) %>% as.data.frame
  colnames(qvalues_df) <- gsub("_p", "_q", pcols)
  data <- cbind(data, qvalues_df)
  #Add signifiance levels if desired
  if (return_sig){
    for (qcol in colnames(qvalues_df)){
      sig_col <- gsub("_q", "_sig", qcol)
      data[, sig_col] <- FALSE
      data[which(data[, qcol] < 0.05), sig_col] <- TRUE
    }
  }
  return(data)
}

###Function: Adjust multiple columns of pvalues by qvalue method. Adjusts columns together
qval_adjustment_mult_groupcols <- function(data, pcols, return_sig = TRUE){
  #Histograms
  for (pcol in pcols){
    if (sum(!is.na(data[, pcol])) > 0){
      hist(data[, pcol], breaks = 20)
    }
  }
  #Adjust
  #Only if not all NA
  if (sum(!is.na(data[, pcols])) > 0){
    qvalues_df <- qvalue(data[, pcols])$qvalues %>% as.data.frame
  } else {
    qvalues_df <- data[, pcols] %>% as.data.frame
  }
  print(qvalues_df)
  colnames(qvalues_df) <- gsub("_p", "_q", pcols)
  data <- cbind(data, qvalues_df)
  #Add signifiance levels if desired
  if (return_sig){
    for (qcol in colnames(qvalues_df)){
      sig_col <- gsub("_q", "_sig", qcol)
      data[, sig_col] <- FALSE
      data[which(data[, qcol] < 0.05), sig_col] <- TRUE
    }
  }
  return(data)
}

###Function: Adjust multiple columns of pvalues together
p.adjust_mult_cols <- function(data, pval_cols, method = "fdr"){
  #Make list
  pval_list <- map(pval_cols, ~data[, .])
  #Unlist and adjust
  pval_vec <- unlist(pval_list)
  adj_pval_vec <- p.adjust(pval_vec, method = method)
  #Convert back to list
  vec_indices <- map(1:length(pval_cols), ~rep(., nrow(data))) %>% unlist
  adj_pval_list <- map(1:length(pval_cols), 
                       ~adj_pval_vec[which(vec_indices == .)])
  data_adj_pvals <- purrr::reduce(adj_pval_list, cbind) %>% as.data.frame
  colnames(data_adj_pvals) <- pval_cols
  rownames(data_adj_pvals) <- rownames(data)
  return(data_adj_pvals)
}



















