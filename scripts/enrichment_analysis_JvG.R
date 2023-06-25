###Last updated:
#20211020

###Background
#Here are all of my functions that perform enrichment
#Fisher's exact
#GSEA
#KSEA

###Packages
library(reshape2)
library(scales)
library(GO.db)
library(org.Mm.eg.db)
library(ksea)
library(rlist)
library(multcomp)
library(ggplot2)
library(tidyverse)
library(purrr)
library(qvalue)






#####Kinase enrichment#####

#Elise's kinase merger function
#kinases_pooled: A vector of kinase names to be merged
##Julian's additions:
#I combine the following kinase groups:
#MSKs seem to be similar and redundant (https://doi.org/10.3389/fcell.2016.00056)
#JNKs seem quite similar (10.1002/mc.20348)
kin_merger <- function(kinases_pooled){
  kinases_pooled[grep("^Akt", kinases_pooled)] <- "Akt"
  
  kinases_pooled[grep("^AMPK", kinases_pooled)] <- "AMPK"
  
  kinases_pooled[grep("^Aur", kinases_pooled)] <- "Aur"
  
  kinases_pooled[grep("^CAMK2", kinases_pooled)] <- "CAMK2"
  
  kinases_pooled[grep("^Chk", kinases_pooled)] <- "Chk"
  
  kinases_pooled[grep("^CK1", kinases_pooled)] <- "CK1"
  
  kinases_pooled[grep("^CK2", kinases_pooled)] <- "CK2"
  
  kinases_pooled[grep("^ERK", kinases_pooled)] <- "ERK"
  
  kinases_pooled[grep("^GSK3", kinases_pooled)] <- "GSK3"
  
  kinases_pooled[grep("^HIPK", kinases_pooled)] <- "HIPK"
  
  kinases_pooled[grep("^IKK", kinases_pooled)] <- "IKK"
  
  kinases_pooled[grep("^JNK", kinases_pooled)] <- "JNK"
  
  kinases_pooled[grep("^KHS", kinases_pooled)] <- "KHS"
  
  kinases_pooled[grep("^LATS", kinases_pooled)] <- "LATS"
  
  kinases_pooled[grep("^LIMK", kinases_pooled)] <- "LIMK"
  
  kinases_pooled[grep("^MARK", kinases_pooled)] <- "MARK"
  
  kinases_pooled[grep("^MKK", kinases_pooled)] <- "MKK"
  
  kinases_pooled[grep("^MSK", kinases_pooled)] <- "MSK"
  
  kinases_pooled[grep("^NDR", kinases_pooled)] <- "NDR"
  
  kinases_pooled[grep("^NEK", kinases_pooled)] <- "NEK"
  
  kinases_pooled[grep("^P38", kinases_pooled)] <- "P38"
  
  kinases_pooled[grep("^p70S6K", kinases_pooled)] <- "p70S6K"
  
  kinases_pooled[grep("^P70S6K", kinases_pooled)] <- "p70S6K"
  
  kinases_pooled[grep("^PAK", kinases_pooled)] <- "PAK"
  
  kinases_pooled[grep("^PDHK", kinases_pooled)] <- "PDHK"
  
  kinases_pooled[grep("^Pim", kinases_pooled)] <- "Pim"
  
  ##PKC subfamilies
  #Classical/conventional
  kinases_pooled[grep("^PKCA", kinases_pooled)] <- "cPKC"
  kinases_pooled[grep("^PKCB", kinases_pooled)] <- "cPKC"
  kinases_pooled[grep("^PKCG", kinases_pooled)] <- "cPKC"
  #Novel
  kinases_pooled[grep("^PKCD", kinases_pooled)] <- "nPKC"
  kinases_pooled[grep("^PKCE", kinases_pooled)] <- "nPKC"
  kinases_pooled[grep("^PKCH", kinases_pooled)] <- "nPKC"
  kinases_pooled[grep("^PKCQ", kinases_pooled)] <- "nPKC"
  kinases_pooled[grep("^PKCT", kinases_pooled)] <- "nPKC"
  #Atypical
  kinases_pooled[grep("^PKCI", kinases_pooled)] <- "aPKC"
  kinases_pooled[grep("^PKCZ", kinases_pooled)] <- "aPKC"
  
  
  kinases_pooled[grep("^PKG", kinases_pooled)] <- "PKG"
  
  kinases_pooled[grep("^PRKD", kinases_pooled)] <- "PRKD"
  
  kinases_pooled[grep("^ROCK", kinases_pooled)] <- "ROCK"
  
  kinases_pooled[grep("^SGK", kinases_pooled)] <- "SGK"
  
  kinases_pooled[grep("^ULK", kinases_pooled)] <- "ULK"
  
  
  return(kinases_pooled)
}


###Function: Pre-processing before KSEA. Removes psites with >= promisc_threshold number of kinases
#Returns boolean vector indicator stats to keep
promisc_psites_trim <- function(stat_names,
                                kin_sub_list,
                                promisc_threshold = 5){
  #Number of kinases for each psites
  psite_kinases <- map(stat_names, function(psite){
    if (is.na(psite)){
      return(NA)
    } else {
      kinase_bools <- map(kin_sub_list, ~(psite %in% .))
      kinases <- names(kin_sub_list)[unlist(kinase_bools)]
      if (length(kinases) == 0){
        kinases <- NA
      }
      return(kinases)
    }
  })
  psite_kin_nums <- map(psite_kinases, 
                        function(x){
                          if (is.na(x[1])){
                            return(0)
                          } else {
                            return(length(x))
                          }
                        }) %>% unlist
  #Filter
  promisc_keep_bool <- psite_kin_nums < promisc_threshold
  return(promisc_keep_bool)
}


###Function: Ochoa KSEA using a list of kinase substrates.
#Adjusts p-values by Benjamini-Hochberg
#kin_sub_list: Named list of kinase substates
#stats: Stats to use for enrichment e.g. log2FC
#stat_names: Names of stats. Must be of the same format as kin_sub_list
#min_sub_num: Minimum number of substrates that will be tested
#trial_n: Number of trials to run in determining significance
#promisc_threshold: An integer "k". All sites with >= k kinases are removed. If NULL, this is ignored
KSEA_from_list <- function(kin_sub_list, stats, stat_names,
                           min_sub_num = 5, trial_n = 1000,
                           promisc_threshold = NULL){
  ##Remove promiscuous psites assigned to >= k kinases
  if (!is.null(promisc_threshold)){
    promisc_keep_bool <- promisc_psites_trim(stat_names = stat_names,
                                             kin_sub_list = kin_sub_list,
                                             promisc_threshold = promisc_threshold)
    #Reassign stats and stat_names
    stats <- stats[promisc_keep_bool]
    stat_names <- stat_names[promisc_keep_bool]
  }
  
  
  ##Give names to stats
  #Give random labels to NAs
  NA_indices <- which(is.na(stat_names))
  if (length(NA_indices) > 0){
    stat_names[NA_indices] <- paste("NA", NA_indices, sep = "_")
  }
  names(stats) <- stat_names
  #Remove NA stats and then order
  stats <- stats[!is.na(stats)] %>%
    sort(decreasing = TRUE)
  #Trim Kinase list
  sub_nums <- map(kin_sub_list, ~length(intersect(., names(stats))))
  kin_sub_list <- kin_sub_list[which(sub_nums >= min_sub_num)]
  #Set up output
  ksea_output_df_temp <- matrix(0,
                                nrow = length(kin_sub_list),
                                ncol = 4)
  rownames(ksea_output_df_temp) <- names(kin_sub_list)
  colnames(ksea_output_df_temp) <- c("ES", "p", "p_adj", "num_subs")
  ##Fill out output
  #Cycle over kinases for enrichment scores and p-values
  for (i in 1:length(kin_sub_list)){
    ksea_temp_kinase_output <- ksea(names(stats),
                                    stats,
                                    kin_sub_list[[i]],
                                    trial = trial_n,
                                    significance = TRUE)
    num_subs <- length(intersect(kin_sub_list[[i]],
                                 names(stats)))
    ksea_output_df_temp[i, c(1, 2, 4)] <- c(ksea_temp_kinase_output$ES,
                                  ksea_temp_kinase_output$p.value,
                                  num_subs)
  }
  #Get adjusted pvalues
  ksea_output_df_temp[, "p_adj"] <- p.adjust(ksea_output_df_temp[, "p"], method = "fdr")
  #Put together
  ksea_output_df_temp <- ksea_output_df_temp[order(ksea_output_df_temp[, 3],
                                                   decreasing = FALSE),]
}

###Function: Ochoa KSEA using a list of kinase substrates, for multiple conditions 
#pvalues are adjusted by Benjamini-Hochberg within each condition
#kin_sub_list: Named list of kinase substates
#num_data: dataframe where columns are stats to use for enrichment e.g. log2FC
#conditions: Names used to label results of each column of num_data
#stat_names: Names of stats. Must be of the same format as kin_sub_list
#min_sub_num: Minimum number of substrates that will be tested
#trial_n: Number of trials to run in determining significance
#merge: if TRUE, merge results from each column into one dataframe
#promisc_threshold: An integer "k". All sites with >= k kinases are removed. If NULL, this is ignored
KSEA_from_list_multconds <- function(num_data, conditions, stat_names,
                                     kin_sub_list, min_sub_num = 5,
                                     trial_n = 1000, merge = FALSE,
                                     promisc_threshold = NULL){
  ##Remove promiscuous psites assigned to >= k kinases
  if (!is.null(promisc_threshold)){
    promisc_keep_bool <- promisc_psites_trim(stat_names = stat_names,
                                             kin_sub_list = kin_sub_list,
                                             promisc_threshold = promisc_threshold)
    #Reassign num_data and stat_names
    num_data <- num_data[promisc_keep_bool, ]
    stat_names <- stat_names[promisc_keep_bool]
  }
  #Loop over conditions to perform enrichment
  KSEA_list <- list()
  for (i in 1:ncol(num_data)){
    #Run KSEA
    #Do not include promisc_threshold, because this has already been performed above
    KSEA <- KSEA_from_list(kin_sub_list = kin_sub_list, 
                           stats = num_data[, i], stat_names = stat_names,
                           min_sub_num = min_sub_num, trial_n = trial_n,
                           promisc_threshold = NULL)
    #Add columns
    KSEA <- as.data.frame(KSEA)
    KSEA$tested <- TRUE
    KSEA_list[[i]] <- KSEA
    names(KSEA_list)[i] <- conditions[i]
  }
  #Expand each data frame to include all tested kinases overall
  tested_kinases <- map(KSEA_list, rownames) %>% unlist %>% unique
  KSEA_list <- map(KSEA_list, function(KSEA){
    untested_kinases <- setdiff(tested_kinases, rownames(KSEA))
    if (length(untested_kinases) > 0){
      KSEA[untested_kinases, ] <- NA
      KSEA[untested_kinases, "tested"] <- FALSE
    }
    KSEA <- KSEA[tested_kinases, ]
    return(KSEA)
  })
  #Merge if desired
  if (merge){
    KSEA_list <- map2(KSEA_list, names(KSEA_list), function(KSEA, name){
      colnames(KSEA) <- paste(name, colnames(KSEA), sep = "_")
      return(KSEA)
    }) %>%
      purrr::reduce(cbind)
  }
  return(KSEA_list)
}



##Heatmap function for KSEA from multiple conditions
#conditions: Prefixes for FC cols and pval cols. E.g. enrichment score cols are of the form condition_ES
#p_suffix: Suffix for p-values (e.g. p, p_adj)
heatmap_KSEA <- function(data, conditions,
                         p_suffix = "p",
                         order_col = NULL, ...){
  #Set up order col
  if (is.null(order_col)){
    order_col <- paste(conditions[1], "_ES", sep = "")
  }
  #Make heatmap
  output_plot <- fc_sig_heatmap(data,
                                fc_cols_w_pvals = paste(conditions, "_ES", sep = ""),
                                pval_cols = paste(conditions, p_suffix, sep = "_"),
                                order_column = order_col,
                                x_axis_names = conditions,
                                yaxis_naming = "keepall",
                                ...)
  return(output_plot)
}



#####Pathway enrichment#####
####Function that enriches from a list of pathways
#Up to date 20210528
#Supply list of relevant pathways, vector of background genes, and vector of DE genes
#Pathway_DE_intersection_threshold: Pathways need this number of DE genes or more to be considered
pathway_enricher_from_list <- function(background_genes,
                                       DE_genes,
                                       pathways_list,
                                       pathway_DE_intersection_threshold = FALSE,
                                       alternative = "greater"){
  
  ##Make pathways relevant 
  #Trim pathway list so it's only pathways with DE genes
  if (pathway_DE_intersection_threshold == FALSE){
    NULL
  } else {
    relevant_pathways_bool <- NULL
    old_pathways_list <- pathways_list
    for (i in 1:length(old_pathways_list)){
      
      relevant_pathways_bool[i] <- length(intersect(DE_genes,
                                                    old_pathways_list[[i]])) >= pathway_DE_intersection_threshold
    }
    
    #Terminate here if no relevant pathways
    if (sum(relevant_pathways_bool) == 0){
      
      return(NULL)
    } else {
      
      pathways_list <- as.list(names(old_pathways_list)[relevant_pathways_bool])
      names(pathways_list) <- names(old_pathways_list)[relevant_pathways_bool]
      for (i in 1:length(pathways_list)){
        
        pathways_list[[i]] <- old_pathways_list[[names(pathways_list)[i]]]
      }
    }
  }
  
  
  
  ##Set up output
  #pval
  #adj_pval
  #Number of DE genes in pathway
  #Number of background genes in pathway
  #Number of total genes in pathway
  output_m <- matrix(NA,
                     nrow = length(pathways_list),
                     ncol = 5)
  rownames(output_m) <- names(pathways_list)
  colnames(output_m) <- c("pval",
                          "adj_pval",
                          "num_DE_genes_in_pathway",
                          "num_background_genes_in_pathway",
                          "num_total_genes_in_pathway")
  
  ##Do stats
  #Do one-sided fisher's exact
  
  
  #Get x, m, n, k
  #m is number of DE genes
  #n is number of background, non-DE genes
  #k is number of background genes in pathway
  #x is number of DE genes in pathway
  
  m <- length(DE_genes)
  n <- length(background_genes) - m
  
  for (i in 1:length(pathways_list)){
    
    temp_pathway <- pathways_list[[i]]
    k <- length(intersect(background_genes,
                          temp_pathway))
    x <- length(intersect(DE_genes,
                          temp_pathway))
    #Run test
    test <- fisher.test(rbind(c(x, k - x),
                              c(m - x, n - k + x)),
                        alternative = alternative)
    
    temp_pval <- test$p.value
    output_m[i, ] <- c(temp_pval,
                       NA,
                       x,
                       k,
                       length(temp_pathway))
  }
  
  ##Adjust pvals
  output_m[, 2] <- p.adjust(output_m[, 1],
                            method = "fdr")
  
  ##Sort
  #By adj_pval. Don't order if length == 1
  if (nrow(output_m) > 1){
    output_m <- output_m[order(output_m[, 2],
                               decreasing = FALSE), ]
  }
  
  ##Add enrichment score
  output_m <- as.data.frame(output_m)
  output_m$ES <- (output_m$num_DE_genes_in_pathway/length(DE_genes))/
    (output_m$num_background_genes_in_pathway/length(background_genes))
  output_m$ES_log2 <- log2(output_m$ES)
  
  ##Return
  return(output_m)
}


####Function: Get GO pathways that intersect with a list of genes
#DE_genes: Genes of interest
#gene_label: Nami gnconvention for the genes, which ahs to be read by AnnotationDBI. Default is SYMBOL (e.g. Akt2, Slc2a4)
#organism_database: The AnnotationDBI organism database e.g. org.Mm.eg.db
#ontology: The desired GO ontology e.g. CC
relevant_GOIDS <- function(DE_genes,
                           gene_label = "SYMBOL",
                           organism_database,
                           ontology){
  relevant_GOIDS_df <- AnnotationDbi::select(organism_database,
                                             keys = DE_genes,
                                             columns = "GO",
                                             keytype = gene_label)
  relevant_GOIDS_df <- relevant_GOIDS_df[which(relevant_GOIDS_df$ONTOLOGY == 
                                                 ontology),]
  relevant_GOIDS <- unique(relevant_GOIDS_df$GO)
  
  ###List of all genes for relevant_GOIDS
  relevant_GOIDS_genes_df <- AnnotationDbi::select(organism_database,
                                                   keys = relevant_GOIDS,
                                                   columns = gene_label,
                                                   keytype = "GO")
  ##Make list
  #Make sure no gene name complete duplicates. Can happen because of different evidence levels
  relevant_GOIDS_genes_list <- as.list(unique(relevant_GOIDS_genes_df$GO))
  names(relevant_GOIDS_genes_list) <- unique(relevant_GOIDS_genes_df$GO)
  for (i in 1:length(relevant_GOIDS_genes_list)){
    
    temp_GOID <- names(relevant_GOIDS_genes_list)[i]
    relevant_GOIDS_genes_list[[i]] <- unique(relevant_GOIDS_genes_df[which(relevant_GOIDS_genes_df$GO == temp_GOID), 
                                                                     gene_label])
  }
  return(relevant_GOIDS_genes_list)
}



####Function that enriches for terms of a given GO ontology
#Up to date 20210528
#Supply vector of background genes, vector of DE genes, gene_label (UNIPROT, SYMBOL, etc), organism database (e.g. org.Mm.eg.db), and ontology (e.g. CC)
GO_enricher <- function(background_genes,
                        DE_genes,
                        gene_label = "SYMBOL",
                        organism_database,
                        ontology,
                        pathway_DE_intersection_threshold,
                        alternative = "greater"){
  
  ###Get relevant GOIDs
  relevant_GOIDS_genes_list <- relevant_GOIDS(DE_genes = DE_genes,
                                              gene_label = gene_label,
                                              organism_database = organism_database,
                                              ontology = ontology)
  
  
  ##Run tests
  output_df <- as.data.frame(pathway_enricher_from_list(background_genes,
                                                        DE_genes,
                                                        relevant_GOIDS_genes_list,
                                                        pathway_DE_intersection_threshold = pathway_DE_intersection_threshold,
                                                        alternative = alternative),
                             stringsasfactors = FALSE)
  
  ##Add GO term and ID
  output_df$GO_term <- AnnotationDbi::select(GO.db,
                                             rownames(output_df),
                                             columns = "TERM",
                                             keytype = "GOID")$TERM
  output_df$GO_id <- rownames(output_df)
  rownames(output_df) <- output_df$GO_term
  
  ##Return
  return(output_df)
}



###Fisher's exact test for DE genes and pathway
#Up to date 20210712
#Supply background genes, DE genes, pathway genes, and specifict alternative
fishers_DE_pathway <- function(background_genes,
                               DE_genes,
                               pathway_genes,
                               alternative = "greater"){
  
  ##Run test
  #m: number of DE genes
  #n: number of background, non-DE genes
  #k: number of background genes in pathway
  #x: number of DE genes in pathway
  m <- length(DE_genes)
  n <- length(background_genes) - m
  k <- length(intersect(background_genes,
                        pathway_genes))
  x <- length(intersect(DE_genes,
                        pathway_genes))
  
  test <- fisher.test(rbind(c(x, k - x),
                            c(m - x, n - k + x)),
                      alternative = alternative)
  pval <- test$p.value
  odds_ratio <- test$estimate
  
  ##Fold enrichment
  fold_enr <- (x/m)/(k/length(background_genes))
  
  ##Output
  output <- c(pval,
              x,
              k,
              length(pathway_genes),
              odds_ratio,
              fold_enr)
  output <- as.list(output)
  names(output) <- c("pval",
                     "num_DE_genes_in_pathway",
                     "num_background_genes_in_pathway",
                     "num_pathway_genes",
                     "odds_ratio",
                     "fold_enr")
  return(output)
}



#####Data processing#####

##Function: Condense list of overrepresentation enrichments to one
condense_overrep_list <- function(data_list){
  ##Get common rows
  common_rows <- map(data_list, rownames) %>% unlist %>% unique
  ##Loop over data
  data_df <- map2(data_list, names(data_list), function(data, name){
    #Add common rows
    data[setdiff(common_rows, rownames(data)), ] <- NA
    data <- data[common_rows, ]
    #Rename cols
    colnames(data) <- paste(name, colnames(data), sep = "_")
    #Return
    return(data)
  }) %>% purrr::reduce(cbind) %>% as.data.frame
}




