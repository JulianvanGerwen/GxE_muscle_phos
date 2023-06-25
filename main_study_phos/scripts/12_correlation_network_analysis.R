###Background
#Here I store code for performing correlation network analysis

###Initialise
library(WGCNA)
library(ComplexHeatmap)
library(circlize)

#####WGCNA functions#####

###Function: Orders data by WGCNA clusters and then performs hierarchical clustering within each cluster
#num_data: Numerical data where columns are variables (to be clustered) and rows are samples
#WGCNA_output: Output of blockwiseModules run on num_data
hclust_within_WGCNAclust <- function(num_data,
                                     WGCNA_output){
  ##Split by clusters
  clust_names <- sort(unique(WGCNA_net$colors))
  data_clust_list <- map(clust_names,
                         ~num_data[, which(WGCNA_output$colors == .)])
  names(data_clust_list) <- clust_names
  ##Loop over clusters and perform hierarchical clustering
  data_clust_list <- map(data_clust_list, function(data){
    dist_m <- dist(t(data))
    hclust <- hclust(dist_m)
    data <- data[, hclust$order]
    return(data)
  })
  ##Merge list
  num_data_ordered <- purrr::reduce(data_clust_list, cbind)
  return(num_data_ordered)
}

#####Cluster interrogation#####
###Function: Enrich in each cluster, combine, adjust pvalues together
#clust_genes must be a list
clust_enricher <- function(clust_genes,
                           bkd_genes,
                           pways_list){
  ##Run enrichment
  enr_compr <- map(clust_genes, function(clust){
    enr <- pathway_enricher_from_list(background_genes = bkd_genes,
                                      DE_genes = clust,
                                      pathways_list = pways_list,
                                      pathway_DE_intersection_threshold = 0) %>% as.data.frame
    return(enr)
  })
  ##Condense to single dfs and adjust pvals together
  enr_compr_df <-  condense_overrep_list(enr_compr)
  #Adjust p-values together
  pval_only_cols <- setdiff(colnames(enr_compr_df) %>% .[grep("_pval$", .)],
                            colnames(enr_compr_df) %>% .[grep("_adj_pval$", .)])
  enr_compr_adjpcomb <- p.adjust_mult_cols(enr_compr_df, 
                                           pval_cols = pval_only_cols)
  colnames(enr_compr_adjpcomb) <- strsplit(colnames(enr_compr_adjpcomb), "_") %>%
    map(~paste(.[1], "_adj_pval_comb", sep = "")) %>% unlist
  enr_compr_df <- cbind(enr_compr_df, enr_compr_adjpcomb)
  return(enr_compr_df)
}
























