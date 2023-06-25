###Background
#Here are functions related to phosphosite plus e.g. loading in data


directory_test <- function(){
  
  setwd(getSrcDirectory(function(x) {x}))
  print(getwd())
}



###Function: General function that loads psp data
#directory: directory where psp data is e.g. "data/psp". If NULL, directory is figured out based on location of this script (not recommended if you are not using this on JvG's computer)
load_psp_data <- function(filename,
                          directory = NULL){
  current_wd <- getwd()
  
  #Load to default location if directory is null
  if (is.null(directory)){
    #Set wd to script
    setwd(getSrcDirectory(function(x) {x}))
    psp_data <- read.delim(paste("..\\..\\data/biol_databases/phosphosite_plus", filename, sep = "/"))
  } else {
    psp_data <- read.delim(paste(directory, filename, sep = "/"))
  }
  setwd(current_wd)
  return(psp_data)
}



###Function: Uniprot site to site group ID
#uniprot_sites: vector of sites to convert of the form e.g. P81122_S604
#psp_phos_data: Dataset of phosphosites to use. If NULL, will load in itself
#restrict_to_phos: Only use phosphorylation modifications (probably redundant)
uniprotsite_to_sitegroupid <- function(uniprot_sites,
                                       psp_phos_data = NULL,
                                       restrict_to_phos = TRUE){
  if (is.null(psp_phos_data)){
    psp_phos_data <- load_psp_data("PSP__phos_dataset_20211007")
  }
  if (restrict_to_phos){
    psp_phos_data <- psp_phos_data[grep("p$", psp_phos_data$MOD_RSD), ]
  }
  #Make uniprot_site in psp_phos_data
  psp_phos_data$uniprot_site <- paste(psp_phos_data$ACC_ID, 
                                      gsub("-.+", "", psp_phos_data$MOD_RSD), 
                                      sep = "_")
  site_group_ids <- lapply(uniprot_sites, function(x) psp_phos_data$SITE_GRP_ID[which(psp_phos_data$uniprot_site == x)])
  names(site_group_ids) <- uniprot_sites
  site_group_ids[vapply(site_group_ids, function(x) {length(x) == 0}, logical(1))] <- NA
  return(site_group_ids)
}



###Function: Search against a psp database, matching by site_group_id
search_psp_by_sitegroupid <- function(sites,
                                      site_format = "uniprot_site",
                                      database_name = NULL,
                                      database_supplied = NULL){
  #Get database
  if (is.null(database_supplied)){
    psp_data <- load_psp_data(database_name)
  } else {
    psp_data <- database_supplied
  }
  
  #Get site group ids if needed
  if (site_format == "uniprot_site"){
    sitegroupids <- uniprotsite_to_sitegroupid(sites)
    sitegroupids <- sitegroupids[vapply(sitegroupids, function(x){length(x) > 0}, logical(1))]
  } else if (site_format == "SITE_GRP_ID") {
    sitegroupids <- sites
  } else {
    return("Incorrect value for site_format")
  }
  
  #Search data
  psp_data_searched <- psp_data[which(psp_data$SITE_GRP_ID %in% sitegroupids), ]
  psp_data_searched$queried_site <- sapply(psp_data_searched$SITE_GRP_ID,
                                           function(x) paste(names(sitegroupids[which(sitegroupids == x)]), collapse = ";"))
  return(psp_data_searched)
}



###Function: Collapse PSP data by site group ID
collapse_psp_by_sitegroupid <- function(psp_data,
                                        desired_cols = NULL){
  collapsed_data <- data.frame(SITE_GRP_ID = unique(psp_data$SITE_GRP_ID))
  
  if (is.null(desired_cols)){
    desired_cols <- colnames(psp_data)[colnames(psp_data) != "SITE_GRP_ID"]
  }
  collapsed_data[, desired_cols] <- NA
  for (col in desired_cols){
    collapsed_data[, col] <- sapply(collapsed_data$SITE_GRP_ID,
                                                       function(x) paste(psp_data[which(psp_data$SITE_GRP_ID == x), col],
                                                                         collapse = ";;;"))
  }
  return(collapsed_data)
}


#####Kinases#####
###Function: Fetch all psites on kinases
#HOW IT WORKS: Gets gene names for each unique kinase in psp_kin_sub_data
#Converts each gene name into caps for human and lowercase for mouse
#Fetches all psites on these genes using psp_data
psp_get_kinasepsites <- function(psp_kin_sub_data,
                             psp_data){
  #Get genes for each kinases
  unique_kinases <- unique(psp_kin_sub_data$KINASE)
  kinase_genes_list <- map(unique_kinases,
                           function(kin){
                             sub_data <- subset(psp_kin_sub_data,
                                                KINASE == kin)
                             return(unique(sub_data$GENE))
                           })
  names(kinase_genes_list) <- unique_kinases
  #Expand to have mouse and human gene names
  kinase_genes_list <- map(kinase_genes_list, function(genes){
    #remove empty or NA genes
    genes <- genes[which(genes != "" &
                           !is.na(genes))]
    caps_genes <- unique(toupper(genes))
    lower_genes <- tolower(caps_genes) %>%
      strsplit("") %>%
      map(function(x){
        x[1] <- toupper(x[1])
        return(paste(x, collapse = ""))
      }) %>% unlist
    return(c(caps_genes, lower_genes))
  })
  
  #Get all phosphosites on each kinase
  psp_kin_autophos_list <- map(kinase_genes_list,
                               function(genes){
                                 psites <- map(genes, function(gene){
                                   data <- subset(psp_data, GENE == gene)
                                   return(unique(data$SITE_GRP_ID))
                                 }) %>% unlist %>% unique
                                 return(psites)
                               })
  return(psp_kin_autophos_list)
}
  

  
  
  
  












