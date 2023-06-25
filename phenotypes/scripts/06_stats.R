###Background
#Here I keep code for performing statistics on phenotype data


###Initialise
library(tidyverse)
library(purrr)

#####t-tests

###Function: Perform multiple t-tests to compare different groups for the same phenotype value
#data: data to use
#Comparions: A list of pairs where each pair contains two conditions to compare
#id_col: A column containing the conditions that are in Comparisons
#value_col: The value to test
#
#outputs a list that contains a list of tests and a summary dataframe with pvalues, adjusted pvalues, and group differences and FCs
ttests_mult <- function(data,
                        comparisons,
                        id_col,
                        value_col){
  #Perform tests
  tests <- map(comparisons,
               function(x){
                 temp_test <- t.test(data[data[, id_col] == x[1], value_col],
                                     data[data[, id_col] == x[2], value_col])
               })
  #Make summary
  pvals <- map(tests, ~.$p.value) %>% unlist()
  adj_pvals <- p.adjust(pvals, method = "fdr")
  diffs <- map(tests, ~(.$estimate[2] - .$estimate[1])) %>% unlist()
  FCs <- map(tests, ~(.$estimate[2]/.$estimate[1])) %>% unlist()
  test_summary <- data.frame(comparisons = map(comparisons, ~paste(., collapse = "vs")) %>% unlist(),
                             p = pvals,
                             adj_p = adj_pvals,
                             diff = diffs,
                             FC = FCs)
  return(list("tests" = tests, "summary" = test_summary))
}



#####ANOVAs
###Function: Perform ANOVA and extract summary
#data: data to use
#value_col: Column of values to test
#A formula for aov where dependent variable is always Value
anova_wsummary <- function(data,
                           value_col,
                           formula = Value ~ Strain*Diet){
  #Set up data
  data$Value <- data[, value_col]
  #Run anova
  temp_aov <- aov(formula, data = data)
  #Get summary
  summary <- summary(temp_aov)[[1]] %>%
    as.data.frame() %>%
    select("F value", "Pr(>F)") %>%
    rename("F" = "F value", "p" = "Pr(>F)")
  rownames(summary) <- gsub(" ", "", rownames(summary))
  summary <- summary[rownames(summary) != "Residuals", ]
  return(list("aov" = temp_aov,
              "summary" = summary))
}

####Specific tests
###Function: Anova of Strain*Diet and t-tests comparing various conditions specified below
#data: Data to use
#value_col: value to test
#diet_compare: Compare CHOW to HFD within each strain
#strain_compare: Compare each strain to C57Bl6J within each diet
tests_strain_diet <- function(data,
                              value_col,
                              diet_compare = TRUE,
                              strain_compare = TRUE){
  #ANOVA
  aov_output <- anova_wsummary(data = data, value_col = value_col,
                               formula = Value ~ Strain*Diet)
  #Set up comparisons for ttests
  #Compare CHOW to HFD within each strain
  diet_comparisons <- map(levels(data$Strain),
                          function(strain){
                            paste(strain, levels(data$Diet), sep = "_")
                          })
  #Compare each strain to B6 within each diet
  strain_comparisons <- levels(pheno_data$Strain) %>%
    .[. != "C57Bl6J"] %>%
    map(~ c("C57Bl6J", .))
  strain_comparisons <- cross(list(strain_comparisons, as.list(levels(pheno_data$Diet)))) %>%
    map(~paste(.[[1]], .[[2]], sep = "_"))
  data <- mutate(data, Strain_Diet = paste(Strain, Diet, sep = "_"))
  #Decide which comparisons to keep
  if (!diet_compare){diet_comparisons <- NULL}
  if (!strain_compare){strain_comparisons <- NULL}
  comparisons <- c(diet_comparisons, strain_comparisons)
  #Run ttests
  ttest_output <- ttests_mult(data = data,
                              comparisons = comparisons,
                              id_col = "Strain_Diet",
                              value_col = value_col)
  return(list("aov" = aov_output,
              "ttest" = ttest_output))
}

###Function: For data separated by Strain, Diet, and Ins, performs Strain*Diet ANOVAs and t-tests on bas or ins data as well as t-tests comparing bas to ins in each condition
#data: Data to use
#value_col: value to test
#diet_compare: Compare CHOW to HFD within each strain
#strain_compare: Compare each strain to C57Bl6J within each diet
tests_strain_diet_ins <- function(data,
                                  value_col){
  ###Strain and Diet effects within bas or ins
  ins_strainvsdiet_tests <- tests_strain_diet(data[data$Ins == "ins", ],
                                              value_col)
  bas_strainvsdiet_tests <- tests_strain_diet(data[data$Ins == "bas", ],
                                              value_col)
  #Modualte comparisons so they contain bas or ins
  ins_strainvsdiet_tests$ttest$summary <- stat_summary_add_basorins(ins_strainvsdiet_tests$ttest$summary,
                                                                    "ins")
  bas_strainvsdiet_tests$ttest$summary <- stat_summary_add_basorins(bas_strainvsdiet_tests$ttest$summary,
                                                                    "bas")
  ###bas vs ins in all conditions
  basvsins_comparisons <- cross(list(levels(data$Diet),
                                     levels(data$Strain))) %>%
    map(~paste(rev(unlist(.)), collapse = "_")) %>%
    map(~paste(., c("bas", "ins"), sep = "_"))
  basvsins_tests <- ttests_mult(pheno_data,
                                comparisons = basvsins_comparisons,
                                id_col = "Condition",
                                value_col = value_col)
  ###Make outputs
  tests_output <- list("ins" = ins_strainvsdiet_tests,
                       "bas" = bas_strainvsdiet_tests,
                       "basvsins" = basvsins_tests)
  return(tests_output)
}



###Function: Big summary table of multiple tests. Tests must have identical composition
#anova_output_list: List of outputs from anova_wsummary or ttests_mult
test_summary_table <- function(test_output_list,
                               pval_col = "p"){
  summary_table <- map2(test_output_list, names(test_output_list),
                        function(test, name){
                          #Add significance codes
                          summary <- test$summary
                          sig_nums <- (summary[, pval_col] < 0.05) + (summary[, pval_col] < 0.01) + (summary[, pval_col] < 0.001)
                          summary$significance <- sapply(sig_nums,
                                                         function(x) paste(rep("*", x), collapse = ""))
                          #Rename
                          colnames(summary) <- paste(name, colnames(summary), sep = "_")
                          return(summary)
                        }) %>%
    purrr::reduce(cbind)
  return(summary_table)
}




#Function: Change comaprison names so that they include ins or bas
#stat_summary: stat_summary as output by ttests_mult
#basorins: self-explanatory
stat_summary_add_basorins <- function(stat_summary,
                                      basorins = "ins"){
  comparisons_split <- strsplit(stat_summary$comparisons, "vs")
  comparisons_w_basorins <- map(comparisons_split,
                                ~paste(.[1], "_", basorins,
                                       "vs",
                                       .[2], "_", basorins,
                                       sep = "")) %>%
    unlist
  stat_summary$comparisons <- comparisons_w_basorins
  return(stat_summary)
}



