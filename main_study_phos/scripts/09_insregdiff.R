####Background
#Store functions for exploring differences in insulin regulation

###initiliase
library(tidyverse)
library(purrr)
library(multcomp)


#####Summaries#####
###Function: Summarise number of proportion of levels in a factor that are changed by Strain or changed in n number of Strains
summary_StrainCHOW_change <- function(data, Var2){
  output_list = list("bool" = summariser_twoway_prop(data,
                                                     Var1 = "insregdiff_StrainCHOW_change_bool",
                                                     Var2 = Var2),
                     "num" = summariser_twoway_prop(data,
                                                    Var1 = "insregdiff_StrainCHOW_change_num",
                                                    Var2 = Var2))
  return(output_list)
}

###Function: Summarise proportions of levels in a factor that are changed by Diet, StrainxDiet, or both
summary_StrainxDiet_change <- function(data, Var2){
  output_list = list("Dietoverall" = summariser_twoway_prop(data,
                                                     Var1 = "insregdiff_Dietoverall_change_bool",
                                                     Var2 = Var2),
                     "Diet" = summariser_twoway_prop(data,
                                                    Var1 = "insregdiff_Diet_change_bool",
                                                    Var2 = Var2),
                     "StrainxDiet" = summariser_twoway_prop(data,
                                                     Var1 = "insregdiff_StrainxDiet_change_bool",
                                                     Var2 = Var2))
  return(output_list)
}



















