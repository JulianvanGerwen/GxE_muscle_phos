###Background
#This scripts contains objects used in this reserach project


#####Levels
Strain_levels <- c("C57Bl6J",
             "NOD",
             "BXH9",
             "CAST",
             "BXD34")
Diet_levels <- c("CHOW", "HFD")
Ins_levels <- c("bas", "ins")
UIDs <- c("C57Bl6J_CHOW_bas_1",
          "C57Bl6J_CHOW_bas_2",
          "C57Bl6J_CHOW_bas_3",
          "C57Bl6J_CHOW_bas_4",
          "C57Bl6J_CHOW_ins_1",
          "C57Bl6J_CHOW_ins_2",
          "C57Bl6J_CHOW_ins_3",
          "C57Bl6J_CHOW_ins_4",
          "C57Bl6J_CHOW_ins_5",
          "C57Bl6J_HFD_bas_1",
          "C57Bl6J_HFD_bas_2",
          "C57Bl6J_HFD_bas_3",
          "C57Bl6J_HFD_bas_4",
          "C57Bl6J_HFD_bas_5",
          "C57Bl6J_HFD_ins_1",
          "C57Bl6J_HFD_ins_2",
          "C57Bl6J_HFD_ins_3",
          "C57Bl6J_HFD_ins_4",
          "C57Bl6J_HFD_ins_5",
          "NOD_CHOW_bas_1",
          "NOD_CHOW_bas_2",
          "NOD_CHOW_bas_3",
          "NOD_CHOW_bas_4",
          "NOD_CHOW_bas_5",
          "NOD_CHOW_ins_1",
          "NOD_CHOW_ins_2",
          "NOD_CHOW_ins_3",
          "NOD_CHOW_ins_4",
          "NOD_CHOW_ins_5",
          "NOD_HFD_bas_1",
          "NOD_HFD_bas_2",
          "NOD_HFD_bas_3",
          "NOD_HFD_bas_4",
          "NOD_HFD_bas_5",
          "NOD_HFD_ins_1",
          "NOD_HFD_ins_2",
          "NOD_HFD_ins_3",
          "NOD_HFD_ins_4",
          "NOD_HFD_ins_5",
          "BXH9_CHOW_bas_1",
          "BXH9_CHOW_bas_2",
          "BXH9_CHOW_bas_3",
          "BXH9_CHOW_bas_4",
          "BXH9_CHOW_ins_1",
          "BXH9_CHOW_ins_2",
          "BXH9_CHOW_ins_3",
          "BXH9_CHOW_ins_4",
          "BXH9_HFD_bas_1",
          "BXH9_HFD_bas_2",
          "BXH9_HFD_bas_3",
          "BXH9_HFD_bas_4",
          "BXH9_HFD_ins_1",
          "BXH9_HFD_ins_2",
          "BXH9_HFD_ins_3",
          "BXH9_HFD_ins_4",
          "BXH9_HFD_ins_5",
          "CAST_CHOW_bas_1",
          "CAST_CHOW_bas_2",
          "CAST_CHOW_bas_3",
          "CAST_CHOW_bas_4",
          "CAST_CHOW_ins_1",
          "CAST_CHOW_ins_2",
          "CAST_CHOW_ins_3",
          "CAST_CHOW_ins_4",
          "CAST_CHOW_ins_5",
          "CAST_HFD_bas_1",
          "CAST_HFD_bas_2",
          "CAST_HFD_bas_3",
          "CAST_HFD_bas_4",
          "CAST_HFD_ins_1",
          "CAST_HFD_ins_2",
          "CAST_HFD_ins_3",
          "CAST_HFD_ins_4",
          "CAST_HFD_ins_5",
          "BXD34_CHOW_bas_1",
          "BXD34_CHOW_bas_2",
          "BXD34_CHOW_bas_3",
          "BXD34_CHOW_bas_4",
          "BXD34_CHOW_bas_5",
          "BXD34_CHOW_ins_1",
          "BXD34_CHOW_ins_2",
          "BXD34_CHOW_ins_3",
          "BXD34_CHOW_ins_4",
          "BXD34_CHOW_ins_5",
          "BXD34_HFD_bas_1",
          "BXD34_HFD_bas_2",
          "BXD34_HFD_bas_3",
          "BXD34_HFD_bas_4",
          "BXD34_HFD_bas_5",
          "BXD34_HFD_ins_1",
          "BXD34_HFD_ins_2",
          "BXD34_HFD_ins_3",
          "BXD34_HFD_ins_4",
          "BXD34_HFD_ins_5",
          "BXD34_HFD_ins_6")
pairwise_levels <- map(list("StrainDiet" = list(Diet_levels, Strain_levels),
                            "StrainIns" = list(Ins_levels, Strain_levels),
                            "DietIns" = list(Ins_levels, Diet_levels)),
                       function(list){
                         cross(list) %>%
                           map(~paste(rev(unlist(.)), collapse = "_")) %>%
                           unlist
                       })
StrainDietIns_levels <- unique(gsub("_\\d+", "", UIDs))
all_levels <- c(pairwise_levels,
                list("Strain" = Strain_levels,
                     "Diet" = Diet_levels,
                     "Ins" = Ins_levels,
                     "StrainDietIns" = StrainDietIns_levels))



#####Colours

###Five strain colours, decided 20220424, updated BXH9 to darker 20220629
fivestraincols_main <- c("C57Bl6J" = "#3780b6",
                         "NOD" = "#71c4b4",
                         "BXH9" = "#f9c45a",
                         "CAST" = "#f89d51",
                         "BXD34" = "#cb493d")

#Colours for insulinr regulation and differences
ins_updown_colours <- c("up" = "#d6604d",
                        "down" = "#4393c3")
updown_enhsupp_colours = c("up_enhanced" = "#cb1a1e", 
                          "up_enhanced;suppressed" = "#d6604d", 
                          "up_suppressed" = "#e0887a",
                          "down_enhanced" = "#1b519c", 
                          "down_enhanced;suppressed" = "#4393c3", 
                          "down_suppressed" = "#72aed2")
updown_enhsupp_labels = c("Up enhanced", "Up both", "Up suppressed",
                          "Down enhanced", "Down both","Down suppressed")
names(updown_enhsupp_labels) <- names(updown_enhsupp_colours)
enhsupp_colours <- c("enhanced" = "#141414",
                     "enhanced;suppressed" = "#9b9b9b",
                     "suppressed" = "#d1d3d3")
enhsupp_labels <- c("Enhanced", "Both", "Suppressed")
names(enhsupp_labels) <- names(enhsupp_colours)


#####Plotting parameters
medium_pointsize <- 1





