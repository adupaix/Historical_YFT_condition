#'#*******************************************************************************************************************
#'@author : Amael DUPAIX
#'@update : 2022-02-20
#'@email : amael.dupaix@ens-lyon.fr
#'#*******************************************************************************************************************
#'@description :  Testing the correlation between Kn and the number of FOBs (2013-2019)
#' The number of FOBs is obtained from the Chapter 4 of Baidai (2020)
#' 
#' Baidai Y (2020) Derivation of direct abundance index for tropical tunas based on their associative behavior with floating objects. Universite de Montpellier, Montpellier, France
#'#*******************************************************************************************************************
#'@revisions:
#'
#'#*******************************************************************************************************************

rm(list = ls())
invisible(gc())

WD <- "/home/adupaix/Documents/These/Axe_3/Historical_YFT_condition/"

DATA_PATH <- file.path(WD, "0.Data/")

FUNC_PATH <- file.path(WD,"1.Functions")
# OUTPUT_PATH <- file.path(WD, "3.Outputs", format(Sys.time()))
OUTPUT_PATH <- file.path(WD, "3.Outputs")

library(plyr)
library(dplyr)
library(lubridate)
library(sf)
library(ggplot2)

#'**************************************************

last_whole_output <- list.dirs(OUTPUT_PATH, recursive = F, full.names = F)
last_whole_output <- grep("whole_dataset", last_whole_output, value = T)
folder_dates <- as.Date(gsub("-whole_dataset", "", last_whole_output))
last_whole_output <- last_whole_output[which(folder_dates == max(folder_dates))]


readRDS(file.path(OUTPUT_PATH, last_whole_output, "df_filtered.rds")) %>%
  dplyr::filter(as.numeric(as.character(fishing_year)) >= 2013) %>%
  plyr::ddply(.variables = c("fishing_year"), summarise, Kn = mean(Kn)) %>%
  dplyr::rename("year" = "fishing_year") -> df


read.csv2(file.path(DATA_PATH, "Baidai2020-estimated_number_of_FOBs.csv")) %>%
  plyr::ddply(.variables = c("year", "timescale"), function(x) sum(x$NFob)) %>%
  plyr::ddply(.variables = c("year"), function(x) mean(x$V1)) %>%
  dplyr::rename("NFob" = "V1") %>%
  dplyr::filter(!is.na(NFob)) -> nFOB

merge(nFOB, df) -> toplot

p = ggplot(toplot, aes(x = NFob, y = Kn))+
  geom_point()+
  geom_label(aes(label = year),
             nudge_x = -100, 
             nudge_y = 0.002)+
  xlab("Estimated number of FOB from Baidai (2020)")+
  ylab("Mean Kn")

ggsave(file.path("/home/adupaix/Documents/These/Axe_3/cours_redaction_article/MEPS/Rev1/Sup_corKnNFob.png"),
       width = 6, height = 6)

cor.test(toplot$NFob, toplot$Kn, method = "spearman")
