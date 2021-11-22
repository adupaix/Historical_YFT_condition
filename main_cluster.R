#'#*******************************************************************************************************************
#'@author : Amael DUPAIX
#'@update : 2021-11-17
#'@email : amael.dupaix@ens-lyon.fr
#'#*******************************************************************************************************************
#'@description :  Second part of the main script to determine if the YFT condition factor (Kn) has been
#' impacted by the introduction of FADs in the Indian Ocean.
#' 
#' This part is the one which runs in the informatic cluster (to do the bootstraping)
#'#*******************************************************************************************************************
#'@revisions
#'#*******************************************************************************************************************

#' launched with the following command
#'    Rscript main_cluster.R [name of the rds object containing the copy of the env] [the number of the line in the command_list.txt]

init <- Sys.time()

Args = commandArgs(trailingOnly = T)

WD <- getwd()
# WD <- "home1/scratch/adupaix/Historical_YFT_condition"

env_name <- as.character(Args[1])
list_of_objects <- readRDS(file.path(WD, env_name))
list2env(list_of_objects, .GlobalEnv)

i <- as.numeric(Args[2])

DATA_PATH <- file.path(WD, "0.Data/")

FUNC_PATH <- file.path(WD,"1.Functions")
# OUTPUT_PATH <- file.path(WD, "3.Outputs", format(Sys.time()))
OUTPUT_PATH <- file.path(WD, "3.Outputs", format(Sys.Date()))
ROUT_PATH <- file.path(WD,"6.Sub-routines")

PLOT_PATH <- file.path(OUTPUT_PATH, "Plots")

seed.i <- seeds[i]

#' ***************
#' Get libraries:
#' ***************
source(file.path(FUNC_PATH, "install_libraries.R"))

srcUsedPackages <- c("plyr", "dplyr","tidyr","sf", "ggplot2", "tibble",
                     "cowplot","RColorBrewer", "MASS","truncnorm", "mgcv", "spdep")

if (cluster == F){srcUsedPackages <- c(srcUsedPackages, "gratia", "parallel", "foreach", "doSNOW")} else {Parallel[1] <- F}

installAndLoad_packages(srcUsedPackages, loadPackages = TRUE)


#' ********************
#' Init in @for loop:
#' ********************
source(file.path(ROUT_PATH, "0.i.Init_in_loop.R"))


#' ********************
#' Get fishing date:
#' ********************
#' at that point if the missing dates are deduced
if(deduce_date == T){
  source(file.path(ROUT_PATH, "2.Calc_fishing_date.R"))
}

#' ********************
#' Get fishing location:
#' ********************
source(file.path(ROUT_PATH, "3.Get_fishing_location.R"))


#' ***************
#' Perform GAM:
#' ***************
source(file.path(ROUT_PATH, "4.GAM.R"))

end <- Sys.time()

cat(format(end-init))
