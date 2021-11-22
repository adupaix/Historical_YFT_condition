#'#*******************************************************************************************************************
#'@author : Amael DUPAIX
#'@update : 2021-11-17
#'@email : amael.dupaix@ens-lyon.fr
#'#*******************************************************************************************************************
#'@description :  Initialization routine
#' load libraries, prepare files names, etc
#'#*******************************************************************************************************************
#'@revisions
#'#*******************************************************************************************************************

list2env(arguments, .GlobalEnv)

#' ***************
#' Get libraries:
#' ***************
source(file.path(FUNC_PATH, "install_libraries.R"))

srcUsedPackages <- c("plyr", "dplyr","tidyr","stringr","lubridate","sf", "ggplot2","tibble",
                     "cowplot","RColorBrewer", "MASS","truncnorm", "mgcv", "spdep")

if (cluster == F){srcUsedPackages <- c(srcUsedPackages, "gratia", "parallel", "foreach", "doSNOW")} else {Parallel[1] <- F}

installAndLoad_packages(srcUsedPackages, loadPackages = TRUE)

#' **********
#' Initialize:
#' **********
# For parallel:
# On Windows, or if don't want to parralelize, set cores number to 1
if (.Platform$OS.type == "windows" | as.logical(Parallel[1]) == F) {
  nb_cores = 1
} else { #use a fraction of the available cores
  nb_cores = trunc(detectCores() * as.numeric(Parallel[2]))
}


# Source functions
source(file.path(FUNC_PATH, "0.Prep_wl_data.R"))
source(file.path(FUNC_PATH, "plot_fig_2.R"))
source(file.path(FUNC_PATH, "subfunctions.R"))


# Initialize names
try(dir.create(PLOT_PATH, showWarnings = F, recursive = T))

fitName <- file.path(PLOT_PATH, "fit_lnorm_dates_goodnessoffit.png")
histFitName <- file.path(PLOT_PATH, "fit_lnorm_dates.png")
fig1Name <- file.path(PLOT_PATH, "Figure1.png")
fig2Name <- file.path(PLOT_PATH, "Figure2.png")

readName <- file.path(OUTPUT_PATH, "README.txt")

summaryName <- file.path(OUTPUT_PATH, "Processing_summary.txt")


# Initialize data processing summary and generate README.txt
sink(summaryName)
cat("Execution time:", format(Sys.time()), "\n\n")
sink()

sink(readName)
cat("Arguments used\n\n")
cat(paste(paste(names(format(arguments)), format(arguments), sep = ": "), collapse = "\n"))
sink()


# Define the seed number
#' @modif
if (nb_of_times_to_run != 1){
  seeds <- round(runif(nb_of_times_to_run, 1, 10^6))
} else {
  seeds = SEED
}

msg <- "\14" ; cat(msg) ; lines.to.cat <- c(msg)
msg <- paste("-------- PREPARING DATA --------") ; cat(msg) ; lines.to.cat <- c(lines.to.cat, msg)
