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

srcUsedPackages <- c("plyr", "dplyr","tidyr","lubridate","sf", "ggplot2","tibble",
                     "cowplot","RColorBrewer", "MASS", "mgcv", "spdep",
                     "doSNOW")

if (cluster == F){srcUsedPackages <- c(srcUsedPackages, "gratia", "parallel")} else {Parallel[1] <- F ; VERBOSE <- F ; generate_plots <- F}

installAndLoad_packages(srcUsedPackages, loadPackages = TRUE, verbose = VERBOSE)

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
if (generate_plots){
  try(dir.create(PLOT_PATH, showWarnings = F, recursive = T))
  
  fitName <- file.path(PLOT_PATH, "fit_lnorm_dates_goodnessoffit.png")
  histFitName <- file.path(PLOT_PATH, "fit_lnorm_dates.png")
  fig1Name <- file.path(PLOT_PATH, "Figure1.png")
  fig1rdsName <- file.path(PLOT_PATH, "Figure1.rds")
  fig2Name <- file.path(PLOT_PATH, "Figure2.png")
}

try(dir.create(OUTPUT_PATH, showWarnings = F, recursive = T))
intermediateDataName <- file.path(OUTPUT_PATH, "df_intermediate.rds")
dfGeneralName <- file.path(OUTPUT_PATH, "df_filtered.rds")
readName <- file.path(OUTPUT_PATH, "README.txt")
argsName <- file.path(OUTPUT_PATH, "README.rds")

summaryName <- file.path(OUTPUT_PATH, "Processing_summary.txt")

# names of the fit of the allometrique law
allomSummaryName <- file.path(OUTPUT_PATH, "Fit_allometric_law_summary.txt")
allomFitName <- file.path(OUTPUT_PATH, "Fit_allometric_law_plot.png")

# Initialize data processing summary and generate README.txt
sink(summaryName)
cat("Execution time:", format(Sys.time()), "\n\n")
sink()

sink(readName)
cat("Arguments used\n\n")
cat(paste(paste(names(format(arguments)), format(arguments), sep = ": "), collapse = "\n"))
sink()

saveRDS(arguments, argsName)


# Define the seed number
#' @modif
if (nb_of_times_to_run != 1){
  seeds <- round(runif(nb_of_times_to_run, 1, 10^8))
} else {
  seeds = SEED
}

if (VERBOSE){
  msg <- "\14" ; cat(msg) ; lines.to.cat <- c(msg)
  msg <- paste("-------- PREPARING DATA --------") ; cat(msg) ; lines.to.cat <- c(lines.to.cat, msg)
}