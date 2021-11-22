#'#*******************************************************************************************************************
#'@author : Amael DUPAIX
#'@update : 2021-11-17
#'@email : amael.dupaix@ens-lyon.fr
#'#*******************************************************************************************************************
#'@description :  Main script to determine if the YFT condition factor (Kn) has been impacted by the introduction of
#' FADs in the Indian Ocean.
#' the @for_study tag is used before each argument, to specify which arguments were used to obtain the study results
#'#*******************************************************************************************************************
#'@revisions
#'#*******************************************************************************************************************

rm(list = ls())
invisible(gc())

WD <- getwd()

DATA_PATH <- file.path(WD, "0.Data/")

FUNC_PATH <- file.path(WD,"1.Functions")
# OUTPUT_PATH <- file.path(WD, "3.Outputs", format(Sys.time()))
OUTPUT_PATH <- file.path(WD, "3.Outputs", format(Sys.Date()))
ROUT_PATH <- file.path(WD,"6.Sub-routines")

PLOT_PATH <- file.path(OUTPUT_PATH, "Plots")

arguments <- list(
#' @reproductibility
nb_of_times_to_run = 10, #' @for_study: set to 1000?
SEED = 123456,

#'@arguments:
#'#**********
#'# Run in parallel ?
#' (in prep_wl_data, to read geometry)
#' First element of the vector:
#'    If F, runs in sequential
#'    if T, runs in parallel
#' Second element of the vector: fraction of the cores to be used
Parallel = c(F, 1/2),

#' If RESET is FALSE, try to read the prepped data
#' else, prep the data in any case
#' @for_study: set to T
RESET = F,

#' When missing, choose to sample fork length (FL) from fish with the same first dorsal length (FDL)
#'  or not
#'  @for_study: set to F
calcfdl = F,

#' Read the geometry column or not
#' need to read from the chr column and change it in geometry format
#' @for_study: set to T
#' @!! if set to F, the second part of the script won't work (it's a rest from before)
getgeom = T,

#' Choose the species of interest
#' one of YFT, BET, or SKJ
#' @for_study: set to "YFT"
species = "YFT",

#' Choose which allometric relationship is going to be used to calculate the condition factor
#' one of "Chassot2015" of "fromData"
#' @for_study: set to "fromData"
allom = "fromData",

#' Either to perform the Moran's tests (T) or not (F)
#' @for_study: set to T
check_spatial_autocorr = F,

#' @Figure_1
#' size classes for filtering data
#' Used to generate Figure 1
#' @for_study: set to c("<75","75-120",">120")
# size_classes = c("40-60","<102",">102"),
# size_classes = c("<75","75-120",">120"),
size_classes_fig1 = c(),

#' @GAM
#' choose if the model is performed on all the individuals
#' or only on a given size class
#' @for_study: set to "all"
size_class_for_model = "<75",

#' choose how to define the size classes used in the model (in the size_class variable)
#' @!! the vector has to be in the format c("<x1","x1-x2",...,"x(n-1)-xn",">xn")
#' @for_study: set to c("<75","75-120",">120")
size_class_levels = c("<75","75-120",">120"),
# size_class_levels = c("<90",">90"),

#' choose if the fishing mode is included in the GAM
#' variables (T) or not (F)
#' @for_study: both are used
fad_fsc = F,

#' choose if the model is performed on all the individuals
#' or only on one of the fishing modes
#' one of "all", "DFAD", "FSC"
#' 
#' if fad_fsc is False, this argument is skipped
#' @for_study: set to "all"
fishing_mode_for_model = "all",

#' Either pick randomly the missing dates from the log-normal
#' distribution fitted on the time between the sampling and the fishing date (T)
#' or discard the data with missing fishing date interval (F)
#' @for_study: set to F
deduce_date = F,

#' Either use the years as single factor levels in the GAM (F)
#' or group them by 5 years intervals (T)
#' @for_study: set to F
year_by_groups = F,

#' Use the Kn in the GAM (F)
#' or transform it using the Geary-Hinkley transformation (T)
#' see https://doi.org/10.1287/mnsc.21.11.1338
#' @for_study: set to T
Kn_transformation = T,

#' Which method is to be used to determine the geomtery
#' either consider the centroid of the multipoint
#' or randomly sample a point from the multipoint
#'    one of "centroid" of "sampling"
#' @for_study: set to "sampling"
geometry_method = "centroid",

#' The limits of the colour scale of
#' the plot representing the lat-lon smooth
#' if NULL, the limits are c(0.9,1.1) if the GH transformation is not applied
#'                     or c(-1,1) if the GH transformation is applied
#' @for_study: set to NULL
# smooth_col_limits = c(0.8,1)
smooth_col_limits = NULL,

#' Choose if the GAMs are generated following this script (F)
#' or if the scripts necessary to generate the GAMs
#' on a cluster are generated (T)
cluster = T

)


#' **********
#' Initialize:
#' **********
source(file.path(ROUT_PATH, "0.Init.R"))


#' **********
#' Read data:
#' **********
data1 <- read.table(file = file.path(DATA_PATH, "Tunabio_OI_20210826_env.csv"),
                    sep = ";", dec=",",stringsAsFactors = F,
                    quote="\"", row.names = NULL,
                    h=T, colClasses = "character")
data2 <- read.table(file = file.path(DATA_PATH, "Tunabio_OI_20210826_specimen.csv"),
                    sep = ";", dec=",",stringsAsFactors = F,
                    quote="\"", row.names = NULL,
                    h=T, colClasses = "character")

data <- base::merge(data1, data2, by = "fish_identifier", all = T)

rm(data1, data2)

#' **********
#' Prep data:
#' **********
data <- prep_wl_data(DATA_PATH, data, calcfdl = calcfdl, read = !RESET, getgeom = getgeom, ncores = nb_cores,
                     summaryName = summaryName)

#' ***************
#' Calculate Kn:
#' ***************

msg <- "\n\n~~~~ Calculating Kn ~~~~\n" ; cat(msg) ; lines.to.cat <- c(lines.to.cat, msg)

data %>% filter(species_code_fao %in% species & gear_code == "PS") -> data

if (allom == "Chassot2015"){
  a = 0.00002459
  b = 2.96670
} else if (allom == "fromData"){
  a = exp(lm(log(whole_fish_weight) ~ log(fork_length), data = data)$coefficients[1])
  b = lm(log(whole_fish_weight) ~ log(fork_length), data = data)$coefficients[2]
}

data %>%
  dplyr::mutate(weight_th = a*fork_length^b,
                Kn = whole_fish_weight / weight_th) -> data

sink(summaryName, append = T)
cat("\n\n\nSpecies and gear filter")
cat("\n------------------------\n")
cat(paste0("\n    - Number of entries after gear filter + species filter (", species, "): ", dim(data)[1]))
sink()

#' ********************
#' Get fishing date:
#' ********************
#' at that point if the entries with missing dates are discarded
if(deduce_date == F){
  source(file.path(ROUT_PATH, "2.Calc_fishing_date.R"))
}

#' ********************
#' Generate Figures:
#' ********************
source(file.path(ROUT_PATH, "1.Generate_figures.R"))



if (cluster == T){
  #' ******************************
  #' Generate script for datarmor
  #' ******************************
  # save all the objects of the environments to a list
  env <- mget(ls())
  fname <- "all_objects.rds"
  # and save that list as an .rds file (which will be loaded again afterwards)
  saveRDS(env, file.path(WD, fname))
  
  #' *************************
  #' Generate the command list
  #' *************************
  commands <- paste0("Rscript ",
                     file.path(WD, "main_cluster.R "),
                     fname, " ",
                     1:nb_of_times_to_run)
  
  file.create(file.path(WD, "commands_bootstrap.txt"))
  cmds <- file(file.path(WD, "commands_bootstrap.txt"), open = "w")
  writeLines(commands, cmds)
  close(cmds)
  
  #source(file.path(ROUT_PATH, "6.Datarmor_prep.R"))
  
} else {
  
  for (i in 1:length(seeds)){
    
    seed.i <- seeds[i]
    
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
    
  }
  
}

