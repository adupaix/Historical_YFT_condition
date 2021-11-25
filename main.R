#'#*******************************************************************************************************************
#'@author : Amael DUPAIX
#'@update : 2021-11-17
#'@email : amael.dupaix@ens-lyon.fr
#'#*******************************************************************************************************************
#'@description :  Main script to determine if the YFT condition factor (Kn) has been impacted by the introduction of
#' FADs in the Indian Ocean.
#'#*******************************************************************************************************************
#' @tags:
#' the @for_study tag is used before each argument, to specify which arguments were used to obtain the study results
#' the @filter tag is used to make explicit where the different filters are applied (not perfectly in agreement with
#'             the names of the sub-routines)
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
nb_of_times_to_run = 1000, #' @for_study: set to 1000?
SEED = 123456,

#'@arguments:
#'#**********
#'# Run in parallel ?
#' (in prep_wl_data, to read geometry)
#' First element of the vector:
#'    If F, runs in sequential
#'    if T, runs in parallel
#' Second element of the vector: fraction of the cores to be used
#' @! if cluster==T, Parallel[1] is automatically set to F
Parallel = c(F, 1/2),

#' If RESET is FALSE, try to read the prepped data
#' else, prep the data in any case
#' @for_study: set to T
RESET = T,

#' If VERBOSE is T, print information at every steps of the script
#' if F, nothing is printed
#' @! if cluster==T, VERBOSE is automatically set to F
VERBOSE = T,

#' When missing, (T) choose to sample fork length (FL) from fish with the same first dorsal length (FDL)
#'  (F) or not
#'  @for_study: set to F
calcfdl = F,

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
check_spatial_autocorr = T,

#' @Figure_1
#' size classes for filtering data
#' Used to generate Figure 1
#' @for_study: set to c("<75","75-120",">120")
# size_classes = c("40-60","<102",">102"),
# size_classes = c("<75","75-120",">120"),
size_classes_fig1 = c("<75","75-120",">120"),

#' @GAM
#' choose if the model is performed on all the individuals
#' or only on a given size class
#' @for_study: set to "all"
size_class_for_model = "all",

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
geometry_method = "sampling",

#' The limits of the colour scale of
#' the plot representing the lat-lon smooth
#' if NULL, the limits are c(0.9,1.1) if the GH transformation is not appliedd
#'                     or c(-1,1) if the GH transformation is applied
#' @for_study: set to NULL
# smooth_col_limits = c(0.8,1)
smooth_col_limits = NULL,

#' Choose if the GAMs are generated following this script (F)
#' or if the scripts necessary to generate the GAMs
#' on a cluster are generated (T)
#' #' @for_study: set to T
cluster = T,

#' Choose if plots are generated or not
#' if (F), the script only performs the model
#' without building figure or coefficient plots
#' However, even if (F) the script generates the
#' diagnostic plots
#' @for_study: set to F
#' @! if cluster==T, generate_plots is automatically set to F
generate_plots = F
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
#' @filter measured FL, measured W, fish caught in the IO
data <- prep_wl_data(DATA_PATH, data, calcfdl = calcfdl, read = !RESET, getgeom = T, ncores = nb_cores,
                     size_class_levels = size_class_levels, summaryName = summaryName, verbose = VERBOSE)

#' ***************
#' Calculate Kn:
#' ***************

if (VERBOSE){
  msg <- "\n\n~~~~ Calculating Kn ~~~~\n" ; cat(msg) ; lines.to.cat <- c(lines.to.cat, msg)
}

#' @filter species of interest, caught by purse seiner
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

#' @filter fishing location not missing
#' Keep only the data entries for which any type of geometry is available (POINT or MULTIPOINT)
data %>% filter(!st_is_empty(data$geometry)) -> data

sink(summaryName, append = T)
cat("\n\n\nMissing location filter")
cat("\n-----------------------\n")
cat(paste0("\n    - Number of entries after missing location filter: ", dim(data)[1]))
sink()

#' ********************
#' Get fishing date:
#' ********************
#' at that point if the entries with missing dates are discarded
#' @filter known fishing date interval, no error in date,
#'         fish caught prior to 2020
if(deduce_date == F){
  source(file.path(ROUT_PATH, "2.Calc_fishing_date.R"))
}

#' ********************
#' Generate Figures:
#' ********************
#' @filter "non duplicated" line:
#'          concerns ~ 400 entries in the YFT data, 156 fish
#'          some rows are duplicated with only the geometry (POINT)
#'          differing
#' @!! in this subroutine, work on "data_fig" not on "data" (only the first of the
#' duplicated rows is kept, as the geographical data is not used)
#' This filter will be applied again on "data" but applying the geometry_method
#'  in sub-routine 3.Get_fishing_location
if (generate_plots){
  source(file.path(ROUT_PATH, "1.Generate_figures.R"))
}



if (cluster == T){
  #' ******************************
  #' Prepare for datarmor script
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
    #' first, read and save data, so the sampling can be performed again
    #'  in the next iteration
    if (i != 1){
      data <- readRDS(intermediateDataName)
    }
    saveRDS(data, file = intermediateDataName)
    
    source(file.path(ROUT_PATH, "0.i.Init_in_loop.R"))
    
    
    #' ********************
    #' Get fishing date:
    #' ********************
    #' at that point if the missing dates are deduced
    #' @filter known fishing date interval, no error in date, fishing location not missing,
    #'         fish caught prior to 2020
    if(deduce_date == T){
      source(file.path(ROUT_PATH, "2.Calc_fishing_date.R"))
    }
    
    #' ********************
    #' Get fishing location:
    #' ********************
    #' @filter "non duplicated" line:
    #'          concerns ~ 400 entries in the YFT data, 156 fish
    #'          some rows are duplicated with only the geometry (POINT)
    #'          differing
    source(file.path(ROUT_PATH, "3.Get_fishing_location.R"))
    
    
    #' ***************
    #' Perform GAM:
    #' ***************
    source(file.path(ROUT_PATH, "4.GAM.R"))
    
  }
  
  file.remove(intermediateDataName)
}

