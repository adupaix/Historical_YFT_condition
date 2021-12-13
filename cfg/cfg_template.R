#'#*******************************************************************************************************************
#'@author : Amael DUPAIX
#'@update : 2021-12-06
#'@email : amael.dupaix@ens-lyon.fr
#'#*******************************************************************************************************************
#'@description :  Configuration file to launch the main script to study the variation of the physiological condition
#' of tuna in the Indian Ocean and the impact of FADs on that condition.
#'#*******************************************************************************************************************
#' @tags:
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

#'#**********
#'@arguments:
#'#**********
arguments <- list(
  #' General arguments
  #'#********************
  
  #' @reproductibility
  #' number of times to perform the geometry sampling and build the GAM
  #' @for_study: set to 1000
  nb_of_times_to_run = 1000,
  #' If we run only once, the seed to use
  #' @for_study: does not matter
  SEED = 123456,
  
  #'# Run in parallel ?
  #' (@! only in prep_wl_data, to read geometry)
  #' First element of the vector:
  #'    If F, runs in sequential
  #'    if T, runs in parallel
  #' Second element of the vector: fraction of the cores to be used
  #' @! if cluster==T, Parallel[1] is automatically set to F
  Parallel = c(F, 1/2),
  
  #' If RESET is FALSE, try to read the prepped data
  #' else, prep the data in any case
  #' @for_study: set to T
  RESET = F,
  
  #' If VERBOSE is T, print information at every steps of the script
  #' if F, nothing is printed
  #' @! if cluster==T, VERBOSE is automatically set to F
  VERBOSE = T,
  
  #' Choose if the GAMs are generated following this script (F)
  #' or if the scripts necessary to generate the GAMs
  #' on a cluster are generated (T)
  #' #' @for_study: set to T
  cluster = F,
  
  #' Arguments for data preparation
  #'#******************************
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
  
  #' @Figure_1
  #' **********
  #' size classes for filtering data
  #' Used to generate Figure 1
  #' @for_study: set to c("<75","75-120",">120")
  # size_classes = c("40-60","<102",">102"),
  # size_classes = c("<75","75-120",">120"),
  size_classes_fig1 = c("<75","75-120",">120"),
  
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
  
  
  #' @GAM
  #' **********
  #' choose if the model is performed on all the individuals
  #' or only on a given size class
  #' @for_study: set to "all"
  size_class_for_model = "all",
  
  #' choose how to define the size classes used in the model (in the size_class variable)
  #' @!! the vector has to be in the format c("<x1","x1-x2",...,"x(n-1)-xn",">xn")
  #' @for_study: set to c("<75","75-120",">120")
  size_class_levels = c("<75","75-120",">120"),
  # size_class_levels = c("<90",">90"),
  
  #' Either to perform the Moran's tests (T) or not (F)
  #' @for_study: set to T
  check_spatial_autocorr = T,
  
  #' choose if the fishing mode is included in the GAM
  #' variables (T) or not (F)
  #' @for_study: both are used
  #'             F for main study
  fad_fsc = F,
  
  #' choose if the model is performed on all the individuals
  #' or only on one of the fishing modes
  #' one of "all", "DFAD", "FSC"
  #' 
  #' if fad_fsc is False, this argument is skipped
  #' @for_study: skipped for the main study
  #'             DFAD and FSC (with fad_fsc T) for supplementary
  fishing_mode_for_model = "all",
  
  #' Values of the explanatory variables chosen as reference levels
  #'  also used to make the spatial prediction
  #' @for_study: ref_var_values = list(fishing_quarter = "1",
  #'                                       fishing_year = 2019,
  #'                                       size_class = "<75")
  #'                                       
  ref_var_values = list(fishing_quarter = "1",
                        fishing_year = 2019,
                        size_class = "<75",
                        fishing_mode = "DFAD"),
  
  
  #' @Plots
  #' **********
  #' The limits of the colour scale of
  #' the plot representing the lat-lon smooth
  #' if NULL, the limits are c(0.9,1.1) if the GH transformation is not appliedd
  #'                     or c(-1,1) if the GH transformation is applied
  #' @for_study: set to NULL
  # smooth_col_limits = c(0.8,1)
  smooth_col_limits = NULL,
  
  #' Choose if plots are generated or not
  #' if (F), the script only performs the model
  #' without building figure or coefficient plots
  #' However, even if (F) the script generates the
  #' diagnostic plots
  #' @for_study: set to F
  #' @! if cluster==T, generate_plots is automatically set to F
  generate_plots = T,
  
  #' Geographical limits of the maps plotted during the study
  #'  also used to generate the prediction and to calculate the density
  #'  of sampling points per cell
  #' @for_study: set to xlm = c(30, 85),
  #'                    ylm = c(-25,20)
  xlm = c(30, 85),
  ylm = c(-25,20)
  
  
)


source(file.path(WD, "main.R"))