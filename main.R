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
OUTPUT_PATH <- file.path(WD, "3.Outputs", format(Sys.time()))
ROUT_PATH <- file.path(WD,"6.Sub-routines")

PLOT_PATH <- file.path(OUTPUT_PATH, "Plots")

arguments <- list(
#' @reproductibility
nb_of_times_to_run = 1, #' @for_study: set to 1000?
SEED = 123456,

#'@arguments:
#'#**********
#'# Run in parallel ?
#' (in prep_wl_data, to read geometry)
#' First element of the vector:
#'    If F, runs in sequential
#'    if T, runs in parallel
#' Second element of the vector: fraction of the cores to be used
Parallel = c(T, 1/2),

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
Kn_transformation = F,

#' Which method is to be used to determine the geomtery
#' either consider the centroid of the multipoint
#' or randomly sample a point from the multipoint
#'    one of "centroid" of "sampling"
#' @for_study: set to "sampling"
geometry_method = "sampling",

#' The limits of the colour scale of
#' the plot representing the lat-lon smooth
#' if NULL, the limits are c(0.9,1.1) if the GH transformation is not applied
#'                     or c(-1,1) if the GH transformation is applied
#' @for_study: set to NULL
# smooth_col_limits = c(0.8,1)
smooth_col_limits = NULL

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
  source(file.path(ROUT_PATH, "1.Calc_fishing_date.R"))
}

#' ********************
#' Generate Figures:
#' ********************
source(file.path(ROUT_PATH, "1.Generate_figures.R"))


for (i in 1:length(seeds)){
  
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
  
  #' If size_class_for_model != 'all', subsample the data to keep only one size class
  if (size_class_for_model != 'all'){
    
    if (size_class_for_model %in% size_classes_fig1){
      
      data <- data_byclass[[grep(size_class_for_model, size_classes_fig1)+1]]
      
    } else {
      
      if (grepl("-", size_class_for_model)){
        l1 <- as.numeric(sub("-.*", "", size_class_for_model))
        l2 <- as.numeric(sub(".*-", "", size_class_for_model))
      } else if (grepl(">", size_class_for_model)){
        l1 <- as.numeric(sub(">", "", size_class_for_model))
        l2 <- Inf
      } else if (grepl("<", size_class_for_model)){
        l1 <- 0
        l2 <- as.numeric(sub("<", "", size_class_for_model))
      }
      
      data %>% filter(fork_length > l1 & fork_length <= l2 ) -> data
      
    }
  }
  
  #' @1. Scale quantitative variables + scale and center geographical variables
  #' + filter data according to the provided arguments
  if (year_by_groups){
    ref <- (min(data$fishing_year)-5):(max(data$fishing_year)+5)
    ref <- ref[which(ref %% 5 == 0)]
    
    f <- function(x, ref) ref[max(which(ref<=x))]
    
    data$fishing_year <- mapply(f, data$fishing_year, MoreArgs = list(ref = ref))
  }
  
  data %>% mutate(scaled_FL = scale(fork_length, scale = T, center = F),
                  scaled_lon = scale(lon, scale = T, center = T),
                  scaled_lat = scale(lat, scale = T, center = T),
                  fishing_month = as.factor(fishing_month),
                  fishing_quarter = as.factor(fishing_quarter),
                  fishing_year = as.factor(fishing_year)) -> data
  
  data$size_class <- NA
  for (k in 1:length(size_class_levels)){
    if (k == 1){
      data$size_class[which(data$fork_length <= as.numeric(sub("<", "", size_class_levels[k])))] <- size_class_levels[k]
    } else if (k == length(size_class_levels)){
      data$size_class[which(data$fork_length > as.numeric(sub(">", "", size_class_levels[k])))] <- size_class_levels[k]
    } else {
      data$size_class[which(data$fork_length > as.numeric(sub("-.*", "", size_class_levels[k])) &
                              data$fork_length <= as.numeric(sub(".*-", "", size_class_levels[k])))] <- size_class_levels[k]
    }
  }

  if (fad_fsc == T){
    data %>% dplyr::filter(fishing_mode %in% c("DFAD","FSC")) %>%
      mutate(fishing_mode = as.factor(fishing_mode)) -> data
    
    if (fishing_mode_for_model != 'all'){
      data %>% filter(fishing_mode == fishing_mode_for_model) -> data
    }
    
  }
  
  if(Kn_transformation == T){
    
    data %>% mutate(Kn = geary.hinkley.transform(Kn,
                                                 weight_th,
                                                 whole_fish_weight,
                                                 cor.method = "pearson")) -> data
  }
  
  #' @2. Testing for colinearity among variables
  # Considered variables : lon, lat, fishing_quarter, size_class, year
  car::vif(lm(Kn ~ scaled_lon + scaled_lat + scaled_FL + fishing_year + fishing_quarter + fishing_mode, data = data))
  
  
  #' @3. Testing for spatial autocorrelation
  #'  @!!! the script was very long if size_class_for_model == "all" (26,000 entries in the df)
  #'  hence, the spatial autocorrelation is tested on a set of points only
  if (check_spatial_autocorr){
    part = 1 ; source(file.path(ROUT_PATH, "4.Spatial_autocorr.R"))
  }
  
  
  #' @4. Visualisation données
  part = 1 ; source(file.path(ROUT_PATH, "5.Generate_gam_plots.R"))
  
  #' @5. Build the model
  if (size_class_for_model == "all"){
    if (fad_fsc == F){
    gam1 <- mgcv::gam(Kn ~ fishing_quarter + fishing_year + size_class + te(scaled_lon, scaled_lat),
                        data = data)
    
    gam2 <- mgcv::gam(Kn ~ fishing_quarter + fishing_year + size_class + s(scaled_lon, scaled_lat),
                      data = data)
    
    } else {
      if (fishing_mode_for_model == "all"){
        gam1 <- mgcv::gam(Kn ~ fishing_quarter + fishing_year + size_class + fishing_mode + te(scaled_lon, scaled_lat),
                          data = data)
        gam2 <- mgcv::gam(Kn ~ fishing_quarter + fishing_year + size_class + fishing_mode + s(scaled_lon, scaled_lat),
                          data = data)
      } else {
        gam1 <- mgcv::gam(Kn ~ fishing_quarter + fishing_year + size_class + te(scaled_lon, scaled_lat),
                          data = data)
        gam2 <- mgcv::gam(Kn ~ fishing_quarter + fishing_year + size_class + s(scaled_lon, scaled_lat),
                          data = data)
      }
    }
  } else {
    if (fad_fsc == F){
      gam1 <- mgcv::gam(Kn ~ fishing_quarter + fishing_year + te(scaled_lon, scaled_lat),
                        data = data)
      
      gam2 <- mgcv::gam(Kn ~ fishing_quarter + fishing_year + s(scaled_lon, scaled_lat),
                        data = data)
      
    } else {
      if (fishing_mode_for_model == "all"){
        gam1 <- mgcv::gam(Kn ~ fishing_quarter + fishing_year + fishing_mode + te(scaled_lon, scaled_lat),
                          data = data)
        gam2 <- mgcv::gam(Kn ~ fishing_quarter + fishing_year + fishing_mode + s(scaled_lon, scaled_lat),
                          data = data)
      } else {
        gam1 <- mgcv::gam(Kn ~ fishing_quarter + fishing_year + te(scaled_lon, scaled_lat),
                          data = data)
        gam2 <- mgcv::gam(Kn ~ fishing_quarter + fishing_year + s(scaled_lon, scaled_lat),
                          data = data)
      }
    }
  }
  
  #' @6. Test spatial autocorrelation on the residuals
  if(check_spatial_autocorr){
    part = 2 ; source(file.path(ROUT_PATH, "4.Spatial_autocorr.R"))
  }
  
  #' @7. Generate the plots associated with the GAMs
  part = 2 ; source(file.path(ROUT_PATH, "5.Generate_gam_plots.R"))
  
  
}



