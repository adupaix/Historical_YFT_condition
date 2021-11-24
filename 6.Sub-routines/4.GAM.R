#'#*******************************************************************************************************************
#'@author : Amael DUPAIX
#'@update : 2021-11-17
#'@email : amael.dupaix@ens-lyon.fr
#'#*******************************************************************************************************************
#'@description :  Sub-routine which performs the GAMs (also get associated plots, etc)
#'#*******************************************************************************************************************
#'@revisions
#'#*******************************************************************************************************************

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

data %>% mutate(scaled_lon = scale(lon, scale = T, center = T),
                scaled_lat = scale(lat, scale = T, center = T),
                fishing_month = as.factor(fishing_month),
                fishing_quarter = as.factor(fishing_quarter),
                fishing_year = as.factor(fishing_year),
                size_class = as.factor(size_class)) -> data

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

#' @2. Testing for spatial autocorrelation
#'  @!!! the script was very long if size_class_for_model == "all" (26,000 entries in the df)
#'  hence, the spatial autocorrelation is tested on a set of points only
if (check_spatial_autocorr){
  part = 1 ; source(file.path(ROUT_PATH, "4.1.Spatial_autocorr.R"))
}


#' @3. Visualisation données
if (generate_plots){
  part = 1 ; source(file.path(ROUT_PATH, "4.2.Generate_gam_plots.R"))
}

#' @4. Build the model
if (size_class_for_model == "all"){
  if (fad_fsc == F){
    # gam1 <- mgcv::gam(Kn ~ fishing_quarter + fishing_year + size_class + te(scaled_lon, scaled_lat),
    #                   data = data)
    
    gam1 <- mgcv::gam(Kn ~ fishing_quarter + fishing_year + size_class + s(scaled_lon, scaled_lat),
                      data = data)
    
  } else {
    if (fishing_mode_for_model == "all"){
      # gam1 <- mgcv::gam(Kn ~ fishing_quarter + fishing_year + size_class + fishing_mode + te(scaled_lon, scaled_lat),
      #                   data = data)
      gam1 <- mgcv::gam(Kn ~ fishing_quarter + fishing_year + size_class + fishing_mode + s(scaled_lon, scaled_lat),
                        data = data)
    } else {
      # gam1 <- mgcv::gam(Kn ~ fishing_quarter + fishing_year + size_class + te(scaled_lon, scaled_lat),
      #                   data = data)
      gam1 <- mgcv::gam(Kn ~ fishing_quarter + fishing_year + size_class + s(scaled_lon, scaled_lat),
                        data = data)
    }
  }
} else {
  if (fad_fsc == F){
    # gam1 <- mgcv::gam(Kn ~ fishing_quarter + fishing_year + te(scaled_lon, scaled_lat),
    #                   data = data)
    
    gam1 <- mgcv::gam(Kn ~ fishing_quarter + fishing_year + s(scaled_lon, scaled_lat),
                      data = data)
    
  } else {
    if (fishing_mode_for_model == "all"){
      # gam1 <- mgcv::gam(Kn ~ fishing_quarter + fishing_year + fishing_mode + te(scaled_lon, scaled_lat),
      #                   data = data)
      gam1 <- mgcv::gam(Kn ~ fishing_quarter + fishing_year + fishing_mode + s(scaled_lon, scaled_lat),
                        data = data)
    } else {
      # gam1 <- mgcv::gam(Kn ~ fishing_quarter + fishing_year + te(scaled_lon, scaled_lat),
      #                   data = data)
      gam1 <- mgcv::gam(Kn ~ fishing_quarter + fishing_year + s(scaled_lon, scaled_lat),
                        data = data)
    }
  }
}

#' @5. Apply a stepAIC to select the most parsimonious model
gam1 <- stepAIC.gam(gam1, verbose = VERBOSE)
# gam2 <- stepAIC.gam(gam2, verbose = VERBOSE)


#' @6. Test spatial autocorrelation on the residuals
if(check_spatial_autocorr){
  part = 2 ; source(file.path(ROUT_PATH, "4.1.Spatial_autocorr.R"))
}

#' @7. Generate the plots associated with the GAMs
if (generate_plots){
  part = 2 ; source(file.path(ROUT_PATH, "4.2.Generate_gam_plots.R"))
}
