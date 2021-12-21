#'#*******************************************************************************************************************
#'@author : Amael DUPAIX
#'@update : 2021-11-17
#'@email : amael.dupaix@ens-lyon.fr
#'#*******************************************************************************************************************
#'@description :  Sub-routine which performs the GAMs (also get associated plots, etc)
#'#*******************************************************************************************************************
#'@revisions
#'#*******************************************************************************************************************

#' @1. If size_class_for_model != 'all', subsample the data to keep only one size class
if (size_class_for_model != 'all'){
  
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
    
  data %>% dplyr::filter(fork_length > l1 & fork_length <= l2 ) -> data
    
}

#' @2. Scale quantitative variables + scale and center geographical variables
#' + filter data according to the provided arguments
if (year_by_groups){
  ref <- (min(data$fishing_year)-5):(max(data$fishing_year)+5)
  ref <- ref[which(ref %% 5 == 0)]
  
  f <- function(x, ref) ref[max(which(ref<=x))]
  
  data$fishing_year <- mapply(f, data$fishing_year, MoreArgs = list(ref = ref))
}

if (fad_fsc == T){
  data %>% dplyr::filter(fishing_mode %in% c("DFAD","FSC")) %>%
    mutate(fishing_mode = as.factor(fishing_mode)) -> data
}
if (fishing_mode_for_model != 'all'){
  data %>% filter(fishing_mode == fishing_mode_for_model) -> data
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
#' 
#' relevel the variables, to have the chosen level as reference
data$fishing_year <- relevel(data$fishing_year, ref = ref_var_values$fishing_year)
data$fishing_quarter <- relevel(data$fishing_quarter, ref = ref_var_values$fishing_quarter)

if (size_class_for_model == "all"){
  
  data$size_class <- relevel(data$size_class, ref = ref_var_values$size_class)
  
  if (fad_fsc == F | fishing_mode_for_model != "all"){
    gam1 <- mgcv::gam(Kn ~ fishing_quarter +
                        fishing_year +
                        size_class +
                        s(lon, lat),
                      data = data)
    
  } else if (fishing_mode_for_model == "all"){
      
      data$fishing_mode <- relevel(data$fishing_mode, ref = ref_var_values$fishing_mode)
      
      gam1 <- mgcv::gam(Kn ~ fishing_quarter +
                          fishing_year +
                          size_class +
                          fishing_mode +
                          s(lon, lat),
                        data = data)
    
  }
} else {
  if (fad_fsc == F | fishing_mode_for_model != "all"){
    gam1 <- mgcv::gam(Kn ~ fishing_quarter +
                        fishing_year +
                        s(lon, lat),
                      data = data)
    
  } else if (fishing_mode_for_model == "all"){
      
      data$fishing_mode <- relevel(data$fishing_mode, ref = ref_var_values$fishing_mode)
      
      gam1 <- mgcv::gam(Kn ~ fishing_quarter +
                          fishing_year +
                          fishing_mode +
                          s(lon, lat),
                        data = data)
    
  }
}

#' @5. Apply a stepAIC to select the most parsimonious model
gam1 <- stepAIC.gam(gam1, verbose = VERBOSE)
# gam2 <- stepAIC.gam(gam2, verbose = VERBOSE)


#' @6. Test spatial autocorrelation on the residuals
if(check_spatial_autocorr){
  part = 2 ; source(file.path(ROUT_PATH, "4.1.Spatial_autocorr.R"))
}

#' @7. Save the GAM
saveRDS(gam1, gamSummary)

#' @8. Generate the plots associated with the GAMs
part = 2 ; source(file.path(ROUT_PATH, "4.2.Generate_gam_plots.R"))

