#'#*******************************************************************************************************************
#'@author : Amael DUPAIX
#'@update : 2021-11-17
#'@email : amael.dupaix@ens-lyon.fr
#'#*******************************************************************************************************************
#'@description :  Generate the plots associated with the GAMs performed
#'#*******************************************************************************************************************
#'@revisions
#'#*******************************************************************************************************************


if (part == 1){
  #' @months (pas utilise dans le model)
  p <- ggplot() +
    geom_boxplot(data=data, aes(x=fishing_month, y = Kn))+
    ylab("Kn")+xlab("Fishing month")
  
  l=1
  ggsave(plotsNames[l], p) ; l= l+1
  
  ggplot() +
    geom_point(data=data %>% plyr::ddply(.variables = "fishing_month", summarise, Kn = mean(Kn)),
               aes(x=as.factor(fishing_month), y = Kn))+
    ylab("Mean Kn")+xlab("Fishing month")
  
  
  #' @quarter (used in the model)
  p <- ggplot() +
    geom_boxplot(data=data, aes(x=fishing_quarter, y = Kn))+
    ylab("Kn")+xlab("Fishing quarter")
  
  ggsave(plotsNames[l], p) ; l= l+1
  
  ggplot() +
    geom_point(data=data %>% plyr::ddply(.variables = "fishing_quarter", summarise, Kn = mean(Kn)),
               aes(x=fishing_quarter, y = Kn))+
    ylab("Mean Kn")+xlab("Fishing quarter")
  
  #' @long_lat (used in the model)
  p <- ggplot() +
    geom_point(data= data, aes(x=lon, y=lat, color = Kn)) +
    xlab("Longitude")+ylab("Latitude")
  
  ggsave(plotsNames[l], p) ; l= l+1
  
  p <- ggplot() + geom_point(data= data, aes(x=lon, y = Kn))+xlab("Longitude")
  ggsave(plotsNames[l], p) ; l= l+1
  p <- ggplot() + geom_point(data= data, aes(x=lat, y = Kn))+xlab("Latitude")
  ggsave(plotsNames[l], p) ; l= l+1
  
  
} else if (part ==2){
  
  if (cluster == F){
    #' @save diagnostic plots of the gams
    pc1 <- gratia::appraise(gam1)
    pc2 <- gratia::appraise(gam2)
    ggsave(diagnoPlotNames[1], pc1, width = 100, height = 100, units = "mm", dpi = "retina")
    ggsave(diagnoPlotNames[2], pc2, width = 100, height = 100, units = "mm", dpi = "retina")
    
    saveRDS(gam1, gamteSummary)
    saveRDS(gam2, gamsSummary)
  } else {
    #' ajouter qqnorm et sauvegarde
  }
  
  #' @save plots of the coefficients of the gams
  countries <- map_data("world")
  
  xlm <- c(min(data$lon), max(data$lon))
  ylm <- c(min(data$lat), max(data$lat))
  
  if (is.null(smooth_col_limits)){
    if(Kn_transformation){
      lims <- c(-1,1)
      bwidth <- 0.25
    } else {
      lims <- c(0.9,1.1)
      bwidth <- 0.025
    }
  } else {
    lims <- smooth_col_limits
    bwidth <- (lims[2]-lims[1])/8
  }
  
  for (i in c(1,3)){
    ct_int_pred <- expand_grid(
      scaled_lon = seq(from=min(data$scaled_lon), 
                       to=max(data$scaled_lon), 
                       length.out = 100),
      scaled_lat = seq(from=min(data$scaled_lat), 
                       to=max(data$scaled_lat), 
                       length.out = 100),
      fishing_quarter = "1",
      fishing_year = min(as.numeric(as.character(data$fishing_year))),
      size_class = size_class_levels[1]
    )
    
    if (fad_fsc){
      ct_int_pred$fishing_mode <- "DFAD"
    }
    
    sc_lon <- attr(data$scaled_lon, "scaled:scale")
    ce_lon <- attr(data$scaled_lon, "scaled:center")
    sc_lat <- attr(data$scaled_lat, "scaled:scale")
    ce_lat <- attr(data$scaled_lat, "scaled:center")
    
    if(i==1){g <- gam1}else{g <- gam2}
    
    ct_int_pred <- predict(g, newdata = ct_int_pred, 
                           se.fit = TRUE) %>%  
      as_tibble() %>% 
      cbind(ct_int_pred) %>%
      mutate(lon = scaled_lon * sc_lon + ce_lon,
             lat = scaled_lat * sc_lat + ce_lat)
    
    predict_latlon <- ggplot(ct_int_pred, aes(x=lon, y=lat)) + 
      coord_sf(xlim = xlm, ylim = ylm, expand = FALSE, crs = st_crs(4326))
    
    predict_latlon2 <- predict_latlon +
      geom_point(data=data, aes(x=lon, y=lat), size = 0.01, color = "black")
    
    predict_latlon <- predict_latlon +
      geom_tile(aes(fill = fit), alpha = 0.8) +
      scale_fill_gradientn("Predicted\nKn", limits = lims, colors = c("blue","grey","red"))+
      geom_contour(aes(z = fit), binwidth = bwidth, colour = "grey40")+
      theme(panel.background = element_blank(),
            panel.border = element_rect(fill=NA))+
      geom_polygon(data=countries, aes(x=long, y=lat, group = group))+
      xlab("Longitude")+ylab("Latitude")
    
    predict_latlon2 <- predict_latlon2 +
      geom_tile(aes(fill = fit), alpha = 0.8) +
      scale_fill_gradientn("Predicted\nKn", limits = lims, colors = c("blue","grey","red"))+
      geom_contour(aes(z = fit), binwidth = bwidth, colour = "grey40")+
      theme(panel.background = element_blank(),
            panel.border = element_rect(fill=NA))+
      geom_polygon(data=countries, aes(x=long, y=lat, group = group))+
      xlab("Longitude")+ylab("Latitude")
    
    ggsave(smoothPlotNames[i], predict_latlon)
    ggsave(smoothPlotNames[i+1], predict_latlon2)
  }
  
  
  #'******************************
  #' Plots of the GAM coefficients
  #'******************************
  p1 <- plot.coeff(gam2, "fishing_year", labelx = "Year")
  p2 <- plot.coeff(gam2, "fishing_quarter", labelx = "Quarter")
  
  ggsave(gamCoeffPlotNames[1], p1)
  ggsave(gamCoeffPlotNames[2], p2)
  
  if (size_class_for_model == 'all'){
    p3 <- plot.coeff(gam2, "size_class", levels_order = size_class_levels[-1], labelx = "Size Class")
    ggsave(gamCoeffPlotNames[3], p3)
  }
  
  if (fad_fsc == T & fishing_mode_for_model == "all"){
    p4 <- plot.coeff(gam2, "fishing_mode", labelx = "School Type")
    ggsave(gamCoeffPlotNames[4], p4)
  }
  
}