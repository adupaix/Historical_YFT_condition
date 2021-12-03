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
  
  #'********************
  #' Plots of Kn values
  #'********************
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
    #' @save diagnostic plots of the gam
    pc1 <- gratia::appraise(gam1)
    
  } else {
    
    pc1 <- qqplot.gam(residuals(gam1))
  }
  
  ggsave(diagnoPlotNames, pc1, width = 100, height = 100, units = "mm", dpi = "retina")
  
  #'****************************************
  #' Plots of the GAM latlong smooth term
  #'****************************************
  #' @save plots of the coefficients of the gams
  if (generate_plots){
    countries <- map_data("world")

    
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
    
    
    ct_int_pred <- expand_grid(
      lon = seq(from=min(data$lon), 
                       to=max(data$lon), 
                       length.out = 100),
      lat = seq(from=min(data$lat), 
                       to=max(data$lat), 
                       length.out = 100),
      fishing_quarter = "1",
      fishing_year = min(as.numeric(as.character(data$fishing_year))),
      size_class = size_class_levels[1]
    )
      
    if (fad_fsc){
      ct_int_pred$fishing_mode <- "DFAD"
    }
      
    ct_int_pred <- predict(gam1, newdata = ct_int_pred, 
                           se.fit = TRUE) %>%  
      as_tibble() %>% 
      cbind(ct_int_pred)
      
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
      
    ggsave(smoothPlotNames[1], predict_latlon)
    ggsave(smoothPlotNames[2], predict_latlon2)
    

    #'******************************
    #' Plots of the GAM coefficients
    #'******************************
    p1 <- plot.coeff(gam1, "fishing_year", labelx = "Year")
    p2 <- plot.coeff(gam1, "fishing_quarter", labelx = "Quarter")
    
    ggsave(gamCoeffPlotNames[1], p1)
    ggsave(gamCoeffPlotNames[2], p2)
    
    if (size_class_for_model == 'all'){
      p3 <- plot.coeff(gam1, "size_class", levels_order = size_class_levels[-1], labelx = "Size Class")
      ggsave(gamCoeffPlotNames[3], p3)
    }
    
    if (fad_fsc == T & fishing_mode_for_model == "all"){
      p4 <- plot.coeff(gam1, "fishing_mode", labelx = "School Type")
      ggsave(gamCoeffPlotNames[4], p4)
    }
  }
}