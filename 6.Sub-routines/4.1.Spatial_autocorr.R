#'#*******************************************************************************************************************
#'@author : Amael DUPAIX
#'@update : 2021-11-17
#'@email : amael.dupaix@ens-lyon.fr
#'#*******************************************************************************************************************
#'@description :  Routine to test the spatial autocorrelation of the data and of the model residuals
#'#*******************************************************************************************************************
#'@revisions
#'#*******************************************************************************************************************


if (part == 1){
  
  #' Test spatial autocorrelation on the variable
  
  msg <- "\n\n~~~~ Performing Moran's I test on Kn ~~~~" ; cat(msg) ; lines.to.cat <- c(lines.to.cat, msg)
  
  if (size_class_for_model == "all"){
    to_sample <- sample.int(dim(data)[1], size = 2000)
  } else {
    to_sample <- 1:(dim(data)[1])
  }
  subdata <- data[to_sample,]
  
  c <- seq(.5,5,.5)
  listw <- list() ; mc <- list() ; test <- list()
  
  for (k in 1:length(c)){
    cat(lines.to.cat)
    cat(paste("\n  - Consider Neighboors if d <=", c[k]*100, "km"))
    cat("\n      - Getting nearest neighboors")
    nb_list <- spdep::dnearneigh(x = st_coordinates(subdata),
                                 d1 = 0,
                                 d2 = c[k]*100,
                                 longlat = T)
    
    msg <- "\n      - Getting weight list from nearest neighboors" ; cat(msg)
    listw[[k]] <- spdep::nb2listw(nb_list, zero.policy = T)
    
    msg <- "\n      - Performing Moran's I tests" ; cat(msg)
    mc[[k]] <- moran.mc(subdata$Kn, listw = listw[[k]], nsim = 99, zero.policy = T)
    test[[k]] <- moran.test(subdata$Kn, listw = listw[[k]], zero.policy = T) 
    
  }
  
  morans.mc <- unlist(lapply(mc, function(x) x$statistic))
  morans.test <- unlist(lapply(test, function(x) x$estimate[1]))
  
  toplot <- data.frame(d = rep(c, 2),
                       I = c(morans.mc, morans.test),
                       moran_type = c(rep("mc",length(c)),
                                      rep("test",length(c))))
  
} else if (part == 2){
  
  #' Test spatial autocorrelation on the residuals
  cat(lines.to.cat)
  msg <- "\n\n~~~~ Performing Moran's I tests on model residuals ~~~~" ; cat(msg) ; lines.to.cat <- c(lines.to.cat, msg)
  
  c <- seq(.5,5,.5)
  mc_gam1 <- list() ; test_gam1 <- list() ; mc_gam2 <- list() ; test_gam2 <- list()
  
  for (k in 1:length(c)){
    cat(lines.to.cat)
    cat(paste("\n  - Consider Neighboors if d <=", c[k]*100, "km"))
    # cat("\n      - Getting nearest neighboors")
    # nb_list <- spdep::dnearneigh(x = st_coordinates(subdata),
    #                              d1 = 0,
    #                              d2 = c[k]*100,
    #                              longlat = T)
    # 
    # msg <- "\n      - Getting weight list from nearest neighboors" ; cat(msg)
    # listw[[k]] <- nb2listw(nb_list, zero.policy = T)
    
    msg <- "\n      - Performing Moran's I tests" ; cat(msg)
    mc_gam1[[k]] <- moran.mc(resid(gam1)[to_sample], listw = listw[[k]], nsim = 99, zero.policy = T)
    test_gam1[[k]] <- moran.test(resid(gam1)[to_sample], listw = listw[[k]], zero.policy = T) 
    mc_gam2[[k]] <- moran.mc(resid(gam2)[to_sample], listw = listw[[k]], nsim = 99, zero.policy = T)
    test_gam2[[k]] <- moran.test(resid(gam2)[to_sample], listw = listw[[k]], zero.policy = T) 
    
  }
  
  morans.mc.res.gam1 <- unlist(lapply(mc_gam1, function(x) x$statistic))
  morans.test.res.gam1 <- unlist(lapply(test_gam1, function(x) x$estimate[1]))
  morans.mc.res.gam2 <- unlist(lapply(mc_gam2, function(x) x$statistic))
  morans.test.res.gam2 <- unlist(lapply(test_gam2, function(x) x$estimate[1]))
  
  toplot2 <- data.frame(d = rep(c, 4),
                        I = c(morans.mc.res.gam1, morans.test.res.gam1, morans.mc.res.gam2, morans.test.res.gam2),
                        moran_type = c(rep("mc_res_gam1",length(c)),
                                       rep("test_res_gam1",length(c)),
                                       rep("mc_res_gam2",length(c)),
                                       rep("test_res_gam2",length(c))))
  
  bind_rows(toplot, toplot2) -> toplot
  
  
  p_moran_residuals <- ggplot(data = toplot, aes(x=d*100, y=I, group = moran_type, color = moran_type)) + 
    geom_point()+ geom_line()+
    scale_color_discrete("Moran's\nfunction")+
    xlab("Distance (km)")+
    ylab("Moran's I")
  
  toplot <- toplot[grep("test", toplot$moran_type),]
  toplot$moran_type <- gsub(pattern = "test$", "Kn", toplot$moran_type)
  toplot$moran_type <- gsub(pattern = "test_res_gam1", "GAM1 residuals", toplot$moran_type)
  toplot$moran_type <- gsub(pattern = "test_res_gam2", "GAM2 residuals", toplot$moran_type)
  
  p_moran_residuals2 <- ggplot(data = toplot, aes(x=d*100, y=I, group = moran_type, color = moran_type)) + 
    geom_point()+ geom_line()+
    scale_color_discrete("Data used to\n calculate Moran's I")+
    xlab("Distance (km)")+
    ylab("Moran's I")
  
  ggsave(moranPlotName, p_moran_residuals2)
  
  rm(mc, test, nb_list, listw, part) ; invisible(gc())
  
  
}

