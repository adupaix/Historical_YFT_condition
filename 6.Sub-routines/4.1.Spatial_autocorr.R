#'#*******************************************************************************************************************
#'@author : Amael DUPAIX
#'@update : 2021-11-17
#'@email : amael.dupaix@ens-lyon.fr
#'#*******************************************************************************************************************
#'@description :  Routine to test the spatial autocorrelation of the data and of the model residuals
#'#*******************************************************************************************************************
#'@revisions 2021-11-23: tested distances to consider 2 points as "linked" changed from 100-1000km to 100-500km
#'#*******************************************************************************************************************


if (part == 1){
  
  #' @1. Test spatial autocorrelation on the variable
  
  if (VERBOSE){
    msg <- "\n\n~~~~ Performing Moran's I test on Kn ~~~~" ; cat(msg) ; lines.to.cat <- c(lines.to.cat, msg)
  }
  
  #' The number are a bit arbitrary, but preliminary tests showed that
  #' subsampling the dataset to look at the spatial autocorrelation
  #' did not change the obtained results
  if (dim(data)[1]>3000){
    to_sample <- sample.int(dim(data)[1], size = 2000)
  } else {
    to_sample <- 1:(dim(data)[1])
  }
  subdata <- data[to_sample,]
  
  c <- seq(1,5,.5)
  listw <- list() ; mc <- list() ; test <- list()
  
  for (k in 1:length(c)){
    if (VERBOSE){
      cat(lines.to.cat)
      cat(paste("\n  - Consider Neighboors if d <=", c[k]*100, "km"))
      cat("\n      - Getting nearest neighboors")
    }
    nb_list <- spdep::dnearneigh(x = st_coordinates(subdata),
                                 d1 = 0,
                                 d2 = c[k]*100,
                                 longlat = T)
    if (VERBOSE){
      msg <- "\n      - Getting weight list from nearest neighboors" ; cat(msg)
    }
    listw[[k]] <- spdep::nb2listw(nb_list, zero.policy = T)
    
    if (VERBOSE){
      msg <- "\n      - Performing Moran's I tests" ; cat(msg)
    }
    mc[[k]] <- moran.mc(subdata$Kn, listw = listw[[k]], nsim = 99, zero.policy = T)
    # test[[k]] <- moran.test(subdata$Kn, listw = listw[[k]], zero.policy = T) 
    
  }
  
  morans.mc <- unlist(lapply(mc, function(x) x$statistic))
  # morans.test <- unlist(lapply(test, function(x) x$estimate[1]))
  
  # toplot <- data.frame(d = rep(c, 2),
  #                      I = c(morans.mc, morans.test),
  #                      moran_type = c(rep("Kn mc",length(c)),
  #                                     rep("Kn test",length(c))))
  toplot <- data.frame(d = c,
                       I = morans.mc,
                       moran_type = rep("Kn",length(c)))
  
} else if (part == 2){
  
  #' @2.1 Test spatial autocorrelation on the residuals
  if (VERBOSE){
    cat(lines.to.cat)
    msg <- "\n\n~~~~ Performing Moran's I tests on model residuals ~~~~" ; cat(msg) ; lines.to.cat <- c(lines.to.cat, msg)
  }
  
  c <- seq(1,5,.5)
  mc_gam1 <- list() ; test_gam1 <- list()
  
  for (k in 1:length(c)){
    if (VERBOSE){
      cat(lines.to.cat)
      cat(paste("\n  - Consider Neighboors if d <=", c[k]*100, "km"))
      msg <- "\n      - Performing Moran's I tests" ; cat(msg)
    }
    mc_gam1[[k]] <- moran.mc(resid(gam1)[to_sample], listw = listw[[k]], nsim = 99, zero.policy = T)
    # test_gam1[[k]] <- moran.test(resid(gam1)[to_sample], listw = listw[[k]], zero.policy = T)
    
  }
  
  morans.mc.res <- unlist(lapply(mc_gam1, function(x) x$statistic))
  # morans.test.res <- unlist(lapply(test_gam1, function(x) x$estimate[1]))
  
  #' @2.2 Generate plot
  # toplot2 <- data.frame(d = rep(c, 2),
  #                       I = c(morans.mc.res, morans.test.res),
  #                       moran_type = c(rep("GAM residuals mc",length(c)),
  #                                      rep("GAM residuals test",length(c))))
  toplot2 <- data.frame(d = c,
                        I = morans.mc.res,
                        moran_type = rep("GAM residuals",length(c)))
  
  bind_rows(toplot, toplot2) -> toplot
  
  p_moran_residuals <- ggplot(data = toplot, aes(x=d*100, y=I, group = moran_type, color = moran_type)) + 
    geom_point()+ geom_line()+
    scale_color_discrete("Data used to\n calculate Moran's I")+
    xlab("Distance (km)")+
    ylab("Moran's I")
  
  saveRDS(p_moran_residuals, moranPlotName)
  
  # rm(mc, test, nb_list, listw, part) ; invisible(gc())
  rm(mc, mc_gam1, nb_list, listw, part) ; invisible(gc())
  
}

