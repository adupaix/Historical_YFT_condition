
points.to.multipoints <- function(data_sf, fact){
  geom <- split(x = data_sf,
                f = fact)
  geom <- Map(st_combine, geom)
  geom <- do.call(c, geom)
  
  data_sf <- data_sf[!duplicated(fact),]
  st_geometry(data_sf) <- geom
  
  return(data_sf)
}


geary.hinkley.transform <- function(Z,X,Y, cor.method = "pearson"){
  (mean(Y)*Z - mean(X))/(sqrt(sd(Y)^2*Z^2-2*cor(X,Y, method = cor.method)*sd(X)*sd(Y)*Z + sd(X)^2))
}


get.coeff <- function(gam, var, type = c("coeff","se","p.value")){
  if(type == "coeff"){
    coeff <- coefficients(gam)
  } else if (type == "se"){
    coeff <- summary(gam)$se
  } else if (type == "p.value"){
    coeff <- summary(gam)$p.pv
  }
  
  coeff <- coeff[grep(var, names(coeff))]
  names(coeff) <- gsub(pattern = paste0(".*",substr(var, nchar(var)-1, nchar(var))),
                       "", names(coeff))
  return(coeff)
}


plot.coeff <- function(gam, var, levels_order = NULL,
                       labelx = var){
  
  coeff <- get.coeff(gam, var, type = "coeff")
  se <- get.coeff(gam, var, type = "se")
  p.value <- get.coeff(gam, var, type = "p.value")
  
  if (!is.null(levels_order)){
    if (length(levels_order) != length(coeff)){
      stop("Error: wrong levels provided")
    } else {
      coeff <- coeff[levels_order]
      se <- se[levels_order]
      p.value <- p.value[levels_order]
    }
  } else {
    levels_order <- names(coeff)
  }
  
  my_colors <- ifelse(p.value < 0.05 / length(p.value), "black", "grey40")
  
  data.frame(coeff) %>% tibble::rownames_to_column(var = var) %>%
    mutate(se = se, col = my_colors) -> df
  
  p <- ggplot()+
    geom_point(data = df, aes(x=factor(!!rlang::sym(var), levels = levels_order),
                              y = coeff), color = df$col)+
    # geom_errorbar(data = df, aes(x=factor(!!rlang::sym(var), levels = levels_order),
    #                              ymin = coeff - se,
    #                              ymax = coeff + se),
    #               color = df$col, width = 0.25/length(levels_order), size = 0.5)+
    # scale_y_continuous(limits = c(-0.15, 0.15))+
    geom_hline(yintercept = 0)+
    xlab(labelx)+ylab("GAM coefficient")+
    theme(axis.text.x = element_text(angle = 90),
          panel.border = element_rect(color = "black", fill = NA),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18))
  
  return(p)
}


qqplot.gam <- function(residus){
  ggplot2::ggplot() +
    ggplot2::stat_qq(ggplot2::aes(sample = residus))+
    ggplot2::stat_qq_line(ggplot2::aes(sample = residus), col = 2)+
    ggplot2::labs(x = "Theoretical quantiles", y = "Sample quantiles")
}


stepAIC.gam <- function(gam, verbose = F){
  
  x <- unlist(strsplit(as.character(gam$formula)[3], " \\+ "))
  k=1
  l=length(x)
  
  while(k <= l){
    
    y <- as.character(gam$formula)[2]
    x <- unlist(strsplit(as.character(gam$formula)[3], " \\+ "))
    df <- as.character(gam$call)[3]
    
    AIC_init <- AIC(gam)
    
    formula <- paste("gam(",y,"~",paste(x, collapse = " + "),", data =",df,")")
    if (verbose){
      cat("~~~ Iteration", k, "~~~\n\nInitial model:", formula)
      cat("\nAIC:", AIC_init, "\n\n")
    }
    
    remove.x <- rep(F, length(x))
    
    if (length(x)>1){
      for (i in 1:length(x)){
        new_formula <- paste("gam(",y,"~",paste(x[-i], collapse = " + "),", data =",df,")")
        AIC.i <- eval(parse(text = paste("AIC(",new_formula,")")))
        
        if(verbose){
          cat(new_formula)
          cat("\nAIC:", AIC.i, "\n\n")
        }
        
        if (AIC.i < AIC_init){
          remove.x[i] <- T
        }
      }
    } else {
      new_formula <- paste("gam(",y,"~ 1, data =",df,")")
      AIC.i <- eval(parse(text = paste("AIC(",new_formula,")")))
      
      if(verbose){
        cat(new_formula)
        cat("\nAIC:", AIC.i, "\n\n")
      }
      
      if (AIC.i < AIC_init){
        remove.x <- T
      }
    }
    
    
    if(any(remove.x)){
      if(length(x)>1){x <- x[!remove.x]} else {x <- 1}
      new_formula <- paste("gam(",y,"~",paste(x, collapse = " + "),", data =",df,")")
      gam <- eval(parse(text = new_formula))
      k <- k+1
    } else {
      k=l+1
    }
  }
  
  return(gam)
}

