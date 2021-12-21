
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



#' @function | Built 2 plots:
#'     - one with value of the coefficient as a function of the levels of the explanatory variable of interest
#'     - one with a violin plot of the coefficients
#' @arguments:
#'  coef.df (list) of dataframes: Each of the data frame contains in columns the coefficients of each level of a variable
#'           and in line each model iteration
#'  p.val.df (list) of dataframes: Each of the data frame contains in columns the p-values of each level of a variable
#'           and in line each model iteration
#'  var (chr) : the name of the variable of interest
#'  levels (vector): levels of the variable of interest, if we want to reorder them on the plot
#'  level_ref (chr): level of reference the variable of interest. The plot is represented with the coefficient of this level set to 0
#'  labelx (chr): legend to print of the x axis of the plot
#'  output_path (chr): path to the folder where the plot will be saved

plot.coef.df <- function(coef.df, p.val.df, var, levels = NULL, level_ref = NULL,
                         labelx = var, output_path){
  
  # get the dataframe of interest in the list
  coef.df <- as.data.frame(coef.df[[var]])
  if (dim(coef.df)[2]==1){
    names(coef.df) <- levels[which(levels != level_ref)]
  } else {
    names(coef.df) <- gsub(pattern = paste0(".*",substr(var, nchar(var)-1, nchar(var))),
                           "", names(coef.df))
  }
  
  # get the columns of interest in the df with coefficients
  p.val.df <- as.data.frame(p.val.df[[var]])
  if (dim(p.val.df)[2]==1){
    names(p.val.df) <- levels[which(levels != level_ref)]
  } else {
    names(p.val.df) <- gsub(pattern = paste0(".*",substr(var, nchar(var)-1, nchar(var))),
                            "", names(p.val.df))
  }
  
  if (!is.null(levels)){
    if (length(levels) != dim(coef.df)[2]+1){
      stop("Error: wrong levels provided. Please note that the reference level also needs to be provided.")
    } else {
      if (is.null(level_ref)){
        stop("Error: please specify the level of reference")
      } else {
        # keep the names of the levels with coefficients different from 0
        nm <- names(coef.df)
        
        # add a column with the reference level
        coef.df <- data.frame(cbind(rep(0,dim(coef.df)[1]),
                                    coef.df))
        p.val.df <- data.frame(cbind(rep(1,dim(p.val.df)[1]),
                                     p.val.df))
        # rename coef.df properly, getting the missing name from the provided variable levels
        names(coef.df) <- c(levels[which(!levels %in% nm)],
                            nm)
        names(p.val.df) <- c(levels[which(!levels %in% nm)],
                            nm)
        
        # reorder the columns according to the provided order in levels
        coef.df <- coef.df[,levels]
        p.val.df <- p.val.df[,levels]
        
        coef.df <- as.data.frame(apply(coef.df, 2, function(x) x-coef.df[,level_ref]))
        
        write.csv(coef.df, file.path(output_path, paste0("coef_figure_",var,".csv")))
      }
    }
  } else {
    levels <- names(coef.df)
  }
  
  data.frame(var = names(coef.df),
             coeff = apply(coef.df, 2, mean),
             sd = apply(coef.df, 2, sd),
             n = dim(coef.df)[1]) %>%
    mutate(se = sd/sqrt(n)) -> df
  
  
  
  #' significance
  signif.df <- as.data.frame(apply(p.val.df, c(1,2), function(x) x<=0.05))
  signi <- apply(signif.df, 2, function(x) 100*sum(x)/length(x))
  l = which(names(signi) == level_ref)
  signi <- as.character(format(signi, digits = 3))
  signi[l] <- ""
  
  rang <- range(coef.df)[2]-range(coef.df)[1]
  
  annotation <- data.frame(x = factor(df$var, levels = levels),
                           y = max(coef.df)+rang/10,
                           signi = signi)
  
  p1 <- ggplot()+
    geom_point(data = df, aes(x=factor(var, levels = levels),
                              y = coeff))+
    geom_errorbar(data = df, aes(x=factor(var, levels = levels),
                                 ymin = coeff - sd,
                                 ymax = coeff + sd),
                  width = 0.25/length(levels), size = 0.5)+
    geom_text(data = annotation,
              aes(x = x, y = y, label = signi),
              angle = 45, color = "grey40")+
    # scale_y_continuous(limits = c(-0.15, 0.15))+
    geom_hline(yintercept = 0)+
    xlab(labelx)+ylab("GAM coefficient")+
    theme(axis.text.x = element_text(angle = 90),
          panel.border = element_rect(color = "black", fill = NA),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18))
  
  coef.df %>% pivot_longer(cols = names(coef.df), names_to = "var", values_to = "coeff") -> df2 
  
  p2 <- ggplot()+
    geom_hline(yintercept = 0)+
    geom_violin(data = df2, aes(x=factor(var, levels = levels),
                                y = coeff))+
    geom_point(aes(x=level_ref, y=0))+
    geom_text(data = annotation,
              aes(x = x, y = y, label = signi),
              angle = 45, color = "grey40")+
    xlab(labelx)+ylab("GAM coefficient")+
    theme(axis.text.x = element_text(angle = 90),
          panel.border = element_rect(color = "black", fill = NA),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18))
  
  return(list(p1,p2))
}
