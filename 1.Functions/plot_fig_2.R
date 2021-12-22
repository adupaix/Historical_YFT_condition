#'#*******************************************************************************************************************
#'@author : Amael DUPAIX
#'@update : 2021-10-15
#'@email : amael.dupaix@ens-lyon.fr
#'#*******************************************************************************************************************
#'@description :  Function to generate figure 2
#' data (data frame): contains all the data to plot
#' var.to.compare (chr): name of the column which will be plotted in the y axis
#' var.grp (chr): name of the column used to group the values from var.to.compare. The column needs to be as factor
#' levels.var.grp (vector of chr): levels of interest of var.grp (values which are not in this levels are discarded)
#' var.x (chr): name of the column containing the data which will be on the x axis
#' 
#' scale.color.title (chr): for the plot, title of the colour scale
#' xlabel (chr): for the plot, x axis title
#' vline (num): for the plot, give the xintersect of the vlines added to the plot
#'#*******************************************************************************************************************
#'@revisions
#'#*******************************************************************************************************************

figure2 <- function(data, var.to.compare, var.grp, levels.var.grp, var.x,
                    scale.color.title = "School\n type",
                    xlabel = "Year",
                    vline = 0
){
  
  data %>% dplyr::filter((!!rlang::sym(var.grp)) %in% levels.var.grp) %>%
    plyr::ddply(.variables = c(var.x,var.grp), summarise, n=n()) %>%
    tidyr::spread((!!rlang::sym(var.grp)), n) %>%
    dplyr::filter(!is.na((!!rlang::sym(levels.var.grp[1]))) & !is.na((!!rlang::sym(levels.var.grp[2])))) -> spl_sizes
  
  var.x_of_int <- spl_sizes[,var.x]
  
  data %>% 
    dplyr::filter((!!rlang::sym(var.x)) %in% var.x_of_int & (!!rlang::sym(var.grp)) %in% levels.var.grp) -> dat_
  
  dat_ %>% plyr::ddply(.variables = c(var.x,var.grp), summarise, m = mean(!!rlang::sym(var.to.compare)),
                 sd = sd(!!rlang::sym(var.to.compare)), n = n()) %>%
    dplyr::mutate(se = sd / sqrt(n)) -> toplot
  
  p <- c()
  
  for (i in 1:length(var.x_of_int)){
    
    data %>% dplyr::filter(!!rlang::sym(var.x) == var.x_of_int[i] &
                             !!rlang::sym(var.grp) == levels.var.grp[1]) -> dat_.1
    data %>% dplyr::filter(!!rlang::sym(var.x) == var.x_of_int[i] &
                             !!rlang::sym(var.grp) == levels.var.grp[2]) -> dat_.2
    
    p <- c(p, wilcox.test(dat_.1[[var.to.compare]],dat_.2[[var.to.compare]])$p.value)
    # summary_tests$Significant_diff[i] <- F
    # if (wilcox.test(yft.i$whole_fish_weight, yft.i$weight_th)$p.value <= 0.05/ntests){
    #   summary_tests$Significant_diff[i] <- T
    # }
    # 
    # summary_tests$Mean_residual[i] <- mean(yft.i$weight_residuals)
    # 
    
  }
  
  df <- data.frame(var.x_of_int = as.numeric(as.character(var.x_of_int)), 
                   p)
  ntests <- dim(df)[1]
  
  signif_y <- df$var.x_of_int[which(df$p <= 0.05/ntests)]
  
  toplot %>% dplyr::filter(!!rlang::sym(var.x) %in% signif_y) -> signi
  
  axis_col <- ifelse(unique(toplot[[var.x]]) %in% signi[[var.x]],"red","black") 
  axis_face <- ifelse(unique(toplot[[var.x]]) %in% signi[[var.x]],"bold","plain")
  
  p3.1 <- ggplot(toplot, aes(x = as.factor(!!rlang::sym(var.x)), y = m, color = !!rlang::sym(var.grp)))+
    scale_color_brewer(scale.color.title, palette = "Set1")+
    geom_point()+
    geom_errorbar(aes(ymin = m - se, ymax = m + se), width = 0.25, size = 0.5)+
    theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          panel.border = element_rect(color = "black", fill = NA),
          legend.justification = c(1,0),
          legend.position = c(0.9,0.1),
          legend.background = element_rect(colour = "black"))+
    ylab(paste("Mean", var.to.compare))
  
  toplot %>% plyr::ddply(.variables = var.x, .fun = function(x){
    (x %>% dplyr::filter(!!rlang::sym(var.grp) == levels.var.grp[1]))$m -
      (x %>% dplyr::filter(!!rlang::sym(var.grp) == levels.var.grp[2]))$m
  }) -> data_hist
  
  p3.2 <- ggplot(data_hist, aes(x = as.factor(!!rlang::sym(var.x)), y = V1))+
    geom_bar(stat = "identity")+
    theme(axis.text.x = element_text(angle = 90, colour = axis_col, face = axis_face),
          panel.border = element_rect(color = "black", fill = NA))+
    ylab(paste0(var.to.compare," (",levels.var.grp[1],") - ", var.to.compare," (",levels.var.grp[2],")"))+
    xlab(xlabel)
  
  if(vline[1] != 0){
    p3.1 <- p3.1 + geom_vline(xintercept = vline)
    p3.2 <- p3.2 + geom_vline(xintercept = vline)
  }
  
  p3 <- ggdraw()+
    draw_plot(p3.1, 0, 1/3, 1, 2/3)+
    draw_plot(p3.2, 0, 0, 1, 1/3)+
    draw_plot_label(c("A","B"), c(0,0), c(1,1/3))
  
  return(p3)
  
}