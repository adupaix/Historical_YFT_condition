#'#*******************************************************************************************************************
#'@author : Amael DUPAIX
#'@update : 2021-11-17
#'@email : amael.dupaix@ens-lyon.fr
#'#*******************************************************************************************************************
#'@description :  Script generating the 2 first figures
#'#*******************************************************************************************************************
#'@revisions
#'#*******************************************************************************************************************


#' ***************
#' Generate Fig1:
#' ***************
#' Evolution of the condition factor (Kn)
#' according to time

if (VERBOSE){
  msg <- "\n\n~~~~ Generating Figure 1 ~~~~\n" ; lines.to.cat <- c(lines.to.cat, msg) ; cat(lines.to.cat)
}

data_byclass <- list()
spl_sizes <- list()
toplot <- list()

data_byclass[[1]] <- data

data %>% ddply(.variables = "fishing_year", summarise, n=n()) %>%
  filter(n > 50) -> spl_sizes[[1]]

ddply(data_byclass[[1]], .variables = "fishing_year", summarise, sd = sd(Kn), m = mean(Kn), n = n()) %>%
  mutate(se = sd / sqrt(n)) %>%
  mutate(group = "all") -> toplot[[1]]

if (length(size_classes_fig1 != 0)){
  for (k in 1:length(size_classes_fig1)){
    
    if (k %in% grep("-", size_classes_fig1)){
      l1 <- as.numeric(sub("-.*", "", size_classes_fig1[k]))
      l2 <- as.numeric(sub(".*-", "", size_classes_fig1[k]))
    } else if (k %in% grep(">", size_classes_fig1)){
      l1 <- as.numeric(sub(">", "", size_classes_fig1[k]))
      l2 <- Inf
    } else if (k %in% grep("<", size_classes_fig1)){
      l1 <- 0
      l2 <- as.numeric(sub("<", "", size_classes_fig1[k]))
    }
    
    data %>% filter(fork_length >= l1 & fork_length <= l2 ) %>%
      ddply(.variables = "fishing_year", summarise, n=n()) %>%
      filter(n > 50) -> spl_sizes[[k+1]]
    
    y_of_int <- spl_sizes[[k+1]]$fishing_year
    
    data %>%
      filter(fishing_year %in% y_of_int & fork_length >= l1 & fork_length <= l2) -> data_byclass[[k+1]]
    
    ddply(data_byclass[[k+1]], .variables = "fishing_year", summarise, sd = sd(Kn), m = mean(Kn), n = n()) %>%
      mutate(se = sd / sqrt(n)) %>%
      mutate(group = size_classes_fig1[k]) -> toplot[[k+1]]
  }
  
  toplot <- bind_rows(toplot)
  toplot$group <- factor(toplot$group, levels = c(size_classes_fig1, "all"))
} else {
  toplot <- toplot[[1]]
}

# get holes in data, to plot a vertical line
diff_following_years <- lead(as.numeric(levels(as.factor(toplot$fishing_year)))) - as.numeric(levels(as.factor(toplot$fishing_year)))
line_pos = which(diff_following_years != 1)+0.5

if (length(size_classes_fig1) != 0){
  fig1 <- ggplot(toplot, aes(x = as.factor(fishing_year), y = m, color = group, group = group))+
    geom_point(size = 0.75, position = position_dodge(0.25))+
    scale_color_brewer("Size class", palette = "Set1")+
    geom_vline(aes(xintercept = line_pos))+
    geom_errorbar(aes(ymin = m - se, ymax = m + se), width = 0.2, size = 0.75, position = position_dodge(0.25))+
    theme(axis.text.x = element_text(angle = 90),
          panel.border = element_rect(color = "black", fill = NA))+
    ylab("Mean Kn")+
    xlab("Fishing year")
} else {
  fig1 <- ggplot(toplot, aes(x = as.factor(fishing_year), y = m))+
    geom_point(size = 0.75, position = position_dodge(0.25))+
    scale_color_brewer("Size class", palette = "Set1")+
    geom_vline(aes(xintercept = line_pos))+
    geom_errorbar(aes(ymin = m - se, ymax = m + se), width = 0.2, size = 0.75, position = position_dodge(0.25))+
    theme(axis.text.x = element_text(angle = 90),
          panel.border = element_rect(color = "black", fill = NA),
          axis.text = element_text(size = 16),
          axis.title = element_text(size = 18))+
    ylab("Mean Kn")+
    xlab("Fishing year")
}


ggsave(fig1Name, fig1, width = 8, height = 6)

# rm(data_by_class) ; invisible(gc())

#' ***************
#' Generate Fig2:
#' ***************
#' Comparison of the condition factor (Kn) between
#' FSC and FAD associated schools

if (fishing_mode_for_model == "all"){
  
  if (VERBOSE){
    msg <- "\n\n~~~~ Generating Figure 2 ~~~~\n" ; cat(msg) ; lines.to.cat <- c(lines.to.cat, msg)
  }
  
  p3 <- figure2(data = data,
                var.to.compare = "Kn",
                var.grp = "fishing_mode",
                levels.var.grp = c("DFAD","FSC"),
                var.x = "fishing_year",
                scale.color.title = "School\n type",
                xlabel = "Year",
                vline = c(1.5,2.5))
  
  ggsave(fig2Name, p3, width = 7, height = 10)
}