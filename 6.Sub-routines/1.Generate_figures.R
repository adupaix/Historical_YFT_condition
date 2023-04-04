#'#*******************************************************************************************************************
#'@author : Amael DUPAIX
#'@update : 2021-11-17
#'@email : amael.dupaix@ens-lyon.fr
#'#*******************************************************************************************************************
#'@description :  Script generating the 2 first figures
#'#*******************************************************************************************************************
#'@revisions :
#' 2022/08/22: for submission as a Rapid communication in the Canadian Journal of Fisheries and Aquatic Science
#'             the ggplot noted fig1 becomes the panel A of Figure 2
#'             the ggplot noted fig2 is not used anymore
#'#*******************************************************************************************************************


#' ***************
#' Generate Fig1:
#' ***************
#' Evolution of the condition factor (Kn)
#' according to time

if (VERBOSE){
  msg <- "\n\n~~~~ Generating Figure 1 ~~~~\n" ; lines.to.cat <- c(lines.to.cat, msg) ; cat(lines.to.cat)
}

#' @filter "non duplicated" line:
#'          concerns ~ 400 entries in the YFT data, 156 fish
#'          some rows are duplicated with only the geometry (POINT)
#'          differing
#' @!! in this subroutine, work on "data_fig" not on "data" (only the first of the
#' duplicated rows is kept, as the geographical data is not used)
#' This filter will be applied again on "data" but applying the geometry_method
#'  in sub-routine 3.Get_fishing_location
data %>% dplyr::filter(!duplicated(fish_identifier)) -> data_fig

#' @save the data.frame used for the figure 1 (general data.frame)
saveRDS(data_fig, dfGeneralName)


data_byclass <- list()
spl_sizes <- list()
toplot <- list()

data_byclass[[1]] <- data_fig

data_fig %>% plyr::ddply(.variables = "fishing_year", summarise, n=n()) %>%
  dplyr::filter(n > 50) -> spl_sizes[[1]]

plyr::ddply(data_byclass[[1]], .variables = "fishing_year", summarise, sd = sd(Kn), m = mean(Kn), n = n()) %>%
  dplyr::mutate(se = sd / sqrt(n)) %>%
  dplyr::mutate(group = "all") -> toplot[[1]]

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
    
    data_fig %>% dplyr::filter(fork_length >= l1 & fork_length <= l2 ) %>%
      plyr::ddply(.variables = "fishing_year", summarise, n=n()) %>%
      dplyr::filter(n > 50) -> spl_sizes[[k+1]]
    
    y_of_int <- spl_sizes[[k+1]]$fishing_year
    
    data_fig %>%
      dplyr::filter(fishing_year %in% y_of_int & fork_length >= l1 & fork_length <= l2) -> data_byclass[[k+1]]
    
    plyr::ddply(data_byclass[[k+1]], .variables = "fishing_year", summarise, sd = sd(Kn), m = mean(Kn), n = n()) %>%
      dplyr::mutate(se = sd / sqrt(n)) %>%
      dplyr::mutate(group = size_classes_fig1[k]) -> toplot[[k+1]]
  }
  
  toplot <- bind_rows(toplot)
  toplot$group <- factor(toplot$group, levels = c(size_classes_fig1, "all"))
} else {
  toplot <- toplot[[1]]
}

# get holes in data, to plot a vertical line
diff_following_years <- lead(as.numeric(levels(as.factor(toplot$fishing_year)))) - as.numeric(levels(as.factor(toplot$fishing_year)))
line_pos = which(diff_following_years != 1)+0.5

# prepare the facet wrap + specify lintype
grp_break = levels(toplot$fishing_year)[which(diff_following_years != 1)]
if (length(grp_break) == 1){
  toplot %>% dplyr::mutate(facet_grp = case_when(as.numeric(as.character(fishing_year))<=grp_break ~ 1,
                                                 as.numeric(as.character(fishing_year))>grp_break ~ 2)) %>%
    dplyr::mutate(facet_grp = paste(group, facet_grp, sep = "_"),
                  linetypes = case_when(group=="all"~1,
                                        group!="all"~2))-> toplot
}

if (length(size_classes_fig1) != 0){
  fig1 <- ggplot(toplot, aes(x = as.factor(fishing_year), y = m, shape = group, color = group, group = group))+
    geom_errorbar(aes(ymin = m - se, ymax = m + se), width = 0.2, size = 0.5, position = position_dodge(0.25))+
    geom_line(aes(group = facet_grp, linetype = group))+
    geom_point(position = position_dodge(0.25), size = 1.5)+
    # facet_wrap(~facet_grp, strip.position = NULL, scales = "free_x")+
    scale_shape_manual(values = 16:19)+
    scale_color_manual(values = c(RColorBrewer::brewer.pal(3, "Set1"), "grey20"))+
    scale_linetype_manual(values = c(2,2,2,1))+
    geom_vline(aes(xintercept = line_pos))+
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust=1),
          panel.border = element_rect(color = "black", fill = NA))+
    ylab("Mean Kn")+
    xlab("Fishing year") +
    labs(shape = "Size class", linetype = "Size class", color = "Size class")
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

saveRDS(fig1, fig1rdsName)
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
  
  p3 <- figure2(data = data_fig,
                var.to.compare = "Kn",
                var.grp = "fishing_mode",
                levels.var.grp = c("FOB","FSC"),
                var.x = "fishing_year",
                scale.color.title = "School\n type",
                xlabel = "Year",
                vline = c(1.5,2.5),
                with.second.panel = F)
  
  ggsave(fig2Name, p3, width = 10, height = 6)
}