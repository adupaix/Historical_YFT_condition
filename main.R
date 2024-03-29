#'#*******************************************************************************************************************
#'@author : Amael DUPAIX
#'@update : 2021-11-17
#'@email : amael.dupaix@ens-lyon.fr
#'#*******************************************************************************************************************
#'@description :  Main script to determine if the YFT condition factor (Kn) has been impacted by the introduction of
#' FADs in the Indian Ocean.
#' Launched by configuration files which are in the cfg folder
#'#*******************************************************************************************************************
#' @tags:
#' the @for_study tag is used before each argument, to specify which arguments were used to obtain the study results
#' the @filter tag is used to make explicit where the different filters are applied (not exactly in agreement with
#'             the names of the sub-routines)
#'#*******************************************************************************************************************
#'@revisions
#'#*******************************************************************************************************************


#' **********
#' Initialize:
#' **********
source(file.path(ROUT_PATH, "0.Init.R"))


#' **********
#' Read data:
#' **********
data1 <- read.table(file = file.path(DATA_PATH, "Tunabio_OI_20210826_env.csv"),
                    sep = ",", dec=".",stringsAsFactors = F,
                    quote="\"", row.names = NULL,
                    h=T, colClasses = "character")
data2 <- read.table(file = file.path(DATA_PATH, "Tunabio_OI_20210826_specimen.csv"),
                    sep = ";", dec=",",stringsAsFactors = F,
                    quote="\"", row.names = NULL,
                    h=T, colClasses = "character")

data <- base::merge(data1, data2, by = "fish_identifier", all = T)

rm(data1, data2)

#' **********
#' Prep data:
#' **********
#' @filter measured FL, measured W, fish caught in the IO
data <- prep_wl_data(DATA_PATH, data, calcfdl = calcfdl, read = !RESET, getgeom = T, ncores = nb_cores,
                     size_class_levels = size_class_levels, summaryName = summaryName, verbose = VERBOSE)

#' ***************
#' Calculate Kn:
#' ***************

if (VERBOSE){
  msg <- "\n\n~~~~ Calculating Kn ~~~~\n" ; cat(msg) ; lines.to.cat <- c(lines.to.cat, msg)
}

# replace DFAD by FOB in fishing_mode
data$fishing_mode[which(data$fishing_mode == "DFAD")] <- "FOB"

#' @filter species of interest, caught by purse seiner
data %>% dplyr::filter(species_code_fao %in% species & gear_code == "PS") -> data

if (allom == "Chassot2015"){
  a = 0.00002459
  b = 2.96670
} else if (allom == "fromData"){
  a = exp(lm(log(whole_fish_weight) ~ log(fork_length), data = data)$coefficients[1])
  b = lm(log(whole_fish_weight) ~ log(fork_length), data = data)$coefficients[2]
}

if(allom == "fromData"){
  lmodel <- lm(log(whole_fish_weight) ~ log(fork_length), data = data)
  sink(allomSummaryName)
  cat("\n\nSummary of the fit of the allometric relationship")
  cat("\n---------------------------------------------------\n")
  cat("\na =", a)
  cat("\nb =", b)
  cat("\nr-square :", summary(lmodel)$r.squared)
  cat("\nFish weight ranged from", min(data$whole_fish_weight),
      "kg for a", min(data$fork_length[which(data$whole_fish_weight == min(data$whole_fish_weight))]),
      "cm fish, to", max(data$whole_fish_weight),
      "kg for a", data$fork_length[which(data$whole_fish_weight == max(data$whole_fish_weight))],
      "cm fish.\nThe average size was", mean(data$fork_length), "cm (SD =", sd(data$fork_length),
      ") and the average weight was", mean(data$whole_fish_weight), "kg (SD =", sd(data$whole_fish_weight), ").")
  sink()
  
  x <- seq(min(data$fork_length),max(data$fork_length), 0.1)
  y <- a*x^b
  
  p <- ggplot()+
    geom_point(data=data, aes(x=fork_length, y=whole_fish_weight), col = "grey40", alpha = 0.5)+
    geom_line(aes(x=x,y=y), col = "red")+
    labs(x = "Fork length (cm)",
         y = "Fish weight (kg)")
  
  ggsave(filename = allomFitName, p)
}

data %>%
  dplyr::mutate(weight_th = a*fork_length^b,
                Kn = whole_fish_weight / weight_th) -> data

sink(summaryName, append = T)
cat("\n\n\nSpecies and gear filter")
cat("\n------------------------\n")
cat(paste0("\n    - Number of entries after gear filter + species filter (", species, "): ", dim(data)[1]))
sink()

#' @filter fishing location not missing
#' Keep only the data entries for which any type of geometry is available (POINT or MULTIPOINT)
data %>% dplyr::filter(!st_is_empty(data$geometry)) -> data

sink(summaryName, append = T)
cat("\n\n\nMissing location filter")
cat("\n-----------------------\n")
cat(paste0("\n    - Number of entries after missing location filter: ", dim(data)[1]))
sink()

#' @filter some fish by there identifier. Explanations:
#'     - A368: among the potential fishing locations, one is in the Atlantic
#'             either it's an error, or, if it's not, we cannot be certain that the
#'             fish comes from the Indian Ocean
#'             @update: corrected on the xlsx and csv files on 12 dec 2021
data %>% dplyr::filter(!(fish_identifier %in% c())) -> data

sink(summaryName, append = T)
cat("\n\n\nChoosen fish filter")
cat("\n-----------------------\n")
cat(paste0("\n    - Number of entries after choosen fish filter: ", dim(data)[1]))
sink()


#' ********************
#' Get fishing date:
#' ********************
#' at that point if the entries with missing dates are discarded
#' @filter known fishing date interval, no error in date,
#'         fish caught prior to 2020
if(deduce_date == F){
  source(file.path(ROUT_PATH, "2.Calc_fishing_date.R"))
}

#' ********************
#' Generate Figures:
#' ********************
#' @filter "non duplicated" line:
#'          concerns ~ 400 entries in the YFT data, 156 fish
#'          some rows are duplicated with only the geometry (POINT)
#'          differing
#' @!! in this subroutine, work on "data_fig" not on "data" (only the first of the
#' duplicated rows is kept, as the geographical data is not used)
#' This filter will be applied again on "data" but applying the geometry_method
#'  in sub-routine 3.Get_fishing_location
if (generate_plots){
  source(file.path(ROUT_PATH, "1.Generate_figures.R"))
} else {
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
}



if (cluster == T){
  #' ******************************
  #' Prepare for datarmor script
  #' ******************************
  # save all the objects of the environments to a list
  env <- mget(ls())
  fname <- file.path(OUTPUT_PATH, "all_objects.rds")
  # and save that list as an .rds file (which will be loaded again afterwards)
  saveRDS(env, fname)
  
  #' *************************
  #' Generate the command list
  #' *************************
  commands <- paste0("Rscript ",
                     file.path(WD, "main_cluster.R "),
                     fname, " ",
                     1:nb_of_times_to_run)
  
  file.create(file.path(WD, paste0("commands_bootstrap",output_suffixe,".txt")))
  cmds <- file(file.path(WD, paste0("commands_bootstrap",output_suffixe,".txt")), open = "w")
  writeLines(commands, cmds)
  close(cmds)
  
  #source(file.path(ROUT_PATH, "6.Datarmor_prep.R"))
  
} else {
  
  for (i in 1:length(seeds)){
    
    seed.i <- seeds[i]
    
    #' ********************
    #' Init in @for loop:
    #' ********************
    #' first, read and save data, so the sampling can be performed again
    #'  in the next iteration
    if (i != 1){
      data <- readRDS(intermediateDataName)
    }
    saveRDS(data, file = intermediateDataName)
    
    source(file.path(ROUT_PATH, "0.i.Init_in_loop.R"))
    
    
    #' ********************
    #' Get fishing date:
    #' ********************
    #' at that point if the missing dates are deduced
    #' @filter known fishing date interval, no error in date, fishing location not missing,
    #'         fish caught prior to 2020
    if(deduce_date == T){
      source(file.path(ROUT_PATH, "2.Calc_fishing_date.R"))
    }
    
    #' ********************
    #' Get fishing location:
    #' ********************
    #' @filter "non duplicated" line:
    #'          concerns ~ 400 entries in the YFT data, 156 fish
    #'          some rows are duplicated with only the geometry (POINT)
    #'          differing
    source(file.path(ROUT_PATH, "3.Get_fishing_location.R"))
    
    
    #' ***************
    #' Perform GAM:
    #' ***************
    source(file.path(ROUT_PATH, "4.GAM.R"))
    
  }
  
  file.remove(intermediateDataName)
}

