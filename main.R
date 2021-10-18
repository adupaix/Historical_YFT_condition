#'#*******************************************************************************************************************
#'@author : Amael DUPAIX
#'@update : 2021-09-27
#'@email : amael.dupaix@ens-lyon.fr
#'#*******************************************************************************************************************
#'@description :  Main script to determine if the YFT condition factor (Kn) has been impacted by the introduction of
#' FADs in the Indian Ocean.
#'#*******************************************************************************************************************
#'@revisions
#'#*******************************************************************************************************************

rm(list = ls())
invisible(gc())

WD <- getwd()

DATA_PATH <- file.path(WD, "0.Data/")

FUNC_PATH <- file.path(WD,"1.Functions")
OUTPUT_PATH <- file.path(WD, "3.Outputs", format(Sys.time()))

PLOT_PATH <- file.path(OUTPUT_PATH, "Plots")

#' @reproductibility
nb_of_times_to_run <- 1
SEED <- 123456

#'@arguments:
#'#**********
#'# Run in parallel ?
#' (in prep_wl_data, to read geometry)
#' First element of the vector:
#'    If F, runs in sequential
#'    if T, runs in parallel
#' Second element of the vector: fraction of the cores to be used
Parallel = c(T, 1/2)

#' If RESET is FALSE, try to read the prepped data
#' else, prep the data in any case
RESET = F

#' When missing, choose to sample fork length (FL) from fish with the same first dorsal length (FDL)
#'  or not
calcfdl = F

#' Read the geometry column or not
#' need to read from the chr column and change it in geometry format
getgeom = T

#' Choose the species of interest
#' one of YFT, BET, or SKJ
species = "YFT"

#' Choose which allometric relationship is going to be used to calculate the condition factor
#' one of "Chassot2015" of "fromData"
allom = "Chassot2015"

#' @Figure_1
#' size classes for filtering data
#' Used to generate Figure 1
size_classes = c("40-60","<102",">102")

size_class_for_model <- "all"

fad_fsc <- T

#' ***************
#' Get libraries:
#' ***************
source(file.path(FUNC_PATH, "install_libraries.R"))

srcUsedPackages <- c("plyr", "dplyr","tidyr","foreach","doSNOW","stringr","lubridate","sf",
                     "parallel","ggplot2","tibble","cowplot","RColorBrewer", "MASS","truncnorm",
                     "mgcv", "spdep", "gratia")

installAndLoad_packages(srcUsedPackages, loadPackages = TRUE)


#' **********
#' Initialize:
#' **********
# For parallel:
# On Windows, or if don't want to parralelize, set cores number to 1
if (.Platform$OS.type == "windows" | as.logical(Parallel[1]) == F) {
  nb_cores = 1
} else { #use a fraction of the available cores
  nb_cores = trunc(detectCores() * as.numeric(Parallel[2]))
}
# Source functions
source(file.path(FUNC_PATH, "0.Prep_wl_data.R"))
source(file.path(FUNC_PATH, "plot_fig_2.R"))

if (nb_of_times_to_run != 1){
  seeds <- round(runif(nb_of_times_to_run, 1, 10^6))
} else {
  seeds = SEED
}


# glm_summaries <- data.frame(matrix(ncol = 3*9+1, nrow = length(seeds), data = NA))
# names(glm_summaries) <- c("seed",
#                           paste(rep(c("(Intercept)","fishing_quarter2","fishing_quarter3",
#                                       "fishing_quarter4", "fishing_year", "scaled_FL",
#                                       "scaled_lon","scaled_lat", "scaled_lon:scaled_lat"), each = 3),
#                                 rep(c("Estimate","Std. Error","Pr(>|t|)"), times = 9),
#                                 sep = "__"))

# Initialize names
try(dir.create(PLOT_PATH, showWarnings = F, recursive = T))

#Initialize data processing summary name
summaryName <- file.path(OUTPUT_PATH, paste0("Data_processing_summary_RESET-",RESET,".txt"))

fitName <- file.path(PLOT_PATH, "fit_lnorm_dates_goodnessoffit.png")
histFitName <- file.path(PLOT_PATH, "fit_lnorm_dates.png")
fig1Name <- file.path(PLOT_PATH, "Figure1.png")
fig2Name <- file.path(PLOT_PATH, "Figure2.png")

vars <- c("month", "quarter", "lonlat", "lon", "lat")
plotsNames <- file.path(PLOT_PATH, paste0("Kn_f-", vars, ".png"))


for (i in 1:length(seeds)){
  
  set.seed(seeds[i])
  
  # glm_summaries[i,1] <- seeds[i]
  
  msg <- "\14" ; cat(msg) ; lines.to.cat <- c(msg)
  msg <- paste("-------- ITERATION", i, "/", length(seeds), " - SEED NB =", seeds[i], "--------") ; cat(msg) ; lines.to.cat <- c(lines.to.cat, msg)
  
  # Initialize names
  diagnoPlotName <- file.path(PLOT_PATH, paste0("diagn_glm_", seeds[i], ".png"))
  moranPlotName <- file.path(PLOT_PATH, paste0("Moran_I_plot_",seeds[i],".png"))
  gamSummary <- file.path(OUTPUT_PATH, paste0("gam_summary_",seeds[i],".rds"))
  glmSummary <- file.path(OUTPUT_PATH, paste0("glm_summary_",seeds[i],".rds"))
  
  # Initialize data processing summary
  sink(summaryName)
  cat("Execution time:", format(Sys.time()), "\n\n")
  sink()
  
  #' **********
  #' Read data:
  #' **********
  data1 <- read.table(file = file.path(DATA_PATH, "Tunabio_OI_20210826_env.csv"),
                      sep = ";", dec=",",stringsAsFactors = F,
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
  prep_output <- prep_wl_data(DATA_PATH, data, calcfdl = calcfdl, read = !RESET, getgeom = getgeom, ncores = nb_cores,
                              summaryName = summaryName)
  data <- prep_output[[1]]
  lines.to.cat <- prep_output[[2]]
  
  #' ***************
  #' Calculate Kn:
  #' ***************
  
  msg <- "\n\n~~~~ Calculating Kn ~~~~\n" ; cat(msg) ; lines.to.cat <- c(lines.to.cat, msg)
  
  data %>% filter(species_code_fao %in% species & gear_code == "PS") -> data
  
  if (allom == "Chassot2015"){
    a = 0.00002459
    b = 2.96670
  } else if (allom == "fromData"){
    a = exp(lm(log(whole_fish_weight) ~ log(fork_length), data = data)$coefficients[1])
    b = lm(log(whole_fish_weight) ~ log(fork_length), data = data)$coefficients[2]
  }
  
  data %>%
    dplyr::mutate(weight_th = a*fork_length^b,
                  Kn = whole_fish_weight / weight_th) -> data
  
  sink(summaryName, append = T)
  cat("\n\n\nSpecies and gear filter")
  cat("\n------------------------\n")
  cat(paste0("\n    - Number of entries after gear filter + species filter (", species, "): ", dim(data)[1]))
  sink()
  
  
  #' *************************
  #' Calculate fishing dates:
  #' *************************
  
  msg <- "\n~~~~ Calculating fishing date from interval and/or sampling dates ~~~~\n" ; cat(msg) ; lines.to.cat <- c(lines.to.cat, msg)
  
  # keep only the entries for which the exact fishing date is known
  data %>% filter(!is.na(fishing_date)) %>%
    # keep in a variable how the fishing date was obtained
    mutate(fishing_date_origin = "exact") -> data_exact_dates
  
  N1 <- dim(data_exact_dates)[1]
  
  # save another part of the data for which a fishing date interval is known
  data %>% filter(is.na(fishing_date) & !is.na(fishing_date_min) & !is.na(fishing_date_max)) %>%
    mutate(fishing_date_origin = "interval") -> data_interval_dates
  
  N2 <- dim(data_interval_dates)[1]
  
  # save the other part of the data frame (where the fishing date is to be deduced)
  data %>% filter(is.na(fishing_date) & is.na(fishing_date_min) & is.na(fishing_date_max)) %>%
    # keep in a variable if the fishing date was the real one or not
    mutate(fishing_date_origin = "deduced") -> data_no_dates
  
  #' sometimes, the fish can be sampled onboard directly
  #' hence, we consider that there s an error only if
  #'       1. if the beginning of the fishing dates interval is after the sampling date
  #'       2. if there is more than 1000 days between the one of the fishing date boundaries and
  #'          the sampling date (often an error in the year, e.g. 2008 entered instead of 2018)
  #'       3. if an exact fishing date is provided and it's after the sampling date
  data_interval_dates %>%
    filter(fishing_date_min > fish_sampling_date | 
             # fishing_date_max >= fish_sampling_date & is.na(fishing_date) |
             # landing_date > fish_sampling_date |
             fish_sampling_date - fishing_date_max > 1000 |
             fish_sampling_date - fishing_date_min > 1000 |
             fish_sampling_date - fishing_date > 1000
    ) %>%
    st_drop_geometry() -> error_interval # errors 1 and 2
  
  data_exact_dates %>%
    filter(fishing_date > fish_sampling_date) %>%
    st_drop_geometry() -> error_exact # errors 3
  
  bind_rows(error_interval, error_exact) %>%
    write.csv(file.path(OUTPUT_PATH,"error_fishing_dates.csv"))
  
  # filter data to remove the above errors
  data_interval_dates %>%
    filter(!(fish_identifier %in% error_interval$fish_identifier)) -> data_interval_dates
  
  data_exact_dates %>%
    filter(!(fish_identifier %in% error_exact$fish_identifier)) -> data_exact_dates
  
  msg <- paste0("    - Number of removed entries (date errors): ", dim(error_interval)[1]+dim(error_exact)[1], "\n") ; cat(msg) ; lines.to.cat <- c(lines.to.cat, msg)
  
  
  # For data with a date interval, take the middle of the interval
  data_interval_dates %>%
    mutate(date1 = fishing_date_min,
           date2 = pmin(fishing_date_max, fish_sampling_date)) %>%
    mutate(fishing_date = date1 + (date2 - date1)/2) %>%
    dplyr::select(-date1, -date2) -> data_interval_dates
  
  
  # bind data with date interval and with exact date
  bind_rows(data_interval_dates, data_exact_dates) %>%
    arrange(fish_identifier) %>%
    mutate(t_fishing_sampling = as.numeric(difftime(fish_sampling_date, fishing_date, units = "days")))-> data_dates
  
  # fit a log normal law to the time between fishing and sampling
  fit_lnorm = fitdistrplus::fitdist(data_dates$t_fishing_sampling, "lnorm")
  
  # plot informations on the fit and save the plots
  png(fitName)
  plot(fit_lnorm)
  dev.off()
  
  fit_plot <- ggplot()+
    geom_histogram(data = data_dates,
                   aes(x=t_fishing_sampling,
                       y=..density..),
                   color = NA,
                   fill = "red",
                   alpha = 0.5,
                   binwidth = 1)+
    stat_function(fun = dlnorm, args = fit_lnorm$estimate, color = "red")+
    xlab("Time between fishing and sampling (days)")+
    ylab("Fraction of entries in the dataset")
  
  ggsave(histFitName, fit_plot)
  
  # for data with missing fishing dates, deduce a fishing date using the above fit
  data_no_dates %>%
    mutate(t_fishing_sampling = rlnorm(dim(data_no_dates)[1],
                                       meanlog = fit_lnorm$estimate[1],
                                       sdlog = fit_lnorm$estimate[2])) %>%
    mutate(fishing_date = fish_sampling_date - as.difftime(t_fishing_sampling, units = "days")) -> data_no_dates
  
  
  #bind data with dates and data with randomly assigned dates
  bind_rows(data_dates, data_no_dates) %>%
    mutate(fishing_year = year(fishing_date),
           fishing_month = month(fishing_date),
           fishing_quarter = quarter(fishing_date, fiscal_start = 12)) -> data
  N3 <- dim(data)[1]
  
  #' @!! keep only data before 2020
  #' sampling on YFT in 2020 was performed in January tuna fished in Nov or Dec 2019
  data %>%
    filter(fishing_year < 2020) -> data
  
  sink(summaryName, append = T)
  cat("\n\n\nDates filter")
  cat("\n------------\n")
  cat("\n  Number of entries with exact fishing date available:", N1)
  cat("\n     Of which",dim(error_exact)[1],"contain an error (removed)")
  cat("\n  Number of entries with fishing date interval:", N2)
  cat("\n     Of which",dim(error_interval)[1],"contain an error (removed)")
  cat("\n  Number of entries with no fishing date available (deduced using the lnorm fit):", dim(data_no_dates)[1])
  cat("\n  Number of entries after 2019 (removed):", N3 - dim(data)[1])
  cat("\n\n    - Number of entries after date filter:",dim(data)[1])
  sink()
  
  rm(data_dates, data_no_dates, data_exact_dates, data_interval_dates, N1, N2, N3) ; invisible(gc())
  
  
  #' ********************
  #' Get fishing location:
  #' ********************
  
  msg <- "\n~~~~ Getting fishing location ~~~~\n" ; cat(msg) ; lines.to.cat <- c(lines.to.cat, msg)
  
  msg <- "    - Sampling unique points from MULTIPOINT geometries:" ; cat(msg) ; lines.to.cat <- c(lines.to.cat, msg)
  
  data %>% filter(geom_type != "MULTIPOINT") %>%
    mutate(geom_sampled = F) -> data_na_point
  
  # marche pas, rempli la RAM
  ## FAIRE TOURNER PAR MORCEAUX
  
  data %>% filter(geom_type == "MULTIPOINT") %>%
    mutate(geom_sampled = T) -> data_multi
  
  sample_size = 1000
  niter <- floor(dim(data_multi)[1]/sample_size) + 1
  
  data_multi_list <- list()
  
  for (k in 1:niter){
    
    cat(lines.to.cat)
    cat("    Sub-sample ",k,"/",niter, "\n")
    
    if(k == niter){
      indexes <- ((k-1)*sample_size+1):(dim(data)[1])
    } else {
      indexes <- ((k-1)*sample_size+1):(k*sample_size)
    }
    
    data_multi[indexes,] %>%
      # change the MULTIPOINT to several lines with a POINT geometry
      st_cast('POINT') %>%
      # keep only one line, randomly
      ddply("fish_identifier", function(x) x[sample.int(x$n_points, 1),],
            .progress = "text") %>%
      # ddply changes the sf data.frame back to a simple data.frame so we change the df back to sf
      st_as_sf() -> data_multi_list[[k]]
  }
  
  
  bind_rows(data_multi_list, data_na_point) %>%
    arrange("fish_identifier") -> data
  
  coordinates <- data.frame(st_coordinates(data))
  
  N1 <- dim(data)[1]
  
  data %>% mutate(lon = coordinates$X,
                  lat = coordinates$Y) %>%
    filter(!is.na(lon) & !is.na(lat)) -> data
  
  N2 <- dim(data)[1]
  Nna <- N1 - N2 ; rm(N1, N2)
  
  sink(summaryName, append = T)
  cat('\n\nPositions filter')
  cat('\n----------------\n')
  cat("\n  Number of entries with exact position available:", dim(data_na_point)[1]-Nna)
  cat("\n  Number of entries with multiple positions available (randomly sampled):", dim(data_multi)[1])
  cat("\n  Number of entries with no available position (removed):", Nna)
  cat("\n\n    - Number of entries after position filter:",dim(data)[1])
  sink()
  
  rm(coordinates, data_multi, data_multi_list, data_na_point, Nna) ; invisible(gc())
  
  
  if (i == 1){
    #' ***************
    #' Generate Fig1:
    #' ***************
    #' Evolution of the condition factor (Kn)
    #' of different size classes according to time
    
    msg <- "\n\n~~~~ Generating Figure 1 ~~~~\n" ; cat(msg) ; lines.to.cat <- c(lines.to.cat, msg)
    
    data_byclass <- list()
    spl_sizes <- list()
    toplot <- list()
    
    data_byclass[[1]] <- data
    
    data %>% ddply(.variables = "fishing_year", summarise, n=n()) %>%
      filter(n > 50) -> spl_sizes[[1]]
    
    ddply(data_byclass[[1]], .variables = "fishing_year", summarise, sd = sd(Kn), m = mean(Kn), n = n()) %>%
      mutate(se = sd / sqrt(n)) %>%
      mutate(group = "all") -> toplot[[1]]
    
    for (k in 1:length(size_classes)){
      
      if (k %in% grep("-", size_classes)){
        l1 <- as.numeric(sub("-.*", "", size_classes[k]))
        l2 <- as.numeric(sub(".*-", "", size_classes[k]))
      } else if (k %in% grep(">", size_classes)){
        l1 <- as.numeric(sub(">", "", size_classes[k]))
        l2 <- Inf
      } else if (k %in% grep("<", size_classes)){
        l1 <- 0
        l2 <- as.numeric(sub("<", "", size_classes[k]))
      }
      
      data %>% filter(fork_length >= l1 & fork_length <= l2 ) %>%
        ddply(.variables = "fishing_year", summarise, n=n()) %>%
        filter(n > 50) -> spl_sizes[[k+1]]
      
      y_of_int <- spl_sizes[[k+1]]$fishing_year
      
      data %>%
        filter(fishing_year %in% y_of_int & fork_length >= l1 & fork_length <= l2) -> data_byclass[[k+1]]
      
      ddply(data_byclass[[k+1]], .variables = "fishing_year", summarise, sd = sd(Kn), m = mean(Kn), n = n()) %>%
        mutate(se = sd / sqrt(n)) %>%
        mutate(group = size_classes[k]) -> toplot[[k+1]]
    }
    
    toplot <- bind_rows(toplot)
    toplot$group <- factor(toplot$group, levels = c(size_classes, "all"))
    
    # get holes in data, to plot a vertical line
    diff_following_years <- lead(as.numeric(levels(as.factor(toplot$fishing_year)))) - as.numeric(levels(as.factor(toplot$fishing_year)))
    line_pos = which(diff_following_years != 1)+0.5
    
    fig1 <- ggplot(toplot, aes(x = as.factor(fishing_year), y = m, color = group, group = group))+
      geom_point(size = 0.75, position = position_dodge(0.25))+
      scale_color_brewer("Size class", palette = "Set1")+
      geom_vline(aes(xintercept = line_pos))+
      geom_errorbar(aes(ymin = m - se, ymax = m + se), width = 0.2, size = 0.75, position = position_dodge(0.25))+
      theme(axis.text.x = element_text(angle = 90),
            panel.border = element_rect(color = "black", fill = NA))+
      ylab("Mean Kn")+
      xlab("Fishing year")
    
    ggsave(fig1Name, fig1)
    
    # rm(data_by_class) ; invisible(gc())
    
    #' ***************
    #' Generate Fig2:
    #' ***************
    #' Comparison of the condition factor (Kn) between
    #' FSC and FAD associated schools
    
    msg <- "\n\n~~~~ Generating Figure 2 ~~~~\n" ; cat(msg) ; lines.to.cat <- c(lines.to.cat, msg)
    
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
  
  #' ***************
  #' Perform GAM:
  #' ***************
  
  #' If size_class_for_model != 'all', subsample the data to keep only one size class
  if (size_class_for_model != 'all'){
    
    if (size_class_for_model %in% size_classes){
      
      data <- data_byclass[[grep(size_class_for_model, size_classes)+1]]
      
    } else {
      
      if (grepl("-", size_class_for_model)){
        l1 <- as.numeric(sub("-.*", "", size_class_for_model))
        l2 <- as.numeric(sub(".*-", "", size_class_for_model))
      } else if (grepl(">", size_class_for_model)){
        l1 <- as.numeric(sub(">", "", size_class_for_model))
        l2 <- Inf
      } else if (grepl("<", size_class_for_model)){
        l1 <- 0
        l2 <- as.numeric(sub("<", "", size_class_for_model))
      }
      
      data %>% filter(fork_length > l1 & fork_length <= l2 ) -> data
      
    }
  }
  
  #' Data saved at this point to study spatial autocorrelation
  #' with seed = 123456
  #' saveRDS(data, file = "/home/adupaix/Documents/These/Axe_3/Historical_YFT_condition/3.Outputs/data_for_spatial_autocorr.rds")
  #'      (dataset of 40-60cm, used for the model on 21.10.14)
  #'      
  saveRDS(data, file = "/home/adupaix/Documents/These/Axe_3/Historical_YFT_condition/3.Outputs/data_for_model.rds")
  #'      (dataset of all sizes, used for the model on 21.10.14) seed = 123456
  
  #' @1. Scale quantitative variables + scale and center geographical variables
  data %>% mutate(scaled_FL = scale(fork_length, scale = T, center = F),
                  scaled_lon = scale(lon, scale = T, center = T),
                  scaled_lat = scale(lat, scale = T, center = T),
                  fishing_month = as.factor(fishing_month),
                  fishing_quarter = as.factor(fishing_quarter),
                  fishing_year = as.factor(fishing_year),
                  size_class = case_when(fork_length <= 60 ~ "<60",
                                         fork_length > 60 & fork_length <= 80 ~ "60-80",
                                         fork_length > 80 & fork_length <= 100 ~ "80-100",
                                         fork_length > 100 & fork_length <= 120 ~ "100-120",
                                         fork_length > 120 ~ ">120")) -> data
  
  if (fad_fsc == T){
    data %>% dplyr::filter(fishing_mode %in% c("DFAD","FSC")) %>%
      mutate(fishing_mode = as.factor(fishing_mode)) -> data
  }
  
  #' @2. Testing for correlation among variables
  # Considered variables : lon, lat, fishing_quarter, size, year
  
  car::vif(lm(Kn ~ scaled_lon + scaled_lat + scaled_FL + fishing_year + fishing_quarter, data = data))
  
  #'                         GVIF Df GVIF^(1/(2*Df))
  #'scaled_lon          1.331296  1        1.153818
  #'scaled_lat          1.582895  1        1.258132
  #'scaled_FL           1.085781  1        1.042008
  #'scaled_fishing_year 1.041272  1        1.020427
  #'fishing_quarter     1.764107  3        1.099227
  #'
  #' @Results: no VIF above 5, no strong correlation between explanatory variables, we keep all of them
  
  #' @3. Testing for spatial autocorrelation
  #'  @!!! the script is very long if size_class_for_model == "all" (26,000 entries in the df)
  
  msg <- "\n\n~~~~ Performing Moran's I test on Kn ~~~~" ; cat(msg) ; lines.to.cat <- c(lines.to.cat, msg)
  
  if (size_class_for_model == "all"){
    to_sample <- sample.int(dim(data)[1], size = 2000)
    subdata <- data[to_sample,]
  } else {
    subdata = data
  }
  
  c <- seq(1,10,1)
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
    listw[[k]] <- nb2listw(nb_list, zero.policy = T)
    
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
  
  #' @4. Choose link function for the model
  
  ggplot(data = data, aes(x=Kn)) + geom_histogram()
  #' @Results: we consider that Kn follows a normal distribution and use the identity function as link function
  #'  hence glm(..., family = gaussian) 
  #'  
  
  #' @5. Visualisation données
  #' 
  #' @months (pas utilise dans le model)
  p <- ggplot() +
    geom_boxplot(data=data, aes(x=fishing_month, y = Kn))+
    ylab("Kn")+xlab("Fishing month")
  
  l=1
  ggsave(plotsNames[l], p) ; l= l+1
  
  ggplot() +
    geom_point(data=data %>% plyr::ddply(.variables = "fishing_month", summarise, Kn = mean(Kn)),
               aes(x=fishing_month, y = Kn))+
    ylab("Mean Kn")+xlab("Fishing month")
  
  
  #' @quarter (utilise dans le model)
  p <- ggplot() +
    geom_boxplot(data=data, aes(x=fishing_quarter, y = Kn))+
    ylab("Kn")+xlab("Fishing quarter")
  
  ggsave(plotsNames[l], p) ; l= l+1
  
  ggplot() +
    geom_point(data=data %>% plyr::ddply(.variables = "fishing_quarter", summarise, Kn = mean(Kn)),
               aes(x=fishing_quarter, y = Kn))+
    ylab("Mean Kn")+xlab("Fishing quarter")
  
  #' @long_lat (utilise dans le model)
  p <- ggplot() +
    geom_point(data= data, aes(x=lon, y=lat, color = Kn)) +
    xlab("Longitude")+ylab("Latitude")
  
  ggsave(plotsNames[l], p) ; l= l+1
  
  p <- ggplot() + geom_point(data= data, aes(x=lon, y = Kn))+xlab("Longitude")
  ggsave(plotsNames[l], p) ; l= l+1
  p <- ggplot() + geom_point(data= data, aes(x=lat, y = Kn))+xlab("Latitude")
  ggsave(plotsNames[l], p) ; l= l+1
  
  #' @6. Build the model
  if (size_class_for_model == "all"){
    if (fad_fsc == F){
    gam_Kn <- mgcv::gam(Kn ~ fishing_quarter + fishing_year + size_class + te(scaled_lon, scaled_lat),
                        data = data)
    glm_Kn <- glm(Kn ~ fishing_quarter + fishing_year + size_class + scaled_lon + scaled_lat + scaled_lon*scaled_lat,
                  data = data)
    } else {
      gam_Kn <- mgcv::gam(Kn ~ fishing_quarter + fishing_year + size_class + fishing_mode + te(scaled_lon, scaled_lat),
                          data = data)
      glm_Kn <- glm(Kn ~ fishing_quarter + fishing_year + size_class + fishing_mode + scaled_lon + scaled_lat + scaled_lon*scaled_lat,
                    data = data)
    }
  } else {
    gam_Kn <- mgcv::gam(Kn ~ fishing_quarter + fishing_year + te(scaled_lon, scaled_lat),
                        data = data)
    glm_Kn <- glm(Kn ~ fishing_quarter + fishing_year + scaled_lon + scaled_lat + scaled_lon*scaled_lat,
                  data = data)
  }
  
  png(diagnoPlotName, width = 1080, height = 1080, units = "px")
  gratia::appraise(gam_Kn)
  dev.off()
  
  saveRDS(gam_Kn, gamSummary)
  saveRDS(glm_Kn, glmSummary)
  
  # for (k in 2:(dim(glm_summaries)[2])){
  #   glm_summaries[i,k] <- summ[sub("__.*", "", names(glm_summaries)[k]),sub(".*__", "", names(glm_summaries)[k])]
  # }
  
  #' @7. Test spatial autocorrelation on the residuals
  cat(lines.to.cat)
  msg <- "\n\n~~~~ Performing Moran's I tests on model residuals ~~~~" ; cat(msg) ; lines.to.cat <- c(lines.to.cat, msg)
  
  c <- seq(1,10,1)
  mc_glm <- list() ; test_glm <- list() ; mc_gam <- list() ; test_gam <- list()
  
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
    mc_glm[[k]] <- moran.mc(resid(glm_Kn)[to_sample], listw = listw[[k]], nsim = 99, zero.policy = T)
    test_glm[[k]] <- moran.test(resid(glm_Kn)[to_sample], listw = listw[[k]], zero.policy = T) 
    mc_gam[[k]] <- moran.mc(resid(gam_Kn)[to_sample], listw = listw[[k]], nsim = 99, zero.policy = T)
    test_gam[[k]] <- moran.test(resid(gam_Kn)[to_sample], listw = listw[[k]], zero.policy = T) 
    
  }
  
  morans.mc.res <- unlist(lapply(mc_glm, function(x) x$statistic))
  morans.test.res <- unlist(lapply(test_glm, function(x) x$estimate[1]))
  morans.mc.res.gam <- unlist(lapply(mc_gam, function(x) x$statistic))
  morans.test.res.gam <- unlist(lapply(test_gam, function(x) x$estimate[1]))
  
  toplot2 <- data.frame(d = rep(c, 4),
                       I = c(morans.mc.res, morans.test.res, morans.mc.res.gam, morans.test.res.gam),
                       moran_type = c(rep("mc_res_glm",length(c)),
                                      rep("test_res_glm",length(c)),
                                      rep("mc_res_gam",length(c)),
                                      rep("test_res_gam",length(c))))
  
  bind_rows(toplot, toplot2) -> toplot
  
  
  p_moran_residuals <- ggplot(data = toplot, aes(x=d*100, y=I, group = moran_type, color = moran_type)) + 
    geom_point()+ geom_line()+
    scale_color_discrete("Moran's\nfunction")+
    xlab("Distance (km)")+
    ylab("Moran's I")
  
  toplot <- toplot[grep("test", toplot$moran_type),]
  toplot$moran_type <- gsub(pattern = "test$", "Kn", toplot$moran_type)
  toplot$moran_type <- gsub(pattern = "test_res_glm", "GLM residuals", toplot$moran_type)
  toplot$moran_type <- gsub(pattern = "test_res_gam", "GAM residuals", toplot$moran_type)
  
  p_moran_residuals2 <- ggplot(data = toplot, aes(x=d*100, y=I, group = moran_type, color = moran_type)) + 
    geom_point()+ geom_line()+
    scale_color_discrete("Data used to\n calculate Moran's I")+
    xlab("Distance (km)")+
    ylab("Moran's I")
  
  ggsave(moranPlotName, p_moran_residuals2)
  
  rm(mc, test, nb_list, listw) ; invisible(gc())
  
  
}
