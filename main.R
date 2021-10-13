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
OUTPUT_PATH <- file.path(WD, "3.Outputs")

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


#' ***************
#' Get libraries:
#' ***************
source(file.path(FUNC_PATH, "install_libraries.R"))

srcUsedPackages <- c("plyr", "dplyr","tidyr","foreach","doSNOW","stringr","lubridate","sf",
                     "parallel","ggplot2","tibble","cowplot","RColorBrewer", "MASS","truncnorm",
                     "ape", "mgcv", "spdep")

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

if (nb_of_times_to_run != 1){
  seeds <- round(runif(nb_of_times_to_run, 1, 10^6))
} else {
  seeds = SEED
}


glm_summaries <- data.frame(matrix(ncol = 3*9+4, nrow = length(seeds), data = NA))
names(glm_summaries) <- c("seed", "moran.mc.Kn.p.value", "moran.test.Kn.p.value",
                          paste(rep(c("(Intercept)","fishing_quarter2","fishing_quarter3",
                                      "fishing_quarter4", "fishing_year", "scaled_FL",
                                      "scaled_lon","scaled_lat", "scaled_lon:scaled_lat"), each = 3),
                                rep(c("Estimate","Std. Error","Pr(>|t|)"), times = 9),
                                sep = "__"),
                          "percent.moran.res.not.signif")


for (i in 1:length(seeds)){
  
  set.seed(seeds[i])
  
  glm_summaries[i,1] <- seeds[i]
  
  msg <- "\14" ; cat(msg) ; lines.to.cat <- c(msg)
  msg <- paste("-------- ITERATION", i, "/", length(seeds), " - SEED NB =", seeds[i], "--------") ; cat(msg) ; lines.to.cat <- c(lines.to.cat, msg)
  
  # Initialize names
  #Initialize data processing summary name
  summaryName <- file.path(OUTPUT_PATH, paste0("Data_processing_summary_RESET-",RESET,".txt"))
  
  try(dir.create(PLOT_PATH, showWarnings = F))
  
  fitName <- file.path(PLOT_PATH, "fit_lnorm_dates_goodnessoffit.png")
  histFitName <- file.path(PLOT_PATH, "fit_lnorm_dates.png")
  fig1Name <- file.path(PLOT_PATH, "Figure1.png")
  fig2Name <- file.path(PLOT_PATH, "Figure2.png")
  
  vars <- c("month", "quarter", "lonlat", "lon", "lat")
  plotsNames <- file.path(PLOT_PATH, paste0("Kn_f-", vars, ".png"))
  diagnoPlotName <- file.path(PLOT_PATH, paste0("diagn_glm_", seeds[i], ".png"))
  
  glmSummary <- file.path(OUTPUT_PATH, "glm_summary.csv")
  
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
  
  
  
  #' ***************
  #' Generate Fig1:
  #' ***************
  #' Evolution of the condition factor (Kn)
  #' of different size classes according to time
  
  msg <- "\n~~~~ Generating Figure 1 ~~~~\n\n" ; cat(msg) ; lines.to.cat <- c(lines.to.cat, msg)
  
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
  
  figure2 <- function(data, var.to.compare, var.grp, levels.var.grp, var.x,
                      scale.color.title = "School\n type",
                      xlabel = "Year",
                      line = 0
                      ){
    
    data %>% dplyr::filter((!!rlang::sym(var.grp)) %in% levels.var.grp) %>%
      plyr::ddply(.variables = c(var.x,var.grp), summarise, n=n()) %>%
      spread((!!rlang::sym(var.grp)), n) %>%
      filter(!is.na((!!rlang::sym(levels.var.grp[1]))) & !is.na((!!rlang::sym(levels.var.grp[2])))) -> spl_sizes
    
    var.x_of_int <- spl_sizes[,var.x]
    
    data %>% 
      filter((!!rlang::sym(var.x)) %in% var.x_of_int & (!!rlang::sym(var.grp)) %in% levels.var.grp) -> dat_
    
    dat_ %>% ddply(.variables = c(var.x,var.grp), summarise, m = mean(!!rlang::sym(var.to.compare)),
                   sd = sd(!!rlang::sym(var.to.compare)), n = n()) %>%
      mutate(se = sd / sqrt(n)) -> toplot
    
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
    
    df <- data.frame(cbind(var.x_of_int, p))
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
    
    toplot %>% ddply(.variables = var.x, .fun = function(x){
      (x %>% filter(!!rlang::sym(var.grp) == levels.var.grp[1]))$m -
        (x %>% filter(!!rlang::sym(var.grp) == levels.var.grp[2]))$m
    }) -> data_hist
    
    p3.2 <- ggplot(data_hist, aes(x = as.factor(!!rlang::sym(var.x)), y = V1))+
      geom_bar(stat = "identity")+
      theme(axis.text.x = element_text(angle = 90, colour = axis_col, face = axis_face),
            panel.border = element_rect(color = "black", fill = NA))+
      ylab(paste0(var.to.compare," (",levels.var.grp[1],") - ", var.to.compare," (",levels.var.grp[2],")"))+
      xlab(xlabel)
    
    if(line != 0){
      p3.1 <- p3.1 + geom_vline(aes(xintercept = line))
      p3.2 <- p3.2 + geom_vline(aes(xintercept = line))
    }
    
    p3 <- ggdraw()+
      draw_plot(p3.1, 0, 1/3, 1, 2/3)+
      draw_plot(p3.2, 0, 0, 1, 1/3)+
      draw_plot_label(c("A","B"), c(0,0), c(1,1/3))
    
    return(p3)
    
  }
  
  p3 <- figure2(data = data,
                var.to.compare = "Kn",
                var.grp = "fishing_mode",
                levels.var.grp = c("DFAD","FSC"),
                var.x = "fishing_year",
                scale.color.title = "School\n type",
                xlabel = "Year",
                line = 1.5)
  
  ggsave(fig2Name, p3, width = 7, height = 10)
  
  
  #' ***************
  #' Perform GLM:
  #' ***************
  
  #' @1. Scale quantitative variables + scale and center geographical variables
  data %>% mutate(scaled_FL = scale(fork_length, scale = T, center = F),
                  scaled_lon = scale(lon, scale = T, center = T),
                  scaled_lat = scale(lat, scale = T, center = T),
                  fishing_month = as.factor(fishing_month),
                  fishing_quarter = as.factor(fishing_quarter),
                  scaled_fishing_year = scale(fishing_year, scale = T, center = F)) -> data
  
  #' @2. Testing for correlation among variables
  # Considered variables : lon, lat, fishing_quarter, size, year
  
  car::vif(lm(Kn ~ scaled_lon + scaled_lat + scaled_FL + scaled_fishing_year + fishing_quarter, data = data))
  
  #'                         GVIF Df GVIF^(1/(2*Df))
  #'scaled_lon          1.331296  1        1.153818
  #'scaled_lat          1.582895  1        1.258132
  #'scaled_FL           1.085781  1        1.042008
  #'scaled_fishing_year 1.041272  1        1.020427
  #'fishing_quarter     1.764107  3        1.099227
  #'
  #' @Results: no VIF above 5, no strong correlation between explanatory variables, we keep all of them
  
  #' @3. Testing for spatial autocorrelation
  #' Calculate distance matrix between points
  #' 
  #' Because we have too much entries in the data.frame (~26,000, hence a distance matrix of 26,000^2)
  #' we test for spatial autocorrelation on a sub-sample, and repeat the operation 1000 times
  
  msg <- "\n\n~~~~ Performing Moran's I test on Kn ~~~~" ; cat(msg) ; lines.to.cat <- c(lines.to.cat, msg)
  
  c <- seq(50, 450, 100)
  mc <- list() ; test <- list() ; nb_list <- list() ; listw <- list()
  
  for (k in 1:length(c)){
    cat(lines.to.cat)
    cat(paste("\n  - Consider Neighboors if d <=", c[k], "km"))
    cat("\n      - Getting nearest neighboors")
    nb_list[[k]] <- spdep::dnearneigh(x = st_coordinates(data),
                                 d1 = 0,
                                 d2 = c[k],
                                 longlat = T)
    
    msg <- "\n      - Getting weight list from nearest neighboors" ; cat(msg)
    listw[[k]] <- nb2listw(nb_list[[k]], zero.policy = T)
    
    msg <- "\n      - Performing Moran's I tests" ; cat(msg)
    mc[[k]] <- moran.mc(data$Kn, listw = listw[[k]], nsim = 99, zero.policy = T)
    test[[k]] <- moran.test(data$Kn, listw = listw[[k]], zero.policy = T)
  }
  
  morans.mc <- unlist(lapply(mc, function(x) x$statistic))
  morans.test <- unlist(lapply(test, function(x) x$statistic))
  
  toplot <- data.frame(d = rep(c, 2),
                       I = c(morans.mc, morans.test),
                       moran_type = c(rep("mc",length(c)),
                                      rep("test",length(c))))
  
  p_moran <- ggplot(data = toplot, aes(x=d, y=I, group = moran_type)) + 
    geom_point()+ geom_line()+
    scale_color_discrete("Moran's\nfunction")+
    xlab("Distance (km)")+
    ylab("Moran's I")
  
  # glm_summaries[i,2] <- mc$p.value
  # glm_summaries[i,3] <- test$p.value
  
  rm(mc, test, nb_list, listw) ; invisible(gc())
  
  # sample_size = 5000
  # niter = 100
  # moran_p_values <- vector(length = niter)
  # 
  # pb <- txtProgressBar(min = 0, max = niter, style = 3)
  # 
  # for (k in 1:niter){
  #   
  #   to_sample <- sample.int(dim(data)[1], sample_size)
  #   
  #   lon_line <- matrix(data = data$lon[to_sample],
  #                      nrow = sample_size,
  #                      ncol = sample_size,
  #                      byrow = T)
  #   
  #   lon_col <- matrix(data = data$lon[to_sample],
  #                     nrow = sample_size,
  #                     ncol = sample_size,
  #                     byrow = F)
  #   
  #   lat_line <- matrix(data = data$lat[to_sample],
  #                      nrow = sample_size,
  #                      ncol = sample_size,
  #                      byrow = T)
  #   
  #   lat_col <- matrix(data = data$lat[to_sample],
  #                     nrow = sample_size,
  #                     ncol = sample_size,
  #                     byrow = F)
  #   
  #   dist_mat <- sqrt((lon_line - lon_col) ^ 2 + (lat_line -lat_col) ^ 2)
  #   
  #   dist_mat_inv <- 1/dist_mat
  #   dist_mat_inv[which(is.infinite(dist_mat_inv))] <- 0
  #   
  #   moran_p_values[k] <- ape::Moran.I(data$Kn[to_sample], dist_mat_inv)$p.value
  #   
  #   setTxtProgressBar(pb, k)
  # }
  # 
  # close(pb)
  # 
  # glm_summaries[i,2] <- length(which(moran_p_values > .05)) / length(moran_p_values)
  # 
  # rm(dist_mat, dist_mat_inv) ; invisible(gc())
  
  # sink(glmSummary)
  # 
  # cat("\n\n ~~~~ P-values of Moran's I tests (spatial autocorrelation) ~~~~\n")
  # 
  # print(quantile(moran_p_values, seq(0,1,.1)))
  # 
  # sink()
  
  #' @Results:
  #'          0%           5%          10%          15%          20%          25%          30%          35%          40%          45%          50%          55% 
  #'0.000000e+00 1.110223e-15 1.002753e-13 1.544831e-12 1.630562e-11 1.470827e-10 4.496569e-10 1.668565e-09 6.738344e-09 2.067790e-08 4.658868e-08 1.279228e-07 
  #'         60%          65%          70%          75%          80%          85%          90%          95%         100% 
  #'3.630626e-07 1.111732e-06 3.171815e-06 9.950312e-06 3.112144e-05 9.731789e-05 4.439412e-04 2.603984e-03 6.791270e-01
  #'
  #'@Results: Moran's I test significant for more than 95% of the 1000 sub-samples
  #'          Hence we can reject the null hypothesis that there is zero spatial autocorrelation in the Kn
  #'       
  #' Hence, we consider latitude, longitude and latitude x longitude interaction as explanatory variables in the model
  
  
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
  glm_Kn <- glm(Kn ~ fishing_quarter + fishing_year + scaled_FL + scaled_lon + scaled_lat + scaled_lon*scaled_lat,
                data = data, family = gaussian)
  
  gam_Kn <- mgcv::gam(Kn ~ fishing_quarter + fishing_year + scaled_FL + te(scaled_lon + scaled_lat))
  
  png(diagnoPlotName, width = 1080, height = 1080, units = "px")
  par(mfrow = c(2, 2), oma = c(0, 0, 2, 0))
  plot(glm_Kn)
  dev.off()
  
  # sink(glmSummary, append = T)
  # cat("\n\n ~~~~ Step AIC process ~~~~\n")
  model <- stepAIC(glm_Kn, trace = T, direction = "both")
  summ <- summary(model)$coefficients
  
  for (k in 3:(dim(glm_summaries)[2]-1)){
    glm_summaries[i,k] <- summ[sub("__.*", "", names(glm_summaries)[k]),sub(".*__", "", names(glm_summaries)[k])]
  }
  
  # Start:  AIC=-56157.99
  # Kn ~ fishing_quarter + fishing_year + scaled_FL + scaled_lon + 
  #   scaled_lat + scaled_lon * scaled_lat
  # 
  # Df Deviance    AIC
  # <none>                       178.99 -56158
  # - scaled_FL              1   179.21 -56129
  # - scaled_lon:scaled_lat  1   179.28 -56118
  # - fishing_quarter        3   181.72 -55768
  # - fishing_year           1   193.52 -54119
  # 
  # Call:  glm(formula = Kn ~ fishing_quarter + fishing_year + scaled_FL + 
  #              scaled_lon + scaled_lat + scaled_lon * scaled_lat, family = gaussian, 
  #            data = data)
  # 
  # Coefficients:
  #   (Intercept)       fishing_quarter2       fishing_quarter3       fishing_quarter4           fishing_year              scaled_FL             scaled_lon  
  # -5.647592              -0.027108              -0.022727              -0.020973               0.003314               0.011633              -0.004301  
  # scaled_lat  scaled_lon:scaled_lat  
  # 0.003067              -0.003294  
  # 
  # Degrees of Freedom: 26165 Total (i.e. Null);  26157 Residual
  # Null Deviance:	    198.9 
  # Residual Deviance: 179 	AIC: -56160
  
  #' Based on the AIC algorithm, the complete model is the best model...Hence we shall keep all the variables
  
  # cat("\n\n ~~~~ GLM summary ~~~~\n")
  # print(summary(model))
  
  # Call:
  #   glm(formula = Kn ~ fishing_quarter + fishing_year + scaled_FL + 
  #         scaled_lon + scaled_lat + scaled_lon * scaled_lat, family = gaussian, 
  #       data = data)
  # 
  # Deviance Residuals: 
  #   Min        1Q    Median        3Q       Max  
  # -0.45398  -0.05460  -0.00677   0.04762   0.68674  
  # 
  # Coefficients:
  #   Estimate Std. Error t value Pr(>|t|)    
  # (Intercept)           -5.648e+00  1.448e-01 -38.989  < 2e-16 ***
  #   fishing_quarter2      -2.711e-02  1.520e-03 -17.830  < 2e-16 ***
  #   fishing_quarter3      -2.273e-02  1.495e-03 -15.198  < 2e-16 ***
  #   fishing_quarter4      -2.097e-02  1.713e-03 -12.246  < 2e-16 ***
  #   fishing_year           3.314e-03  7.195e-05  46.067  < 2e-16 ***
  #   scaled_FL              1.163e-02  2.093e-03   5.559 2.74e-08 ***
  #   scaled_lon            -4.301e-03  5.948e-04  -7.231 4.94e-13 ***
  #   scaled_lat             3.067e-03  6.414e-04   4.782 1.74e-06 ***
  #   scaled_lon:scaled_lat -3.294e-03  5.099e-04  -6.461 1.06e-10 ***
  #   ---
  #   Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
  # 
  # (Dispersion parameter for gaussian family taken to be 0.00684308)
  # 
  # Null deviance: 198.87  on 26165  degrees of freedom
  # Residual deviance: 178.99  on 26157  degrees of freedom
  # AIC: -56158
  # 
  # Number of Fisher Scoring iterations: 2
  
  # sink()
  
  
  #' @7. Test spatial autocorrelation on the residuals
  cat(lines.to.cat)
  msg <- "\n\n~~~~ Performing Moran's I tests on model residuals ~~~~" ; cat(msg) ; lines.to.cat <- c(lines.to.cat, msg)
  
  c <- seq(50, 450, 100)
  mc <- list() ; test <- list()
  
  for (k in 1:length(c)){
    cat(lines.to.cat)
    cat(paste("\n  - Consider Neighboors if d <=", c[k], "km"))

    msg <- "\n      - Performing Moran's I tests" ; cat(msg)
    mc[[k]] <- moran.mc(glm_Kn$residuals, listw = listw[[k]], nsim = 99, zero.policy = T)
    test[[k]] <- moran.test(glm_Kn$residuals, listw = listw[[k]], zero.policy = T)
  }
  
  morans.mc <- unlist(lapply(mc, function(x) x$statistic))
  morans.test <- unlist(lapply(test, function(x) x$statistic))
  
  toplot <- data.frame(d = rep(c, 2),
                       I = c(morans.mc, morans.test),
                       moran_type = c(rep("mc",length(c)),
                                      rep("test",length(c))))
  
  p_moran_residuals <- ggplot(data = toplot, aes(x=d, y=I, group = moran_type)) + 
    geom_point()+ geom_line()+
    scale_color_discrete("Moran's\nfunction")+
    xlab("Distance (km)")+
    ylab("Moran's I")
  
  glm_summaries[i,2] <- mc$p.value
  glm_summaries[i,3] <- test$p.value
  
  rm(mc, test, nb_list, listw) ; invisible(gc())
  
  sample_size = 5000
  niter = 100
  moran_p_values <- vector(length = niter)
  
  pb <- txtProgressBar(min = 0, max = niter, style = 3)
  
  for (k in 1:niter){
    
    to_sample <- sample.int(dim(data)[1], sample_size)
    
    lon_line <- matrix(data = glm_Kn$data$lon[to_sample],
                       nrow = sample_size,
                       ncol = sample_size,
                       byrow = T)
    
    lon_col <- matrix(data = glm_Kn$data$lon[to_sample],
                      nrow = sample_size,
                      ncol = sample_size,
                      byrow = F)
    
    lat_line <- matrix(data = glm_Kn$data$lat[to_sample],
                       nrow = sample_size,
                       ncol = sample_size,
                       byrow = T)
    
    lat_col <- matrix(data = glm_Kn$data$lat[to_sample],
                      nrow = sample_size,
                      ncol = sample_size,
                      byrow = F)
    
    dist_mat <- sqrt((lon_line - lon_col) ^ 2 + (lat_line -lat_col) ^ 2)
    
    dist_mat_inv <- 1/dist_mat
    dist_mat_inv[which(is.infinite(dist_mat_inv))] <- 1
    
    moran_p_values[k] <- ape::Moran.I(glm_Kn$residuals[to_sample], dist_mat_inv)$p.value
    
    setTxtProgressBar(pb, k)
  }
  
  close(pb)
  
  glm_summaries[i,dim(glm_summaries)[2]] <- length(which(moran_p_values > .05)) / length(moran_p_values)
  
  rm(dist_mat, dist_mat_inv) ; invisible(gc())
  
  
}



write.csv(glm_summaries, file = glmSummary)