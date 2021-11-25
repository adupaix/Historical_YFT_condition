#'#*******************************************************************************************************************
#'@author : Amael DUPAIX
#'@update : 2021-11-17
#'@email : amael.dupaix@ens-lyon.fr
#'#*******************************************************************************************************************
#'@description :  Calculate the fishing date from the fishing dates interval or from the landing date
#'#*******************************************************************************************************************
#'@revisions
#'#*******************************************************************************************************************

if (VERBOSE){
  msg <- "\n~~~~ Calculating fishing date from interval and/or sampling dates ~~~~\n" ; cat(msg) ; lines.to.cat <- c(lines.to.cat, msg)
}

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
#' @remark: this data.frame (data_no_dates) should be empty, because when the fishing date interval
#'          is absent, it means the fishing trip could not be deduced, hence the geometry is also
#'          empty (neither a POINT nor a MULTIPOINT)
#'          But just in case, I leave it...


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

if (VERBOSE){
  msg <- paste0("    - Number of removed entries (date errors): ", dim(error_interval)[1]+dim(error_exact)[1], "\n") ; cat(msg) ; lines.to.cat <- c(lines.to.cat, msg)
}

# For data with a date interval, take the middle of the interval
data_interval_dates %>%
  mutate(date1 = fishing_date_min,
         date2 = pmin(fishing_date_max, fish_sampling_date)) %>%
  mutate(fishing_date = date1 + (date2 - date1)/2) %>%
  dplyr::select(-date1, -date2) -> data_interval_dates


# bind data with date interval and with exact date
bind_rows(data_interval_dates, data_exact_dates) %>%
  dplyr::arrange("fish_identifier") %>%
  mutate(t_fishing_sampling = as.numeric(difftime(fish_sampling_date, fishing_date, units = "days")))-> data_dates

if (deduce_date & cluster == F){
  # fit a log normal law to the time between fishing and sampling
  fit_lnorm = fitdistrplus::fitdist(data_dates$t_fishing_sampling, "lnorm")
  
  # plot informations on the fit and save the plots
  if (generate_plots){
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
  }
  
  # for data with missing fishing dates, deduce a fishing date using the above fit
  data_no_dates %>%
    mutate(t_fishing_sampling = rlnorm(dim(data_no_dates)[1],
                                       meanlog = fit_lnorm$estimate[1],
                                       sdlog = fit_lnorm$estimate[2])) %>%
    mutate(fishing_date = fish_sampling_date - as.difftime(t_fishing_sampling, units = "days")) -> data_no_dates
  
  
  #bind data with dates and data with randomly assigned dates
  bind_rows(data_dates, data_no_dates) -> data
  
  if (VERBOSE){
    msg <- paste0("    - Number of deduced entries (missing date): ", dim(data_no_dates)[1], "\n") ; cat(msg) ; lines.to.cat <- c(lines.to.cat, msg)
  }
  
} else {
  
  if (VERBOSE){
    msg <- paste0("    - Number of removed entries (missing date): ", dim(data_no_dates)[1], "\n") ; cat(msg) ; lines.to.cat <- c(lines.to.cat, msg)
  }
  
  data <- data_dates
}

data %>%
  mutate(fishing_year = year(fishing_date),
         fishing_month = as.factor(month(fishing_date)),
         fishing_quarter = as.factor(quarter(fishing_date, fiscal_start = 12))) -> data

N3 <- dim(data)[1]

#' Last @filter
#' Keep only data before 2020
#'     sampling on YFT in 2020 was mainly performed in January on tuna fished in Nov or Dec 2019
data %>%
  filter(fishing_year < 2020) %>%
  mutate(fishing_year = as.factor(fishing_year)) -> data

sink(summaryName, append = T)
cat("\n\n\nDates filter")
cat("\n------------\n")
cat("\n  Number of entries with exact fishing date available:", N1)
cat("\n     Of which",dim(error_exact)[1],"contain an error (removed)")
cat("\n  Number of entries with fishing date interval:", N2)
cat("\n     Of which",dim(error_interval)[1],"contain an error (removed)")
if (deduce_date){
  cat("\n  Number of entries with no fishing date available (deduced using the lnorm fit):", dim(data_no_dates)[1])
} else {
  cat("\n  Number of entries with no fishing date available (removed):", dim(data_no_dates)[1])
}
cat("\n  Number of entries after 2019 (removed):", N3 - dim(data)[1])
cat("\n\n    - Number of entries after date filter:",dim(data)[1])
sink()

rm(data_dates, data_no_dates, data_exact_dates, data_interval_dates, N1, N2, N3) ; invisible(gc())

