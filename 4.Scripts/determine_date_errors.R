#' Script pour tester les erreurs de dates dans la base de données BIO
#' 

ocean = "OA"

WD <- getwd()

DATA_PATH <- file.path(WD, "0.Data/")

FUNC_PATH <- file.path(WD,"1.Functions")
OUTPUT_PATH <- file.path(WD, "3.Outputs")

#' ***************
#' Get libraries:
#' ***************
source(file.path(FUNC_PATH, "install_libraries.R"))

srcUsedPackages <- c("plyr", "dplyr", "tidyr")

installAndLoad_packages(srcUsedPackages, loadPackages = TRUE)

# Open data

data_files <- list.files(path = DATA_PATH)
env_file <- grep(pattern = paste0(ocean, ".*env.csv"), data_files, value = T)
specimen_file <- grep(pattern = paste0(ocean, ".*specimen.csv"), data_files, value = T)

data1 <- read.table(file = file.path(DATA_PATH, env_file),
                    sep = ";", dec=",",stringsAsFactors = F,
                    quote="\"", row.names = NULL,
                    h=T, colClasses = "character")
data2 <- read.table(file = file.path(DATA_PATH, specimen_file),
                    sep = ";", dec=",",stringsAsFactors = F,
                    quote="\"", row.names = NULL,
                    h=T, colClasses = "character")

data <- base::merge(data1, data2, by = "fish_identifier", all = T)

data %>%
  mutate(fish_sampling_date = as.Date(fish_sampling_date, format = "%Y-%m-%d"),
         fishing_date = as.Date(fishing_date, format = "%Y-%m-%d"),
         fishing_date_min = as.Date(fishing_date_min, format = "%Y-%m-%d"),
         fishing_date_max = as.Date(fishing_date_max, format = "%Y-%m-%d"),
         landing_date = as.Date(landing_date, format = "%Y-%m-%d")) -> data

# keep only the entries for which the exact fishing date is known
data %>% filter(!is.na(fishing_date)) %>%
  # keep in a variable how the fishing date was obtained
  mutate(fishing_date_origin = "exact") -> data_exact_dates

# save another part of the data for which a fishing date interval is known
data %>% filter(is.na(fishing_date) & !is.na(fishing_date_min) & !is.na(fishing_date_max)) %>%
  mutate(fishing_date_origin = "interval") -> data_interval_dates

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
  ) -> error_interval # errors 1 and 2

data_exact_dates %>%
  filter(fishing_date > fish_sampling_date) -> error_exact # errors 3

bind_rows(error_interval, error_exact) %>%
  dplyr::select(fish_identifier, ocean_code, species_code_fao,
                fishing_date_min, fishing_date_max, fishing_date, landing_date, fish_sampling_date) %>%
  write.csv(file.path(OUTPUT_PATH, paste0("error_fishing_dates_",ocean,".csv")))
