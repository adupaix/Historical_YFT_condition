#~~~~ Data exploration ~~
library(tidyr)
library(plyr)
library(dplyr)
library(lubridate)
library(ggplot2)

data <- read.csv(file = "ddd_dupaix.csv")

dim(data)

## select the variables of interest
data %>% dplyr::select(fish_identifier, fish_sampling_date, species_code_fao,
                       total_length, fork_length, first_dorsal_length,
                       whole_fish_weight, sex, macro_maturity_stage,
                       gonad_total_weight, liver_weight, ocean_code,
                       gear_code, landing_site, landing_date,
                       fishing_date, fishing_date_min, fishing_date_max,
                       geometry) %>%
  # add a column year
  # change the date column in class Date
  mutate(year = year(as.Date(fish_sampling_date, format = "%m/%d/%Y")),
         fish_sampling_date = as.Date(fish_sampling_date, format = "%m/%d/%Y")) %>%
  # change the lengths to numeric
  mutate(fork_length = as.numeric(as.character(fork_length)),
         total_length = as.numeric(as.character(total_length)),
         first_dorsal_length = as.numeric(as.character(first_dorsal_length))) %>%
  # change the weights to numeric
  mutate(whole_fish_weight = as.numeric(as.character(whole_fish_weight)),
         gonad_total_weight = as.numeric(as.character(gonad_total_weight)),
         liver_weight = as.numeric(as.character(liver_weight))) %>%
  # keep only data for which a length was measured
  dplyr::filter(!is.na(fork_length) | !is.na(total_length) | !is.na(first_dorsal_length)) %>%
  # add a column to follow where the fork length comes from
  mutate(fl_origin = as.factor(ifelse(!is.na(fork_length), "measured", "deduced"))) -> data

dim(data) ## -> 30 lines deleted


data %>% dplyr::filter(!is.na(whole_fish_weight)) -> data

dim(data) ## -> 7408 lines deleted


# Deduction of missing fork lengths (FL) from first dorsal length (FDL)


## data frame with the a first dorsal length but no fork length measured
data %>% dplyr::filter(is.na(fork_length)) -> fdl_no_fl

## number of entries per species with a measured FDL and no FL
ddply(fdl_no_fl, .variables = "species_code_fao", summarise, n())


# df with both FL and FDL measured
data %>% dplyr::filter(!is.na(first_dorsal_length) & !is.na(fork_length)) -> both_l


ggplot(both_l, aes(x = first_dorsal_length, y = fork_length, color = sex))+
  geom_point(alpha = 0.3)+
  facet_wrap(~species_code_fao)+
  geom_smooth(method = "lm", se = FALSE, color = "black")+
  xlab("First Dorsal Length")+ylab("Fork Length")

# R2 of the linear regression
ddply(both_l, .variables = "species_code_fao", function(x) summary(lm(fork_length~first_dorsal_length, x))$r.squared)



# delete SKJ from fdl_no_fl
fdl_no_fl %>% dplyr::filter(species_code_fao != "SKJ") -> fdl_no_fl

sd_max = 5
spl_size_min = 100

sds <- c()
spl_sizes <- c()

for (i in 1:dim(fdl_no_fl)[1]){
  fdl.i <- fdl_no_fl[i,"first_dorsal_length"]
  sp.i <- fdl_no_fl[i,"species_code_fao"]
  id.i <- fdl_no_fl[i, "fish_identifier"]
  
  fl_to_draw_from <- (data %>% dplyr::filter(first_dorsal_length == fdl.i &
                                               species_code_fao == sp.i &
                                               !is.na(fork_length)
  )
  )$fork_length
  
  sds <- c(sds, sd(fl_to_draw_from))
  spl_sizes <- c(spl_sizes, length(fl_to_draw_from))
  # if there are enough FL measured for the same FDL and if their standard deviation is not too big
  if (length(fl_to_draw_from) > spl_size_min & sd(fl_to_draw_from) < sd_max){
    new_fl <- sample(fl_to_draw_from, 1)
    data[which(data$fish_identifier==id.i), "fork_length"] <- new_fl
  }
  
}

data %>% dplyr::filter(fl_origin == "deduced") -> res

summary(res$fork_length) ## -> 529 values of FL not replaced (sd_max = 5, spl_size_min = 100)

data %>% dplyr::filter(!is.na(fork_length)) -> data

dim(data)

