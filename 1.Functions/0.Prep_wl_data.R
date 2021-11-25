#'#*******************************************************************************************************************
#'@author : Amael DUPAIX
#'@update : 2021-11-17
#'@email : amael.dupaix@ens-lyon.fr
#'#*******************************************************************************************************************
#'@description :  Function to prepare the dataset (apply first filters, etc)
#'#*******************************************************************************************************************
#'@revisions
#'#*******************************************************************************************************************

prep_wl_data <- function(DATA_PATH,
                         data,
                         calcfdl = TRUE,
                         calcfdl.sd_max = 5, calcfdl.spl_size = 100,
                         read = TRUE,
                         getgeom = TRUE, # rest from before, in this study getgeom is always set to T
                         ncores = 1,
                         size_class_levels,
                         summaryName,
                         verbose = F){
  
  if (verbose){
    msg <- "\n\n~~~~ Preparing weight-length data ~~~~\n\n" ; cat(msg) ; lines.to.cat <<- c(lines.to.cat, msg)
  }
  
  preped_file <- file.path(DATA_PATH,"preped_wl_data")
  
  if (calcfdl){
    preped_file <- paste0(preped_file, "_calcfdl")
  }
  if (getgeom){
    preped_file <- paste0(preped_file, "_geom")
  }
  preped_file <- paste0(preped_file, "_sc", paste(size_class_levels, collapse = "."))
  
  preped_file <- paste0(preped_file, ".rds")
  
  
  if (file.exists(preped_file) & read){
    if (verbose){
      msg <- paste("  ~~~ Reading prepared data file from ",preped_file,"\n") ; cat(msg) ; lines.to.cat <<- c(lines.to.cat, msg)
    }
    data <- readRDS(preped_file)
    
    sink(summaryName, append = T)
    cat("\n\n- Number of entries in read dataset:", dim(data)[1])
    sink()
    
  } else {
    
    N1 <- dim(data)[1]
    sink(summaryName, append = T)
    cat("- Total number of entries:", N1)
    sink()
    
    if (verbose){
      msg <- "  ~~~ Applying first filters\n" ; cat(msg) ; lines.to.cat <<- c(lines.to.cat, msg)
    }
    
    ## select the variables of interest
    data %>% dplyr::select(fish_identifier, ocean_code, quadrant, gear_code,
                           fishing_mode, landing_site, landing_date,
                           fishing_date, fishing_date_min, fishing_date_max,
                           fishing_hour, sea_surface_temp, geometry,
                           fish_identifier_origin, fish_sampling_date, species_code_fao,
                           total_length, fork_length, low_jaw_fork_length, disc_width,
                           first_dorsal_length, body_height, body_width,
                           straight_carapace_length, dorsal_mantle_length,pectoral_length,
                           curved_fork_length, middle_thorax_girth,
                           first_thorax_girth, second_thorax_girth, head_length,
                           mouth_length, opening_mouth_height,
                           opening_mouth_width, whole_fish_weight, gutted_fish_weight, sex, macro_maturity_stage,
                           gonad_total_weight, gonad_1_weight, gonad_2_weight, liver_weight,
                           empty_stomach_weight, stomach_length, stomach_fullness_stage, stomach_prey_group,
                           rest_viscera_weight, full_stomach_weight, first_tag_number, second_tag_number
                           ) %>%
      # add a column year
      # change the date column in class Date
      mutate(year = year(as.Date(fish_sampling_date, format = "%Y-%m-%d")),
             fish_sampling_date = as.Date(fish_sampling_date, format = "%Y-%m-%d"),
             fishing_date = as.Date(fishing_date, format = "%Y-%m-%d"),
             fishing_date_min = as.Date(fishing_date_min, format = "%Y-%m-%d"),
             fishing_date_max = as.Date(fishing_date_max, format = "%Y-%m-%d"),
             landing_date = as.Date(landing_date, format = "%Y-%m-%d")) %>%
      # change the lengths to numeric
      mutate(fork_length = as.numeric(gsub(",",".", fork_length)),
             total_length = as.numeric(gsub(",",".", total_length)),
             first_dorsal_length = as.numeric(gsub(",",".", first_dorsal_length))) %>%
      # change the weights to numeric
      mutate(whole_fish_weight = as.numeric(gsub(",",".", whole_fish_weight)),
             gonad_total_weight = as.numeric(gsub(",",".", gonad_total_weight)),
             liver_weight = as.numeric(gsub(",",".", liver_weight))) %>%
      # keep only data for which a length was measured
      dplyr::filter(!is.na(fork_length) | !is.na(total_length) | !is.na(first_dorsal_length)) %>%
      # add a column to follow where the fork length comes from
      mutate(fl_origin = as.factor(ifelse(!is.na(fork_length), "measured", "deduced"))) -> data
    
    # add a column containing the size_class
    data$size_class <- NA
    for (k in 1:length(size_class_levels)){
      if (k == 1){
        data$size_class[which(data$fork_length <= as.numeric(sub("<", "", size_class_levels[k])))] <- size_class_levels[k]
      } else if (k == length(size_class_levels)){
        data$size_class[which(data$fork_length > as.numeric(sub(">", "", size_class_levels[k])))] <- size_class_levels[k]
      } else {
        data$size_class[which(data$fork_length > as.numeric(sub("-.*", "", size_class_levels[k])) &
                                data$fork_length <= as.numeric(sub(".*-", "", size_class_levels[k])))] <- size_class_levels[k]
      }
    }
    data$size_class <- as.factor(data$size_class)
    
    N2 <- dim(data)[1]
    if (verbose){
      msg <- "    - Deleted data with no length measurement\n" ; cat(msg) ; lines.to.cat <<- c(lines.to.cat, msg)
    }
    
    
    data %>% dplyr::filter(!is.na(whole_fish_weight)) -> data
    N3 <- dim(data)[1]
    if (verbose){
      msg <- "    - Deleted data with no weight measurement\n" ; cat(msg) ; lines.to.cat <<- c(lines.to.cat, msg)
    }
    
    data %>% dplyr::filter(ocean_code == "IO") -> data
    N4 <- dim(data)[1]
    if (verbose){
      msg <- "    - Deleted data of fish not from the IO\n\n" ; cat(msg) ; lines.to.cat <<- c(lines.to.cat, msg)
    
      msg <- paste("~~~ Sampling missing FL values among other individuals with the same FDL: ",calcfdl,"\n") ; cat(msg) ; lines.to.cat <<- c(lines.to.cat, msg)
    }
    
    if (calcfdl){
      data <- calc_FL_from_FDL(data, sd_max = calcfdl.sd_max, spl_size_min = calcfdl.spl_size, verbose = verbose)
    } else{
      if (verbose){
        msg <- paste("    - Number of deleted entries (no fork length available):",sum(is.na(data$fork_length)), "\n\n") ; cat(msg) ; lines.to.cat <<- c(lines.to.cat, msg)
      }
      data %>% dplyr::filter(!is.na(fork_length)) -> data
    }
    
    N5 <- dim(data)[1]
    
    if (verbose){
      msg <- "  ~~~ Getting geometry information from str: " ; cat(msg) ; lines.to.cat <<- c(lines.to.cat, msg)
    }
    
    if (getgeom){
      
      data %>% mutate(geom_type = gsub("[^[:alpha:]]", "", as.character(geometry)),
                      geom_coord = gsub("[[:alpha:]]", "", as.character(geometry))) -> data
      
      sample_size = 1000
      niter <- floor(dim(data)[1]/sample_size) + 1
      
      data_list <- list()
      
      #' run on subsamples of the data, because, for a reason I ignore, the foreach
      #' does not succeed in binding rows when there are too much iterations (in the foreach loop)
      #'      @maybe? because of the multicombine option (when we use .combine = rbind, .multicombine is set to T by default)
      #'      hence it tries to combine a lot of rows together and hence needs to maintain the connection open with
      #'      several "cores" (the number of R processes running at the same time is > 4...)
      for (k in 1:niter){
        
        if (verbose){
          cat(lines.to.cat)
          cat("    Sub-sample ",k,"/",niter, "\n")
        }
        
        if(k == niter){
          indexes <- ((k-1)*sample_size+1):(dim(data)[1])
        } else {
          indexes <- ((k-1)*sample_size+1):(k*sample_size)
        }
        
        sub_data <- data[indexes,]
        
        cl <- makeCluster(ncores)
        doSNOW::registerDoSNOW(cl)
        
        if (VERBOSE){
          pb <- txtProgressBar(min = 1, max = dim(sub_data)[1], style = 3)
          progress <- function(n) setTxtProgressBar(pb, n)
          opts <- list(progress = progress)
        } else {
          opts <- list()
        }
        
        data_list[[k]] <- foreach(i = 1:dim(sub_data)[1],
                                  .combine = dplyr::bind_rows,
                                  .options.snow = opts,
                                  .multicombine = F,
                                  .packages = c("sf", "dplyr")) %dopar% {
                          
                          # extract one line of the data
                          data.i <- sub_data[i,]
                          
                          # from the coordinates which have this format: "((x1 y1),(x2 y2))"
                          # get a vector which is as follow: c(x1,y1,x2,y2)
                          str_coords <- as.numeric(gsub("[^0-9.-]", "", unlist(strsplit(gsub(",", " ", trimws(data.i$geom_coord)), " "))))
                          
                          # create a data.frame with the extracted coordinates:
                          #  x  y
                          # x1 y1
                          # x2 y2
                          coords <- data.frame(cbind(x = str_coords[seq(length(str_coords)) %% 2 == 1], y = str_coords[seq(length(str_coords)) %% 2 == 0]))
                          
                          n_pts <- 0
                          
                          # create a geometry from the coordinates and the geom_type
                          if (data.i$geom_type == "POINT" | data.i$geom_type == "na"){
                            geom <- st_sfc(st_point(as.matrix(coords)), crs = 4326)
                            n_pts <- 1
                          } else if (data.i$geom_type == "MULTIPOINT"){
                            geom <- st_sfc(st_multipoint(as.matrix(coords)), crs = 4326)
                            n_pts <- dim(coords)[1]
                          }
                          
                          # create an sf, data.frame from the coordinates and the rest of the values
                          st_sf(data.i, geometry = geom) %>% dplyr::select(-geom_coord) %>%
                            mutate(n_points = n_pts)
                          
                          
                          
                        }
        
        if (VERBOSE){
          close(pb)
        }
        
        foreach::registerDoSEQ()
        
      }
      
      data <- bind_rows(data_list)
      
    }
    
    saveRDS(data, preped_file)
    
    sink(summaryName, append = T)
    cat("\n\nFirst filters")
    cat("\n-------------\n")
    cat("\n  Number of entries with no fork length (FL) and no fish dorsal length (FDL) measured:", N1-N2)
    cat("\n  Number of entries with no fish weight measured:", N2-N3)
    cat("\n  Number of entries of fish not from the IO:", N3-N4)
    cat("\n\n  Fork length deduced from fish dorsal length:", calcfdl)
    if (calcfdl == F){
      cat("\n    Number of entries with no FL (but with FDL) measured (removed):", N4-N5)
    }
    cat("\n\n      - Number of entries after first filters:", dim(data)[1])
    sink()
    
    
  }
  
  return(data)
  
}

#%############

calc_FL_from_FDL <- function(data, sd_max = 5, spl_size_min = 100, verbose = F){
  
  data %>%
    dplyr::filter(is.na(fork_length)) %>%
    dplyr::filter(species_code_fao != "SKJ") -> fdl_no_fl
  
  if (verbose){
    cat(paste("    Number of values to deduce:", dim(fdl_no_fl)[1], "\n"))
  }
  
  # df with both FL and FDL measured
  data %>% dplyr::filter(!is.na(first_dorsal_length) & !is.na(fork_length)) -> both_l
  
  for (i in 1:dim(fdl_no_fl)[1]){
    
    
    fdl.i <- fdl_no_fl[i,"first_dorsal_length"]
    sp.i <- fdl_no_fl[i,"species_code_fao"]
    id.i <- fdl_no_fl[i, "fish_identifier"]
    
    fl_to_draw_from <- (data %>% dplyr::filter(first_dorsal_length == fdl.i &
                                                 species_code_fao == sp.i &
                                                 !is.na(fork_length)
    )
    )$fork_length
    # if there are enough FL measured for the same FDL and if their standard deviation is not too big
    if (length(fl_to_draw_from) > spl_size_min &  sd(fl_to_draw_from)< sd_max){
      new_fl <- sample(fl_to_draw_from, 1)
      data[which(data$fish_identifier==id.i), "fork_length"] <- new_fl
    }
    
  }
  
  data %>% dplyr::filter(fl_origin == "deduced") -> res
  
  if (verbose){
    cat(paste("    Number of deleted entries (no replacement found):",sum(is.na(res$fork_length)), "\n")) ## -> 661 values of FL not replaced (sd_max = 5, spl_size_min = 100)
  }
  
  data %>% dplyr::filter(!is.na(fork_length)) -> data
  
  return(data)
  
}
