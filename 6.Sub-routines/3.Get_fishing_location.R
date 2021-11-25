#'#*******************************************************************************************************************
#'@author : Amael DUPAIX
#'@update : 2021-11-17
#'@email : amael.dupaix@ens-lyon.fr
#'#*******************************************************************************************************************
#'@description :  Subroutine to get the fishing locations associated with each line
#'#*******************************************************************************************************************
#'@revisions
#'#*******************************************************************************************************************

#' The geometry is in several forms in the data
#'     1. 1 line per fish with a MULTIPOINT (the whole trip is contained)
#'     2. Several lines per fish, each line with a POINT (the well allowed to narrow down the fishing sets but not to a unique set)
#'     3. 1 line per fish with a POINT (the exact fishing set could be deduced)
#'     4. 1 line with an empty POINT (no information present on the fishing location)
#'
#' The (1) is contained in data_multi
#' The (2) in data_hidden_multi (changed to MULTIPOINTS and binded to data_multi afterwards)
#' The (3) and (4) are in data_na_point, and the (4) will be discarded afterwards in the script

if (VERBOSE){
  msg <- "\n~~~~ Getting fishing location ~~~~\n" ; cat(msg) ; lines.to.cat <- c(lines.to.cat, msg)
}

data %>% filter(geom_type != "MULTIPOINT") %>%
  mutate(geom_sampled = F) -> data_na_point

data %>% filter(geom_type == "MULTIPOINT") %>%
  mutate(geom_sampled = T) -> data_multi

duplicated_fish <- unique(data_na_point$fish_identifier[duplicated(data_na_point$fish_identifier)])
data_na_point %>%
  filter(fish_identifier %in% duplicated_fish) -> data_hidden_multi

data_hidden_multi <- points.to.multipoints(data_hidden_multi, data_hidden_multi$fish_identifier)

data_hidden_multi %>% mutate(geom_sampled = T) -> data_hidden_multi

data_na_point %>%
  filter(!fish_identifier %in% duplicated_fish) -> data_na_point

bind_rows(data_multi, data_hidden_multi) %>%
  dplyr::arrange("fish_identifier") -> data_multi

if (geometry_method == "sampling"){
  if (VERBOSE){
    msg <- "    - Sampling unique points from MULTIPOINT geometries:" ; cat(msg) ; lines.to.cat <- c(lines.to.cat, msg)
  }
  
  sample_size = 1000
  niter <- floor(dim(data_multi)[1]/sample_size) + 1
  
  data_multi_list <- list()
  
  for (k in 1:niter){
    
    if (VERBOSE){
      cat(lines.to.cat)
      cat("    Sub-sample ",k,"/",niter, "\n")
    }
    
    if(k == niter){
      indexes <- ((k-1)*sample_size+1):(dim(data_multi)[1])
    } else {
      indexes <- ((k-1)*sample_size+1):(k*sample_size)
    }
    
    data_multi[indexes,] %>%
      # change the MULTIPOINT to several lines with a POINT geometry
      st_cast('POINT') %>%
      # keep only one line, randomly
      ddply("fish_identifier", function(x) x[sample.int(x$n_points, 1),],
            .progress = ifelse(VERBOSE, "text", "none")) %>%
      # ddply changes the sf data.frame back to a simple data.frame so we change the df back to sf
      st_as_sf() -> data_multi_list[[k]]
  }
  
  txt <- "randomly sampled"
  
} else if (geometry_method == "centroid"){
  
  if (VERBOSE){
    msg <- "    - Considering the centroid from MULTIPOINT geometries\n" ; cat(msg) ; lines.to.cat <- c(lines.to.cat, msg)
  }
  
  data_multi %>% st_transform(3857) %>%
    st_centroid() %>%
    st_transform(4326) -> data_multi_list
  
  txt <- "centroid considered"
}

bind_rows(data_multi_list, data_na_point) %>%
  dplyr::arrange("fish_identifier") -> data

coordinates <- data.frame(st_coordinates(data))

N1 <- dim(data)[1]

data %>% mutate(lon = coordinates$X,
                lat = coordinates$Y) %>%
  mutate(scaled_lon = scale(lon, scale = T, center = T), 
         scaled_lat = scale(lat, scale = T, center = T)) -> data

N2 <- dim(data)[1]
Nna <- N1 - N2 ; rm(N1, N2)

sink(summaryName, append = T)
cat('\n\nPositions filter')
cat('\n----------------\n')
cat("\n  Number of entries with exact position available:", dim(data_na_point)[1]-Nna)
cat("\n  Number of entries with multiple positions available (", txt, "):", dim(data_multi)[1])
cat("\n  Number of entries with no available position (removed):", Nna)
cat("\n\n    - Number of entries after position filter:",dim(data)[1])
sink()

rm(coordinates, data_multi, data_multi_list, data_na_point, Nna) ; invisible(gc())