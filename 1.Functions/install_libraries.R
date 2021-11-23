

#' Title
#'
#' @param srcUsedPackages : require packages
#' @param loadPackages : load or not after installation
installAndLoad_packages <-function(srcUsedPackages = srcUsedPackages, loadPackages =TRUE,
                                   verbose = F)
{
  #---Installed packages
  userPackages <- as.data.frame(installed.packages())
  
  #---Packages to install
  neededPackages <- srcUsedPackages[!srcUsedPackages %in% userPackages$Package]
  
  #---Download and installing packages
  if (verbose){
    cat("Installing  and Loading of required packages....\n")
  }
  if(length(neededPackages) > 0)
  {
    for(i in 1:length(neededPackages))
    {
      install.packages(neededPackages[i], quiet = TRUE, verbose = TRUE,
                       repos = "http://cran.us.r-project.org")
    }
  }
  correclty_installed <- srcUsedPackages[srcUsedPackages %in% userPackages$Package]
  not_installed <- srcUsedPackages[!srcUsedPackages %in% userPackages$Package]
  if(length(not_installed) >0)
  {
    warning(paste0("Error at installation of the following package(s) : \n", not_installed, "\n"))
  }
  
  #---Loading packages
  if(loadPackages == TRUE)
  {
    notLoaded <- NA
    for(i in 1: length(correclty_installed))
    {
      res <- library(package = eval(correclty_installed[i]), character.only = TRUE, logical.return = TRUE, quietly = TRUE, warn.conflicts = FALSE)
      if(res == FALSE)
      {
        notLoaded <- na.omit(c(notLoaded, correclty_installed[i]))
      }
    }
    if(!is.na(notLoaded))
    {
      warning(paste0("Error at the loading of the The following package(s) : \n", notLoaded))
    }
  }
}

