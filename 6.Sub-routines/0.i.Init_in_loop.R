#'#*******************************************************************************************************************
#'@author : Amael DUPAIX
#'@update : 2021-11-17
#'@email : amael.dupaix@ens-lyon.fr
#'#*******************************************************************************************************************
#'@description :  Initialization routine INSIDE THE FOR LOOP
#' load libraries, prepare files names, etc
#'#*******************************************************************************************************************
#'@revisions
#'#*******************************************************************************************************************

set.seed(seed.i)

if (VERBOSE){
  msg <- "\14" ; cat(msg) ; lines.to.cat <- c(msg)
  msg <- paste("-------- ITERATION", i, "/", length(seeds), " - SEED NB =", seed.i, "--------") ; cat(msg) ; lines.to.cat <- c(lines.to.cat, msg)
}

# Initialize names
if (generate_plots){
  try(dir.create(file.path(OUTPUT_PATH, seed.i, "Plots"), showWarnings = F, recursive = T))
  smoothPlotNames <- c(file.path(OUTPUT_PATH, seed.i, "Smooth_latlon_s.png"),
                       file.path(OUTPUT_PATH, seed.i, "Smooth_latlon_s_withpoints.png"),
                       file.path(OUTPUT_PATH, seed.i, "Smooth_latlon_te.png"),
                       file.path(OUTPUT_PATH, seed.i, "Smooth_latlon_te_withpoints.png"))
  gamCoeffPlotNames <- c(file.path(OUTPUT_PATH, seed.i, "gam_coeff_year.png"),
                         file.path(OUTPUT_PATH, seed.i, "gam_coeff_quarter.png"),
                         file.path(OUTPUT_PATH, seed.i, "gam_coeff_size_class.png"),
                         file.path(OUTPUT_PATH, seed.i, "gam_coeff_fishing_mode.png"))
  vars <- c("month", "quarter", "lonlat", "lon", "lat")
  plotsNames <- file.path(OUTPUT_PATH, seed.i, "Plots", paste0("Kn_f-", vars, ".png"))
}

diagnoPlotNames <- c(file.path(OUTPUT_PATH, seed.i, "diagn_gam_te.png"),
                     file.path(OUTPUT_PATH, seed.i, "diagn_gam_s.png"))
moranPlotName <- file.path(OUTPUT_PATH, seed.i, "Moran_I_plot.png")

gamteSummary <- file.path(OUTPUT_PATH, seed.i, "gam_te.rds")
gamsSummary <- file.path(OUTPUT_PATH, seed.i, "gam_s.rds")
summaryName <- file.path(OUTPUT_PATH, seed.i, "Processing_summary.txt")


# Initialize data processing summary and generate README.txt
sink(summaryName)
cat("Execution time:", format(Sys.time()), "\n\n")
sink()

# Verify arguments
if (fad_fsc == F){
  fishing_mode_for_model <- "all"
}