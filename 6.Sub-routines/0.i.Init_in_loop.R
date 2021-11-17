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

set.seed(seeds[i])

# glm_summaries[i,1] <- seeds[i]

msg <- "\14" ; cat(msg) ; lines.to.cat <- c(msg)
msg <- paste("-------- ITERATION", i, "/", length(seeds), " - SEED NB =", seeds[i], "--------") ; cat(msg) ; lines.to.cat <- c(lines.to.cat, msg)

# Initialize names
try(dir.create(file.path(OUTPUT_PATH, seeds[i], "Plots"), showWarnings = F, recursive = T))
diagnoPlotNames <- c(file.path(OUTPUT_PATH, seeds[i], "diagn_gam_te.png"),
                     file.path(OUTPUT_PATH, seeds[i], "diagn_gam_s.png"))
moranPlotName <- file.path(OUTPUT_PATH, seeds[i], "Moran_I_plot.png")
gamteSummary <- file.path(OUTPUT_PATH, seeds[i], "gam_te.rds")
gamsSummary <- file.path(OUTPUT_PATH, seeds[i], "gam_s.rds")
smoothPlotNames <- c(file.path(OUTPUT_PATH, seeds[i], "Smooth_latlon_te.png"),
                     file.path(OUTPUT_PATH, seeds[i], "Smooth_latlon_te_withpoints.png"),
                     file.path(OUTPUT_PATH, seeds[i], "Smooth_latlon_s.png"),
                     file.path(OUTPUT_PATH, seeds[i], "Smooth_latlon_s_withpoints.png"))
gamCoeffPlotNames <- c(file.path(OUTPUT_PATH, seeds[i], "gam_coeff_year.png"),
                       file.path(OUTPUT_PATH, seeds[i], "gam_coeff_quarter.png"),
                       file.path(OUTPUT_PATH, seeds[i], "gam_coeff_size_class.png"),
                       file.path(OUTPUT_PATH, seeds[i], "gam_coeff_fishing_mode.png"))
summaryName <- file.path(OUTPUT_PATH, seeds[i], "Processing_summary.txt")

vars <- c("month", "quarter", "lonlat", "lon", "lat")
plotsNames <- file.path(OUTPUT_PATH, seeds[i], "Plots", paste0("Kn_f-", vars, ".png"))

# Initialize data processing summary and generate README.txt
sink(summaryName)
cat("Execution time:", format(Sys.time()), "\n\n")
sink()

# Verify arguments
if (fad_fsc == F){
  fishing_mode_for_model <- "all"
}