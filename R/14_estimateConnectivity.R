
library(MigConnectivity)
library(terra)
library(tidyterra)
library(data.table)

source("R/00_Setup.R")

### targetDist -----

mydata_transformed <- readRDS( file.path(wd$bin, "mydata_transformed.rds") )

samplingLocations <- mydata_transformed %>%
  dplyr::select(lon, lat) %>%
  distinct %>%
  dplyr::mutate(locationCode = LETTERS[row_number()])

mydata_samples <- left_join(mydata_transformed, samplingLocations, by = join_by(lon, lat))

targetDist <- samplingLocations %>%
  st_as_sf(coords = c("lon", "lat"), crs = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") %>%
  st_distance()

targetDist_numeric <- matrix(data = as.numeric(targetDist), nrow = dim(targetDist)[1], ncol = dim(targetDist)[2])

subsampledLocationsPath <- file.path(wd$bin, "subsampledLocations.rds")
if(!file.exists(subsampledLocationsPath)) {
  n <- 10000
  wd$tmp_dist <- file.path(wd$bin,"tmp_dist")
  whichFiles <- list.files(wd$tmp_dist, pattern = "distDir_.*.csv", full.names = T)

  locations <- pbapply::pblapply(whichFiles, function(x){
    whichName <- gsub(".csv", "", gsub("distDir_", "", basename(x)))
    df <- fread(x)
    set.seed(42)
    df %>%
      sample_n(weight = value, size = n, replace = T) %>%
      mutate(rn = row_number())
  }) %>%
    bind_rows()
  saveRDS(locations, subsampledLocationsPath)
} else {
  locations <- readRDS(subsampledLocationsPath)
}


# # Estimate MC ---------
#
# Code estimates MC, not in parallel.
#
# tmpMCPath <- file.path(wd$bin, "tmp_mc")
# if(!dir.exists(tmpMCPath)) dir.create(tmpMCPath)
#
# pbapply::pblapply(1:1000, function(i){
#
#   savePath <- file.path(tmpMCPath, paste0(i, ".rds"))
#   if(file.exists(savePath)) return(NULL)
#
#   simLocations1 <- locations[locations$rn == i, ]
#   originDist <- simLocations1 %>%
#     st_as_sf(coords = c("x", "y"), crs = myCRS) %>%
#     st_distance()
#   origDist_numeric <- matrix(data = as.numeric(originDist), nrow = dim(originDist)[1], ncol = dim(originDist)[2])
#
#
#   ### relative abundance ---------
#   # One individual per origin site (centroid)
#   originRelAbund <- rep(1/dim(origDist_numeric)[1],dim(origDist_numeric)[1])
#
#   ### Transition probabilities ------
#
#   # Rows = number of origins (1/individual), sums to 1
#   # columns = number of targets
#
#   nrow <- dim(origDist_numeric)[1]
#   ncol <- dim(targetDist_numeric)[1]
#
#   psi <- matrix(data = rep(0, nrow*ncol), nrow = nrow, ncol = ncol)
#   # Add names for clearer indexing.
#   rownames(psi) <- simLocations1$ID
#   colnames(psi) <- samplingLocations$locationCode
#
#   for(j in 1:nrow(psi)) {
#     # For a given row, find which column needs to have a 1 in it.
#     whichCol <- unlist(mydata_samples[mydata_samples$SampleName == rownames(psi)[j], "locationCode"])
#     psi[j, which(colnames(psi) == whichCol)] <- 1
#   }
#
#   ### Calculate ------
#   MC_ests <- calcMC(origDist_numeric, targetDist_numeric, originRelAbund, psi, sampleSize = dim(origDist_numeric)[1])
#
#   # Save.
#   saveRDS(MC_ests, savePath)
#
# })
#


# Estimate mc in parallel -------------------------------------------------------


n_cores <- detectCores()
cluster <- makeCluster(n_cores - 1)
registerDoParallel(cluster)
foreach(i = 1:10000, .packages = c("tidyverse", "sf", "MigConnectivity")) %dopar% {

  savePath <- file.path(tmpMCPath, paste0(i, ".rds"))
  if(!file.exists(savePath)) {


    simLocations1 <- locations[locations$rn == i, ]
    originDist <- simLocations1 %>%
      st_as_sf(coords = c("x", "y"), crs = myCRS) %>%
      st_distance()
    origDist_numeric <- matrix(data = as.numeric(originDist), nrow = dim(originDist)[1], ncol = dim(originDist)[2])


    ### relative abundance ---------
    # One individual per origin site (centroid)
    originRelAbund <- rep(1/dim(origDist_numeric)[1],dim(origDist_numeric)[1])

    ### Transition probabilities ------

    # Rows = number of origins (1/individual), sums to 1
    # columns = number of targets

    nrow <- dim(origDist_numeric)[1]
    ncol <- dim(targetDist_numeric)[1]

    psi <- matrix(data = rep(0, nrow*ncol), nrow = nrow, ncol = ncol)
    # Add names for clearer indexing.
    rownames(psi) <- simLocations1$ID
    colnames(psi) <- samplingLocations$locationCode

    for(j in 1:nrow(psi)) {
      # For a given row, find which column needs to have a 1 in it.
      whichCol <- unlist(mydata_samples[mydata_samples$SampleName == rownames(psi)[j], "locationCode"])
      psi[j, which(colnames(psi) == whichCol)] <- 1
    }

    ### Calculate ------
    MC_ests <- calcMC(origDist_numeric, targetDist_numeric, originRelAbund, psi, sampleSize = dim(origDist_numeric)[1])

    # Save.
    saveRDS(MC_ests, savePath)
  }
}
stopCluster(cluster)


# Load and summarise ------------------------------------------------------

n_cap <- 1000

MC_files <- list.files(tmpMCPath, full.names = T)[1:n_cap]
MC_ests_all <- lapply(MC_files, function(x) {
  out <- readRDS(x)
  return(out)
  }) %>%
  unlist()

quantile(MC_ests_all, c(0.025, 0.5, 0.975)) %>%
  signif(1)
