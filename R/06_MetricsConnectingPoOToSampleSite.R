# Load objects -------------------------------------------------------------------
library(data.table)
rerun <- FALSE

mydata_transformed_ID <- readRDS( file.path(wd$bin, "mydata_transformed.rds") ) %>%
  group_by(lon, lat) %>%
  dplyr::mutate(sampleSiteID = cur_group_id())

# Assemble lat/lon coordinates --------------------------------------------

# Because of how distance is most accurately measured (e.g., by distGeo), I'm
# going to convert my surface coordinates (currently in equal area projections)
# to lat/lon.

# All maps are projected across the same extent, so just load a random one.
coords_aea <- list.files(wd$tmp_df, pattern = "df_list.*rds$", full.names = T) %>%
  sample(1) %>%
  lapply(readRDS) %>%
  bind_rows() %>%
  dplyr::select(x,y) %>%
  distinct()
# Get unique coordinates to convert
coords_dd <- coords_aea %>%
  SpatialPointsDataFrame(coords = .[, 1:2], proj4string = CRS(myCRS)) %>%
  st_as_sf() %>%
  st_transform(crs = 4326) %>%
  st_coordinates() %>%
  as.data.frame %>%
  dplyr::rename(x_dd = X, y_dd = Y) %>%
  data.frame(., coords_aea)


# Functions ----------------------------------------------------------------
# Function that extracts distance and bearing relating each potential cell of
# geographic origin with individual's sample site.

getDistanceDirection <- function(
  rowNumber, dataframe, fromLat, toLat, fromLon, toLon, getDistance = TRUE,
  getDirection = TRUE, roundTo = 2
){

  p1 <- c( dataframe[ rowNumber, fromLon ], dataframe[ rowNumber, fromLat ] )
  p2 <- c( dataframe[ rowNumber, toLon ],   dataframe[ rowNumber, toLat ] )
  myResults <- list()
  if( getDistance == TRUE ){
    dist_km   <- round( geosphere::distGeo(p1, p2) / 1000 , roundTo) #Convert to km, round.
    myResults <- cbind(myResults, dist_km)
  }
  if( getDirection == TRUE ){
    theta_from_site   <- round( geosphere::bearing(p2, p1), roundTo)
    theta_from_origin <- round( geosphere::bearing(p1, p2), roundTo)
    myResults         <- cbind(myResults, theta_from_site, theta_from_origin)
  }
  return(myResults)

}

# Apply -------------------------------------------------------------------

# Measure distance and direction from every cell to each sampling location.
wd$tmp_dist <- file.path(wd$bin,"tmp_dist")
if(!dir.exists(wd$tmp_dist) ) dir.create(wd$tmp_dist)
if(length(list.files(wd$tmp_dist)) > 1 ) stop("Files already exist in this directory!")

# For each unique sample site, find distance and direction from each potential
# origin to that sample site.

toCombos <- mydata_transformed_ID %>%
  dplyr::select(sampleSiteID, lat,lon) %>%
  distinct() %>%
  arrange(sampleSiteID)

for(i in 1:nrow(toCombos)) {
  myfileLocation <- file.path(wd$tmp_dist, paste0("distDir_site_",i, ".csv"))
  if(rerun | !file.exists(myfileLocation)) {

    print(paste("Working on", i, "of", nrow(toCombos)))

    mydata_FromTo <- data.frame(coords_dd, toCombos[i,])

    mdf <- pbmcapply::pbmclapply( # mclapply, alternatively
      FUN = getDistanceDirection, mc.cores = parallel::detectCores() - 1,
      1:nrow(mydata_FromTo),
      dataframe =  mydata_FromTo,
      fromLat = "lat", toLat = "y_dd", fromLon = "lon", toLon = "x_dd",
      getDistance = TRUE, getDirection = TRUE
    ) %>%
      lapply(as.data.frame) %>%
      plyr::ldply() %>%
      data.frame(mydata_FromTo, .)

    fwrite(mdf, file = myfileLocation )

  }
}


