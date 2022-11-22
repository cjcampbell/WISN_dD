# Load objects -------------------------------------------------------------------
library(data.table)
rerun <- FALSE

wd$tmp_df <- file.path(wd$bin,"tmp_df")

mydata_transformed <- readRDS( file.path(wd$bin, "mydata_transformed.rds") )

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

# Assemble lat/lon coordinates --------------------------------------------

# Because of how distance is most accurately measured (e.g., by distGeo), I'm
# going to convert my surface coordinates (currently in equal area projections)
# to lat/lon.

rerun <- FALSE

whichFiles <- list.files(wd$tmp_df, pattern = "df_list.*rds$", full.names = T)

wd$tmp_dist <- file.path(wd$bin,"tmp_dist")
if(!dir.exists(wd$tmp_dist) ) dir.create(wd$tmp_dist)
if(length(list.files(wd$tmp_dist)) > 1 ) stop("Files already exist in this directory!")

for( x in whichFiles ) {

  s1 <- x %>%
    lapply(readRDS) %>%
    bind_rows() %>%
    dplyr::filter(method == "raw")

  myfileLocation <- file.path(wd$tmp_dist, paste0("distDir_",s1$ID[1], ".csv"))

  if( rerun | !file.exists(myfileLocation) ) {

    print(paste("Working on", s1$ID[1]))

    coords_dd <- s1 %>%
      SpatialPointsDataFrame(coords = .[, 1:2], proj4string = CRS(myCRS)) %>%
      st_as_sf() %>%
      st_transform(crs = 4326) %>%
      st_coordinates() %>%
      as.data.frame %>%
      dplyr::rename(x_dd = X, y_dd = Y) %>%
      data.frame(., s1)

    mydata_FromTo <- left_join(
      coords_dd,
      dplyr::select(mydata_transformed, SampleName, lat, lon),
      by = c("ID" = "SampleName")
      )

    mdf <- pbmcapply::pbmclapply( # mclapply, alternatively
      FUN = getDistanceDirection, mc.cores = parallel::detectCores() - 1,
      1:nrow(mydata_FromTo),
      dataframe =  mydata_FromTo,
      fromLat = "y_dd", toLat = "lat", fromLon = "x_dd", toLon = "lon",
      getDistance = TRUE, getDirection = TRUE
    ) %>%
      lapply(as.data.frame) %>%
      plyr::ldply() %>%
      data.frame(mydata_FromTo, .)

    fwrite(mdf, file = myfileLocation )
  } else {
    print(paste("Already done with", s1$ID[1]))
  }
}



# Extract statistics for each individual. ---------------------------------

csumThreshold <- 0.25

statsOut <- list.files(wd$tmp_dist, pattern = "distDir_.*.csv", full.names = T) %>%
  lapply(., function(x){
    df <- fread(x)
    # Find minimum distance ----
    ## Convert to cumulative sum
    df <- df %>%
      arrange(value) %>%
      mutate(
        csum = cumsum(value),
        orig = case_when(csum >= csumThreshold ~ 1, TRUE ~ 0)
        )
    ## Find minimum distance traveled from nearest point over cumulative sum threshold.
    minDist <- df %>%
      dplyr::filter(orig == 1) %>%
      dplyr::arrange(dist_km) %>%
      slice(1)
    ## Extract details ready for plotting.
    sampleSiteLocation <-
      sf::sf_project(from = st_crs(4326), to = myCRS, minDist[,c("lon", "lat")])
    minDistDeets <- data.frame(
      minDist,
      lon_aea = sampleSiteLocation[1],
      lat_aea = sampleSiteLocation[2]
    )

    # Define potential region of origin ----
    set.seed(42)
    fromWest <- df %>%
      sample_n(weight = value, size = 1e6, replace = T) %>%
      group_by(x,y, x_dd, y_dd) %>%
      dplyr::summarise(n=n()) %>%
      dplyr::ungroup() %>%
      dplyr::mutate(isEast = case_when(x_dd > -100 ~ 1, TRUE ~ 0)) %>%
      dplyr::group_by(isEast) %>%
      dplyr::summarise(propSim = sum(n) / 1e6) %>%
      dplyr::filter(isEast == 0) %>%
      dplyr::select(propSim) %>%
      unlist

    minDistDeets <- data.frame(minDistDeets, probWestOrigin = fromWest)
    return(minDistDeets)

  }) %>%
  bind_rows()

fwrite(statsOut, file = file.path(wd$bin, "distDirStats.csv"), row.names = F)
