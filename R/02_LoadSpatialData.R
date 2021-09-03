
# Setup -------------------------------------------------------------------
download_GADM <- FALSE

# This script assumes that there are candidate isoscapes downloaded from isoMAP
# somewhere in the wd$data directory.
reload_isoscapes <- FALSE

# This script assumes that there are IUCN rangemaps somewhere in the wd$iucn
# directory.
reload_IUCN_rangemaps <- FALSE

# Load GADM data ----------------------------------------------------------
if(download_GADM == TRUE){

  library(rmapshaper)
  message("Loading GADM data...")

  # First, make NoAm_sf (an "sf" object) fit extent object (my_extent).
  # my_extent_aea_st <- st_bbox(st_transform(st_as_sfc(st_bbox(my_extent, crs = 4326)), myCRS))
  # saveRDS(my_extent_aea_st, file = file.path(wd$bin, "my_extent_aea_st.rds"))

  # Get GADM data to state level.
  USA <- raster::getData('GADM', path = wd$bin, country='USA', level=0)
  MEX <- raster::getData('GADM', path = wd$bin, country='MEX', level=0)
  CAN <- raster::getData('GADM', path = wd$bin, country='CAN', level=0)
  GTM <- raster::getData('GADM', path = wd$bin, country='GTM', level=0)

  # Prepare to remove areas outside of desired extent.
  # ## Remove Hawaii
  USA_1 <- raster::getData('GADM', path = wd$bin, country='USA', level=1)
  hawaii <- USA_1[USA_1@data$NAME_1 == "Hawaii",] # Select Hawaii
  hawaii_simpl <- hawaii %>%
    sf::st_as_sf() %>%
    st_transform(crs = myCRS) %>%
    st_simplify(preserveTopology = FALSE, dTolerance = 5000) %>%
    st_buffer(dist = 1e6)
  # ## Remove water bodies.
  USA_2 <- raster::getData('GADM', path = wd$bin, country='USA', level=2)
  MEX_2 <- raster::getData('GADM', path = wd$bin, country='MEX', level=2)
  CAN_2 <- raster::getData('GADM', path = wd$bin, country='CAN', level=2)
  GTM_2 <- raster::getData('GADM', path = wd$bin, country='GTM', level=2)
  waterbodies <- lapply(list(USA_2,MEX_2,CAN_2,GTM_2), function(x){
    x[(x$ENGTYPE_2) == "Water body",]
  }) %>%
    do.call(rbind, .) %>%
    sf::st_as_sf() %>%
    st_transform(crs = myCRS) %>%
    st_simplify(preserveTopology = FALSE, dTolerance = 5000)

  # Combine into one polygon, convert to sf object.
  NoAm <- raster::bind(
    MEX, USA, CAN, GTM#, BLZ, SLV, HND, NIC
  ) %>%
    sf::st_as_sf(.) %>%
    st_transform(crs = myCRS) %>%
    st_simplify(preserveTopology = FALSE, dTolerance = 5e3) %>%
    st_difference(., hawaii_simpl) # Remove Hawaii

  saveRDS(NoAm, file = file.path(wd$bin, "NoAm.rds"))

  NoAm_boundary_aea <- NoAm %>%
    st_buffer(dist = 5e4) %>%
    rmapshaper::ms_erase(., waterbodies)  # Remove water bodies

  saveRDS(NoAm_boundary_aea, file = file.path(wd$bin, "NoAm_boundary_aea.rds"))

} else message("Not redownloading GADM Data...")


# Load isoscapes ----------------------------------------------------------

if(reload_isoscapes == TRUE){
  message("reloading isoscapes...")

  # Location of my local isoscape and NoAm boundary directory is stored in object 'bigDataStorage'
  source("~/WISN_dD/.Rprofile")

  NoAm_boundary_aea <- readRDS(
    file.path(bigDataStorage, "NoAm_maps", "NoAm_boundary_aea.rds")
    ) %>%
    st_as_sf() %>%
    st_transform(crs = myCRS)

  # Function to extend/crop/mask by above.
  ECM <- function(rasterLayer){
    rasterLayer %>%
      raster::projectRaster(., crs = myCRS) %>%
      raster::extend( ., my_extent_aea ) %>%
      raster::crop(   ., my_extent_aea ) %>%
      # Following line not needed-- no need to resample for one surface.
      #raster::resample(., refIsoscape) %>%
      raster::mask(   ., NoAm_boundary_aea  )
  }

  myisoscape <- list()
  myisoscape$directory          <- file.path(bigDataStorage, "isoscapes")
  myisoscape$path_pattern       <- "66100"
  myisoscape$isoscape_pattern   <- "predkrig.tiff$"
  myisoscape$sd_pattern         <- "stdkrig.tiff$"

  myisoscape$isoscape <- list.files(
    myisoscape$directory, pattern =  myisoscape$isoscape_pattern,
    recursive = TRUE, full.names = TRUE
  ) %>%
    grep(myisoscape$path_pattern, ., value = TRUE) %>%
    raster::raster(.) %>%
    ECM(.)

  myisoscape$sd <- list.files(
    myisoscape$directory, pattern = myisoscape$sd_pattern,
    recursive = TRUE, full.names = TRUE
  ) %>%
    grep(myisoscape$path_pattern, ., value = TRUE) %>%
    raster::raster(.) %>%
    ECM(.)

  # Save.
  save(myisoscape, file = file.path(wd$bin, "my_isoscapes.RData"))

} else message("Not reloading isoscapes, loading saved version...")


# Load and buffer IUCN Rangemaps -----------------------------------------------------
if(reload_IUCN_rangemaps == TRUE){

  load(file.path(wd$bin, "my_isoscapes.RData"), verbose = TRUE)
  if(!exists("NoAm_boundary_aea")) {
    NoAm_boundary_aea <- readRDS( file.path(wd$bin, "NoAm_boundary_aea.rds") )
  }

  locationIUCNData <- dir(wd$data, pattern = "redlist_species_data", full.names = T)

  # Load IUCN rangemap.
  # Convert to simple features object and reproject to myCRS.
  IUCNmaps <-
    locationIUCNData %>%
    st_read(layer = "data_0") %>%
    st_as_sf(crs = 4326) %>%
    st_transform(crs = myCRS) %>%
    st_simplify(preserveTopology = TRUE, dTolerance = 5000) %>%
    st_make_valid()

  range_tidy <- IUCNmaps %>%
    st_buffer(dist = 10e3) %>%
    st_combine()

  # Convert buffered rangemaps to rasters with appropriate
  ex_rast <- myisoscape$isoscape
  ex_rast[] <- 1
  range_raster <- raster::mask(
    ex_rast,
    mask = as_Spatial(range_tidy),
    updatevalue = NA
    ) %>%
    raster::mask(., NoAm_boundary_aea) %>%
    raster::crop(., my_extent_aea)

  saveRDS(range_raster, file = file.path(wd$bin, "range_raster.rds" ))

} else message("Not reloading IUCN rangemaps...")
