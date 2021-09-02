
# Load libraries required for analyses.
library(isocat)
library(dplyr)
library(sf)

# Load libraries relied on for a few analyses.
library(assignR)
library(chron)
library(geosphere)
library(ggmap)
library(ggplot2)
library(lubridate)
library(lwgeom)
library(measurements)
library(purrr)
library(readxl)
library(rgdal)
library(rmapshaper)
library(smatr)
library(stringr)
library(tidyr)

# Make an object to help navigate the subdirectories.
# Be sure to fill this in if you're reproducing these analyses!
my_dir_path <- "/Users/cjcampbell/WISN_dD"
if(!exists("wd")) { wd <- list() }
wd$R       <- file.path( my_dir_path, "R" )
wd$bin     <- file.path( my_dir_path, "bin" )
wd$data    <- file.path( my_dir_path, "data" )
wd$figs    <- file.path( my_dir_path, "figs" )

# Check for presence of subdirectories. Create if needed.
invisible({
  lapply(wd, function(i) if( dir.exists(i) != 1 ) dir.create(i) )
})


# Define extent of spatial analysis.
# Units == meters
my_extent_aea <- raster::extent(
  -60e5, 51e5,
  -52e5, 60e5
)

# CRS for aea projection:
myCRS <- "+proj=aea +lat_0=40 +lon_0=-96 +lat_1=20 +lat_2=60 +x_0=0 +y_0=0
+ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
