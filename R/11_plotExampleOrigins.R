source("~/WISN_dD/R/00_Setup.R")
library(ggstar)
library(data.table)
library(ggpubr)
wd$tmp_df <- file.path(wd$bin,"tmp_df")
wd$tmp_dist <- file.path(wd$bin,"tmp_dist")
csumThreshold <- 0.25
countries <- rnaturalearth::countries110 %>%
  st_as_sf() %>%
  st_transform(myCRS)

# Load sample locations.
mydata_transformed <- readRDS( file.path(wd$bin, "mydata_transformed.rds") )
pts <- as.data.frame(sf::sf_project(pts = as.matrix(mydata_transformed[,c("lon", "lat")]), from = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs ", to = myCRS))
names(pts) <- c("lon_aea", "lat_aea")
mydata_transformed <- data.frame(mydata_transformed, pts )

# Load IUCN rangemap.
range_raster <- readRDS( file.path(wd$bin, "range_raster.rds" ) )
locationIUCNData <- dir(wd$data, pattern = "redlist_species_data", full.names = T)
IUCNmaps <-
  locationIUCNData %>%
  st_read(layer = "data_0") %>%
  st_as_sf(crs = 4326) %>%
  st_transform(crs = myCRS) %>%
  st_simplify(preserveTopology = TRUE, dTolerance = 5000) %>%
  st_make_valid()
breeding <- IUCNmaps %>% filter(SEASONAL %in% c(1,2))
other <- IUCNmaps %>% filter(SEASONAL %in% c(3,4))

# Load world map
wrld <- rnaturalearthhires::countries10 %>%
  st_as_sf() %>%
  st_transform(crs = myCRS) %>%
  st_crop(my_extent_aea)

# Load maps for names
mapNames <- list.files(
  file.path( wd$bin, "maps"), pattern = "Combined.*.grd$", full.names = TRUE, recursive = T) %>%
  raster::stack() %>%
  names()

# Specify plot structures.
whichOnes <- c(323, 3, 107)


# Make indiv plots --------------------------------------------------------

myplots <- lapply(whichOnes, function(x){

  # Setup
  writeLines(paste("Working on number", x))
  thisFile <- list.files(wd$tmp_dist, pattern = "distDir_.*.csv", full.names = T)[x]
  df <- fread(thisFile)
  ID <- gsub(".csv","", gsub("distDir_","", basename(thisFile)))
  probSurface <- list.files(file.path(wd$bin, "maps"), pattern = paste0("df_list_", ID, ".csv"), full.names = T) %>%
    fread() %>%
    dplyr::filter(method == "raw") %>%
    arrange(value) %>%
    mutate(
      csum = cumsum(value),
      orig = case_when(csum >= csumThreshold ~ 1, TRUE ~ 0)
    )

  stopifnot(probSurface$ID[1] == ID)

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
  deets <- df %>%
    sample_n(weight = value, size = 1e6, replace = T) %>%
    group_by(x,y, x_dd, y_dd) %>%
    dplyr::summarise(n=n()) %>%
    as.data.frame()

  myline <- matrix(
    data = c(-100, -100, 0, 90),
    byrow = F, nrow = 2) %>%
    st_linestring() %>%
    st_sfc(crs = 4326) %>%
    st_transform(crs = myCRS)

  allotmentDeets <- deets %>%
    dplyr::mutate(isEast = case_when(x_dd > -100 ~ 1, TRUE ~ 0)) %>%
    dplyr::group_by(isEast) %>%
    dplyr::summarise(propSim = sum(n) / 1e6) %>%
    dplyr::mutate(
      propSim = signif(propSim, 2),
      xposition = case_when(isEast == 1 ~ 8e5, TRUE ~ -8e5)
    )

  p <- ggplot() +
    geom_sf(
      wrld, mapping=aes(),
      fill = "white", color = "grey20", size = 0.25) +
    geom_tile(probSurface, mapping=aes(x=x,y=y,fill=csum, color = csum)) +
    scale_fill_viridis_c(
      "Cumulative sum of probability of origin",
      option = "turbo", direction = 1,
      guide = guide_colourbar(title.position = "top")
    ) +
    scale_color_viridis_c(
      "Cumulative sum of probability of origin",
      option = "turbo", direction = 1,
      guide = guide_colourbar(title.position = "top")
    ) +
    geom_sf(breeding, mapping = aes(),
            fill = NA, color = "black", size = 1) +
    geom_sf(other, mapping = aes(),
            fill = "grey70", alpha = 0.5, color = NA, size = 1) +
    geom_star(
      data = dplyr::filter(mydata_transformed, SampleName == ID),
      mapping = aes(x=lon_aea, y = lat_aea),
      fill = "black", color = "black", size = 3.5
    ) +
    geom_sf(myline, mapping = aes(), color = "black", linetype = "dotted", size = 2) +
    scale_x_continuous(
      breaks = seq(-200,200,by=20)
    ) +
    scale_y_continuous(
      breaks = seq(0,90,by=10)
    ) +
    coord_sf(
      xlim = c(extent(range_raster)[1]+1.5e6, extent(range_raster)[2]-2.4e6),
      ylim = c(extent(range_raster)[3]+3.5e6, extent(range_raster)[4]-1.8e6)
    ) +
    theme_bw() +
    theme(
      plot.background = element_blank(),
      panel.background = element_blank(),
      legend.position = c(0.2,0.2),
      legend.title = element_text(vjust = 0.5),
      strip.background = element_blank(),
      axis.title = element_blank()
    )
  return(p)
})

# Combine -----------------------------------------------------------------

appendCenteredLegend <- function(listOfPlots, legendArgs = NULL) {
  # Function that takes list of plots (that could share a legend),
  # returns a list of plots with legends removed and a centered legend horizontal legend appended.

  # Extract and format legend.
  legendTarget <- listOfPlots[[1]] +
    theme(
      legend.position = "bottom",
      legend.title.align = 0.5,
      legend.title = element_text(hjust = 0.5)
    )
  if(!is.null(legendArgs)) {
    legendTarget <- legendTarget + legendArgs
  }
  mylegend0 <- ggpubr::get_legend( legendTarget )

  # Remove legends from plots.
  listOfPlots2 <- lapply(listOfPlots, function(x) {
    x <- x +
      theme(legend.position = "none")
    return(x)
  })

  # Append legend to plotlist.
  listOfPlots2[[ length(listOfPlots2) + 1 ]] <- mylegend0

  # Return new plotlist.
  return(listOfPlots2)

}

myplotList <- appendCenteredLegend(
  myplots,
  legendArgs = list(
    theme(legend.key.width = unit(1.8, "cm"))
    )
  )
bigplot <- ggarrange(
  plotlist = myplotList,
  ncol = 1, heights = c(1,1,1,0.2),
  labels = c(LETTERS[1:3], as.character(NA))
  )
ggsave(bigplot, filename = file.path(wd$figs, "ExampleOrigins.png"), width = 4, height = 10)
