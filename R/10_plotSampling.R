library(data.table)
library(terra)

allResults <- fread(file.path(wd$bin, "allResults.csv")) %>%
  dplyr::mutate(
    westOrigin = case_when(
      probWestOrigin >= 0.8 ~ 1,
      probWestOrigin <= 0.2 ~ 0,
      TRUE ~ as.numeric(NA)
    )
  )

pts <- allResults %>%
  group_by(lat_aea, lon_aea) %>%
  slice(1) %>%
  dplyr::select(lat_aea, lon_aea) %>%
  na.omit()

# if(!exists("NoAm")) {
#   NoAm <- readRDS( file.path(wd$bin, pattern = "NoAm_sf.rds") ) %>%
#     st_simplify(dTolerance = 100) %>%
#     st_transform(myCRS)
# }

countries <- rnaturalearthdata::countries110 %>%
  st_as_sf() %>%
  st_transform(myCRS) %>%
  st_crop(my_extent_aea) %>%
  dplyr::filter(continent %in% c("North America", "South America"))

states <- rnaturalearthhires::states10 %>%
  st_as_sf() %>%
  st_transform(myCRS) %>%
  st_crop(my_extent_aea)

# Load rangemaps
locationIUCNData <- dir(wd$data, pattern = "redlist_species_data", full.names = T)
IUCNmaps <-
  locationIUCNData %>%
  st_read(layer = "data_0") %>%
  st_as_sf(crs = 4326) %>%
  st_transform(crs = myCRS) %>%
  st_simplify(preserveTopology = TRUE, dTolerance = 5000) %>%
  st_make_valid()


# inset sampling plot -----------------------------------------------------

# Specify bounds of inset plot
inset_xlim <- c(0.95e6, 1.55e6)
inset_ylim <- c(-1.65e6, -1.0e6)

## Load wetland data -----
if(!file.exists( file.path(wd$bin, "wetlandRast.tif") )) {
  # https://www.fws.gov/node/264847
  wet_ga <- sf::st_read(file.path(wd$data, "GA_shapefile_wetlands"), "GA_Wetlands")
  wet_fl <- sf::st_read(file.path(wd$data, "FL_shapefile_wetlands"), "FL_Wetlands") %>%
    dplyr::select(names(wet_ga))
  trans <- lapply(list(wet_fl, wet_ga), st_transform, crs = myCRS )
  #trans[[1]] <- trans[[1]][,names(trans[[1]]) %in% names(trans[[2]])]
  wet_all <- do.call(rbind, trans) %>%
    st_simplify(preserveTopology = TRUE, dTolerance = 5000) %>%
    st_transform(crs = myCRS) %>%
    st_crop(xmin = inset_xlim[1] - 5e5, xmax = inset_xlim[2] + 5e5, ymin = inset_ylim[1] - 5e5, ymax = inset_ylim[2] + 5e5)
  saveRDS(wet_all, file  = file.path(wd$bin, "wet_all.rds"))

  # Rasterize wetland data.
  library(terra)
  r <- rast(nrows = 3000, ncols = 3000, crs = myCRS , xmin = inset_xlim[1], xmax = inset_xlim[2], ymin = inset_ylim[1], ymax = inset_ylim[2])
  rasterize(wet_all, r, field = "WETLAND_TY", filename = file.path(wd$bin, "wetlandRast.tif"), overwrite = T)

} else {
 if(!exists("wet_rast")) {
   r2 <- rast( file.path(wd$bin, "wetlandRast.tif") )
   wet_rast <- r2 %>%
     as.data.frame(xy = TRUE, na.rm =T) %>%
     dplyr::mutate(
       value = case_when(
         WETLAND_TY %in% c("Lake", "Freshwater Pond", "Riverine") ~ "Open water",
         WETLAND_TY %in% c( "Freshwater Emergent Wetland", "Freshwater Forested/Shrub Wetland") ~ "Wetland",
         TRUE ~ as.character(NA)
       )
     ) %>%
     na.omit()
   }
}

library(ggspatial)

insetMap <- ggplot() +
  geom_sf(states, mapping = aes(), fill = "white") +
  # Plot wetlands
  geom_tile(wet_rast, mapping = aes(x=x,y=y,color=value, fill = value)) +
  scale_fill_manual(
    breaks = c("Wetland", "Open water"),
    values = c("#98D7D5","#6C7FB2"),
    labels = c("Freshwater wetland", "Open water")
    ) +
  scale_color_manual(
    breaks = c("Wetland", "Open water"),
    values = c("#98D7D5","#6C7FB2"),
    labels = c("Freshwater wetland", "Open water")
    ) +
  # Plot state lines
  geom_sf(states, mapping = aes(), fill = NA) +
  # Plot points.
  geom_point(
    pts, mapping = aes(x=lon_aea, y = lat_aea),
    shape = 13, size = 3, color= "black") +

  # Add scale bar and arrow.
  annotation_scale(location = "tr", width_hint = 0.25) +
  annotation_north_arrow(location = "tr", which_north = "true",
                         pad_x = unit(0.0, "in"), pad_y = unit(0.2, "in"),
                         style = north_arrow_minimal) +
  # Other details
  coord_sf(
    xlim = c(inset_xlim[1] + 3e4, inset_xlim[2]),
    ylim = c(inset_ylim[1], inset_ylim[2] - 4e4),
    ) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    axis.title = element_blank(),
    axis.text = element_blank(),
    legend.position = c(0.2,0.65),
    legend.title = element_blank()
    )


# continent map -----------------------------------------------------------
breeding <- IUCNmaps %>%
  dplyr::filter(SEASONAL %in% c(1,2)) %>%
  st_union()
nonbreeding <- IUCNmaps %>%
  dplyr::filter(SEASONAL %in% c(1,3)) %>%
  st_union()

seasoncolors <- c(
  "#688E26",
  "#FAA613",
  "#550527"
)

main_map <- ggplot() +
  geom_sf(countries,   mapping = aes(), fill = NA) +
  geom_sf(dplyr::filter(IUCNmaps, SEASONAL %in% c(2)), mapping = aes(), alpha = 0.5, color = seasoncolors[1], fill = seasoncolors[1]) +
  geom_sf(dplyr::filter(IUCNmaps, SEASONAL %in% c(3)), mapping = aes(), alpha = 0.5, color = seasoncolors[2], fill = seasoncolors[2]) +
  geom_sf(dplyr::filter(IUCNmaps, SEASONAL %in% c(1)), mapping = aes(), alpha = 0.5, color = seasoncolors[3], fill = seasoncolors[3]) +
  geom_text(mapping = aes(x=0, y = 1.5e6),      label = "Breeding",    color = seasoncolors[1], hjust = 0.5, size = 4) +
  geom_text(mapping = aes(x=0, y = -0.5e6),     label = "Nonbreeding", color = seasoncolors[2], hjust = 0.5, size = 4) +
  geom_text(mapping = aes(x=-1.48e6, y = 0.5e6), label = "Year-round",  color = seasoncolors[3], hjust = 0.5, size = 4) +
  coord_sf(
    xlim = c(-4.5e6, 4.5e6),
    ylim = c(-3.8e6, 4e6)
  ) +
  theme_void() +
  theme(
    panel.background = element_rect(fill = "white", color = "black", linewidth = 0.5)
  ) +
  # Draw line around inset
   geom_rect(
     color = "black",
     linewidth = 0.75,
     fill = NA,
     mapping = aes(
      xmin = inset_xlim[[1]] - 1e5,
      xmax = inset_xlim[[2]] + 1e5,
      ymin = inset_ylim[[1]] - 1e5,
      ymax = inset_ylim[[2]] + 1e5
     )
  )


# Load image --------------------------------------------------------------

library(jpeg)
snipePhoto <- magick::image_read("data/Photo 278588771, (c) Andrew Thomas, some rights reserved (CC BY-NC), uploaded by Andrew Thomas_cropped.jpeg")
photo1 <- draw_image(snipePhoto)

# Combine, no photo-----------------------------------------------------------------

library(cowplot)

togetherMap <- ggdraw(insetMap) +
  draw_plot(
    main_map,
    x = 0, y = 0,
    width  = 0.5, height = 0.55
    ) +
  draw_plot_label(
    c("B", "A"),
    c(0.05, 0.2),
    c(0.55, 0.99),
    size = 12
  )

ggsave(togetherMap, file= file.path(wd$figs, "samplingmap.png"), width = 6, height = 4, dpi = 400)


# Add photo ---------------------------------------------------------------

togetherMap2 <- plot_grid(
  {ggdraw() + photo1},
  {
    ggdraw(insetMap) +
      draw_plot(
        main_map,
        x = -0.05, y = -0.05,
        width  = 0.5, height = 0.55
      )
  }
  ) +
  draw_plot_label(
    c("A",  "C", "B"),
    c(0.05, 0.48, 0.50),
    c(0.99, 0.44, 0.99),
    size = 12
  )

ggsave(togetherMap2, file= file.path(wd$figs, "samplingmap2.png"), width = 8, height = 4, dpi = 400)



# Extra analysis: identify countries in nonbreeding range -----------------

main_map +
  geom_sf_label(data = countries, aes(label = name))


st_intersection(countries, IUCNmaps) %>%
  group_by(sovereignt) %>%
  slice(1) %>%
  View()

# How many Canadian provinces?
can <- geodata::gadm("CAN", level = 1, path = "tmp") %>%
  st_as_sf() %>%
  st_transform(crs = myCRS) %>%
  st_simplify(preserveTopology = TRUE, dTolerance = 5000) %>%
  st_make_valid()

st_intersection(can, IUCNmaps) %>%
  group_by(NAME_1) %>%
  slice(1)

main_map +
  geom_sf(can, mapping = aes(), fill = NA) +
  geom_sf_label(can, mapping = aes(label = NAME_1))
