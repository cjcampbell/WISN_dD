library(terra)
library(tidyterra)
library(ggpubr)

load( file.path(wd$bin, "my_isoscapes.RData") )
range_raster <- readRDS( file.path(wd$bin, "range_raster.rds" ) )
iso_cropped <- mask(myisoscape$isoscape, range_raster) %>% rast()

iso_feather <- iso_cropped*1.16+23.57

myline <- matrix(
  data = c(-100, -100, 0, 90),
  byrow = F, nrow = 2) %>%
  st_linestring() %>%
  st_sfc(crs = 4326) %>%
  st_transform(crs = myCRS)
mybreaks <- seq(-200,200,by=10)

map <- ggplot() +
  tidyterra::geom_spatraster_contour_filled(iso_feather, mapping = aes(), breaks = mybreaks) +
  scale_fill_viridis_d("Feather isoscape values", option = "plasma") +
  geom_sf(myline, mapping = aes(), color = "black", linetype = "dotted", size = 2) +
  coord_sf(
    xlim = c(extent(range_raster)[1]+1.5e6, extent(range_raster)[2]-2.4e6),
    ylim = c(extent(range_raster)[3]+5e6, extent(range_raster)[4]-1.8e6)
  ) +
  ggtitle("Wilson's snipe breeding habitat feather isoscape")


# Extract percentages -----------------------------------------------------

precip_df <- iso_cropped %>%
  terra::project("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") %>%
  as.data.frame(xy=T) %>%
  dplyr::mutate(
    under100 = factor(case_when(x < -100 ~ 1, TRUE ~0)),
    lat_bin = cut(y, breaks = mybreaks),
    precip_bin = cut(d2h_GS, breaks = seq(-200,200,by=20))
  )

feather_df <- iso_feather %>%
  terra::project("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs") %>%
  as.data.frame(xy=T) %>%
  dplyr::mutate(
    under100 = factor(case_when(x < -100 ~ 1, TRUE ~0)),
    lat_bin = cut(y, breaks = mybreaks)
    )

p_lon <- feather_df %>%
  ggplot() +
  aes(y = d2h_GS, x = under100, group = under100) +
  geom_jitter(aes(color = d2h_GS), height = 0, width = 0.3) +
  scale_color_viridis_c("Feather isoscape values", option = "plasma", limits = c(-200,-10)) +
  geom_violin(fill = NA) +
  geom_boxplot(fill = NA) +
  scale_x_discrete("Longitude", labels = c(" < -100", " >= -100"))


p_lat <- feather_df %>%
  ggplot() +
  aes(y = d2h_GS, x = lat_bin, group = lat_bin) +
  geom_jitter(aes(color = d2h_GS), height = 0, width = 0.3) +
  scale_color_viridis_c("Feather isoscape values", option = "plasma", limits = c(-200,-10)) +
  geom_violin(fill = NA) +
  geom_boxplot(fill = NA) +
  scale_x_discrete("Latitude")

out <-ggarrange(map, ggarrange(p_lon, p_lat, ncol = 1, common.legend = T), ncol = 2 , widths = c(2,1))
ggsave(out, filename = file.path("figs", "featherIsoscape_details.png"))
