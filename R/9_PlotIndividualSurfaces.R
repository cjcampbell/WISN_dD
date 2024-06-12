
library(dplyr)
library(sf)
library(ggpubr)
library(ggstar)

if(!exists("maps_df")) {
  maps_df <- list.files(wd$tmp_df, pattern = "df_list.*rds$", full.names = T) %>%
    lapply(readRDS) %>%
    bind_rows()
  }

# This is specific to my local computer.
source("~/WISN_dD/.Rprofile")
NoAm <- readRDS( file.path(wd$bin, "NoAm.rds")) %>%
  st_transform(myCRS)
countries <- rnaturalearth::countries110 %>%
  st_as_sf() %>%
  st_transform(myCRS)

range_raster <- readRDS( file.path(wd$bin, "range_raster.rds" ) )

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

breeding <- IUCNmaps %>% filter(SEASONAL %in% c(1,2))
other <- IUCNmaps %>% filter(SEASONAL %in% c(3,4))

mydata_transformed <- readRDS( file.path(wd$bin, "mydata_transformed.rds") )
pts <- as.data.frame(sf::sf_project(pts = as.matrix(mydata_transformed[,c("lon", "lat")]), from = "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs ", to = myCRS))
names(pts) <- c("lon_aea", "lat_aea")
mydata_transformed <- data.frame(mydata_transformed, pts )


# Make individual plots.
df_OR <- maps_df %>% dplyr::filter(method == "OR")
df_cs <- maps_df %>% dplyr::filter(method == "raw") %>%
  dplyr::group_by(ID) %>%
  arrange(value) %>%
  dplyr::mutate(cumsum = cumsum(value))


plotID <- function(i, save = TRUE) {
  p <- dplyr::filter(df_OR, ID == i) %>%
    ggplot() +
    # geom_sf(
    #   countries, mapping=aes(),
    #   fill = "white", color = "grey20", size = 0.25) +
    geom_sf(
      NoAm, mapping=aes(),
      fill = "grey90", color = "grey20", size = 0.25) +
    geom_tile(mapping=aes(x=x,y=y,fill=value, color = value)) +
    scale_fill_viridis_c(
      "Odds of origin",
      option = "turbo", direction = 1
    ) +
    scale_color_viridis_c(
      "Odds of origin",
      option = "turbo", direction = 1
    ) +
    geom_sf(breeding, mapping = aes(),
            fill = NA, color = "black", size = 1) +
    geom_sf(other, mapping = aes(),
            fill = "grey40", alpha = 0.5, color = NA, size = 1) +
    geom_star(
      data = dplyr::filter(mydata_transformed, SampleName == i),
      mapping = aes(x=lon_aea, y = lat_aea),
      fill = "#E9D539FF", color = "black", size = 3
      ) +
    coord_sf(
      xlim = c(extent(range_raster)[1]+1e6, extent(range_raster)[2]-2.4e6),
      ylim = c(extent(range_raster)[3]+3.5e6, extent(range_raster)[4]-1.8e6)
    ) +
    theme_bw() +
    theme(
      legend.position = c(0.2,0.2),
      strip.background = element_blank(),
      axis.title = element_blank()
    ) +
    ggtitle(i)

  if(save == TRUE){
    ggsave(p, filename = file.path(wd$figs, paste0("map_", i, ".png")))
  } else {
    plot(p)
  }

}

# Make plots for all.
lapply(SamplesToRerun, plotID , save = T)

mydata_transformed %>%
  filter(dDprecip > -50) %>%
  arrange(desc(dDprecip)) %>%
  dplyr::select(SampleName) %>%
  unlist %>%
  lapply(., plotID , save = F)



# Make aggregate surfaces -------------------------------------------------
load(file.path(wd$bin, "mydata_clustered.Rdata"))

mdf <- full_join(mydata_transformed, mydata_clustered)

aggSurfaces <- lapply(unique(mdf$OriginCluster), function(cl){

  indivs <- filter(mdf, OriginCluster == cl) %>%
    dplyr::select(SampleName) %>% unlist

  cells <- filter(df_OR, ID %in% indivs) %>%
    dplyr::mutate(over = case_when(value >= 2/3 ~ 1, TRUE ~ 0))

  aggSurgface <- cells %>%
    group_by(x,y) %>%
    dplyr::summarise(n = sum(over)) %>%
    ungroup %>%
    mutate(OriginCluster = cl)

  return(aggSurgface)

}) %>%
  bind_rows()

fwrite(aggSurfaces, file = file.path(wd$bin, "aggSurfaces.csv"), row.names = F)

p_aggSurfaces <- aggSurfaces %>%
  ggplot() +
  geom_sf(
    NoAm, mapping=aes(),
    fill = "white", color = "grey20", size = 0.25) +
  geom_tile(mapping=aes(x=x,y=y,fill=log10(n), color = log10(n))) +
  scale_color_gradient(high = "#41B7C4", low = "#FEFFD9") +
  scale_fill_gradient( high = "#41B7C4", low = "#FEFFD9") +
  geom_sf(breeding, mapping = aes(),
          fill = NA, color = "black", size = 1) +
  geom_sf(other, mapping = aes(),
          fill = "grey50", alpha = 0.5, color = NA, size = 1) +
  geom_text(
    data = mdf %>% group_by(OriginCluster) %>% dplyr::summarise(n=n()),
    aes(label = paste("n = ", n)),
    x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.1
  ) +
  coord_sf(
    xlim = c(extent(range_raster)[1], extent(range_raster)[2]),
    ylim = c(extent(range_raster)[3], extent(range_raster)[4])
  ) +
  theme_bw() +
  theme(
    strip.background = element_blank(),
    axis.title = element_blank()
  ) +
  facet_wrap(~OriginCluster)

ggsave(plot = p_aggSurfaces, filename = file.path(wd$figs, "aggSurfaces.png"))




# Make min dist plot ------------------------------------------------------
wd$tmp_df <- file.path(wd$bin,"tmp_df")
df <- list.files(wd$tmp_df, full.names = T)[100] %>%
  readRDS() %>%
  dplyr::filter(method == "OR")

allResults <- data.table::fread(file.path(wd$bin, "allResults.csv"))
md <- allResults %>%
  dplyr::filter(SampleName == df$ID[1])


p <- df %>%  ggplot() +
  # geom_sf(
  #   countries, mapping=aes(),
  #   fill = "white", color = "grey20", size = 0.25) +
  geom_sf(
    NoAm, mapping=aes(),
    fill = "grey90", color = "grey20", size = 0.25) +
  geom_tile(mapping=aes(x=x,y=y,fill=value, color = value)) +
  scale_fill_viridis_c(
    "Odds of origin",
    option = "turbo", direction = 1
  ) +
  scale_color_viridis_c(
    "Odds of origin",
    option = "turbo", direction = 1
  ) +
  geom_sf(breeding, mapping = aes(),
          fill = NA, color = "black", size = 1) +
  geom_sf(other, mapping = aes(),
          fill = "grey40", alpha = 0.5, color = NA, size = 1) +
  geom_segment(
    data = md,
    mapping = aes(x = x, y= y, xend = lon_aea, yend = lat_aea),
    arrow = arrow(length = unit(0.5, "cm"))
  ) +
  geom_star(
    data = md,
    mapping = aes(x=lon_aea, y = lat_aea),
    fill = "#E9D539FF", color = "black", size = 3
  ) +
  coord_sf(
    xlim = c(extent(range_raster)[1]+1e6, extent(range_raster)[2]-2.4e6),
    ylim = c(extent(range_raster)[3]+3.5e6, extent(range_raster)[4]-1.8e6)
  ) +
  theme_bw() +
  theme(
    legend.position = c(0.2,0.2),
    strip.background = element_blank(),
    axis.title = element_blank()
  ) +
  ggtitle(df$ID[1])
ggsave(p, file = file.path(wd$figs, paste0("map_line_", df$ID[1], ".png")))


