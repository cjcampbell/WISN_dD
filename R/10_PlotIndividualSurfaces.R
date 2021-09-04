
library(dplyr)
library(sf)
library(ggpubr)

if(!exists("maps_df")) maps_df <- readRDS(file.path(wd$bin, "maps_df.rds"))

# This is specific to my local computer.
source("~/WISN_dD/.Rprofile")
NoAm <- readRDS( file.path(bigDataStorage, "NoAm_maps", "NoAm.rds"))




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


# Make individual plots.
df <- maps_df %>% dplyr::filter(method == "OR")

plotID <- function(i, save = TRUE) {
  p <- dplyr::filter(df, ID == i) %>%
    ggplot() +
    geom_sf(
      NoAm, mapping=aes(),
      fill = "white", color = "grey20", size = 0.25) +
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
            fill = "grey50", alpha = 0.5, color = NA, size = 1) +
    coord_sf(
      xlim = c(extent(range_raster)[1], extent(range_raster)[2]),
      ylim = c(extent(range_raster)[3], extent(range_raster)[4])
    ) +
    theme_bw() +
    theme(
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
lapply(unique(df$ID), plotID , save = T)

mydata_transformed %>%
  filter(dDprecip > -50) %>%
  arrange(desc(dDprecip)) %>%
  dplyr::select(SampleName) %>%
  unlist %>%
  lapply(., plotID , save = F)



# Make aggregate surfaces -------------------------------------------------

mdf <- full_join(mydata_transformed, mydata_clustered)

aggSurfaces <- lapply(unique(mdf$OriginCluster), function(cl){

  indivs <- filter(mdf, OriginCluster == cl) %>%
    dplyr::select(SampleName) %>% unlist

  cells <- filter(df, ID %in% indivs) %>%
    dplyr::mutate(over = case_when(value >= 2/3 ~ 1, TRUE ~ 0))

  aggSurgface <- cells %>%
    group_by(x,y) %>%
    dplyr::summarise(n = sum(over)) %>%
    ungroup %>%
    mutate(OriginCluster = cl)

  return(aggSurgface)

}) %>%
  bind_rows()


aggSurfaces %>%
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
