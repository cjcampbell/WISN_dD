
library(dplyr)
library(sf)
library(ggpubr)

if(!exists("maps_df")) maps_df <- readRDS(file.path(wd$bin, "maps_df.rds"))

# This is specific to my local computer.
source("~/WISN_dD/.Rprofile")
NoAm <- readRDS( file.path(bigDataStorage, "NoAm_maps", "NoAm.rds"))

# Make individual plots.
df <- maps_df %>% dplyr::filter(method == "OR")

plotID <- function(ID) {
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
    coord_sf() +
    theme_bw() +
    theme(
      strip.background = element_blank(),
      axis.title = element_blank()
    ) +
    ggtitle(i)
  ggsave(p, filename = file.path(wd$figs, paste0("map", i, ".png")))
}

lapply(unique(df$ID), plotID )



