
# Setup -------------------------------------------------------------------
theme_set(
  theme(plot.margin=unit(c(0,0,0,0), "mm"))
)


load( file.path(wd$bin, "mydata_clustered.Rdata"), verbose = T )
namesByGroup <- mydata_clustered %>%
  split(.$OriginCluster) %>%
  map(~dplyr::select(., SampleName) %>% unlist)

thesefiles <- list.files(file.path(wd$bin, "maps"), pattern = "df_list_", full.names = T)

csumThreshold <- 0.25


if(!file.exists(file.path(wd$bin, "clusterAggSurfaces.csv"))) {

  clusterAggSurfaces <- lapply(1:length(namesByGroup), function(i) {

    y <- namesByGroup[[i]]

    allOrigByClust <- lapply(y, function(x) {
      grep(thesefiles, pattern = x, value = T) %>%
        fread() %>%
        dplyr::filter(method == "raw") %>%
        arrange(value) %>%
        mutate(
          csum = cumsum(value),
          orig = case_when(csum >= csumThreshold ~ 1, TRUE ~ 0)
        )
    }) %>%
      bind_rows()

    aggSurface <- allOrigByClust %>%
      group_by(x,y) %>%
      dplyr::summarise(sum = sum(orig)) %>%
      dplyr::mutate(Cluster = i)
    return(aggSurface)

  }) %>%
    bind_rows()

  fwrite(clusterAggSurfaces, file = file.path(wd$bin, "clusterAggSurfaces.csv"), row.names = F)

}

if(!exists("clusterAggSurfaces")) {
  clusterAggSurfaces <- fread(file.path(wd$bin, "clusterAggSurfaces.csv"))
  }


# Inset plot stuff --------------------------------------------------------

# Specify bounds of inset plot
inset_xlim <- c(0.95e6, 1.55e6)
inset_ylim <- c(-1.65e6, -1.0e6)

states <- rnaturalearthhires::states10 %>%
  st_as_sf() %>%
  st_transform(myCRS) %>%
  st_crop(my_extent_aea)

# load spatial stuff ------------------------------------------------------

# Load world map
wrld <- rnaturalearthhires::countries10 %>%
  st_as_sf() %>%
  st_transform(crs = myCRS) %>%
  st_crop(my_extent_aea)

# plot cluster maps --------------------------------------------------------------------

numByClust <- mydata_clustered %>%
  dplyr::group_by(OriginCluster) %>%
  dplyr::summarise(n_indivs = n()) %>%
  dplyr::rename(Cluster = OriginCluster)


myClustMaps <- lapply(1:4, function(i) {

  out <- clusterAggSurfaces %>%
    dplyr::filter(Cluster == i) %>%
    left_join(numByClust) %>%
    dplyr::mutate(
      prop = sum/n_indivs,
      prop = case_when(prop == 0 ~ as.numeric(NA), TRUE ~ prop)) %>%
    ggplot() +
    geom_sf(
      wrld, mapping=aes(),
      fill = "white", color = "grey20", size = 0.25) +
    geom_tile(
      mapping=aes(
        x=x,y=y,
        #fill=log10(sum), color = log10(sum)
        fill=prop, color = prop
        )
      ) +
    # scale_color_gradient(high = "#41B7C4", low = "#FEFFD9", limits = c(0,1)) +
    # scale_fill_gradient( high = "#41B7C4", low = "#FEFFD9", limits = c(0,1)) +
    scale_color_gradient(NULL, high = "#35A0AC", low = "#FEFFD9", limits = c(0,1)) +
    scale_fill_gradient( NULL, high = "#35A0AC", low = "#FEFFD9", limits = c(0,1)) +

    geom_sf(breeding, mapping = aes(),
            fill = NA, color = "black", size = 1) +
    geom_sf(other, mapping = aes(),
            fill = "grey70", alpha = 0.5, color = NA, size = 1) +
    # geom_text(
    #   data = dplyr::filter(numByClust, Cluster == i),
    #   aes(label = paste("n = ", n_indivs)),
    #   x = -Inf, y = -Inf, hjust = -0.1, vjust = -0.1
    # ) +
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
      axis.line = element_line(color = "grey50"),
      strip.background = element_blank(),
      axis.title = element_blank()
    ) +
    # Draw line around inset.
    geom_rect(
      color = "black",
      linewidth = 0.7,
      fill = NA,
      mapping = aes(
        xmin = inset_xlim[[1]] - 1e5,
        xmax = inset_xlim[[2]] + 1e5,
        ymin = inset_ylim[[1]] - 1e5,
        ymax = inset_ylim[[2]] + 1e5
      )
    )

  return(out)

})

## Get legend. -------------

myClustMaps[[1]] %>%
  get_legend(position = "bottom") %>%
  ggsave(filename = file.path(wd$figs, "legend_aggMap.png"), width = 3, height = 1)

# Create distance traveled histograms -----------------------------------

allResults <- fread(file.path(wd$bin, "allResults.csv")) %>%
  dplyr::mutate(
    westOrigin = case_when(
      probWestOrigin >= 0.8 ~ 1,
      probWestOrigin <= 0.2 ~ 0,
      TRUE ~ as.numeric(NA)
    )
  )


distByDirAndCluster <- lapply(1:4, function(i) {
  allResults %>%
    dplyr::filter(OriginCluster == i) %>%
    dplyr::mutate(westOrigin = as.character(westOrigin)) %>%
    ggplot() +
    geom_histogram(aes(dist_km, group = westOrigin, fill = westOrigin), binwidth = 500) +
    scale_x_continuous(
      NULL,
      breaks = seq(0,5000,by = 1000)-500/2,
      labels = c(0, paste0(seq(1,5,by=1), "k")),
      limits = c(0,5000)
      ) +
    scale_y_continuous(NULL, expand = expansion(add = c(0,0)), limits = c(0,60)) +
    scale_fill_manual(
      "Region of origin",
      values = c("#D11149", "#AE7709", "black"),
      breaks = c("0","1", NA),
      labels = c("East (> -100°W)", "West (<= -100°W)", "Region unspecified"),
      na.value = "#29335C"
    ) +
    theme_minimal() +
    theme(
      axis.line = element_line(color = "grey50"),
      legend.position = c(0.1,0.8),
      axis.ticks = element_line(),
      panel.grid = element_blank()
    )
})

commonDistDirLegend <- {
    allResults %>%
    dplyr::mutate(westOrigin = as.character(westOrigin)) %>%
    ggplot() +
    geom_histogram(aes(dist_km, group = westOrigin, fill = westOrigin), binwidth = 250) +
    scale_x_continuous("Minimum distance traveled (km)", breaks = seq(0,4000,by = 500), limits = c(0,4300)) +
    scale_y_continuous(expand = c(0,1)) +
    scale_fill_manual(
      "Region of origin",
      values = c("#D11149", "#AE7709", "black"),
      breaks = c("0","1", NA),
      labels = c("East (> -100°W)", "West (<= -100°W)", "Region unspecified"),
      na.value = "#29335C"
    ) +
    theme_minimal() +
    theme(
      plot.margin = unit(c(0,0,0,0), "cm"),
      legend.position = c(0.5,0.5),
      axis.ticks = element_line(),
      panel.grid = element_blank()
    )
  } %>%
  ggpubr::get_legend()



# Create maps of sampled regions ----------------------------------------

sampledRegionMaps <- lapply(1:4, function(i) {

  allResults %>%
    dplyr::filter(OriginCluster == i) %>%
    dplyr::mutate(westOrigin = as.character(westOrigin)) %>%
    dplyr::mutate(westOrigin = as.character(westOrigin)) %>%
    group_by(westOrigin, lon_aea, lat_aea) %>%
    dplyr::summarise(n=n()) %>%
    ggplot() +
    geom_sf(states, mapping = aes()) +
    geom_point(aes(size = n, x=lon_aea, y=lat_aea, color = westOrigin), shape = 21) +
    scale_color_manual(
      "Region of origin",
      values = c("#D11149", "#AE7709", "black"),
      breaks = c("0","1", NA),
      labels = c("East (> -100°W)", "West (<= -100°W)", "Region unspecified"),
      na.value = "#29335C"
    ) +
    scale_x_continuous(
      breaks = seq(-200,200,by=20), expand = c(0,0)
    ) +
    scale_y_continuous(
      breaks = seq(0,90,by=10), expand = c(0,0)
    ) +
    scale_size_continuous(range = c(1,8), limits = c(1,40), breaks = seq(0,100,by=10)) +
    coord_sf(xlim = inset_xlim, ylim = inset_ylim, clip = "on") +
    theme_bw() +
    theme(
      axis.text = element_blank(),
      axis.title = element_blank(),
      legend.position = "none",
      plot.background = element_blank(),
      panel.background = element_blank()
    ) +
    # Draw line around inset.
    geom_rect(
      color = "black",
      linewidth = 0.75,
      fill = NA,
      mapping = aes(
        xmin = inset_xlim[[1]],
        xmax = inset_xlim[[2]],
        ymin = inset_ylim[[1]],
        ymax = inset_ylim[[2]]
      )
    )
})

sampleLegend <- { allResults %>%
    dplyr::mutate(westOrigin = as.character(westOrigin)) %>%
    dplyr::mutate(westOrigin = as.character(westOrigin)) %>%
    group_by(westOrigin, lon_aea, lat_aea) %>%
    dplyr::summarise(n=n()) %>%
    ggplot() +
    geom_sf(states, mapping = aes()) +
    geom_point(aes(size = n, x=lon_aea, y=lat_aea, color = westOrigin), shape = 21) +
    scale_color_manual(
      "Region of origin",
      values = c("#D11149", "#AE7709", "black"),
      breaks = c("0","1", NA),
      labels = c("East (> -100°W)", "West (<= -100°W)", "Region unspecified"),
      na.value = "#29335C"
    ) +
    scale_x_continuous(
      breaks = seq(-200,200,by=20)
    ) +
    scale_y_continuous(
      breaks = seq(0,90,by=10)
    ) +
    scale_size_continuous(range = c(1,8), limits = c(1,40), breaks = c(1,10,25)) +
    theme_minimal() +
    theme(axis.title = element_blank())+
    coord_sf(xlim = inset_xlim, ylim = inset_ylim)
  }

get_legend(sampleLegend, position = "bottom") %>%
  ggsave(filename = file.path(wd$figs, "legend-sampledregions.png"), height = 3, width = 10)


# Create individual maps of combined stuff -----------------------------------------

myClustMaps2 <- appendCenteredLegend(myClustMaps)

ggarrange(plotlist = myClustMaps2, ncol = 1) %>%
  ggsave(filename = file.path(wd$figs, "clustermaps.png"), height = 8, width = 6)


appendCenteredLegend(sampledRegionMaps) %>%
  ggarrange(plotlist = c(.), ncol = 1) %>%
  ggsave(filename = file.path(wd$figs, "sampledRegionsByCluster.png"), height = 8, width = 6)



ggsave(sampleLegend, filename = file.path(wd$figs, "sampledRegionsColors.png"))


distByDirAndCluster2 <- appendCenteredLegend(distByDirAndCluster)


ggarrange(plotlist = distByDirAndCluster2[1:4], ncol = 1) %>%
  ggarrange(., commonDistDirLegend, ncol = 1, heights = c(1,0.4)) %>%
  ggsave(filename = file.path(wd$figs, "distancesByCluster.png"), height = 8, width = 6)



# Create cluster-specific rows --------------------------------------------

library(patchwork)
p1 <-
  {myClustMaps[[1]]         + theme(legend.position = "none") } +
  {sampledRegionMaps[[1]]   + theme(legend.position = "none") + ylab(NULL)} +
  {distByDirAndCluster[[1]] + theme(legend.position = "none") + xlab(NULL) } +

  {myClustMaps[[2]]         + theme(legend.position = "none") } +
  {sampledRegionMaps[[2]]   + theme(legend.position = "none") + ylab(NULL)} +
  {distByDirAndCluster[[2]] + theme(legend.position = "none") + xlab(NULL) } +

  {myClustMaps[[3]]         + theme(legend.position = "none") } +
  {sampledRegionMaps[[3]]   + theme(legend.position = "none") + ylab(NULL)} +
  {distByDirAndCluster[[3]] + theme(legend.position = "none") + xlab(NULL) } +

  {myClustMaps[[4]]         + theme(legend.position = "none") } +
  {sampledRegionMaps[[4]]   + theme(legend.position = "none") + ylab(NULL)} +
  {distByDirAndCluster[[4]] + theme(legend.position = "none") + xlab(NULL) } +


  plot_layout(
    widths = c(1,1,1,0.7),
    design =c(

      area(t=1,b=1,l=1,r=1),
      area(t=1,b=1,l=2,r=2),
      area(t=1,b=1,l=3,r=4),

      area(t=2,b=2,l=1,r=1),
      area(t=2,b=2,l=2,r=2),
      area(t=2,b=2,l=3,r=4),

      area(t=3,b=3,l=1,r=1),
      area(t=3,b=3,l=2,r=2),
      area(t=3,b=3,l=3,r=4),

      area(t=4,b=4,l=1,r=1),
      area(t=4,b=4,l=2,r=2),
      area(t=4,b=4,l=3,r=4)
  )
  )

ggsave(p1, filename = file.path(wd$figs, "alignmentTest.png"), width = 7, height = 7)
