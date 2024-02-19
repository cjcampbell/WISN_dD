library(data.table)
library(sf)
library(ggstatsplot)

# Exploratory direction analyses ------------------------------------------

allResults <- fread(file.path(wd$bin, "allResults.csv")) %>%
  dplyr::mutate(
    westOrigin = case_when(
      probWestOrigin >= 0.8 ~ 1,
      probWestOrigin <= 0.2 ~ 0,
      TRUE ~ as.numeric(NA)
      ),
    yday = as.numeric(yday(Date)),
    yday2 = case_when(yday <= 200 ~ (yday + 200), TRUE ~ yday),
    month = as.numeric(month(Date))
  )

## Do definitive east/west origins have a clear relationship w sampling location?
library(cowplot)
library(sf)

p_count <- allResults %>%
  dplyr::mutate(westOrigin = as.character(westOrigin)) %>%
  ggplot() +
  aes(x=westOrigin) +
  geom_bar(aes(fill = westOrigin)) +
  geom_text(stat='count', aes(label=..count..), vjust=-.05) +
  scale_y_continuous(NULL, limits = c(0,225), expand = c(0,0)) +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    legend.position = "none"
    ) +
  scale_x_discrete( "Region of origin", breaks = c("0", "1", NA), labels =  c("East\n(>- 100°W)", "West\n(<= -100°W)", "Region\nunspecified")) +
  scale_fill_manual("Region of origin", breaks = c("0", "1", NA), values = c("#D11149", "#AE7709", "black"), na.value = "#29335C")

states <- rnaturalearth::ne_states(country = "United States of America") %>%
  st_as_sf() %>%
  sf::st_transform(myCRS) %>%
  dplyr::filter(name == "Florida")

p_map <- allResults %>%
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
  scale_size_continuous(range = c(1,8)) +
  theme_minimal() +
  theme(axis.title = element_blank())+
  coord_sf()

p_together <- ggdraw(p_map) +
  draw_plot(p_count, x = 0.15, y = 0.1, width = 0.3, height = 0.4)
ggsave(p_together, file = file.path(wd$figs, "originMap.png"))


# Summarize distances traveled ---------
allResults %>%
  dplyr::mutate(westOrigin = as.character(westOrigin)) %>%
  ggplot() +
  geom_histogram(aes(dist_km, group = westOrigin, fill = westOrigin), binwidth = 250) +
  scale_x_continuous("Minimum distance traveled (km)", breaks = seq(0,4000,by = 500), limits = c(0,5000)) +
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
    legend.position = c(0.1,0.8),
    axis.ticks = element_line(),
    panel.grid = element_blank()
    )
ggsave(file.path(wd$figs, "histogramMinDist.png"))

set.seed(42)
allResults %>%
  dplyr::mutate(westOrigin = as.character(westOrigin)) %>%
  ggplot() +
    aes(x=lat, y = dist_km) +
    geom_point(aes(color= westOrigin), position = position_jitter(w = 0.1, h = 10)) +
    stat_smooth(method = "lm", color = "black") +
  scale_y_continuous("Minimum distance traveled (km)") +
  scale_x_continuous("Latitude of sample site") +
  scale_color_manual(
    "Region of origin",
    values = c("#D11149", "#AE7709", "black"),
    breaks = c("0","1", NA),
    labels = c("East (> -100°W)", "West (<= -100°W)", "Region unspecified"),
    na.value = "#29335C"
  ) +
  theme_minimal()
ggsave(file.path(wd$figs, "lat-v-MinDist.png"))


allResults %>%
  ggstatsplot::ggbetweenstats(y = dist_km, x = westOrigin)
mean(allResults$dist_km, na.rm = T)
sd(allResults$dist_km, na.rm = T)
summary(allResults$dist_km, na.rm = T)

# boxplots ----------------------------------------------------------------

mydata_transformed <- readRDS( file.path(wd$bin, "mydata_transformed.rds") )
load(file.path(wd$bin, "mydata_clustered.Rdata"))
df <- left_join(mydata_transformed, mydata_clustered) %>%
  left_join(allResults) %>%
  dplyr::filter(!is.na(d2H)) %>%
  dplyr::mutate(
    year = lubridate::year(Date),
    winter = case_when(Date <= as.Date("2021-06-15") ~ "2020-1", Date > as.Date("2021-06-15") ~ "2021-2")
    )

ggstatsplot::ggbetweenstats(
  data = df, y = d2H, x = winter,
  ggplot.component	= list(
    ggplot2::theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ),
    ggplot2::ylab(expression(bold(paste(delta^2 ~ H[feather], " (", "\u2030", ", VSMOW)") ))),
    ggplot2::xlab("Collection year")
    )
)


ggstatsplot::ggbetweenstats(
  data = df, y = d2H, x = OriginCluster,
  ggplot.component	= list(
    ggplot2::theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ),
    ggplot2::ylab(expression(bold(paste(delta^2 ~ H[feather], " (", "\u2030", ", VSMOW)") )))
  )
  )

# ggplot(df) +
#   geom_boxplot(aes(y=d2H, x = OriginCluster, group = OriginCluster))


# Dist by cluster ---------------------------------------------------------

set.seed(42); allResults %>%
  dplyr::mutate(
    westOrigin = case_when(
    probWestOrigin >= 0.8 ~ "W",
    probWestOrigin <= 0.2 ~ "E",
    TRUE ~ "Unknown"
  )) %>%
  ggplot() +
  geom_point(aes(y = dist_km, color = westOrigin, x = OriginCluster), position = position_jitter(width = 0.33)) +
  scale_y_continuous("Minimum distance traveled (km)", breaks = seq(0,4000,by = 500), limits = c(0,5000)) +
  scale_color_manual(
    "Region of origin",
    values = c("#D11149", "#AE7709", "black"),
    breaks = c("E","W", "Unknown"),
    labels = c("East (> -100°W)", "West (<= -100°W)", "Region unspecified"),
    na.value = "#29335C"
  ) +
  theme_minimal() +
  theme(
    legend.position = c(0.1,0.8),
    axis.ticks = element_line(),
    panel.grid = element_blank()
  )

ggstatsplot::ggbetweenstats(
  data = allResults, y = dist_km, x = OriginCluster,
  ggplot.component	= list(
    ggplot2::theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank()
    ),
    ggplot2::ylab(expression(bold(paste(delta^2 ~ H[feather], " (", "\u2030", ", VSMOW)") )))
  )
)

allResults %>%
  ggplot() +
  geom_histogram(aes(dist_km, group = factor(OriginCluster), fill = factor(OriginCluster)), binwidth = 250) +
  scale_x_continuous("Minimum distance traveled (km)", breaks = seq(0,4000,by = 500), limits = c(0,5000)) +
  scale_y_continuous(expand = c(0,1)) +
  scale_fill_viridis_d("OriginCluster") +
  theme_minimal() +
  theme(
    legend.position = c(0.1,0.8),
    axis.ticks = element_line(),
    panel.grid = element_blank()
  )




# distance over time --------------------------------------------------------------------

allResults %>%
  ggplot() +
  geom_point(aes(x=Date, y = dist_km))

  geom_histogram(aes(dist_km, group = factor(OriginCluster), fill = factor(OriginCluster)), binwidth = 250) +
  scale_x_continuous("Minimum distance traveled (km)", breaks = seq(0,4000,by = 500), limits = c(0,5000)) +
  scale_y_continuous(expand = c(0,1)) +
  scale_fill_viridis_d("OriginCluster") +
  theme_minimal() +
  theme(
    legend.position = c(0.1,0.8),
    axis.ticks = element_line(),
    panel.grid = element_blank()
  )


# Clusters by year --------------------------------------------------------

allResults %>%
    dplyr::mutate(
      year = lubridate::year(Date),
      winter = case_when(Date <= as.Date("2021-06-15") ~ "2020-1", Date > as.Date("2021-06-15") ~ "2021-2")
    ) %>%
    ggplot() +
    geom_bar(aes( fill= OriginCluster,x= winter, group = OriginCluster), position = "dodge")

allResults %>%
  dplyr::mutate(
    year = lubridate::year(Date),
    winter = case_when(Date <= as.Date("2021-06-15") ~ "2020-1", Date > as.Date("2021-06-15") ~ "2021-2")
  ) %>%
  dplyr::group_by(winter, OriginCluster) %>%
  dplyr::summarise(n=n()) %>%
  pivot_wider(names_from = winter, values_from = n) %>%
  dplyr::select(-OriginCluster) %>%
  chisq.test()


# d2H v sample Lat --------------------------------------------------------

# d2H v lat
m1 <- lm(data = allResults, d2H ~ lat)
summary(m1)
confint(m1, "lat")
library(sjPlot)
plot_model(m1, type = "pred")

ggstatsplot::ggcoefstats(m1)

# distance v. latitude.
# This is expected (~5% of variation in dist_km is explained by sampling latitude)
# and frankly, circular (since dist_km is calculated with respect to sampling location).
# So, good to know, but probably not worth reporting.
m2 <- lm(data = allResults, dist_km ~ lat)
summary(m2)
plot_model(m2, type = "pred")


ggplot(allResults) +
  aes(x=lat, y= dist_km) +
  geom_point() +
  stat_smooth(method = "lm")



# d2H v month of sampling -------------------------------------------------

m3 <- lm(data = allResults, d2H ~ yday2)
summary(m3)
confint(m3, "yday2")
plot_model(m3, type = "pred")
ggstatsplot::ggcoefstats(m3)

# Quick check of month
ggbetweenstats(data = allResults, y = d2H, x = month)





