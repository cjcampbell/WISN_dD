
# Exploratory direction analyses ------------------------------------------


allResults <- fread(file.path(wd$bin, "allResults.csv")) %>%
  dplyr::mutate(
    westOrigin = case_when(
      probWestOrigin >= 0.8 ~ 1,
      probWestOrigin <= 0.2 ~ 0,
      TRUE ~ as.numeric(NA)
      )
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
