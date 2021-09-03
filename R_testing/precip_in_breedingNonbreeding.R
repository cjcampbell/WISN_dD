# what precip values are present in the breeding range?

load(file.path(wd$bin, "my_isoscapes.RData"), verbose = TRUE)
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

isoscape_breeding <- mask(myisoscape$isoscape, breeding)
isoscape_breeding[] %>% hist

isoscape_breeding_feathers <- (isoscape_breeding*0.95)-23
isoscape_breeding_feathers[] %>% hist
mydata_transformed$d2H %>% hist


source("~/WISN_dD/R/03_LoadApplyTransferFunctions.R")
mydata_transformed %>%
  dplyr::mutate(
    d2H_1 = (d2H ) / 0.95,
    d2H_Sullins_grow = (d2H - 23.57) / 1.16,
    d2H_Sullins_ann = (d2H - 23.51) / 1.06,
    d2H_2 = (d2H + 20),
    d2H_3 = (d2H + 20) / 0.95,
    d2H_4 = (d2H + 20) / 1.36,
    d2H_5 = (d2H - 5) / 1.9,
    d2H_6 = (d2H + 74) / 0.57
  ) %>%
  dplyr::select_at(vars(starts_with("d2H"))) %>%
  pivot_longer(everything()) %>%
  rbind(., data.frame(name = "isoscape", value = isoscape_breeding[])) %>%
  ggplot() +
  geom_histogram(aes(value)) +
  facet_wrap(~name, ncol = 1, scales = "free_y") +
  geom_vline(xintercept = -50, color = "red") +
  theme_minimal()
