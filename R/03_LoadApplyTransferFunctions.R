# Load  ----------------------------------------------------------

if(!exists("mydata")) source(file.path(wd$R, "01_loadIsotopedata.R"))

# Apply transfer functions, define molt status  ---------------------------------

# Sullins et al 2016 woodcock, growing season

mydata_transformed <- mydata %>%
  dplyr::mutate_if(is.factor, as.character) %>%
  dplyr::mutate(
    dDprecip =  (d2H - 23.57) / 1.16,
    sdResid = 12.6
  )

saveRDS(mydata_transformed, file = file.path(wd$bin, "mydata_transformed.rds"))

