
# Setup -------------------------------------------------------------------

mydata_transformed <- readRDS( file.path(wd$bin, "mydata_transformed.rds") )

load( file.path(wd$bin, "my_isoscapes.RData") )

mapPath <- file.path( wd$bin, "maps")
if(!dir.exists(mapPath)) dir.create(mapPath)

# Create maps. -----------------------------------------------------------

isocat::isotopeAssignmentModel(
  ID               = mydata_transformed$SampleName,
  isotopeValue     = mydata_transformed$dDprecip,
  SD_indv          = mydata_transformed$sdResid,
  precip_raster    = myisoscape$isoscape,
  precip_SD_raster = myisoscape$sd,
  savePath         = mapPath,
  additionalModels = range_raster,
  nClusters = FALSE
)


# Convert to other formats (for future analyses) --------------------------

# Load maps.
maps_cropped <- list.files(
  mapPath, pattern = "Combined.*.grd$", full.names = TRUE, recursive = T) %>%
  raster::stack()

# Also calculate probability quantiles.
maps_quantile_stack <-
  lapply(1:nlayers(maps_cropped), function(i) {
    print(i)
    isocat::makeQuantileSurfaces(maps_cropped[[i]])
    }) %>%
  stack()
names(maps_quantile_stack) <- paste0(names(maps_cropped), "_quantile")
writeRaster(maps_quantile_stack, filename = file.path(mapPath, "quantileProbabilityMaps.grd"), overwrite = TRUE)

# And odds ratios.
maps_odds_stack <-
  lapply(1:nlayers(maps_cropped), function(i){
    isocat::makeOddsSurfaces(maps_cropped[[i]])
    }) %>%
  stack()
names(maps_odds_stack) <- paste0(names(maps_cropped), "_OR")
writeRaster(maps_odds_stack, filename = file.path(mapPath, "ORProbabilityMaps.grd"), overwrite = TRUE)

# Combine, make dataframe. -----------
# Do this in batches b/c it's pretty resource-intensive.
names(maps_cropped) <- paste0(names(maps_cropped), "_raw")

wd$tmp_df <- file.path(wd$bin,"tmp_df")
if(!dir.exists(wd$tmp_df) ) dir.create(wd$tmp_df)

maps_cropped_df_list <- pbmcapply::pbmclapply(
  1:nlayers(maps_cropped), mc.cores = 4, function(i){
    mdf <- stack(maps_cropped[[i]], maps_quantile_stack[[i]], maps_odds_stack[[i]]) %>%
      # So raster::as.data.frame has given me TWO big troubles here.
      # Some weird bug in raster::data.frame is messing up column names when long = TRUE.
      # AND if na.rm = TRUE, it seems to throw out *any* cell with an NA, not just a particular cell with an NA.
      # Very unhelpful if you have different ranges in your stack!
      raster::as.data.frame(xy = TRUE, long = FALSE, na.rm = FALSE) %>%
      # Hackey fixes:
      tidyr::pivot_longer(-c("x", "y"), names_to = "layer", values_to = "value") %>%
      dplyr::filter(!is.na(value))  %>%
      tidyr::separate(col = layer, into = c("ID", "method"), sep = "_(?=[^_]+$)")
    saveRDS(mdf, file = file.path(wd$tmp_df, paste0("df_list_", i, ".rds")))
})

maps_df <- list.files(wd$tmp_df, pattern = "df_list.*rds$", full.names = T) %>%
  lapply(readRDS) %>%
  bind_rows()

# Save.
saveRDS(maps_df, file = file.path(wd$bin, "maps_df.rds"))
