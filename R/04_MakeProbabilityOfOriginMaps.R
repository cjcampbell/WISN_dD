
# Setup -------------------------------------------------------------------

mydata_transformed <- readRDS( file.path(wd$bin, "mydata_transformed.rds") )

load( file.path(wd$bin, "my_isoscapes.RData") )

mapPath <- file.path( wd$bin, "maps")
if(!dir.exists(mapPath)) dir.create(mapPath)

range_raster <- readRDS( file.path(wd$bin, "range_raster.rds" ) )

# Create maps. -----------------------------------------------------------

pbapply::pblapply(1:nrow(mydata_transformed), function(x){
  comboFileName <-  paste0("Combined_", mydata_transformed$SampleName[x])
  if(!file.exists(file.path(mapPath, paste0(comboFileName, ".grd")))) {
    isocat::isotopeAssignmentModel(
      ID               = mydata_transformed$SampleName[x],
      isotopeValue     = mydata_transformed$dDprecip[x],
      SD_indv          = mydata_transformed$sdResid[x],
      precip_raster    = myisoscape$isoscape,
      precip_SD_raster = myisoscape$sd,
      savePath         = mapPath,
      additionalModels = range_raster,
      additionalModel_name = comboFileName,
      nClusters = FALSE
    )
  }
})



# Convert to other formats (for future analyses) --------------------------

# Load maps.
maps_cropped <- list.files(
  mapPath, pattern = "Combined.*.grd$", full.names = TRUE, recursive = T) %>%
  raster::stack()

# Also calculate probability quantiles.
maps_cumulativeSum_stack <-
  lapply(1:nlayers(maps_cropped), function(i) {
    out <- maps_cropped[[i]] %>%
      as.data.frame(xy=T, na.rm = F) %>%
      arrange(.[[3]]) %>%
      mutate(cumsum = cumsum(.[[3]])) %>%
      dplyr::select(x,y,cumsum) %>%
      rasterFromXYZ()
    names(out) <- names(maps_cropped[[i]])
    crs(out) <- crs(maps_cropped[[i]])
    return(out)
    }) %>%
  stack()
names(maps_cumulativeSum_stack) <- paste0(names(maps_cropped), "_csum")
writeRaster(maps_cumulativeSum_stack, filename = file.path(mapPath, "maps_cumulativeSum_stackProbabilityMaps.grd"), overwrite = TRUE)

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
  1:nlayers(maps_cropped), mc.cores = 2, function(i){
    mdf <- stack(maps_cropped[[i]], maps_cumulativeSum_stack[[i]], maps_odds_stack[[i]]) %>%
      # So raster::as.data.frame has given me TWO big troubles here.
      # Some weird bug in raster::data.frame is messing up column names when long = TRUE.
      # AND if na.rm = TRUE, it seems to throw out *any* cell with an NA, not just a particular cell with an NA.
      # Very unhelpful if you have different ranges in your stack!
      raster::as.data.frame(xy = TRUE, long = FALSE, na.rm = FALSE) %>%
      # Hackey fixes:
      tidyr::pivot_longer(-c("x", "y"), names_to = "layer", values_to = "value") %>%
      dplyr::filter(!is.na(value))  %>%
      tidyr::separate(col = layer, into = c("ID", "method"), sep = "_(?=[^_]+$)")
    ID <- mdf$ID[1]
    data.table::fwrite(mdf, file = file.path(wd$bin, "maps", paste0("df_list_", ID, ".csv")), row.names = F)
})

# maps_df <- list.files(wd$tmp_df, pattern = "df_list.*rds$", full.names = T) %>%
#   lapply(readRDS) %>%
#   bind_rows()
#
# # Save.
# saveRDS(maps_df, file = file.path(wd$bin, "maps_df.rds"))

# wd$tmp_df %>%
#   list.files(full.names = T) -> mf
# lapply(mf, function(x) {
#   out <- readRDS(x)
#   fwrite(out, file = file.path(wd$bin, "maps", paste0("df_list_", out$ID[1],".csv")), row.names = F)
# })
