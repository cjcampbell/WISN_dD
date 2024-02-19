
# Combine -----------------------------------------------------------------
mydata_transformed <- readRDS( file.path(wd$bin, "mydata_transformed.rds") )

distDirStats <- fread(file.path(wd$bin, "distDirStats.csv"))
load(file.path(wd$bin, "mydata_clustered.Rdata"), verbose = T)

allResults <- mydata_transformed %>%
  left_join(., distDirStats, by = c("lon", "lat", "SampleName" = "ID")) %>%
  left_join(., mydata_clustered, by = "SampleName")

fwrite(allResults, file = file.path(wd$bin, "allResults.csv"), row.names = F)
