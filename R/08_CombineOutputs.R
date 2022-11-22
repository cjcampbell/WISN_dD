
# Combine -----------------------------------------------------------------
mydata_transformed <- readRDS( file.path(wd$bin, "mydata_transformed.rds") )

distDirStats <- fread(file.path(wd$bin, "distDirStats.csv"))

allResults <- mydata_transformed %>% left_join(statsOut)

fwrite(allResults, file = file.path(wd$bin, "allResults.csv"), row.names = F)
