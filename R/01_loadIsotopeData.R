
# Load metadata.
meta_path <- list.files(
  path = wd$data, pattern = "SnipeSamples_4.20.21.xlsx",
  full.names = TRUE, recursive = T)
meta <- read_excel(meta_path)

# Load dD results.
res_path <- list.files(
  path = wd$data, pattern = "Carpenter results FWC 210611",
  full.names = TRUE, recursive = T)
res <- read_excel(res_path)

# Check how to join df.
if(sum(meta$SampleName %in% res$`Sample ID`) != nrow(res)) {
  stop("Sample ID do not match up.")
}

mydata <- full_join(meta, res, by = c("SampleName" = "Sample ID"))

# Optionally, check it out and make sure it looks okay.
# View(mydata)

