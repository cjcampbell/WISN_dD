
# Load metadata.
meta_path <- list.files(
  path = wd$data, pattern = "Carpenter_SnipeSamples_9.22.21.xlsx",
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

mydata <- full_join(meta, res, by = c("SampleName" = "Sample ID")) %>%
  filter(!is.na(d2H))

# Optionally, check it out and make sure it looks okay.
# View(mydata)


# Georeferencing ----------------------------------------------------------

# We conducted automatic georefencing followed by manual vallidation.
# Code for automatic georeferencing is below:

rerunAutomaticGeoreferencing <- FALSE
if(rerunAutomaticGeoreferencing) {

  df <- mydata %>%
    select(c(Location, County)) %>%
    distinct %>%
    dplyr::mutate(
      verbatimCounty = County,
      County = case_when(
        County == "Broward" & Location == "Three Forks" ~ Brevard,
        TRUE ~ County
      ) ,

      lookup = case_when(
        is.na(County) ~ paste0(Location, ", FL"),
        TRUE ~  paste0(Location, ", ", County, " County, FL")
      )
    )


  georef <- TRUE

  if(georef == TRUE){
    library(ggmap)
    if(!has_google_key()){
      myAPIkey <- read.table("~/googleAPIkey.txt") %>% unlist %>% paste
      register_google(key = myAPIkey)
      if(!has_google_key()) stop("Google API Key not found. See '?ggmap::register_google`.")
    }

    addresses <- df$lookup

    coded <- lapply(addresses, function(i){
      x <- geocode(
        i, output = "latlona", messaging = F,
        override_limit = T, force = T, source = "google"
      )
      while(is.na(x$lon)){
        Sys.sleep(10)
        x <- geocode(
          i, output = "more", messaging = F,
          override_limit = T, force = T, source = "google"
        )
      }
      Sys.sleep(1)
      tb <- cbind(ADDRESS = i, x)
      return(tb)
    } )

    geo <- plyr::ldply(coded) %>%
      dplyr::rename(lookupAddress = ADDRESS, returnAddress = address) %>%
      dplyr::select(lookupAddress, returnAddress, everything())

    write.csv(geo, file = file.path(wd$bin, "georeferenced.csv"), row.names = F)

  }

} else {

  # Load manually vallidated georeferencing.
  geo_path <- list.files(
    path = wd$data, pattern = "georeferenced_updated9.22.21.xlsx$",
    full.names = TRUE, recursive = T)
  geo <- read_excel(geo_path) %>%
    dplyr::filter(`...1` != 15) %>% # remove typo row
    dplyr::mutate(
      lookupAddress = case_when(
        lookupAddress == "Half Moon WMA, FL" ~ "Half Moon WMA, Sumter County, FL",
        lookupAddress == "Three Forks Conservation Area, Brevard County, FL" ~ "Three Forks, Brevard County, FL",
        TRUE ~ lookupAddress
      )
    ) %>%
    dplyr::select(-`...1`) %>%
    dplyr::rename( localityNotes = `...6` )

  }

## Combine.
df1 <- mydata %>%
  select(c(Location, County)) %>%
  distinct %>%
  dplyr::mutate(
    lookupAddress = case_when(
      is.na(County) ~ paste0(Location, ", FL"),
      TRUE ~  paste0(Location, ", ", County, " County, FL")
    )
  ) %>%
  full_join(., geo) %>%
  dplyr::select(Location, County, lon, lat) %>%
  filter(!is.na(lon)) %>%
  distinct()

mydata <- left_join(mydata, df1, by = c("Location", "County"))
