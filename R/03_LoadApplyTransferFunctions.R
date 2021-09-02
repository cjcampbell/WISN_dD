# Load  ----------------------------------------------------------

if(!exists("mydata")) source(file.path(wd$R, "01_loadIsotopedata.R"))

# Apply transfer functions, define molt status  ---------------------------------

# Placeholder:
# Common name	              n	  Sex	            Life stage	  Equation	                Tissue	  R^2	  P-value	  Precip	Isoscape source	  Continent	  Question	Regression	Equation source 	            Original dataset	  Notes
# Passerines and waterfowl  19  not specified	  not specified	d2Hf = 0.96*d2Hp -23.96	  feather	  0.95	<0.001	  GSD	Bowen et al. 2005	    NoA       	2.1	      RMA	        Hobson et al. 2010  Condor	  Clark et al. 2006

# sdResid is SYNTHETIC atm.

mydata_transformed <- mydata %>%
  dplyr::mutate_if(is.factor, as.character) %>%
  dplyr::mutate(
    dDprecip =  (d2H + 23.96) / 0.96,
    sdResid = 10
  )

saveRDS(mydata_transformed, file = file.path(wd$bin, "mydata_transformed.rds"))

