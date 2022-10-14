


library(tidyverse)
library(rgdal)
library(raster)
library(ncdf4)

# Extract some basic stats from GV phenology image
getAnnualMeanGV <- function(gv_pheno)
{
  mean <- calc(gv_pheno[[1:12]], fun=mean)
  #max <- calc(gv_pheno, fun=max)
  #min <- calc(gv_pheno, fun=min)
  #max_month <- calc(gv_pheno, fun=which.max)
  #min_month <- calc(gv_pheno, fun=which.min)
  
  #return( stack(mean, max, min, max_month, min_month) )
  return(mean)
}


# Choose five scenes for comparison - 2015 to 2019
#   the Cave Fire on the South Base occurred in 2016
pheno_2015 <- stack("D:/landsat_vafb/temp/stacked_scenes/pheno_NDVI_2015_annual.tif")
pheno_2016 <- stack("D:/landsat_vafb/temp/stacked_scenes/pheno_NDVI_2016_annual.tif")
pheno_2017 <- stack("D:/landsat_vafb/temp/stacked_scenes/pheno_NDVI_2017_annual.tif")
pheno_2018 <- stack("D:/landsat_vafb/temp/stacked_scenes/pheno_NDVI_2018_annual.tif")
pheno_2019 <- stack("D:/landsat_vafb/temp/stacked_scenes/pheno_NDVI_2019_annual.tif")
pheno_2020 <- stack("D:/landsat_vafb/temp/stacked_scenes/pheno_NDVI_2020_annual.tif")
pheno_2021 <- stack("D:/landsat_vafb/temp/stacked_scenes/pheno_NDVI_2021_annual.tif")
# Mean GV across growing season
mean_2015 <- getAnnualMeanGV(pheno_2015)
mean_2016 <- getAnnualMeanGV(pheno_2016)
mean_2017 <- getAnnualMeanGV(pheno_2017)
mean_2018 <- getAnnualMeanGV(pheno_2018)
mean_2019 <- getAnnualMeanGV(pheno_2019)
mean_2020 <- getAnnualMeanGV(pheno_2020)
mean_2021 <- getAnnualMeanGV(pheno_2021)

# Riparian Visualization
# Class Mask
classes <- raster("D:/SERDP/Vandenberg/Landsat/classified/classes_2019.tif", overwrite=TRUE)
# for now I'm using the 2019 image for Riparian and the 2015 image for the upland classes (due to fire changes)
riparian_mask <- classes == 1
riparian_mask[riparian_mask == 0] <- NA
riparian_2015 <- mean_2015*riparian_mask
riparian_2016 <- mean_2016*riparian_mask
riparian_2017 <- mean_2017*riparian_mask
riparian_2018 <- mean_2018*riparian_mask
riparian_2019 <- mean_2019*riparian_mask
plot(crop(chaparral_2018-chaparral_2016, extent(715000,735000,3828105,3855000)), zlim=c(-.1,.1))
writeRaster(riparian_2018 - riparian_2015,
            "D:/SERDP/scripts/riparian_groundwater_dependence_vafb/riparian_change.tif", overwrite=TRUE)

# Chaparral Visualization
classes <- raster("D:/SERDP/Vandenberg/Landsat/classified/classes_2015.tif", overwrite=TRUE)
# again using earlier before-fire image for upland classes
chaparral_mask <- classes == 3
chaparral_mask[chaparral_mask == 0] <- NA
chaparral_2015 <- mean_2015*chaparral_mask
chaparral_2016 <- mean_2016*chaparral_mask
chaparral_2017 <- mean_2017*chaparral_mask
chaparral_2018 <- mean_2018*chaparral_mask
chaparral_2019 <- mean_2019*chaparral_mask
plot(crop(chaparral_2018-chaparral_2016, extent(715000,735000,3828105,3845000)), zlim=c(-.1,.1))
writeRaster(chaparral_2018 - chaparral_2015,
            "D:/SERDP/scripts/riparian_groundwater_dependence_vafb/chaparral_change.tif", overwrite=TRUE)

# Grassland Visualization
grassland_mask <- classes == 4
grassland_mask[grassland_mask == 0] <- NA
grassland_2015 <- mean_2015*grassland_mask
grassland_2018 <- mean_2018*grassland_mask
plot(crop(grassland_2018-grassland_2015, extent(710685,739875,3828105,3868000)), zlim=c(-.1,.1))
writeRaster(grassland_2018 - grassland_2015,
            "D:/SERDP/scripts/riparian_groundwater_dependence_vafb/grassland_change.tif", overwrite=TRUE)

# total change fire 
writeRaster(mean_2018 - mean_2016,
            "D:/SERDP/scripts/riparian_groundwater_dependence_vafb/fire_change2.tif", overwrite=TRUE)
writeRaster(mean_2016 - mean_2015,
            "D:/SERDP/scripts/riparian_groundwater_dependence_vafb/fire_change_2015_2016.tif", overwrite=TRUE)
writeRaster(mean_2017 - mean_2015,
            "D:/SERDP/scripts/riparian_groundwater_dependence_vafb/fire_change_2015_2017.tif", overwrite=TRUE)
writeRaster(mean_2018 - mean_2015,
            "D:/SERDP/scripts/riparian_groundwater_dependence_vafb/fire_change_2015_2018.tif", overwrite=TRUE)
writeRaster(mean_2019 - mean_2015,
            "D:/SERDP/scripts/riparian_groundwater_dependence_vafb/fire_change_2015_2019.tif", overwrite=TRUE)
writeRaster(mean_2020 - mean_2015,
            "D:/SERDP/scripts/riparian_groundwater_dependence_vafb/fire_change_2015_2020.tif", overwrite=TRUE)
writeRaster(mean_2021 - mean_2015,
            "D:/SERDP/scripts/riparian_groundwater_dependence_vafb/fire_change_2015_2021.tif", overwrite=TRUE)


