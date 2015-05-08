# 1_BForest
#############################################
## Interaction Abies Alba <-> Picea Abies  ##
##           Black Forest NP               ##
##            Caspar Gross                 ##
#############################################

library(raster)
source("0_functions.R")

#Load Digital Elevation Model for National Park Nördlicher SW
dem <- raster ("Data/dem/n48_e008_1arc_v3.bil")
np_extent <- extent (c(8.15,8.4,48.5,48.7))
dem<-crop(dem, np_extent)
plot(dem)
points(8.276757, 48.578119) #Riesenköpfle UTM

#Convert to GK3
dem_gk <- projectRaster(dem, crs="+init=epsg:31467")  # Project to GK Zone 3
np_gk_extent <- extent (c(3437400, 3455900, 5374500, 5395500))
dem_gk <- crop(dem_gk, np_gk_extent)
plot(dem_gk)
points(3446740, 5382515) #Riesenköpfle GK



