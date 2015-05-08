library(raster)
library(rgdal)

#########################################################################
### Create LandClim maps 
#########################################################################


### LandClim map "slope" and "aspect".
slope <- terrain(dem, filename="Data/slopeAspect.tif", opt='slope', unit="degrees", overwrite = T)

###  LandClim map "soil".
soil <- slope
soil[] <- 8


###  LandClim map "landtype".
landtype <- raster(slope)
landtype[] <- 1  ### Make Normal Landtype.


### Aspect
aspect <- raster(slope)
aspect[] <- 0


###  LandClim map "nitrogen".
nitro <- raster(slope)
nitro[] <- 1  ### Make "City" Landtype.

### Create raster-stack
maps <- stack(dem, slope, aspect, soil, landtype, nitro)
names(maps) <- c("dem", "slope", "aspect", "soil", "landtype", "nitro")
plot(maps)

# maps25 <- resampleLandClimMaps(LandClimRasterStack=maps)
maps25 <- maps
plot(maps25$dem)



