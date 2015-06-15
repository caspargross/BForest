#########################################################################
### Create LandClim maps 
#########################################################################

create_LandClim_Maps <- function(dem, no_aspect=F) {
require(raster)
require(rgdal)
#dem <- dem_gk#´
#ex <- aui
if (no_aspect) {
  ### Without Slope
  slope <- dem 
  slope[] <- 0 
  ### Without Aspect
  aspect <- dem 
  aspect[] <- 1 } else { 
        slope <- terrain(dem, opt='slope', unit="degrees", overwrite = T) 
        aspect <- terrain(dem, opt='aspect', unit="degrees", overwrite = T) 
  }


###  LandClim map "soil".
soil <- raster(slope)
soil[] <- 8


###  LandClim map "landtype".
landtype <- raster(slope)
landtype[] <- 1  ### Make Normal Landtype.


###  LandClim map "nitrogen".
nitro <- raster(slope)
nitro[] <- 1  ### Make "City" Landtype.

### Create raster-stack
maps <- stack(dem, slope, aspect, soil, landtype, nitro)
names(maps) <- c("dem", "slope", "aspect", "soil", "landtype", "nitro")
plot(maps, main="Original Maps")

print(ex)
### Resample LancClimMaps
if (all(c(25,25) == res(maps25))) maps25 <- maps else maps25 <- resample_LandClim_maps(maps)
plot(maps25>0, main="Resampled Maps (25)")
maps25
}









