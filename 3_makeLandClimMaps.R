#########################################################################
### Create LandClim maps 
#########################################################################

create_LandClim_Maps <- function(dem, no_aspect=F) {
require(raster)
require(rgdal)
#dem <- np_dem
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
soil <- slope
soil[!is.na(soil[])] <- 8


###  LandClim map "landtype".
landtype <- slope
landtype[!is.na(landtype[])] <- 1  ### Make Normal Landtype.


###  LandClim map "nitrogen".
nitro <- slope
nitro[!is.na(nitro[])] <- 1  ### no influence, standard value.

### Create raster-stack
maps <- stack(dem, slope, aspect, soil, landtype, nitro)
names(maps) <- c("dem", "slope", "aspect", "soil", "landtype", "nitro")
plot(maps, main="Original Maps")

### Resample LancClimMaps
maps25 <- resample_LandClim_maps(maps)
#plot(maps25>0, main="Resampled Maps (25)")
maps25
}









