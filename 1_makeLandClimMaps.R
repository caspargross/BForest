library(raster)
library(rgdal)

#########################################################################
### Create LandClim maps 
#########################################################################


### Create simplified  Model Landscape 

gk_projection<-CRS("+init=epsg:31468") # GK Zone 4 (Berchtesgaden)
nr <-200  #corresponds to 5000m
nc <- 4   # corresponds to 100m
res <- 25   #25m grid cells
ex <- extent(0, nc*res, 0, nr*res)
dem <- raster(nrows=nr, ncols=nc, ex)
projection(dem) <- gk_projection
dem

hmax <- 2000
hmin <- 500
gradient <- (hmax-hmin)/(nr*res)

dem[,seq(nc)] <- rev(hmin+(seq(nr)*gradient*res))
#dem[] <- rev(1:(nc*nr)) # 
plot(dem)

### LandClim map "slope" and "aspect".
#slope <- terrain(dem, filename="Data/slopeAspect.tif", opt='slope', unit="degrees", overwrite = T)
#plot(slope)
slope <- dem
slope[] <- 0


###  LandClim map "soil".
soil <- slope
soil[] <- 8


###  LandClim map "landtype".
landtype <- raster(slope)
landtype[] <- 1  ### Make "City" Landtype.


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



