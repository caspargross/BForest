###########################
## Dispersal Dynamics    ##
##     Area: Model BF    ##
###########################
setwd("/home/caspar/Bachelor_Thesis/BForest/")
library(ModellingTools)
library(sp)
library(rgdal)
library(raster)


dis_landscape <- function (alt) {
## Flat Gradient Model Landscape
gk_projection<-CRS("+init=epsg:31467") # GK Zone 3 (Black Forest)
nr <-200  #corresponds to 5000m
nc <- 200   # corresponds to 1000m
res <- 25   #25m grid cells
ex <- extent(0, nc*res, 0, nr*res)
dem_dis <- raster(nrows=nr, ncols=nc, ex)
projection(dem_dis) <- gk_projection
dem_dis[] <- alt
maps25_dis <- create_LandClim_Maps(dem_dis, no_aspect=T)
}

### Create map list
list_alt <- as.list (seq (400, 1800, 200))
maps_list <- lapply (list_alt, dis_landscape)


## Create Simulations
for (i in seq_along (maps_list)) { create_inputdir (paste("dis_init_fs-pa_", i, sep=""),
                                                   climpath="Data/DWD/climate_feldberg.dat",
                                                   species=c("fagusilv", "piceabie"), ex=F,
                                                   LandClimRasterStack=maps_list[[i]],
                                                   inputfile=F,
                                                   ctlfile="Data/Landclim/ctl_bforest.xml",
                                                   landtypefile="Data/Landclim/landtype.xml") }

for (i in seq_along (maps_list)) { create_inputdir (paste("dis_init_aa_", i, sep=""),
                                                    climpath="Data/DWD/climate_feldberg.dat",
                                                    species=c("fagusilv", "piceabie"), ex=F,
                                                    LandClimRasterStack=maps_list[[i]],
                                                    inputfile=F,
                                                    ctlfile="Data/Landclim/ctl_bforest.xml",
                                                    landtypefile="Data/Landclim/landtype.xml") }

## Run LandClim Model
lapply(paste("dis_init_aa_", seq_along(maps_list), sep=""), run_landclim_model)



