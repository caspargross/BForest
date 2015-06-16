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
nr <-100  #corresponds to 5000m
nc <- 100   # corresponds to 1000m
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
                                                   ctlfile="Data/Landclim/ctl_bforest_50.xml",
                                                   landtypefile="Data/Landclim/landtype.xml") }

for (i in seq_along (maps_list)) { create_inputdir (paste("dis_init_aa_", i, sep=""),
                                                    climpath="Data/DWD/climate_feldberg.dat",
                                                    species=c("abiealba"), ex=F,
                                                    LandClimRasterStack=maps_list[[i]],
                                                    inputfile=F,
                                                    ctlfile="Data/Landclim/ctl_bforest_50.xml",
                                                    landtypefile="Data/Landclim/landtype.xml") }

## Run LandClim Model
lapply(paste("dis_init_aa_", seq_along(maps_list), sep=""), function(x) run_landclim_model(x, ctl_file="ctl_bforest_50.xml"))
lapply(paste("dis_init_fs-pa_", seq_along(maps_list), sep=""), function(x) run_landclim_model(x, ctl_file="ctl_bforest_50.xml"))

## Create Inputfiles

initlist1 <- as.list( paste("Simulations/dis_init_fs-pa_", seq_along(maps_list), "/Output/fullOut_50.csv", sep=""))
initlist2 <- as.list( paste("Simulations/dis_init_aa_", seq_along(maps_list), "/Output/fullOut_50.csv", sep=""))
outputlist<- as.list( paste("Data/Init_State/dis_1x1_2x100_alt", seq_along(maps_list),  sep=""))

combine_input (initlist1,  # Inputfile1
               initlist2,  # Inputfile2
               outputfile = outputlist,
               npatch_col=1, npatch_row=1, patch_width=2, patch_length=100,  #Form of patches
               aui=extent(maps_list[[1]]) )


## Run Model with Dispersal
for (i in seq_along (maps_list)) { create_inputdir (paste("dis_1x100_", i, sep=""),
                                                    climpath="Data/DWD/climate_feldberg.dat",
                                                    species=c("abiealba"), ex=F,
                                                    LandClimRasterStack=maps_list[[i]],
                                                    inputfile="Data/Init_State/dis_1x1_2x100.csv",
                                                    ctlfile="Data/Landclim/ctl_bforest_dis_2000.xml",
                                                    landtypefile="Data/Landclim/landtype.xml") }
lapply(paste("dis_init_aa_", seq_along(maps_list), sep=""), function(x) run_landclim_model(x, ctl_file="ctl_bforest_dis_2000.xml"))

