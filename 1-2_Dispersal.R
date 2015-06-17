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
outputlist<- as.list( paste("Data/Init_State/dis_1x1_2x100_alt", seq_along(maps_list), ".csv", sep=""))

combine_input (initlist1,  
               initlist2, 
               outputlist,
               1, 1, 2, 100,  
               extent(maps_list[[1]]))

## Check Input Files

input_files <- as.list(paste("Data/Init_State/dis_1x1_2x100_alt", 1:8,".csv", sep=""))
input_files <- lapply(input_files, fread)
input_files <- lapply(input_files, function(x) rev_ycoordsDT(x, extent(maps_list[[1]])))
levelplot(out2rasterDT(input_files[[1]], var="species"))

## Run Model with Dispersal
for (i in seq_along (maps_list)) { create_inputdir (paste("dis_1x100_", i, sep=""),
                                                    climpath="Data/DWD/climate_feldberg.dat",
                                                    species=c("abiealba", "fagusilv", "piceabie"), ex=F,
                                                    LandClimRasterStack=maps_list[[i]],
                                                    inputfile=paste("Data/Init_State/dis_1x1_2x100_alt", i,".csv", sep=""),
                                                    ctlfile="Data/Landclim/ctl_bforest_dis_2000.xml",
                                                    landtypefile="Data/Landclim/landtype.xml") }

lapply(paste("dis_1x100_", seq_along(maps_list), sep=""), function(x) run_landclim_model(x, ctl_file="ctl_bforest_dis_2000.xml"))






############ REDO INPUT FILES!
# error search

err_init1 <- lapply(initlist1, fread)
err_init1<- lapply(err_init1, rev_ycoordsDT)
levelplot(out2rasterDT(err_init1[[6]], var="species"))

err_init2 <- lapply(initlist2, fread)
err_init2<- lapply(err_init2, rev_ycoordsDT)
levelplot(out2rasterDT(err_init2[[6]], var="species"))




#### Load Dispersal Model results
list_results <- as.list (paste("Simulations/dis_1x100_5/Output/fullOut_", seq(5, 200, 5), ".csv", sep=""))
dis_results <- lapply(list_results, fread)
dis_results <- lapply(dis_results, rev_ycoordsDT)
ras_dis_results <- lapply(dis_results, function(x) out2raster(x, var="species"))

levelplot(out2rasterDT(dis_results[[20]], var="species"))


