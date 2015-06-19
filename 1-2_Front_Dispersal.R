###########################
## Dispersal Dynamics    ##
##   Front Propagation   ##
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
nr <-400  #corresponds to 25000m
nc <- 100   # corresponds to 25000m
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
for (i in seq_along (maps_list)) { create_inputdir (paste("dis_front_fs-pa_", i, sep=""),
                                                   climpath="Data/DWD/climate_feldberg.dat",
                                                   species=c("fagusilv", "piceabie"), ex=F,
                                                   LandClimRasterStack=maps_list[[i]],
                                                   inputfile=F,
                                                   ctlfile="Data/Landclim/ctl_bforest_50.xml",
                                                   landtypefile="Data/Landclim/landtype.xml") }

for (i in seq_along (maps_list)) { create_inputdir (paste("dis_front_aa_", i, sep=""),
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

initlist1 <- as.list( paste("Simulations/dis_front_fs-pa_", seq_along(maps_list), "/Output/fullOut_50.csv", sep=""))
initlist2 <- as.list( paste("Simulations/dis_front_aa_", seq_along(maps_list), "/Output/fullOut_50.csv", sep=""))
outputlist<- as.list( paste("Data/Init_State/dis_front_alt", seq_along(maps_list), ".csv", sep=""))

mask <- raster(extent(maps_list[[1]]), resolution=res(maps_list[[1]]))
mask[396:400,] <- 1   ## Starting Front: 100m

combine_input (initlist1,  
               initlist2, 
               outputlist,
               1, 1, 1, 1,  
               mask)

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

plot(extent(maps_list[[1]]))


## Check Dispersal Output files
output_files <- as.list(paste("Simulations/dis_1x100_4/Output/fullOut_", seq(5,200,5), ".csv", sep=""))
output_files <- lapply(output_files, fread)
output_files <- lapply(output_files, function(x) rev_ycoordsDT(x, extent(maps_list[[1]])))
for (i in 1:length(output_files)) {
  png(filename=paste("Animate/dis_1x100_4_dec", i,".png", sep=""))
  print(levelplot(out2rasterDT(output_files[[i]], var="species")))
  dev.off()
}

#### Load Dispersal Model results
list_results <- as.list (paste("Simulations/dis_1x100_1/Output/fullOut_", seq(5, 200, 5), ".csv", sep=""))
dis_results <- lapply(list_results, fread)
dis_results <- lapply(dis_results, function(x) rev_ycoordsDT(x, extent(maps_list[[1]])))
#ras_dis_results <- lapply(dis_results, function(x) out2rasterDT(x, var="species"))

## Plot Species Distribution
for (i in 1:length(output_files)) {
  #png(filename=paste("Animate/dis_1x100_4_dec", i,".png", sep=""))
  print(levelplot(out2rasterDT(output_files[[i]], var="species")))
  Sys.sleep(0.5)
  #dev.off()
}


######### Create Plot of Established Distances (Age>60yr)

dist_quartile<-data.table("Year"=0, "Elevation"=0, "Distance"=0)
elevations <- as.integer(list_alt[])
for (j in seq_along (maps_list)) {
  
list_results <- as.list (paste("Simulations/dis_1x100_",j,"/Output/fullOut_", seq(5, 200, 5), ".csv", sep=""))
dis_results <- lapply(list_results, fread)
dis_results <- lapply(dis_results, function(x) rev_ycoordsDT(x, extent(maps_list[[1]])))
print(paste("Loaded map",j))

distances <- list()

for (i in 1:length(dis_results)) {
    DT<-dis_results[[i]]
    DT[,biomass_cohort:=biomass*stems]
    y_origin <- 1200
    DT <- DT[,.(biomass_cohort=sum(biomass_cohort)),by=list(species, cell, xcoord, ycoord,age)]
    setkey(DT,species)
    DT <- DT["abiealba"]
    setkey(DT, age)
    DT <- DT[age>=60]
    DT[,dist:= abs(ycoord - y_origin),]
    #hist(DT$dist)
    distances[[i]]<-DT
  }

  
   for (i in 1:length(distances)) {
      dist_quartile <- rbind(dist_quartile, list( i*50, list_alt[[j]], quantile(distances[[i]]$dist, 0.95, na.rm=T)))
  }
}

dist_quartile <- dist_quartile[-1,]
dist_quartile$Elevation <- as.factor(dist_quartile$Elevation)

library(ggplot2)
distplot <- ggplot(dist_quartile, aes(x=Year, y=Distance, group=Elevation, col = Elevation))
distplot+geom_line()


######### Create Plot of Established Distances (Age>60yr)

for (j in seq_along (maps_list)) {
  
  list_results <- as.list (paste("Simulations/dis_1x100_",j,"/Output/fullOut_", seq(5, 200, 5), ".csv", sep=""))
  dis_results <- lapply(list_results, fread)
  dis_results <- lapply(dis_results, function(x) rev_ycoordsDT(x, extent(maps_list[[1]])))
  print(paste("Loaded map",j))
  
  distances <- list()
  
  for (i in 1:length(dis_results)) {
    DT<-dis_results[[i]]
    DT[,biomass_cohort:=biomass*stems]
    y_origin <- 1200
    DT <- DT[,.(biomass_cohort=sum(biomass_cohort)),by=list(species, cell, xcoord, ycoord,age)]
    setkey(DT,species)
    DT <- DT["abiealba"]
    setkey(DT, age)
    DT <- DT[age>=60]
    DT[,dist:= abs(ycoord - y_origin),]
    #hist(DT$dist)
    distances[[i]]<-DT
    png(paste("Animate/abiealba_est_alt", j, "dec", i, ".png", sep=""))
    print(plot(extent(maps_list[[1]]), type="n", main=paste("Elevation:",list_alt[[j]],"Year:",i*50), xlab="Latitude", ylab="Longitude"))
    print(plot(out2rasterDT(DT), add=T, legend=F, col="chartreuse4"))
    dev.off()
  }
  
}

dist_quartile <- dist_quartile[-1,]
dist_quartile$Elevation <- as.factor(dist_quartile$Elevation)

library(ggplot2)
distplot <- ggplot(dist_quartile, aes(x=Year, y=Distance, group=Elevation, col = Elevation))
distplot+geom_line()






