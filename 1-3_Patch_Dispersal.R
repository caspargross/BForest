###########################
## Dispersal Dynamics    ##
##   Patch Dispersal     ##
##     Area: Model BF    ##
###########################

setwd("/home/caspar/Bachelor_Thesis/BForest/")
library(ModellingTools)
library(sp)
library(rgdal)
library(raster)


dis_landscape_p <- function (alt) {
  ## Flat Gradient Model Landscape
  gk_projection<-CRS("+init=epsg:31467") # GK Zone 3 (Black Forest)
  nr <-200  #corresponds to 5000m
  nc <- 200   # corresponds to 5000m
  res <- 25   #25m grid cells
  ex <- extent(0, nc*res, 0, nr*res)
  dem_dis <- raster(nrows=nr, ncols=nc, ex)
  projection(dem_dis) <- gk_projection
  dem_dis[] <- alt
  maps25_dis <- create_LandClim_Maps(dem_dis, no_aspect=T)
}

### Create map list

list_alt_p <- as.list (seq (400, 1800, 200))
maps_list_p <- lapply (list_alt_p, dis_landscape_p)


## Create Simulations
for (i in seq_along (maps_list_p)) { create_inputdir (paste("dis_patch_init_fs-pa_", i, sep=""),
                                                    climpath="Data/DWD/climate_feldberg.dat",
                                                    species=c("fagusilv", "piceabie"), ex=F,
                                                    LandClimRasterStack=maps_list_p[[i]],
                                                    inputfile=F,
                                                    ctlfile="Data/Landclim/ctl_bforest_50.xml",
                                                    landtypefile="Data/Landclim/landtype.xml") }

for (i in seq_along (maps_list_p)) { create_inputdir (paste("dis_patch_init_aa_", i, sep=""),
                                                    climpath="Data/DWD/climate_feldberg.dat",
                                                    species=c("abiealba"), ex=F,
                                                    LandClimRasterStack=maps_list_p[[i]],
                                                    inputfile=F,
                                                    ctlfile="Data/Landclim/ctl_bforest_50.xml",
                                                    landtypefile="Data/Landclim/landtype.xml") }

## Run LandClim Model for Background World
mclapply(as.list(paste("dis_patch_init_fs-pa_", seq_along(maps_list_p), sep="")), function(x) run_landclim_model(x, ctl_file="ctl_bforest_50.xml"), mc.cores=3)
mclapply(as.list(paste("dis_patch_init_aa_", seq_along(maps_list_p), sep="")), function(x) run_landclim_model(x, ctl_file="ctl_bforest_50.xml"), mc.cores=3)


## Check Background World
input_files <- as.list(paste("Simulations/dis_patch_init_fs-pa_", 1:8,"/Output/fullOut_50.csv", sep=""))
input_files <- lapply(input_files, fread)
input_files <- lapply(input_files, function(x) rev_ycoordsDT(x, extent(maps_list_p[[1]])))
levelplot(out2rasterDT(input_files[[1]], var="species"))



#### Create Inputfiles
mask <- raster(extent(maps_list_p[[1]]), resolution=res(maps_list_p[[1]]))

m1 <- (create_mask_p(1,1,20,20, mask, p=TRUE)) # 1x400
m2 <- (create_mask_p(1,2,20,10, mask, p=TRUE)) # 2x200
m3 <- (create_mask_p(2,2,10,10, mask, p=TRUE)) # 4x100
m4.1 <- (create_mask_p(3,3,6,7, mask, p=TRUE)) # 8x50  # 6 6x7 patches --> 252 Cells
m4.2 <- (create_mask_p(1,3,7,7, mask, p=TRUE)) # 8x50  # 3 7x7 patches --> 147 Cells | Total: 399 Cells
m4 <- m4.1 | m4.2 # Combine rasters
m5 <- (create_mask_p(4,4,5,5, mask, p=TRUE))   # 16x25
m6 <- (create_mask_p(5,4,4,5, mask, p=TRUE))   # 20x20
m7 <- (create_mask_p(5,5,4,4, mask, p=TRUE))   # 25x16
m8.1 <- (create_mask_p(7,7,3,3, mask, p=TRUE)) # 50x8  # 49 3x3 patches --> 441 Cells
m8.2 <- (create_mask_p(7,1,2,3, mask, p=TRUE, vshift=0)) # remove 7 2*3 = 42 Patches | Total: 399 Cells
m8 <- m8.1 - !is.na(m8.2)
m8[m8[]==0] <- NA
m9 <- (create_mask_p(10,10,2,2, mask, p=TRUE)) # 100x4
m10 <-(create_mask_p(10,20,2,1, mask, p=TRUE)) # 200x2
m11 <-(create_mask_p(20,20,1,1, mask, p=TRUE)) # 400x1


for (i in 1:11) {

initlist1_p <- as.list( paste("Simulations/dis_patch_init_fs-pa_", seq_along(maps_list_p), "/Output/fullOut_50.csv", sep=""))
initlist2_p <- as.list( paste("Simulations/dis_patch_init_aa_", seq_along(maps_list_p), "/Output/fullOut_50.csv", sep=""))
outputlist_p<- as.list( paste("Data/Init_State/dis_patch_init_alt_", seq_along(maps_list_p),"_m", i, ".csv", sep=""))


combine_input (initlist1_p,  
               initlist2_p, 
               outputlist_p,
               1, 1, 2, 100,  
               extent(maps_list_p[[1]]),
               ma=eval(as.name(paste("m",i, sep=""))))

}

## Check Input Files

input_files <- as.list(paste("Data/Init_State/dis_patch_init_alt_", 1:8,"_m11.csv", sep=""))
input_files <- lapply(input_files, fread)
input_files <- lapply(input_files, function(x) rev_ycoordsDT(x, extent(maps_list_p[[1]])))
levelplot(out2rasterDT(input_files[[4]], var="species"))

## Run Model with Dispersal

for (j in 1:1){ 
for (i in seq_along (maps_list)) { create_inputdir (paste("dis_1x100_", i, sep=""),
                                                    climpath="Data/DWD/climate_feldberg.dat",
                                                    species=c("abiealba", "fagusilv", "piceabie"), ex=F,
                                                    LandClimRasterStack=maps_list[[i]],
                                                    paste("Data/Init_State/dis_patch_init_alt_", seq_along(maps_list_p),"_m", j, ".csv", sep=""),
                                                    ctlfile="Data/Landclim/ctl_bforest_dis_2000.xml",
                                                    landtypefile="Data/Landclim/landtype.xml") }

mclapply(as.list(paste("dis_1x100_", seq_along(maps_list), sep="")), function(x) run_landclim_model(x, ctl_file="ctl_bforest_dis_2000.xml"), mc.cores=3)
}



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






