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
  require(parallel)
  
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
  for (i in seq_along (maps_list)) { create_inputdir (paste("dis_front_fs-pa_fast", i,  sep=""),
                                                      climpath="Data/DWD/climate_feldberg.dat",
                                                      species=c("fagusilv", "piceabie"), ex=F,
                                                      LandClimRasterStack=maps_list[[i]],
                                                      inputfile=F,
                                                      ctlfile="Data/Landclim/ctl_bforest_50.xml",
                                                      landtypefile="Data/Landclim/landtype.xml") }
  
  for (i in seq_along (maps_list)) { create_inputdir (paste("dis_front_aa_fast", i,  sep=""),
                                                      climpath="Data/DWD/climate_feldberg.dat",
                                                      species=c("abiealba"), ex=F,
                                                      LandClimRasterStack=maps_list[[i]],
                                                      inputfile=F,
                                                      ctlfile="Data/Landclim/ctl_bforest_50.xml",
                                                      landtypefile="Data/Landclim/landtype.xml") }
  
  
  
  ## Run LandClim Model
  lapply(paste("dis_front_aa_fast", seq_along(maps_list),  sep=""), function(x) run_landclim_model(x, ctl_file="ctl_bforest_50.xml"))
  lapply(paste("dis_front_fs-pa_fast", seq_along(maps_list),  sep=""), function(x) run_landclim_model(x, ctl_file="ctl_bforest_50.xml"))
  
  
  
  ## Create Inputfiles
  initlist1 <- as.list( paste("Simulations/dis_front_fs-pa_fast", seq_along(maps_list), "/Output/fullOut_50.csv", sep=""))
  initlist2 <- as.list( paste("Simulations/dis_front_aa_fast",  seq_along(maps_list), "/Output/fullOut_50.csv", sep=""))
  outputlist<- as.list( paste("Data/Init_State/dis_front_alt_fast",  seq_along(maps_list), ".csv", sep=""))
  
  mask <- raster(extent(maps_list[[1]]), resolution=res(maps_list[[1]]))
  mask[396:400,] <- 1   ## Starting Front: 100m
  
  combine_input (initlist1,  
                 initlist2,
                 outputlist,
                 1, 1, 1, 1,  
                 aui=extent(mask),
                 ma=mask)
  
  ## Check Input Files
  
  #input_files <- as.list(paste("Data/Init_State/dis_1x1_2x100_alt", 1:8,".csv", sep=""))
  #input_files <- lapply(input_files, fread)
  #input_files <- lapply(input_files, function(x) rev_ycoordsDT(x, extent(maps_list[[1]])))
  #levelplot(out2rasterDT(input_files[[1]], var="species"))
  
  
  
  
  for (k in 1:3) {   
    
    ## Create Random Weather files
    random_weather_bf(paste("Data/DWD/climate_feldberg_",k,".dat", sep=""))
    ## Run Model with Dispersal
    for (i in seq_along (maps_list)) { create_inputdir (paste("dis_front_full_fast", i, "rep_", k, sep=""),
                                                        climpath=paste("Data/DWD/climate_feldberg_",k,".dat", sep=""),
                                                        species=c("abiealba", "fagusilv", "piceabie"), ex=F,
                                                        LandClimRasterStack=maps_list[[i]],
                                                        inputfile=paste("Data/Init_State/dis_front_alt_fast",  i, ".csv", sep=""),
                                                        ctlfile="Data/Landclim/ctl_bforest_dis_2000.xml",
                                                        landtypefile="Data/Landclim/landtype.xml",
                                                        clim_rename=T) }
    
  }
  sink("Logs/log_dis_front.txt")
  for (k in 1:3) {
    print(paste("REPETITION NO.", k))
    mclapply(as.list(paste("dis_front_full_fast", seq_along(maps_list), "rep_",k , sep="")), function(x) run_landclim_model(x, ctl_file="ctl_bforest_dis_2000.xml"), mc.cores=3)
  }
  print("SUCCESS")
  sink()


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




d
#### Load Dispersal Model results
list_results <- as.list (paste("Simulations/dis_front_full_fast4rep_2/Output/fullOut_", seq(5, 200, 5), ".csv", sep=""))
dis_results <- lapply(list_results, fread)
dis_results <- lapply(dis_results, function(x) rev_ycoordsDT(x, extent(maps_list[[1]])))
#ras_dis_results <- lapply(dis_results, function(x) out2rasterDT(x, var="species"))

## Plot Species Distribution
for (i in 1:length(output_files)) {
  #png(filename=paste("Animate/dis_1x100_4_dec", i,".png", sep=""))
  print(levelplot(out2rasterDT(dis_results[[30]], var="species")))
  Sys.sleep(0.5)
  #dev.off()
}

  
  ####TEST STUFF
  TT <- dis_results[[20]]
  TT <- TT[species == "abiealba",.(row_bio=sum(biomass*stems)), by=row]
  quantile(TT$row)
  plot(TT$row, TT$row_bio, type="l")
  
  ######### Create Plot of Established Distances (Age>60yr)

dist_quartile<-data.table("Year"=0, "Elevation"=0, "Distance"=0, "Repetition"=0)
elevations <- as.integer(list_alt[])
for (k in 1:3) { 

print(paste("Repetition", k))
for (j in seq_along (maps_list)) {
distances <- list() 
list_results <- as.list (paste("Simulations/dis_front_full_fast",j,"rep_",k,"/Output/fullOut_", seq(5, 200, 5), ".csv", sep=""))
dis_results <- lapply(list_results, fread)
dis_results <- lapply(dis_results, function(x) rev_ycoordsDT(x, extent(maps_list[[1]])))
print(paste("Loaded map",j))

for (i in 1:length(dis_results)) {
    DT<-dis_results[[i]]
    DT[,biomass_cohort:=biomass*stems]
    y_origin <- 9950
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
      dist_quartile <- rbind(dist_quartile, list( i*50, list_alt[[j]], quantile(distances[[i]]$dist, 0.95, na.rm=T), k))
  }
}
}

dist_quartile <- dist_quartile[-1,]
dist_quartile$Elevation <- as.factor(dist_quartile$Elevation)
dist_means <- dist_quartile[,.(mean=mean(Distance, na.rm=T),sd=sd(Distance, na.rm=T)),by=.(Year, Elevation)]
dist_means[, dmean:=(mean-shift(mean, type="lag")), by=.(Elevation)]

library(ggplot2)
distplot <- ggplot(dist_means, aes(x=Year, group=Elevation, col = Elevation))+
geom_line(aes(y=mean)) +
theme_cas() 
distplot


ddistplot <- ggplot(dist_means, aes(x=Year, group=Elevation, col = Elevation))+
    geom_line (aes(y=dmean )) +
    theme_cas() 
ddistplot
  
  
distplot2 <- ggplot(dist_quartile[Repetition==5,,], aes(x=Year, y=Distance, group=Elevation, col = Elevation))
distplot2+geom_line()



######### Create Plot of Established Distances (Age>60yr)

for (j in seq_along (maps_list)) {
  
  list_results <- as.list (paste("Simulations/dis_1x100",j,"/Output/fullOut_", seq(5, 200, 5), ".csv", sep=""))
  dis_results <- lapply(list_results, fread)
  dis_results <- lapply(dis_results, function(x) rev_ycoordsDT(x, extent(maps_list[[1]])))
  print(paste("Loaded map",j))
  
  distances <- list()
  
  for (i in 1:length(dis_results)) {
    DT<-dis_results[[i]]
    DT[,biomass_cohort:=biomass*stems]
    y_origin <- 9950
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

####
# Make nice maps
  ele <- c(2,4,6)
  dec <- c(10,50,100)
  gr <- expand.grid(ele,dec)
 plt_front <- as.list(paste("Simulations/dis_front_full_fast", gr[,1], "rep_2/Output/fullOut_", gr[,2], ".csv" ,sep=""))
 plt_front <- lapply(plt_front, fread)
 plt_front <- lapply(plt_front, function(x) rev_ycoordsDT(x, extent(maps_list[[1]])))
 levelplot(out2rasterDT(plt_front[[2]]))
plt_front[[12]]