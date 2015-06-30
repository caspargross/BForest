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
for (k in 2:8) {
  #input_files <- as.list(paste("Data/Init_State/dis_1x1_2x100_alt", 1:8,".csv", sep=""))
  #input_files <- lapply(input_files, fread)
  #input_files <- lapply(input_files, function(x) rev_ycoordsDT(x, extent(maps_list[[1]])))
  #levelplot(out2rasterDT(input_files[[1]], var="species"))
  
  ## Run Model with Dispersal
  for (i in seq_along (maps_list)) { create_inputdir (paste("dis_front_full_fast", i, "rep_", k, sep=""),
                                                      climpath="Data/DWD/climate_feldberg.dat",
                                                      species=c("abiealba", "fagusilv", "piceabie"), ex=F,
                                                      LandClimRasterStack=maps_list[[i]],
                                                      inputfile=paste("Data/Init_State/dis_front_alt_fast",  i, ".csv", sep=""),
                                                      ctlfile="Data/Landclim/ctl_bforest_dis_2000.xml",
                                                      landtypefile="Data/Landclim/landtype.xml") }
  
  k<-2
  mclapply(paste("dis_front_full_fast", seq_along(maps_list), "rep_",k , sep=""), function(x) run_landclim_model(x, ctl_file="ctl_bforest_dis_2000.xml"), mc.cores=3)
  
}