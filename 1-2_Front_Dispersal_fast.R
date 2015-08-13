require(parallel)
getwd()

dis_landscape_frnew <- function (alt) {
  ## Flat Gradient Model Landscape
  gk_projection<-CRS("+init=epsg:31467") # GK Zone 3 (Black Forest)
  nr <-80  #corresponds to 25000m
  nc <- 20   # corresponds to 25000m
  res <- 25   #25m grid cells
  ex <- extent(0, nc*res, 0, nr*res)
  dem_dis <- raster(nrows=nr, ncols=nc, ex)
  projection(dem_dis) <- gk_projection
  dem_dis[] <- alt
  maps25_frnew <- create_LandClim_Maps(dem_dis, no_aspect=T)
}

### Create map list

list_alt <- as.list (seq (400, 1600, 200))
maps_list_frnew <- lapply (list_alt, dis_landscape_frnew)


## Create Simulations
for (i in seq_along (maps_list_frnew)) { create_inputdir (paste("dis_frnew_fs-pa", i,  sep=""),
                                                    climpath="Data/DWD/climate_feldberg.dat",
                                                    species=c("fagusilv", "piceabie"), ex=F,
                                                    LandClimRasterStack=maps_list_frnew[[i]],
                                                    inputfile=F,
                                                    ctlfile="Data/Landclim/ctl_bforest_50.xml",
                                                    landtypefile="Data/Landclim/landtype.xml") }

for (i in seq_along (maps_list_frnew)) { create_inputdir (paste("dis_frnew_aa", i,  sep=""),
                                                    climpath="Data/DWD/climate_feldberg.dat",
                                                    species=c("abiealba"), ex=F,
                                                    LandClimRasterStack=maps_list_frnew[[i]],
                                                    inputfile=F,
                                                    ctlfile="Data/Landclim/ctl_bforest_50.xml",
                                                    landtypefile="Data/Landclim/landtype.xml") }



## Run LandClim Model
lapply(paste("dis_frnew_aa", seq_along(maps_list_frnew),  sep=""), function(x) run_landclim_model(x, ctl_file="ctl_bforest_50.xml"))
lapply(paste("dis_frnew_fs-pa", seq_along(maps_list_frnew),  sep=""), function(x) run_landclim_model(x, ctl_file="ctl_bforest_50.xml"))



## Create Inputfiles
initlist_fspaa <- as.list( paste("Simulations/dis_frnew_fs-pa", seq_along(maps_list_frnew), "/Output/fullOut_50.csv", sep=""))
initlist_aa <- as.list( paste("Simulations/dis_frnew_aa",  seq_along(maps_list_frnew), "/Output/fullOut_50.csv", sep=""))

# Initialisation for the first cycle of the run
t <- 20 # Timesteps for the calculation
th_y <- 2
initlist1 <- initlist_fspaa
initlist2 <- initlist_aa

# Run-Cycle
for (i in seq(t)){
t <- 1

outputlist<- as.list( paste("Data/Init_State/frnew/dis_frnew_alt",  seq_along(maps_list_frnew), "step_", t, ".csv", sep=""))
mask_frnew <- raster(extent(maps_list_frnew[[1]]), resolution=res(maps_list_frnew[[1]]))
mask_frnew[1:th_y,] <- 1   ## Starting Front: 100m

combine_input (initlist1,  
               initlist2,
               outputlist,
               1, 1, 1, 1,  
               aui=extent(mask_frnew),
               ma=mask_frnew)

create_inputdir (paste("dis_frnew", i, "rep_", k, sep=""),
                 climpath=paste("Data/DWD/climate_feldberg_",k,".dat", sep=""),
                 species=c("abiealba", "fagusilv", "piceabie"), ex=F,
                 LandClimRasterStack=maps_list[[i]],
                 inputfile=paste("Data/Init_State/dis_front_alt_fast",  i, ".csv", sep=""),
                 ctlfile="Data/Landclim/ctl_bforest_dis_2000.xml",
                 landtypefile="Data/Landclim/landtype.xml",
                 clim_rename=T)




for (k in 2:3) {   
  
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

