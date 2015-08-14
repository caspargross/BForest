require(parallel)
getwd()



############## 
#Functions:
##############

# create landscape
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

# create .xml for different timesteps / initfiles
changexml <- function (x) {
  ctl <- readLines("Data/Landclim/ctl_frnew.xml")
  string <- paste("\t\t<decades>",x,"</decades></full>")
  ctl[64] <- string
  writeLines(ctl, "Data/Landclim/ctl_frnew.xml")
  # change line no 64
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


# Run-Cycle
n <- length(maps_list_frnew)
th_y <- list(4) ## Starting Front: 100m
th_y <- rep(th_y, n)
t <- 50  ## Number of Decades 

initlist2 <- initlist_aa
for (t in seq(t)){
outputlist<- as.list( paste("Data/Init_State/frnew/dis_frnew_alt",  1:n , "step_", t, ".csv", sep=""))

mask_frnew <- list(raster(extent(maps_list_frnew[[1]]), resolution=res(maps_list_frnew[[1]])))
mask_frnew <- rep(mask_frnew, n)
for (i in 1:n){ mask_frnew[[i]][1:th_y[[1]],] <- 1 }

## Combine Input +mask
init1 <- lapply(initlist_fspaa, function (x) fread(x, integer64="numeric"))
init2 <- lapply(initlist2, function (x) fread(x, integer64="numeric")) 
init1 <- lapply(init1, function(x) rev_ycoordsDT(x, extent(mask_frnew[[1]])))
init2 <- lapply(init2, function(x) rev_ycoordsDT(x, extent(mask_frnew[[1]])))
out <- list()
out <- mapply(mask_apply, x=init1, y=init2, mask=mask_frnew,  SIMPLIFY=F, USE.NAMES = F )
for (i in 1:length(out)) {write.table(out[[i]], file=outputlist[[i]], row.names=F, col.names=T, sep=",", quote=F)}

for (i in seq_along(maps_list_frnew))
  create_inputdir (paste("dis_frnew/dis_frnew_alt", i, "_step", t, sep=""),
                 climpath=paste("Data/DWD/climate_feldberg_",k,".dat", sep=""),
                 species=c("abiealba", "fagusilv", "piceabie"), ex=F,
                 LandClimRasterStack=maps_list_frnew[[i]],
                 inputfile=outputlist[[i]],
                 ctlfile="Data/Landclim/ctl_frnew.xml",
                 landtypefile="Data/Landclim/landtype.xml",
                 clim_rename=T)

mclapply(as.list(paste("dis_frnew/dis_frnew_alt", seq_along(maps_list_frnew), "_step", t, sep="")), function(x) run_landclim_model(x, ctl_file="ctl_frnew.xml"), mc.cores=3)
frfiles <- as.list(paste("Simulations/dis_frnew/dis_frnew_alt", seq_along(maps_list_frnew), "_step",t,"/Output/fullOut_1.csv",  sep=""))
print(frfiles)
initlist2 <- frfiles
frfiles <- lapply(frfiles, fread)
frfiles <- lapply(frfiles, function(x) rev_ycoordsDT(x, aui=extent(maps_list_frnew[[2]])))
if(plot==T) print(levelplot(out2rasterDT(frfiles[[4]])))
foo <- function(x) {
  DT <- x
  DT[,bio_cohort:=biomass*stems,]
  th_y <- DT[species=="abiealba" & age>70 , max(ycoord)] +160    ## IMPORTANT STEP: Cutoff Distance = Max range of esdtablished trees (>60) + maximal seedling distance
  th_y <- ceiling(th_y/25)  # round for 25m patch cell
  th_y
}

th_y <- lapply(frfiles, foo)
print(th_y[1:n])
}


# LOAD ALL RESULT FILES
stats_frnew <- data.table(timestep=NA, elevation=NA, year=NA, ratio_bio_aa=NA, mass_bio_aa=NA, n_cohorts=NA)
for (t in 1:t)
t <- 1
print(paste("Timestep", t))
res_frnew <- as.list(paste("Simulations/dis_frnew/dis_frnew_alt",1:n, "_step", t,"/Output/fullOut_1.csv", sep=""))
res_frnew <- lapply(res_frnew, fread)
res_frnew <- lapply(res_frnew, rev_ycoordsDT, aui=extent(maps_list_frnew[[2]]))
for (i in 1:n) {
  DT <- res_frnew[[i]]
  DT[,bio_cohort:=biomass*stems,]
  DT[,.( bio_tot=sum(bio_cohort)), by=.(species, cell, xcoord, ycoord, elevation)]
  DT[,.(dom_species = )]
  print(levelplot(out2rasterDT(DT)))
  
  
}

### NICE PRINTS
ele <- c(2,4,6)
dec <- c(10,25,50)
gr <- expand.grid(ele,dec)
frnew_plot <- as.list(paste("Simulations/dis_frnew/dis_frnew_alt", gr[,1], "_step", gr[,2], "/Output/fullOut_1.csv" ,sep=""))
frnew_plot <- lapply(frnew_plot, fread)
frnew_plot <- lapply(frnew_plot, rev_ycoordsDT, aui=extent(maps_list_frnew[[2]]))


for (i in 1:length(frnew_plot))
lapply(frnew_plot, function(x) print(levelplot(out2rasterDT(x))))
