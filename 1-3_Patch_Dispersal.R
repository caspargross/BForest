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


sink("Logs/log_dis_patch.txt")
for (k in 1:1) {
print(paste("Patch pattern no", k, "startedt at:", Sys.time()))
  
initlist1_p <- as.list( paste("Simulations/dis_patch_init_fs-pa_", seq_along(maps_list_p), "/Output/fullOut_50.csv", sep=""))
initlist2_p <- as.list( paste("Simulations/dis_patch_init_aa_", seq_along(maps_list_p), "/Output/fullOut_50.csv", sep=""))
outputlist_p<- as.list( paste("Data/Init_State/dis_patch_init_alt_", seq_along(maps_list_p),"_m", k, ".csv", sep=""))

print("Loading Files [OK]")

combine_input (initlist1_p,  
               initlist2_p, 
               outputlist_p,
               1, 1, 2, 100,  
               extent(maps_list_p[[1]]),
               ma=eval(as.name(paste("m",k, sep=""))))

print("Creating Input File  [OK]")
## Check Input Files

#input_files <- as.list(paste("Data/Init_State/dis_patch_init_alt_", 1:8,"_m11.csv", sep=""))
#input_files <- lapply(input_files, fread)
#input_files <- lapply(input_files, function(x) rev_ycoordsDT(x, extent(maps_list_p[[1]])))
#levelplot(out2rasterDT(input_files[[4]], var="species"))

## Run Model with Dispersal

for (i in seq_along (maps_list)) { create_inputdir (paste("dis_patch_alt_", i,"_m",k, sep=""),
                                                    climpath="Data/DWD/climate_feldberg.dat",
                                                    species=c("abiealba", "fagusilv", "piceabie"), ex=F,
                                                    LandClimRasterStack=maps_list_p[[i]],
                                                    paste("Data/Init_State/dis_patch_init_alt_", i,"_m", k, ".csv", sep=""),
                                                    ctlfile="Data/Landclim/ctl_bforest_dis_2000.xml",
                                                    landtypefile="Data/Landclim/landtype.xml") }
print(paste("Created Input Directory for Pattern", k, "Elevation", i, "[OK]"))
mclapply(as.list(paste("dis_patch_alt_", seq_along(maps_list), "_m", k, sep="")), function(x) run_landclim_model(x, ctl_file="ctl_bforest_dis_2000.xml"), mc.cores=3)
print(paste("Finished Model Calculation Pattern", k, "Elevation", i, "[OK]"))
}
print("SUCCESS!")
sink()



## Check Dispersal Output files
output_files <- as.list(paste("Simulations/dis_patch_alt_1_m1/Output/fullOut_", seq(5,200,5), ".csv", sep=""))
output_files <- lapply(output_files, fread)
output_files <- lapply(output_files, function(x) rev_ycoordsDT(x,  aui_rev=extent(maps_list_p[[1]])))
for (i in 1:length(output_files)) {
  png(filename=paste("Animate/dis_patch_alt1_m1_dec", i,".png", sep=""))
  print(levelplot(out2rasterDT(output_files[[25]], var="species")))
  dev.off()
}


######### Calculation of Output-Analysis 
stats_p <- data.table(mask=NA, altitude=NA, year=NA, ratio_bio_aa=NA, mass_bio_aa=NA, n_cohorts=NA)  #Data.table with stats
plot <- T
for (m in 1:1) {
  
  for (j in seq_along (maps_list_p)) {  ## Different Alitudes!
  
    list_results_p <- as.list (paste("Simulations/dis_patch_alt_", j,"_m",m,"/Output/fullOut_", seq(5, 200, 5), ".csv", sep=""))
    dis_results_p <- lapply(list_results_p, fread)
    dis_results_p <- lapply(dis_results_p, function(x) rev_ycoordsDT(x, extent(maps_list_p[[1]])))
    print(paste("Loaded map",j, "Patch Pattern no:", m))
     
  
    for (i in 1:length(dis_results)) { ## Different Decades !
    
      PDT<-dis_results_p[[i]]
      PDT[,bio_cohort:=biomass*stems] #Calculate cohort biomas
      PDT <- PDT[,.(bio_cohort=sum(bio_cohort), age=max(age)),by=list(species, cell, xcoord, ycoord)]
      setkey(PDT,species)
      PDT[,bio_aa:=as.double(NA),]  ## Create empty columnt (NA) in double type
      PDT["abiealba", bio_aa:=bio_cohort, by=.(cell)]
      PDT2 <- PDT[,.(bio_cell=sum(bio_cohort), bio_aa), by=.(xcoord, ycoord, cell)]
      PDT2[is.na(PDT2)] = 0  # replace NA with 0
      PDT2 <- PDT2[,.(bio_cell=mean(bio_cell), bio_aa=max(bio_aa), bio_quo=max((bio_aa/bio_cell)*100)), by=.(xcoord, ycoord, cell)]
    
      stats_p <- rbind(stats_p, list(m, list_alt_p[[j]], i*50, PDT2[,round((sum(bio_aa)/sum(bio_cell))*100,2)], PDT2[,sum(bio_aa),], nrow(PDT["abiealba", .(age=max(age)) , by="cell"][age>=70,,])))
      print(paste("finished decade", i))
    
      if (plot == T) {
        png(filename=paste("Animate/dis_patch_alt",j,"_m1_dec", i,".png", sep=""))
        print(levelplot(out2rasterDT(PDT2, var="bio_aa")))
        dev.off()
      
      }
    }
  }
}    
    
stats_p <- stats_p[-1,]
stats_p$altitude <- as.factor(stats_p$altitude)
stats_p$mask <- as.factor(stats_p$mask)
write.table(stats_p, "Data/Results/results_patch_dispersal.txt", sep="\t", row.names=F)

stats_p[, max_ratio:=0.7*max(ratio_bio_aa), by=.(mask, altitude)]



result_patch <- stats_p[altitude!=1800,,]

### Plot the biomass ratios
library(ggplot2)
distplot <- ggplot(result_patch, aes(x=year, group=mask, col = mask)) + theme_bw()
distplot + geom_line(aes(y=ratio_bio_aa)) + labs(title = "Mask 1", y="Biomass percentage Abies alba / Total") + facet_wrap(~ altitude, ncol=1)


### Plot the threshhold times
th_time <- stats_p[ratio_bio_aa>=max_ratio, .(min(year)) , by=.(mask, altitude, max_ratio)]
th_time <- th_time[altitude!=1800,,]

timeplot <- ggplot(th_time, aes(x=V1, group=mask, col=mask)) + theme_bw()
timeplot + geom_point(aes(y=max_ratio)) + facet_wrap(~ altitude, ncol=1)


#create subset for logistic regression
logtest <- result_patch[altitude==1200, , ]
logtest <- logtest[mask==4,,]
plot(log(logtest$ratio_bio_aa/100) ~ logtest$year)
logfit <- lm(log(ratio_bio_aa/100) ~ year, data=logtest)
abline(logfit)

stats_p$altitude
stats_p$year
stats_p$ratio_bio_aa

library(rgl)
plot3d(stats_p$altitude,
       stats_p$year,
       stats_p$ratio_bio_aa,
       col=stats_p$mask)

