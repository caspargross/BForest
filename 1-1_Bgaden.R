###########################
## 2D Elevation Gradient ##
##  Area: Berchtesgaden ##
###########################
setwd("/home/caspar/Bachelor_Thesis/BForest/")
library(ModellingTools)
library(sp)
library(rgdal)
library(raster)

## Create Gradient Model Landscape
gk_projection<-CRS("+init=epsg:31468") # GK Zone 4 (Berchtesgaden)
nr <-200  #corresponds to 5000m
nc <- 4   # corresponds to 100m
res <- 25   #25m grid cells
ex <- extent(0, nc*res, 0, nr*res)
dem_bg <- raster(nrows=nr, ncols=nc, ex)
projection(dem_bg) <- gk_projection
dem_bg

hmax <- 2000
hmin <- 500
gradient <- (hmax-hmin)/(nr*res)

dem_bg[,seq(nc)] <- rev(hmin+(seq(nr)*gradient*res))
#dem[] <- rev(1:(nc*nr)) # 
plot(dem_bg)

maps25_bg <- create_LandClim_Maps(dem_bg, no_aspect=T)

##  Create Full Model Input (without browsing)
create_inputdir("full_br0", climpath="Data/DWD/climate_wendelstein.dat", ctlfile="Data/Landclim/ctl_bgaden.xml", species=c("abiealba", "piceabie", "larideci", "fagusilv"), LandClimRasterStack=maps25_bg)
create_inputdir("full_br1", climpath="Data/DWD/climate_wendelstein.dat", ctlfile="Data/Landclim/ctl_bgaden.xml", species=c("abiealba", "piceabie", "larideci", "fagusilv"), LandClimRasterStack=maps25_bg)
create_inputdir("full_br2", climpath="Data/DWD/climate_wendelstein.dat", ctlfile="Data/Landclim/ctl_bgaden.xml", species=c("abiealba", "piceabie", "larideci", "fagusilv"), LandClimRasterStack=maps25_bg)
create_inputdir("full_br3", climpath="Data/DWD/climate_wendelstein.dat", ctlfile="Data/Landclim/ctl_bgaden.xml", species=c("abiealba", "piceabie", "larideci", "fagusilv"), LandClimRasterStack=maps25_bg)
create_inputdir("full_br4", climpath="Data/DWD/climate_wendelstein.dat", ctlfile="Data/Landclim/ctl_bgaden.xml", species=c("abiealba", "piceabie", "larideci", "fagusilv"), LandClimRasterStack=maps25_bg)
create_inputdir("full_br5", climpath="Data/DWD/climate_wendelstein.dat", ctlfile="Data/Landclim/ctl_bgaden.xml", species=c("abiealba", "piceabie", "larideci", "fagusilv"), LandClimRasterStack=maps25_bg)

mclapply(as.list(paste("full_br", 0:5, sep="")), function(x) run_landclim_model(x, ctl_file="ctl_bgaden.xml"), mc.cores=3)


### Load Output of Full Mode
res_bg <- lapply(as.list(paste("Simulations/full_br", 0:5,"/Output/fullOut_50.csv", sep="")), fread)
res_bg <- lapply(res_bg, function(x) rev_ycoordsDT(x, extent(dem_bg)))


sumbio_bg <- data.table(elevation=NA, species=NA, bio_sp=NA, rf_light=NA,  rf_moisture=NA, rf_degreeDay=NA, browsing=NA )  #Data.table with stats

for (i in 1:length(res_bg)) {
  BGDT <- res_bg[[i]]
  BGDT[,bio_cohort:=biomass*stems*4]
  BGDT <- BGDT[, .(bio_sp=sum(bio_cohort), rf_light=min(light_rf),  rf_moisture=min(moisture_rf), rf_degreeDay=min(degreeDay_rf)), by=.(cell, elevation, species)]
  BGDT <- BGDT[, .(bio_sp=sum(bio_sp), rf_light=mean(rf_light),  rf_moisture=mean(rf_moisture), rf_degreeDay=mean(rf_degreeDay)), by=.(elevation, species)]
  BGDT[,browsing:=(i-1)*0.1,]
  sumbio_bg <- rbind(sumbio_bg, BGDT)
}

sumbio_bg <- sumbio_bg[-1,,]
sumbio_rm <- sumbio_bg[, .(bio_sp=runmean(bio_sp, 15), elevation), by=.(species, browsing)]
## Start Plotting

require(ggplot2)
require(caTools)


bg_plot <- ggplot(sumbio_rm, aes(col=species))+ theme_bw()
bg_plot + geom_line (aes(x=elevation, y=bio_sp)) + facet_wrap(~ browsing, ncol=1)

sumbio<- sumbio_bg[,.(bio_sum=sum(bio_sp)), by=.(species, browsing)]
setkey(sumbio, species)
bg_plot <- ggplot(sumbio, aes(fill=species))+ theme_bw()
bg_plot + geom_bar (aes(x=browsing, y=bio_sum), stat = "identity", position="stack")



bg_bio_plot

xy <- sumbio_bg[species=="larideci"]
matplot(xy$elevation, xy$bio_sp, type="l")


a <- eleg$decade==50
matplot(eleg$elevation[a], eleg[a,colnames(eleg) %in% species], type="l", lty=lty,
plot_elevation_gradient(eleg, c("abiealba", "larideci"), 50)


eleg <
