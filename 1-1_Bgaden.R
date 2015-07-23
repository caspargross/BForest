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



### Evaluate and plot Biomass
sumbio_bg <- data.table(elevation=NA, species=NA, bio_sp=NA, browsing=NA )  #Data.table with stats
for (i in 1:length(res_bg)) {
  BGDT <- res_bg[[i]]
  BGDT[,bio_cohort:=biomass*stems*4]
  BGDT <- BGDT[, .(bio_sp=sum(bio_cohort)), by=.(cell, elevation, species)]
  BGDT <- BGDT[, .(bio_sp=sum(bio_sp)), by=.(elevation, species)]
  BGDT[,browsing:=(i-1)*0.1,]
  sumbio_bg <- rbind(sumbio_bg, BGDT)
}

sumbio_bg <- sumbio_bg[-1,,]   
require(caTools) ## Calculate the Running Median
sumbio_rm <- sumbio_bg[, .(bio_sp=runmean(bio_sp, 15), elevation), by=.(species, browsing)]  
## Start Plotting
require(ggplot2)

## Plot the Biomass Elevation Gradient
bg_plot <- ggplot(sumbio_rm, aes(col=species))+ theme_bw()
bg_plot + geom_line (aes(x=elevation, y=bio_sp)) + facet_wrap(~ browsing, ncol=1)

## Plot Total Biomass ~ Browsing Intensity (Stacked Barplot)
sumbio <- sumbio_bg[,.(bio_sum=sum(bio_sp)), by=.(species, browsing)]
sumbio[,sp_sum:=sum(bio_sum), by=.(browsing)]
sumbio[,bio_perc:=(bio_sum/sp_sum)*100]  # Scale up to Percentage
setkey(sumbio, species)

bg_plot <- ggplot(sumbio, aes(fill=species))+ theme_bw() +labs(x="Browsing Intensity", y="Percentage of Total Biomass", fill="Species")
bg_plot + geom_bar (aes(x=as.factor(browsing), y=bio_perc), stat = "identity", position="stack")


### Evaluate and plot Reduction Factors

rf_bg<- data.table()
for (i in 1:length(res_bg)) {
  BGDT <- res_bg[[i]]
  BGDT <- BGDT[, .(idx=seq.int(nrow(BGDT)), cell, species, elevation, slowGrowth, age, stems, biomass, bio_cohort=biomass*stems, rf_light=light_rf,  rf_moisture=moisture_rf, rf_degreeDay=degreeDay_rf),]
  BGDT[,browsing:=(i-1)*0.1,]
  #BGDT[, is_zero := sum(rf_light, rf_moisture, rf_degreeDay)==0, by=.(idx)]
  BGDT[, lim_rf:=which(c(rf_light, rf_moisture, rf_degreeDay)==min(rf_light, rf_moisture, rf_degreeDay)), by=.(idx)]
  rf_bg <- rbind(rf_bg, BGDT)
}

rf_bg_factor <- rf_bg[, .(lim_rf=factor(lim_rf, levels= c(1,2,3), labels = c("Light", "Moisture", "Temperature")), elevation, species,]
rf_bg_factor <- rf_bg_factor[]

## Plot the reduction factors
rf_plot <- ggplot(rf_bg[browsing==0.0], aes(x=elevation, col=species, shape=species)) + theme_bw() + ylim(0, 1)+ labs(x="Elevation (m.a.s.l)",y="Reduction Factor", col="Species", shape="Species")
#RF Light
rf_plot + geom_point(aes( y=rf_light)) + ggtitle("Light")  # + geom_rug(aes(y=rf_light), alpha=0.4)
#RF Moisture
rf_plot + geom_point(aes( y=rf_moisture)) + ggtitle("Moisture")
#RF DegreeDays
rf_plot + geom_point(aes( y=rf_degreeDay)) + ggtitle("Temperature")

## Plot the limiting reduction factor
rf_plot_lim <- ggplot(rf_bg, aes(x=elevation))
rf_plot_lim + geom_point(aes(x=as.factor(elevation), y=lim_rf), stat="identity")                                                                                   

