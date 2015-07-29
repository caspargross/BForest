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

hmax <- 2200
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
sumbio_rm <- sumbio_bg[, .(bio_sp=runmean(bio_sp, 9), elevation), by=.(species, browsing)]  
## Start Plotting
require(ggplot2)

## Plot the Biomass Elevation Gradient (FULL)
bg_ele_full<- ggplot(sumbio_rm, aes(col=species))+
  theme_cas() +
  geom_line (aes(x=elevation, y=bio_sp), size=1.0) +
  facet_wrap(~ browsing, ncol=1) +
  labs(x="Elevation in m.a.s.l", y="Percentage of Total Biomass", col="Species") +
  scale_x_continuous("Elevation in m.a.s.l", breaks=br, labels=br)  
bg_ele_full 


## Plot the Biomass Elevation Gradient (ONE)
bg_ele<- ggplot(sumbio_rm[browsing==0.0,,], aes(col=species))+
  theme_cas() +
  geom_line (aes(x=elevation, y=bio_sp), size=1.0) +
  #facet_wrap(~ browsing, ncol=1) +
  labs(x="Elevation in m.a.s.l", y="Percentage of Total Biomass", col="Species") +
  scale_x_continuous("Elevation in m.a.s.l", breaks=br, labels=br)  
bg_ele 

## Plot Total Biomass ~ Browsing Intensity (Stacked Barplot)
sumbio <- sumbio_bg[,.(bio_sum=sum(bio_sp)), by=.(species, browsing)]
sumbio[,sp_sum:=sum(bio_sum), by=.(browsing)]
sumbio[,bio_perc:=(bio_sum/sp_sum)*100]  # Scale up to Percentage
setkey(sumbio, species)

bg_plot <- ggplot(sumbio, aes(fill=species))+ theme_bw() +labs(x="Browsing Intensity", y="Percentage of Total Biomass", fill="Species") +
geom_bar (aes(x=as.factor(browsing), y=bio_perc), stat = "identity", position="stack")
bg_plot

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


rf_bg_f <- rf_bg[, .(lim_rf=factor(lim_rf, levels= c(1,2,3), labels = c("Light", "Moisture", "Temperature")), elevation, species),]
setkey(rf_bg_f, elevation)
rf_bg_f <- rf_bg_f[, .(count=summary(as.factor(lim_rf)), rf_type=rep(c("Light", "Moisture", "Temperature"), length(unique(elevation)))), by=.(elevation, species)]
#rf_bg_f$count <- range01(rf_bg_f$count)
rf_bg_f[,rm_count:=runmean(count, 10), by=.(rf_type, species)]  



## Plot the reduction factors
rf_plot <- ggplot(rf_bg[browsing==0.0], aes(x=elevation, col=species, shape=species)) + theme_bw() + ylim(0, 1)+ labs(x="Elevation (m.a.s.l)",y="Reduction Factor", col="Species", shape="Species")
#RF Light
rf_plot + geom_point(aes( y=rf_light)) + ggtitle("Light")  # + geom_rug(aes(y=rf_light), alpha=0.4)
#RF Moisture
rf_plot + geom_point(aes( y=rf_moisture)) + ggtitle("Moisture")
#RF DegreeDays
rf_plot + geom_point(aes( y=rf_degreeDay)) + ggtitle("Temperature")

## Plot the limiting reduction factor


rf_plot_lim <- ggplot(rf_bg_f, aes(x=elevation,  colour=species, shape=species, group=interaction(species, rf_type))) +
  theme_cas() +
  #scale_colour_manual(values=cbbPalette) +
  scale_x_continuous("Elevation in m.a.s.l", breaks=br, labels=br)  +
  labs(y = "Occurence as limiting growth factor", colour="Species", shape="Species") +
  scale_shape_manual(values=c(21,22,23,24)) +
  geom_point(aes(y=count)) + 
  geom_line(aes(y=rm_count, k=7), size=1.1) + 
  facet_wrap(~rf_type, ncol=1)

rf_plot_lim


#### EXPORT THE IMAGES AS .pdf Files

#Elevation Gradient FULL
pdf(file ="Figures/2dgradient_full.pdf", width=10, height=10)
print(bg_ele_full)
dev.off()

# Elevation Gradient (Browsing =0)
pdf(file ="Figures/2dgradient.pdf", width=10, height=4)
print(bg_ele)
dev.off()

# Histogramm Biomass
pdf(file ="Figures/bg_plot.pdf", width=4, height=4)
print(bg_plot)
dev.off()

# Reduction Factors Plot
pdf(file ="Figures/rf_plot_lim.pdf", width=10, height=8)
print(rf_plot_lim)
dev.off()

