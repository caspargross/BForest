###########################
## 2D Elevation Gradient ##
##  Area: Berchtesgaden ##
###########################
setwd("/home/caspar/Bachelor_Thesis/BForest/")
library(ModellingTools)
library(sp)
library(rgdal)
library(raster)
library(caTools) ## Calculate the Running Median
library(ggplot2)

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
hmin <- 400
gradient <- (hmax-hmin)/(nr*res)

dem_bg[,seq(nc)] <- rev(hmin+(seq(nr)*gradient*res))
#dem[] <- rev(1:(nc*nr)) # 
plot(dem_bg)

maps25_bg <- create_LandClim_Maps(dem_bg, no_aspect=T)

##  Create Full Model Input (without browsing)

bg_list <- as.list(paste("full_br", 0:10, sep=""))

lapply(bg_list, function (X) create_inputdir(X, 
                                            climpath="Data/DWD/climate_wendelstein.dat", 
                                            ctlfile="Data/Landclim/ctl_bgaden.xml", 
                                            species=c("abiealba", "piceabie", "larideci", "fagusilv"), 
                                            LandClimRasterStack=maps25_bg))


mclapply(bg_list, function(x) run_landclim_model(x, ctl_file="ctl_bgaden.xml"), mc.cores=3)



### Load Output of Full Mode
res_bg <- lapply(as.list(paste("Simulations/full_br", 0:5,"/Output/fullOut_50.csv", sep="")), fread)
res_bg <- lapply(res_bg, function(x) rev_ycoordsDT(x, extent(dem_bg)))



### Analyse Biomass

sumbio_bg <- data.table(elevation=NA, species=NA, bio_sp=NA, browsing=NA )  #Data.table with stats
for (i in 1:length(res_bg)) {
  BGDT <- res_bg[[i]]
  BGDT[,bio_cohort:=biomass*stems*4]
  BGDT <- BGDT[, .(bio_sp=sum(bio_cohort)), by=.(cell, elevation, species)]
  BGDT <- BGDT[, .(bio_sp=sum(bio_sp)), by=.(elevation, species)]
  BGDT[,browsing:=(i-1)*0.1,]
  sumbio_bg <- rbind(sumbio_bg, BGDT)
}

sumbio_bg <- sumbio_bg[-1,,]  # Delete first (empty) line
setkey(sumbio_bg, elevation, species, browsing)

sumbio_rm <- sumbio_bg
p.full <- data.table(expand.grid(elevation=unique(sumbio_rm$elevation), species=unique(sumbio_rm$species), browsing=unique(sumbio_rm$browsing)))
setkey(p.full, elevation, species, browsing)
setkey(sumbio_rm, elevation, species, browsing)
sumbio_rm <- merge(sumbio_rm, p.full, all.y=T)
sumbio_rm$bio_sp[is.na(sumbio_rm$bio_sp)]<-0
sumbio_rm[,bio_rm:=runmean(bio_sp,5), by=.(species,browsing)]

## Plotting ##
##############

# Colour Palettes:
cas_palette <- c("#E55934","#9BC53D", "#5BC0EB", "#FDE74C")
cas_palette2 <- c("#669900", "#0072B2", "#D55E00", "#FA7921"   )
# Breaks
br <- seq(400, 2000, 200)
## Plot the Biomass Elevation Gradient (Stacked Areas)
bg_ele_stack<- ggplot(sumbio_rm, aes(group=species, fill=species))+
  theme_cas() +
  geom_area(aes(x=elevation, y=bio_rm), position="stack") +
  geom_line(aes(x=elevation, y=bio_rm), colour="grey30", position = "stack") +
  facet_wrap(~ browsing, ncol=1) +
  labs(x="Elevation in m.a.s.l", y="Total biomass in t/ha", fill="Species") +
  scale_x_continuous("Elevation in m a.s.l", breaks=br, labels=br) +
  scale_fill_manual(values=cas_palette)
  
bg_ele_stack 


## Plot the Biomass Elevation Gradient (FULL)
bg_ele_full<- ggplot(sumbio_rm, aes(col=species))+
  theme_cas() +
  geom_line (aes(x=elevation, y=bio_rm), size=1.0) +
  facet_wrap(~ browsing, ncol=1) +
  labs(x="Elevation in m.a.s.l", y="Biomass in t/ha", col="Species") +
  scale_x_continuous("Elevation in m.a.s.l", breaks=br, labels=br)+
  scale_colour_manual(values=cas_palette)
bg_ele_full 


## Plot the Biomass Elevation Gradient (ONE)
bg_ele<- ggplot(sumbio_rm[browsing==0.0,,], aes(col=species))+
  theme_cas() +
  scale_color_brewer(palette="Set1")+
  geom_line (aes(x=elevation, y=bio_sp), size=1.0) +
  #facet_wrap(~ browsing, ncol=1) +
  labs(x="Elevation in m.a.s.l", y="Biomass in t/ha", col="Species") +
  scale_x_continuous("Elevation in m.a.s.l", breaks=br, labels=br)  
bg_ele 

## Plot Total Biomass ~ Browsing Intensity (Stacked Barplot)
sumbio <- sumbio_bg[,.(bio_sum=sum(bio_sp)), by=.(species, browsing)]
sumbio[,sp_sum:=sum(bio_sum), by=.(browsing)]
sumbio[,bio_perc:=(bio_sum/sp_sum)*100]  # Scale up to Percentage
setkey(sumbio, species)

bg_plot <- ggplot(sumbio, aes(fill=species))+ theme_bw() +labs(x="Browsing Intensity", y="Percentage of Total Biomass", fill="Species") +
scale_fill_manual(values=cas_palette)+
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


rf_bg_f <- rf_bg[, .(lim_rf=factor(lim_rf, levels= c(1,2,3), labels = c("Light", "Moisture", "Temperature")), elevation, species, stems, biomass, age=as.double(age)),] #convert limits to factor
setkey(rf_bg_f, elevation)
stemtable <- rf_bg_f[, .(stemsum=sum(stems)), by=.(elevation, lim_rf, species)]
names(stemtable)[names(stemtable)[]=="lim_rf"] <- "rf_type"
rf_bg_f <- rf_bg_f[,.(count=summary(lim_rf), med_age=median(age), rf_type=rep(c("Light", "Moisture", "degreeDays"))), by=.(elevation, species)]
# Merge tables
rf_bg_f <- merge(stemtable, rf_bg_f, by=c("species", "elevation", "rf_type" ), all.y=T)
rf_bg_f[is.na(stemsum), stemsum:=0, ]
rf_bg_f[, stemcount:=count*stemsum]


blub <- rf_bg[, .(mlight=median(rf_light), mmoisture=median(rf_moisture), mdegreedays=median(rf_degreeDay)), .(elevation, species)]
rf_bg_f<- merge(rf_bg_f, blub, by=c("elevation", "species"))

#Normalize
rf_bg_f[, count01:=count/(sum(count))*100, by=.(elevation, species)]
rf_bg_f[, stemcount01:=stemcount/(sum(stemcount)), by=.(elevation, species)]
#Calculate Running Mean
rf_bg_f[,rm_count01:=runmean(count01, 10), by=.(rf_type, species)]  
rf_bg_f[,rm_stemcount01:=runmean(stemcount01, 10), by=.(rf_type, species)]  
rf_bg_f[,sum_stemsum:=sum(stemsum), by=.(elevation, species)]
rf_bg_f[,sum_stemsum01:=range01(stemsum), by=.(species)]
## Insert NAs instead of 0 red.factor
rf_bg_f[mlight == 0, mlight:=NA]
rf_bg_f[mdegreedays == 0, mdegreedays:=NA]
rf_bg_f[mmoisture == 0, mmoisture:=NA]
#### Plot the limiting reduction factor


rf_plot_lim <- ggplot(rf_bg_f, aes(x=elevation)) +
  theme_cas_big() +
  scale_colour_manual(values=cas_palette2) +
  scale_x_continuous("Elevation in m a.s.l", breaks=br, labels=br)  +
  labs(y = "Occurence as limiting growth factor (%)", colour="Factor", shape="Factor") +
  scale_shape_manual(values=c(21,22,23,24)) +
  geom_point(aes(y=count01,  colour=rf_type, shape=rf_type, group=interaction(rf_type, species))) + 
  geom_line(aes(y=rm_count01,  colour=rf_type, shape=rf_type, group=interaction(rf_type, species)), size=0.8) +
  facet_wrap(~species, ncol=1)

rf_plot_lim


rf_plot_limvalues <- ggplot(rf_bg_f, aes(x=elevation)) +
  theme_cas_big() +
  scale_colour_manual(values=cas_palette2) +
  scale_x_continuous("Elevation in m a.s.l", breaks=br, labels=br)  +
  labs(y = "Occurence as limiting growth factor (%)", colour="Factor", shape="Factor") +
  scale_shape_manual(values=c(21,22,23,24)) +
  #geom_point(aes(y=count01,  colour=rf_type, shape=rf_type, group=interaction(rf_type, species))) + 
  #geom_line(aes(y=rm_count01,  colour=rf_type, shape=rf_type, group=interaction(rf_type, species)), size=0.8) +
  geom_line(aes(y=mlight), col=cas_palette2[1] )+
  geom_line(aes(y=mdegreedays), col=cas_palette2[3]) +
  geom_line(aes(y=mmoisture), col=cas_palette2[2]) +
  facet_wrap(~species, ncol=1)

rf_plot_limvalues


rf_plot_limvalues <- ggplot(rf_bg_f, aes(x=elevation)) +
  theme_cas_big() +
  scale_colour_manual(values=cas_palette2) +
  scale_x_continuous("Elevation in m a.s.l", breaks=br, labels=br)  +
  labs(y = "Occurence as limiting growth factor (%)", colour="Factor", shape="Factor") +
  scale_shape_manual(values=c(21,22,23,24)) +
  #geom_point(aes(y=count01,  colour=rf_type, shape=rf_type, group=interaction(rf_type, species))) + 
  #geom_line(aes(y=rm_count01,  colour=rf_type, shape=rf_type, group=interaction(rf_type, species)), size=0.8) +
  geom_line(aes(y=med_age), col=cas_palette2[4] )+
  facet_wrap(~species, ncol=1)

rf_plot_limvalues

#### Calculate analytics values
rf_bg[, .(total_sum=sum(stems), med_age=as.integer(median(age)), med_biomass=median(biomass)), by=.(species)] 



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

#############################
# Plot aerial photo of the elevation gradient
############################

### PLot logistic growth functtion
curve((0.08*x)(1-))


