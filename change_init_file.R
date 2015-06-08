##########################
# Manipulate Input Files #
#                        #
###########################


type <- "mono_piceabie_dec40"
init <- read.csv(paste("Data/Init_State/", type, ".csv", sep=""))


raster_init <- raster_from_output(init, spec="piceabie", variable="elevation")
plot(is.na(raster_init))
plot(raster_init)

dat <- init
dat$biomass_cohort <- dat$biomass * dat$stems

ex <- extent(min(dat$xcoord), max(dat$xcoord), min(dat$ycoord), max(dat$ycoord))
r <- raster (ex, res=c(25,25), crs="+init=epsg:31467")
r<- rasterize(cbind(dat$xcoord, dat$ycoord), r, field=dat$elevation)
plot(r)

plot(dem_gk)
plot(extent(bluu), add=T, col="red")
plot(extent(raster_init), add=T, col="blue")




dat
dat<-subset(dat, dat$species==spec, select =  c("xcoord","ycoord", "elevation"))

raster <- rasterFromXYZ(dat, res=c(25,25), crs="+init=epsg:31467")
plot(raster)
plot(bluu)
dat
head(dat)
which(is.na(dat))

plot(dat$ycoord~dat$xcoord)