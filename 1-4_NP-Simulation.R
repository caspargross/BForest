###################################
### Simulation: Black Forest NP ###
###################################

library(rgdal)
setwd("/home/caspar/Bachelor_Thesis/BForest/")


## load shapefiles
# NP Border Shapefile
ogrInfo("Data/np/np_border/np_border.shp", "np_border")
np_shp <- readOGR("Data/np/np_border/np_border.shp", "np_border")

#NP DEM Raster
dem_bawue <- raster("Data/np/dem_bawue_gk.tif")


#Check CRS
crs(np_shp)
crs(np_dem)

np_dem <- mask(dem_bawue, np_shp)
np_dem <- crop(np_dem, np_shp)
np_dem <- resample(np_dem, raster(res=c(25,25), ex=extent(np_dem), crs=crs(np_dem)))
# Conver in Spatial Pixels Dataframe for plotting in ggplot2
df.np_dem <- as(np_dem, "SpatialPixelsDataFrame")
df.np_dem<- as.data.frame(df.np_dem)
head(df.np_dem)

#Create hillshade
np_sl = terrain(np_dem, opt='slope')
np_as = terrain(np_dem, opt='aspect')
hill = hillShade(np_sl, np_as, 40, 270)


# Make LandClim Files
np_map25 <- create_LandClim_Maps(np_dem)

## FULL SPECIES (mixed forest) (500yr)
create_inputdir ("np_full_old",
                 climpath="Data/DWD/climate_feldberg.dat",
                 species=c("fagusilv", "piceabie", "abiealba"), ex=F,
                 LandClimRasterStack=np_map25,
                 inputfile=F,
                 ctlfile="Data/Landclim/ctl_bforest_50.xml",
                 landtypefile="Data/Landclim/landtype.xml")
run_landclim_model("np_full_old", ctl_file="ctl_bforest_50.xml")
np_full_old <- fread("Simulations/np_full_old/Output/fullOut_50.csv")
#np_full_test <- rev_ycoordsDT(np_full_old, extent(np_map25[[2]]))
levelplot(out2rasterDT(np_full_old))

## Full species (new 20 year)
create_inputdir ("np_full_new",
                 climpath="Data/DWD/climate_feldberg.dat",
                 species=c("piceabie", "fagusilv", "abiealba"), ex=F,
                 LandClimRasterStack=np_map25,
                 inputfile=F,
                 ctlfile="Data/Landclim/ctl_bforest_2.xml",
                 landtypefile="Data/Landclim/landtype.xml")
run_landclim_model("np_full_new", ctl_file="ctl_bforest_2.xml")
np_full_new <- fread("Simulations/np_full_new/Output/fullOut_2.csv")
np_full_new[,bio_cohort:=biomass*stems, by=.(cell)]
levelplot(out2rasterDT(np_full_new, var="bio_cohort"))

## Picea abies (monoculture)
create_inputdir ("np_mono_pa",
                 climpath="Data/DWD/climate_feldberg.dat",
                 species=c("piceabie"), ex=F,
                 LandClimRasterStack=np_map25,
                 inputfile=F,
                 ctlfile="Data/Landclim/ctl_bforest_50.xml",
                 landtypefile="Data/Landclim/landtype.xml")
run_landclim_model("np_mono_pa", ctl_file="ctl_bforest_50.xml")
np_mono_pa <- fread("Simulations/np_mono_pa/Output/fullOut_50.csv")
np_mono_pa[,bio_cohort:=biomass*stems, by=.(cell)]
levelplot(out2rasterDT(np_mono_pa, var="bio_cohort"))


#np_full_test <- rev_ycoordsDT(np_full_test, extent(np_map25[[2]]))
levelplot(out2rasterDT(np_full_test))

## Read layer files
np_grinden <- readOGR("Data/np/grinden_gk.shp", "grinden_gk")
np_new <- readOGR("Data/np/new_gk.shp", "new_gk")
np_mixed <- readOGR("Data/np/mixed_gk.shp", "mixed_gk")

#Test plot
plot(np_shp, col="darkgreen", axes=T)
plot(np_grinden, add=T, col="blue")
plot(np_new, add=T, col="lightblue")
plot(np_mixed, add=T, col="lightgreen")
plot(hill, col=grey(0:100/100), alpha=0.35, add=T, legend=F)

### Create masks for different forest types
r1 <- np_dem
r1[!is.na(r1[])] <- 1
m_grinden <- mask(r1, np_grinden)
m_new <- mask(r1, np_new)
m_mixed <-  mask(r1, np_mixed)
m_new <- m_new - !is.na(m_mixed)
m_mixed[m_mixed[]==0] <- NA
m_full_new <- (m_grinden | m_new)
m_mono <- r1 - !is.na (m_new | m_mixed | m_grinden)
m_mono[m_mono[]==0] <- NA

plot(m_mono, col="darkgreen", legend=F)
plot(m_mixed, col="lightgreen", add=T, legend=F)
plot(m_full_new, col="yellow", add=T, legend=F)
 # make combination raster:
np_m_aui <- as.factor(r1)
np_m_aui[!is.na(m_mono[])] <-  1
np_m_aui[!is.na(m_full_new[])] <-  2
np_m_aui[!is.na(m_mixed[])] <-  3
ratify(np_m_aui)
levels(np_m_aui) <- data.frame(ID = 1:3, type=c("Picea abie mono", "young mixed", "old mixed"))

levelplot(np_m_aui, par.settings=PuOrTheme)

## Apply masks to create initial files
## load input files
in_mono <- fread ("Simulations/np_mono_pa/Output/fullOut_50.csv")
in_new <- fread ("Simulations/np_full_new/Output/fullOut_2.csv")
in_mixed <- fread ("Simulations/np_full_old/Output/fullOut_50.csv")

temp1 <- mask_apply(in_mono, in_mixed, m_mixed)
np_input1 <- mask_apply(temp1, in_new, m_full_new)
write.table(np_input1, file="Data/Init_State/np_input1.csv", row.names=F, col.names=T, sep=",", quote=F)
levelplot(out2rasterDT(np_input1))


# B R O W S I N G = 0.%
#====================================================================================================================
### Create Simulation with input file. 
create_inputdir ("np_disp_1",
                 climpath="Data/DWD/climate_feldberg.dat",
                 species=c("piceabie", "abiealba", "fagusilv"), ex=F,
                 LandClimRasterStack=np_map25,
                 inputfile="Data/Init_State/np_input1.csv",
                 ctlfile="Data/Landclim/ctl_bforest_np.xml",
                 landtypefile="Data/Landclim/landtype.xml",
                 )
run_landclim_model("np_disp_1", ctl_file="ctl_bforest_np.xml")

# INPUT Analysis
inp_np <- fread("Data/Init_State/np_input1.csv")
inp_np <- rev_ycoordsDT(inp_np, aui=extent(np_dem))
inp_np <- flipy(inp_np)
inp_np[, bio_cohort:=biomass*stems, ]
inp_np <- inp_np[,.(cell, xcoord, ycoord, elevation, species, age, biomass, bio_cohort)]
inp_np <- inp_np[, .( max_age=max(age), min_age=min(age), sum_bio=sum(bio_cohort) ), by=.(cell, xcoord, ycoord, elevation, species)]
inp_np[, cell_bio:=sum(sum_bio), by=cell]  # Total Biomass of Cells
inp_np[, bio_ratio:=sum_bio/cell_bio, ]
inp_np[species == "abiealba", th_aa:=(ifelse (bio_ratio>0.3, TRUE, FALSE)),  ]
inp_np[, year:=0 ,]

#  RESULT Analysis:
res_np1 <- as.list(paste("Simulations/np_disp_1/Output/fullOut_", 1:50,".csv", sep=""))
res_np1 <- lapply(res_np1, fread)
res_np1 <- lapply(res_np1, function(x) setkey(x, col, row))
res_np1 <- lapply(res_np1, rev_ycoordsDT, aui=extent(np_dem))
res_np1 <- lapply(res_np1, flipy) #### SLOOOW (no data.table operation)


for (i in seq_along(res_np1)) {
  DT <- res_np1[[i]]
  DT[, bio_cohort:=biomass*stems, ]
  DT <- DT[,.(cell, xcoord, ycoord, elevation, species, age, biomass, bio_cohort)]
  DT <- DT[, .( max_age=max(age), min_age=min(age), sum_bio=sum(bio_cohort) ), by=.(cell, xcoord, ycoord, elevation, species)]
  DT[, cell_bio:=sum(sum_bio), by=cell]  # Total Biomass of Cells
  DT[, bio_ratio:=sum_bio/cell_bio, ]
  DT[species == "abiealba", th_aa:=(ifelse (bio_ratio>0.3, TRUE, FALSE)),  ]
  DT[,year:= i*10,]
  res_np1[[i]] <- DT
}


######################
#  PLOT THE RESULTS  #
######################
#Combine all files for plot in single data.table
npyr <- c(3, 5, 10, 15, 20) # decades to plot
DT_NP_plot <- data.table(cell=integer(), xcoord=integer(), ycoord=integer(), elevation=integer(), species=character(), max_age=integer(), min_age=integer(), sum_bio=numeric(), cell_bio=numeric(), bio_ratio=numeric(), th_aa=logical(), year=integer())
DT_NP_plot <- rbind(DT_NP_plot, inp_np)
for (i in npyr) {DT_NP_plot <- rbind(DT_NP_plot, res_np1[[i]])}

# plot
np_plot <- ggplot(DT_NP_plot, aes(x=xcoord, y=ycoord)) +
  theme_cas() +
  theme(legend.position = "bottom")+
  theme(axis.text.y = element_text(angle = 90, hjust=0.5)) +  #Flip ycoordinates
  coord_fixed() +
  geom_tile(aes(fill=species),  data= DT_NP_plot[species == "fagusilv"], fill="grey80")+
  geom_tile(aes(fill=species),  data= DT_NP_plot[species == "piceabie"], fill="wheat1")+
  #geom_tile(aes(fill=bio_ratio), data=DT_NP_plot[species=="abiealba" & bio_ratio > 0.1]) +
  geom_tile(aes(fill=bio_ratio), data=DT_NP_plot[species=="abiealba"]) +
  scale_fill_gradient(low="lightgreen", high="darkgreen", limits=c(0, 1), na.value="transparent")+
  stat_contour(data=df.np_dem, aes(x=x, y=y, z=dem_bawue_gk), col="grey20", size = 0.1) +
  labs( y="Longitude", x="Latitude", fill="Biomass ratio \n A.alba / Total") +
  facet_wrap(~year, ncol=2)
 

np_plot

res_np1[[]][]
#Convert List with Results into raster format.
raster_np1 <- lapply(res_np1, function(x) out2rasterDT(x, var="bio_ratio"))
raster_np2 <- lapply(raster_np1, function(x) focal(x, w=matrix(1,3,3), fun=sum))

levelplot(raster_np1[[3]])
levelplot(raster_np2[[3]])


#### Links for the aerial photographs
http://sg.geodatenzentrum.de/wms_dop40?FORMAT=image%2Fjpeg&TRANSPARENT=FALSE&SERVICE=WMS&VERSION=1.1.1&REQUEST=GetMap&STYLES=&LAYERS=0&SRS=EPSG%3A25832&BBOX=441500.00000000,5382400.00000000,452500.00000000,5392400.00000000&WIDTH=1000&HEIGHT=1000
http://sg.geodatenzentrum.de/wms_dop40?FORMAT=image%2Fjpeg&TRANSPARENT=FALSE&SERVICE=WMS&VERSION=1.1.1&REQUEST=GetMap&STYLES=&LAYERS=0&SRS=EPSG%3A25832&BBOX=443063.87159437,5386418.9716511,445509.8564995,5388864.9565562&WIDTH=256&HEIGHT=256