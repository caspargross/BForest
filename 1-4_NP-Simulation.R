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
#load input files
in_mono <- fread ("Simulations/np_mono_pa/Output/fullOut_50.csv")
in_new <- fread ("Simulations/np_full_new/Output/fullOut_2.csv")
in_mixed <- fread ("Simulations/np_full_old/Output/fullOut_50.csv")

temp1 <- mask_apply(in_mono, in_mixed, m_mixed)
np_input1 <- mask_apply(temp1, in_new, m_full_new)
write.table(np_input1, file="Data/Init_State/np_input1.csv", row.names=F, col.names=T, sep=",", quote=F)
levelplot(out2rasterDT(np_input1))

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

# Analysis:
res_np1 <- as.list(paste("Simulations/np_disp_1/Output/fullOut_", 1:50,".csv", sep=""))
res_np1 <- lapply(res_np1, fread)
res_np2 <- res_np1

levelplot(out2rasterDT(res_np1[[6]]))

for (i in seq_along(res_np1)) {
  DT <- res_np2[[i]]
  DT[, bio_cohort:=biomass*stems, ]
  DT <- DT[,.(cell, xcoord, ycoord, elevation, species, age, biomass, bio_cohort)]
  DT <- DT[, .( max_age=max(age), min_age=min(age), sum_bio=sum(bio_cohort) ), by=.(cell, xcoord, ycoord, elevation, species)]
  DT[, cell_bio:=sum(sum_bio), by=cell]  # Total Biomass of Cells
  DT[, bio_ratio:=sum_bio/cell_bio, ]
  DT[species == "abiealba", th_aa:=(ifelse (bio_ratio>0.3, TRUE, FALSE)),  ]
  res_np1[[i]] <- DT
}


aa_th <- DT[species == "abiealba", .(xcoord, ycoord, th_aa=any(th_aa)), by=cell]
setkey(aa_th, xcoord, ycoord)

x <- unique(aa_th$xcoord)
y <- unique(aa_th$ycoord)
image(unique(aa_th$xcoord), unique(aa_th$ycoord), 1)

plot(x,y)

np_mono_pa <- fread("Simulations/np_mono_pa/Output/fullOut_50.csv")
np_mono_pa[,bio_cohort:=biomass*stems, by=.(cell)]
levelplot(out2rasterDT(np_mono_pa, var="bio_cohort"))


ggraster <- ggplot (DT[species == "abiealba",,], aes(x=xcoord, y=ycoord)) +
    coord_fixed(ratio=1) +
    theme_cas()+
    geom_tile(aes(fill=(bio_ratio)))
ggraster


for (i in 1:50) {
ggraster <- ggplot (res_np1[[i]], aes(x=xcoord, y=ycoord)) +
  coord_fixed(ratio=1) +
  theme_cas()+
  geom_point(aes(col=species, size=bio_ratio), alpha=0.1)
  
print(paste("No", i))
png(file=paste("Animate/np_1_dec", i,".png", sep=""), width=500 , height=400)
print(ggraster)
dev.off()
}


DT <- res_np1[[30]]








ra <- out2rasterDT(np_full_test)
levelplot(ra)
extent(ra)

http://sg.geodatenzentrum.de/wms_dop40?FORMAT=image%2Fjpeg&TRANSPARENT=FALSE&SERVICE=WMS&VERSION=1.1.1&REQUEST=GetMap&STYLES=&LAYERS=0&SRS=EPSG%3A25832&BBOX=441500.00000000,5382400.00000000,452500.00000000,5392400.00000000&WIDTH=1000&HEIGHT=1000

http://sg.geodatenzentrum.de/wms_dop40?FORMAT=image%2Fjpeg&TRANSPARENT=FALSE&SERVICE=WMS&VERSION=1.1.1&REQUEST=GetMap&STYLES=&LAYERS=0&SRS=EPSG%3A25832&BBOX=443063.87159437,5386418.9716511,445509.8564995,5388864.9565562&WIDTH=256&HEIGHT=256
plot(np_shp)

download.file("http://sg.geodatenzentrum.de/wms_dop40?FORMAT=image%2Fjpeg&TRANSPARENT=FALSE&SERVICE=WMS&VERSION=1.1.1&REQUEST=GetMap&STYLES=&LAYERS=0&SRS=EPSG%3A25832&BBOX=443063.87159437,5386418.9716511,445509.8564995,5388864.9565562&WIDTH=256&HEIGHT=256", "Data/np/np_bg.jpg")