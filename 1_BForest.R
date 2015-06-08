# 1_BForest
#############################################
## Interaction Abies Alba <-> Picea Abies  ##
##           Black Forest NP               ##
##            Caspar Gross                 ##
#############################################
setwd("/home/caspar/Bachelor_Thesis/BForest/")
library(raster)
library(devtools)

#install_github("KDolos/LandClimTools")
library(LandClimTools)


#Load Digital Elevation Model for National Park Nördlicher SW
dem_utm <- raster ("Data/dem/n48_e008_1arc_v3.bil")
np_extent <- extent (c(8.15,8.4,48.5,48.7))
dem_utm <-crop(dem_utm, np_extent)
plot(dem_utm)
points(8.276757, 48.578119) #Riesenköpfle UTM

#Convert to GK3
dem_gk <- projectRaster(dem_utm, crs="+init=epsg:31467")  # Project to GK Zone 3
np_gk_extent <- extent (c(3437400, 3455900, 5374500, 5395500))
dem_gk <- crop(dem_gk, np_gk_extent)
plot(dem_gk)
points(3446740, 5382515) #Riesenköpfle GK

writeRaster(dem_gk, file="Data/dem/dem_overview.grd")
dem_gk <- raster("Data/dem/dem_overview.grd")

#define Extent area of interest
aui <- extent( c(3446000, 3447500, 5381500, 5383500))
plot (dem_gk)
plot(aui, add=T)

source("3_makeLandClimMaps.R")
maps25<-create_LandClim_Maps(dem_gk)

#### TEST
create_inputdir("testbf", species=c("abiealba", "fagusilv"), ex=aui)
run_landclim_model("testbf")
bio_out_testbf <- read.csv("Simulations/testbf/Output/elevation_biomass_out.csv")
plot_elevation_gradient(bio_out_testbf, c("abiealba", "fagusilv"))


xytree <- tree_coordinates("Simulations/testbf/Output/fullOut_30.csv", silent=F)
colours <- c(rgb(0,1,0.3,0.3),rgb(0.3,0.2,0.6,0.3))
plot_forest(xytree, c("abiealba","fagusilv"), colours)



#### 3 Species, but without browsing.
create_inputdir("3simple", species=c("abiealba", "fagusilv", "piceabie"), ex=aui)
run_landclim_model("3simple")

ele3simple<-read.csv("Simulations/3simple/Output/elevation_biomass_out.csv")
plot_elevation_gradient(ele3simple, species=c("abiealba", "piceabie","fagusilv"))


r_aa_3simple<-raster_from_output(simname="3simple", decade="40", spec="abiealba", variable="biomass_cohort")
r_fs_3simple<-raster_from_output(simname="3simple", decade="40", spec="fagusilv", variable="biomass_cohort")
r_pa_3simple<-raster_from_output(simname="3simple", decade="40", spec="piceabie", variable="biomass_cohort")

plot(r_aa_3simple)
plot(r_fs_3simple)
plot(r_pa_3simple)

make_3d_plot(r_aa_3simple)
make_3d_plot(r_fs_3simple)
make_3d_plot(r_pa_3simple)

r_aa_3simple_dd<-raster_from_output(simname="3simple", decade="40", spec="abiealba", variable="degreeDay_rf")
r_fs_3simple_dd<-raster_from_output(simname="3simple", decade="40", spec="fagusilv", variable="degreeDay_rf")
r_pa_3simple_dd<-raster_from_output(simname="3simple", decade="40", spec="piceabie", variable="degreeDay_rf")

plot(r_aa_3simple_dd)
plot(r_fs_3simple_dd)
plot(r_pa_3simple_dd)

r_aa_3simple_mo<-raster_from_output(simname="3simple", decade="40", spec="abiealba", variable="light_rf")
r_fs_3simple_mo<-raster_from_output(simname="3simple", decade="40", spec="fagusilv", variable="light_rf")
r_pa_3simple_mo<-raster_from_output(simname="3simple", decade="40", spec="piceabie", variable="light_rf")

plot(r_aa_3simple_mo)
plot(r_fs_3simple_mo)
plot(r_pa_3simple_mo)


#### 3 Species with browsing 0-3
create_inputdir("3br", species=c("abiealba", "fagusilv", "piceabie"), ex=aui)
run_landclim_model("3br")


r_aa_3br_mo<-raster_from_output(simname="3br", decade="40", spec="abiealba", variable="light_rf")
r_fs_3br_mo<-raster_from_output(simname="3br", decade="40", spec="fagusilv", variable="light_rf")
r_pa_3br_mo<-raster_from_output(simname="3br", decade="40", spec="piceabie", variable="light_rf")

ele<-read.csv("Simulations/3br/Output/elevation_biomass_out.csv")
plot_elevation_gradient(ele, species=c("abiealba", "piceabie","fagusilv"))


plot(r_aa_3br_mo)
plot(r_fs_3br_mo)
plot(r_pa_3br_mo)

#### Silver Fir Monoculture
create_inputdir("mono_piceabie", species=c("piceabie"), ex=aui)
run_landclim_model("mono_piceabie")
ele_mp<-read.csv("Simulations/mono_piceabie/Output/elevation_biomass_out.csv")
plot_elevation_gradient(ele_mp, species=c("piceabie"), selection=43)


#### Abie alba Monoculture
create_inputdir("mono_abiealba", species="abiealba", ex=aui)
run_landclim_model("mono_abiealba")
ele_aa <-read.csv("Simulations/mono_abiealba/Output/elevation_biomass_out.csv")
plot_elevation_gradient(ele_aa, species=c("abiealba"), selection=43)



