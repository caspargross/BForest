##########################
# Manipulate Input Files #
#                        #
###########################


type1 <- "mono_piceabie_dec50"
init1 <- read.csv(paste("Data/Init_State/", type, ".csv", sep=""))

type2 <-"mono_abiealba_dec50"
init2 <- read.csv(paste("Data/Init_State/", type, ".csv", sep=""))

### Reverse and fix YCoordinates
init1 <- rev_ycoords(init1)
init2 <- rev_ycoords(init2)

plot(out2raster(init1))
plot(out2raster(init1, var="biomass_cohort"))

### Create Mask for superposition
mask <- out2raster(init1)
mask[] <- NA
plot(mask)
plot()
######### Problemsuche
in1 <- raster("Simulations/3br/Input/dem.asc")
out1 <- read.csv("Simulations/3simple/Output/fullOut_4.csv")
ex_out1 <- extent(c(min(out1$xcoord), max(out1$xcoord), min(out1$ycoord), max(out1$ycoord)))
plot(dem_gk)
plot(extent(in1), add=T, col="red")
plot(r, add=T, col="blue")

extent(in1)
ex_out1


