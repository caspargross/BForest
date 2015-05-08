library(raster)
library(ModellingTools)
setwd("~/StudiProjekte/LandClim_SilverFir/SilverFir/Simulations/present")

list.files(paste(getwd(), "/Input/", sep=""))
inputpath <- paste(getwd(), "/Input/", sep="")
dem <- raster(paste(inputpath, "dem.asc", sep=""))
dem
x11()
plot(dem)

landscapeSize <- prod(dim(dem)) * prod(res(dem)) / (100*100)  ### ncell * cellsize/(100*100) -> in hectare
sp <- read.species.xml(paste(inputpath, "species.xml", sep=""))
sp
sp$name
figurepath <- paste("figures/")
dir.create(figurepath)      ### Create a new folder

### Elev Output: Biomass
result <- read.csv("Output/elevation_biomass_out.csv", strip.white=TRUE)

mainspecies <- sp$name
cols <- rainbow(length(mainspecies))
cols
lty <- 1

x11()
plotElevationGradient(elevationBiomassOut=result, species=mainspecies, selection=10, lty=1,  cols= rainbow(length(species)), plotlegend=TRUE)

### Read in "fullOutGrowth".
full <- read.csv("Output/fullOutGrowth_30.csv", strip.white=TRUE)   ### Forest after 300 Jears succession, Initialization "bare ground" and "global seed rain".
head(full)
summary(full)

###Task: Use function aggregate in order to get gap-biomasses of the species.
full$biomass.cohort <- full$biomass * full$stems
biom <- aggregate(full$biomass.cohort, list(species=full$species, row=full$row, col=full$col), sum)
head(biom)

### Plot maps in a pdf.
pdf(paste(figurepath, "maps.pdf", sep=""))
plot(dem, main="DEM")
for(s in 1:length(sp$name)) {
  biom.species <- biom[biom$species==sp$name[s], ]
  map <- matrix(0, nrow=max(biom$row), ncol=max(biom$col))
  for (i in 1:nrow(biom.species)) map[biom.species$row[i], biom.species$col[i]] <- biom.species$x[i]
  map <- raster(map)
  plot(map, col=rev(heat.colors(100)), main=sp$name[s])
}
dev.off()

### Running mean um den Höhengradienten zu glätten.
#install.packages("caTools")
library("caTools")

dec <- 40  ### Dekade 40 auswählen.
result <- read.csv("Output/elevation_biomass_out.csv", strip.white=TRUE)  ### Daten einlesen.
elev <- result[result$decade==dec,]  ### Nur die Daten einer Dekade nehmen.
elev <- apply(elev, 2, function(x) round(runmean(x, k=20),2))   ### Für alle Spalten gleichzeitig das running mean berechnen.
matplot(elev[,"elevation"], elev[,2:10], type="l")  ### Plotten und noch aufhübshcen (s.o.).

