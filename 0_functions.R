########################
## Functions LandClim ##
########################

run_landclim_model<-function(sim_name, ctl_file="ctl_bforest.xml", lcpath="/Data/Landclim/LandClim"){
  oldwd <- getwd()
  dir.create(paste("Simulations/",sim_name,"/Output",sep=""))
  dir.create(paste("Simulations/",sim_name,"/Output",sep=""))   ### Erstellt das Verzeichnis, in das die Ergebnisse geschrieben werden. Es muss so heißen, wie im LandClim-Controll-File angegeben.
  setwd(paste("Simulations/",sim_name,"/Input",sep="")) ### Damit der "system"-Befehl funktioniert, muss in R das working directory der Ordner "Input" der entsprechenden Simulation sein.
  print(getwd())
  print("Start Model")
  print(paste(lcpath,ctl_file, sep=" "))
  if (Sys.info()[1]=="Linux") system(paste(oldwd,lcpath," ",ctl_file, sep=""))   ### FOR LINUX Die Simulation wird ausgeführt und schreibt die Ergebnisse ins Verzeichnis ../Output.
  setwd(paste(oldwd,"/Simulations/",sim_name,sep="")) ### Zurücksetzen des Working Directories.
  move_Output(newname = "Output")
  setwd(oldwd) ### Zurücksetzen WD
}

create_inputdir <- function (sim_name,
                             climpath="Data/DWD/climate_feldberg.dat",
                             ctlpath="Data/Landclim/",
                             species=c("abiealba", "fagusilv", "larideci"), ex ) {
  simdir<-paste("Simulations/",sim_name,"/Input" ,sep="")
  dir.create(simdir, recursive=TRUE)
  file.copy(climpath, simdir)
  file.copy(paste(ctlpath, "ctl_bforest.xml", sep=""), simdir)
  file.copy(paste(ctlpath, "landtype.xml", sep=""), simdir)
  specieslist<-read_species_xml("Data/Species_Full/species.xml")
  specieslist <- specieslist[specieslist$name %in% species,]
  write_species_xml(specieslist, file=paste(simdir,"/species.xml", sep=""))
  oldwd<-getwd()
  setwd(simdir)
  write_LandClim_maps(LandClimRasterStack=maps25, nodata_value="-9999", lcResolution=25, ex=ex)
  setwd(oldwd)                  
}

move_Output <- function(newname="Output1"){  
  ### Copy LandClim-Output to the standart location "Output".
  fi <-  list.files()[grep("Output", list.files())][-1]
  file.copy(fi, paste("Output/", fi, sep=""), overwrite=T)
  file.remove(fi)
  file.rename("Output", newname)
  
  ### Rename files
  fi <- list.files(newname)
  fi_new <- substring(fi, 8)
  file.rename(paste(newname, "/", fi, sep=""),paste(newname, "/", fi_new, sep=""))
}

raster_from_output <- function (simname, decade=NA , spec , variable, res=c(25,25), crs="+init=epsg:31467") {
  require(raster)
  
  ifelse (is.na(decade) ,
    dat <- simname,
    dat <- read.csv(paste("Simulations/",simname,"/Output/fullOut_",decade,".csv", sep=""), strip.white=T)
  )
  dat$biomass_cohort <- dat$biomass * dat$stems
  dat<-subset(dat, dat$species==spec, select =  c("xcoord","ycoord", variable))
  #dat$ycoord<-rev(dat$ycoord)
  raster <- rasterFromXYZ(dat, res=res, crs=crs)
  raster
}

make_3d_plot <- function (raster, w=matrix(1/9,nrow=3,ncol=3)) {
  require(rasterVis)
  require(rgl)
  raster[is.na(raster)] <- 0
  raster <- focal(raster, w, fun=mean)
  plot3D(raster)
}

plot_succession <- function(elevationBiomassOut, species, elevation, elev_var=20, lty=1, cols= rainbow(length(species)), plotlegend=TRUE){
  a <-elevationBiomassOut$elevation==c( 800 - 0:elev_var, elevation + 0:elev_var)
  matplot(elevationBiomassOut$decade[a], elevationBiomassOut[a,colnames(elevationBiomassOut) %in% species], type="l", lty=lty, col=cols, xlab="Decade", ylab="Biomass (t/ha)")
  if(plotlegend) legend("topright", legend=species, lty=lty, col=cols, bg="white")
}


