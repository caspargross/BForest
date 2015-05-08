explained.deviance <- function(model, null.model=NULL){
  if(class(model)== "glm")  dev <-   (model$null.deviance - model$deviance) / model$null.deviance
  if(class(model)[1] == "multinom")  {
    dev <-(null.model$deviance - model$deviance) / null.model$deviance
  }
  dev
}



### LANDCLIM ####
resampleLandClimMaps <- function(LandClimRasterStack, targetResolution=25){
  r <- raster(LandClimRasterStack, layer=1)
  rt <- raster(extent(r), crs=projection(r))
  res(rt) <- targetResolution
  
  foo <- function(x){   
    rre <- resample(x, rt)
    rre[is.na(rre)] <- -9999
    rre 
  }
  res <- lapply(unstack(LandClimRasterStack), foo)
  stack(res)
}

writeLandClimMaps <- function(LandClimRasterStack, nodata_value="-9999", lcResolution=25, ex){
    
  LandClimRasterStack_list <- lapply(unstack(LandClimRasterStack), function(x) crop(x, ex))
  rs <- stack(LandClimRasterStack_list)
  names(rs) <- names(LandClimRasterStack)
  writeRaster(rs, "landClimMaps.tif", overwrite=T)
  rm(rs)
  foo <- function(x){
    
    sink(paste(names(x), ".asc", sep=""))
    writeLines(c(paste("ncols", ncol(x)), paste("nrows", nrow(x)), paste("xllcorner", xmin(x)), paste("yllcorner", ymin(x)), paste("cellsize", lcResolution), paste("nodata_value ", nodata_value)))
    sink()    
    
    write.table(matrix(round(x[]), nrow=nrow(x), ncol=ncol(x), byrow=T), file=paste(names(x), ".asc", sep=""), append=T, quote = FALSE, row.names = FALSE, col.names = FALSE)
    
  }  
  
  lapply(LandClimRasterStack_list, function(x) foo(x))
}

makeClimateChange <- function(inputPath, outputPath, dt=0, dn=0){
  header<- readLines(inputPath)
  header <- header[1:14]
  clim<- read.table(inputPath, skip=11)
  clim[,2:13] <- clim[,2:13] + dt
  clim[,14:25] <- clim[,14:25] + dn
  clim[,14:25][clim[,14:25] < 0] <- 0
  writeLines(header, outputPath)  
  write.table(clim, outputPath, append=T, row.names=F, col.names=F)  
}

moveOutput <- function(newname="Output1"){  
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

biomass.to.dbh <- function(biomass, leafHabit, allometry="SCHUMACHER") {
  if(allometry=="SCHUMACHER") {
    if(!all(leafHabit %in% c("EVERGREEN", "BROADLEAFEVERGREEN", "DECIDUOUS"))) stop ("unknown leaf habit")
    dbh <- ifelse(leafHabit == "EVERGREEN", exp(3.800 + 0.451 * log(biomass)),
                  exp(3.708 + 0.475 * log(biomass)))
  }
  
  if(allometry=="POWER") {
    carbon_kg <- 500 * biomass # assumption: 50% of dry weight is carbon
    dbh <- exp(1.3481648 + 0.3977240 * log(carbon_kg))
    
  }
  return(dbh)
  
}


dbh.to.biomass <- function(dbh, leafHabit, allometry="SCHUMACHER") {
  if(allometry=="SCHUMACHER") {
    if(!all(leafHabit %in% c("EVERGREEN", "BROADLEAFEVERGREEN", "DECIDUOUS"))) stop ("unknown leaf habit")
    biomass <- ifelse(leafHabit == "EVERGREEN", exp((log(dbh) - 3.800)/0.451),
                      exp((log(dbh)- 3.708)/0.475))
    
  }
  
  if(allometry=="POWER") {
    carbon_kg <- exp((log(dbh) - 1.3481648) / 0.3977240)
    biomass <- carbon_kg/500   # assumption: 50% of dry weight is carbon
    
  }
  return(biomass) 
}


calculateLandscapeSize <- function(dem){
  prod(dim(dem)) * prod(res(dem)) / (100*100)
}


plotElevation <- function(elevationBiomassOut, species, decade,   lty=1,  cols= rainbow(length(species))){
  
  a <- elevationBiomassOut$elevation==decade
  matplot(elevationBiomassOut$decade[a], elevationBiomassOut[a,colnames(elevationBiomassOut) %in% species], type="l", lty=lty, col=cols, xlab="Decade", ylab="Biomass (t/ha)", main=paste("Elevation = ",elev, "m.a.s." ,sep=""))
  if(legend==TRUE) legend("topright", legend=species, lty=lty, col=cols, bg="white")
}

plotSuccession <- function(elevationBiomassOut, species, elevation,  lty=1,  cols= rainbow(length(species))){ 
  a <- elevationBiomassOut$elevation==elevation
  matplot(elevationBiomassOut$decade[a], elevationBiomassOut[a,colnames(elevationBiomassOut) %in% species], type="l", ltycalculateLandscapeSize <- function(dem){
  prod(dim(dem)) * prod(res(dem)) / (100*100)
}=lty, col=cols, xlab="Decade", ylab="Biomass (t/ha)", main=paste("Elevation = ", elev, " m.a.s." ,sep=""))
  if(legend==T) legend("topright", legend=species, lty=lty, col=cols, bg="white")
}


