########################
## Functions LandClim ##
########################

range01 <- function(x){(x-min(x))/(max(x)-min(x))}


theme_cas <- function (base_size = 12, base_family = "") {
  theme_grey(base_size = base_size, base_family = base_family) %+replace% 
  theme(axis.text = element_text(size = rel(0.8)), axis.ticks = element_line(colour = "black"), 
        legend.key = element_rect(colour = "grey80"),
        panel.background = element_rect(fill = "white", colour = NA),
        panel.border = element_rect(fill = NA, colour = "grey50"),
        panel.grid.major = element_line(colour = "grey90", size = 0.2),
        panel.grid.minor = element_line(colour = "grey98"),
        strip.background = element_rect(colour = "grey50", size = 0.5))
}

run_landclim_model<-function(sim_name, ctl_file="ctl_bforest.xml", lcpath="/Data/Landclim/LandClim"){
  oldwd <- getwd()
  if (file.exists(paste("Simulations/",sim_name,"/Output",sep=""))) unlink(paste("Simulations/",sim_name,"/Output",sep=""), recursive=TRUE)
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
                             species=c("abiealba", "piceabie"), ex=F,
                             LandClimRasterStack=maps25,
                             inputfile=F,
                             ctlfile="Data/Landclim/ctl_bforest.xml",
                             landtypefile="Data/Landclim/landtype.xml",
                             clim_rename=F){
  simdir<-paste("Simulations/",sim_name,"/Input" ,sep="")
  if (file.exists(simdir)) unlink(simdir, recursive=TRUE)
  dir.create(simdir, recursive=TRUE)
  ifelse(clim_rename==T,  file.copy(climpath, paste(simdir,"/climate_feldberg.dat", sep="")), file.copy(climpath, simdir))
  file.copy(landtypefile, simdir)
  if (inputfile!=F) {
    file.copy(inputfile, paste(simdir,"/tree_init.csv", sep=""))
   }
  file.copy(ctlfile, simdir) 
  specieslist<-read_species_xml("Data/Species_Full/species.xml")
  specieslist <- specieslist[specieslist$name %in% species,]
  write_species_xml(specieslist, file=paste(simdir,"/species.xml", sep=""))
  oldwd<-getwd()
  setwd(simdir)
  write_LandClim_maps(LandClimRasterStack, nodata_value="-9999", lcResolution=25, ex)
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



rev_ycoords <-function(out, aui_rev=aui, res_rev=c(25,25)) 
{
  yc <- out$ycoord
  foo <- function (x) (aui_rev@ymax+res_rev[2])-(which(unique(yc)==x)*res_rev[2])
  out$ycoord <- sapply(yc, foo)
  out <- out[order(-out$ycoord, out$xcoord),]
  out
}

rev_ycoordsDT <-function(out, aui_rev=aui, res_rev=c(25,25)) 
{ require(data.table)
  yc <- seq (from=aui_rev@ymin, to=aui_rev@ymax, by=res_rev[2])
  out$ycoord <- 0
  out$ycoord <- as.integer(out$ycoord)
  for (i in 1:length(unique(out$row)))  {out[row==i, ycoord:=yc[i]]}
  out[order(-rank(ycoord), xcoord)]
}

out2raster <- function (dat, var="elevation")
{
  ex <- extent(min(dat$xcoord), max(dat$xcoord), min(dat$ycoord), max(dat$ycoord)) 
  r <- raster (ex=ex, res=c(25,25), crs="+init=epsg:31467")
  if (var=="species" | var=="dom_species") {   ## This function creates a factorized raster for species information. Plot with levelplot (rasterVis)
      r<- rasterize(cbind(dat$xcoord, dat$ycoord), r, field=as.numeric(dat[[var]]))
      r <- ratify(r) 
      rat <- levels(r)[[1]]
      rat$legend <- levels(dat$species)
      levels(r) <- rat
      } else {r<- rasterize(cbind(dat$xcoord, dat$ycoord), r, field=dat[[var]])}
  r
}

## For data.table class
out2rasterDT <- function (dat, var="species"){ 
  ex <- extent(min(dat$xcoord), max(dat$xcoord), min(dat$ycoord), max(dat$ycoord)) 
  r <- raster (ex=ex, res=c(25,25), crs="+init=epsg:31467")
  if (var=="species" | var=="dom_species") {   ## This function creates a factorized raster for species information. Plot with levelplot (rasterVis)
    r<- rasterize(cbind(dat$xcoord, dat$ycoord), r, field=as.numeric(as.factor(dat[[var]])))
    r <- ratify(r) 
    rat <- levels(r)[[1]]
    rat$legend <- levels(as.factor(dat$species))
    levels(r) <- rat
  } else {r<- rasterize(cbind(dat$xcoord, dat$ycoord), r, field=dat[[var]])}
  r
}


out2rasterDT2 <- function (dat, var){ 
  ex <- extent(min(dat$xcoord), max(dat$xcoord), min(dat$ycoord), max(dat$ycoord)) 
  r <- raster (ex=ex, res=c(25,25), crs="+init=epsg:31467")
  r<- rasterize(cbind(dat$xcoord, dat$ycoord), r, field=dat[[var]])
  r
}

### Create Mask for superposition

create_mask <- function (npatch_row, npatch_col, patch_width, patch_length, outputfile, p=F)
{
  #outputfile <- init1
  centerrow <- NULL
  centercol <- NULL
  if (class(outputfile)[1]=="RasterLayer") {mask <- outputfile 
  } else {mask <- out2rasterDT(outputfile)}
  #print(plot(mask))
  mask[] <- NA
  print("empty mask loaded")  
  for (i in seq(npatch_row)) {centerrow[i]<-i*(ceiling(nrow(mask)/(npatch_row+1)))}
  for (i in seq(npatch_col)) {centercol[i]<-i*(ceiling(ncol(mask)/(npatch_col+1)))}

  xpatch <- sapply(centerrow, function (x) floor((x-(patch_width/2)):(x+(patch_width/2))))
  ypatch <- sapply(centercol, function (x) floor((x-(patch_length/2)):(x+(patch_length/2))))

  mask[c(xpatch), c(ypatch)]<-1
  print("mask created")
  if (p==T) print(plot(mask, col="red", add=T, legend=F))
  mask
}

create_mask <- function (npatch_row, npatch_col, patch_width, patch_length, outputfile, p=F)
{
  #outputfile <- init1
  
  centerrow <- NULL
  centercol <- NULL
  if (class(outputfile)[1]=="RasterLayer") {mask <- outputfile 
  } else {mask <- out2rasterDT(outputfile)}
  #print(plot(mask))
  mask[] <- NA
  print("empty mask loaded")
  
  for (i in seq(npatch_row)) {centerrow[i]<-i*(ceiling(nrow(mask)/(npatch_row+1)))}
  for (i in seq(npatch_col)) {centercol[i]<-i*(ceiling(ncol(mask)/(npatch_col+1)))}
  
  xpatch <- sapply(centerrow, function (x) floor((x-(patch_width/2)):(x+(patch_width/2))))
  ypatch <- sapply(centercol, function (x) floor((x-(patch_length/2)):(x+(patch_length/2))))
  
  mask[c(xpatch), c(ypatch)]<-1
  print("mask created")
  if (p==T) print(plot(mask, col="red", legend=F))
  mask
}





create_mask_p <- function (npatch_row, npatch_col, patch_width, patch_length, outputfile, p=F, vshift=0){
  
  xc_plus <- (nrow(mask)/npatch_row)*1:npatch_row
  xc_minus<- sort(nrow(mask)-(nrow(mask)/npatch_row)*1:npatch_row)
  xcorners <- mapply(function(x,y) mean(c(x,y))+vshift, x=xc_plus, y=xc_minus)
  xcorners <- round(xcorners - (patch_width/2))
  
  
  yc_plus <- (ncol(mask)/npatch_col)*1:npatch_col
  yc_minus<- sort(ncol(mask)-(ncol(mask)/npatch_col)*1:npatch_col)
  ycorners <- mapply(function(x,y) mean(c(x,y)), x=yc_plus, y=yc_minus)
  ycorners <- round(ycorners - (patch_length/2))
  
  xpatch <- sapply(xcorners, function (x) x:(x+(patch_width-1)))
  ypatch <- sapply(ycorners, function (x) x:(x+(patch_length-1)))
  
  mask[c(xpatch), c(ypatch)]<-1
  print("mask created")
  if (p==T) print(plot(mask, col="red", legend=F))
  mask
}


## Create percentage Biomass DT from Outputfile
biomass_dt <- function (file) {
  DT <- fread(file)
  DT <- rev_ycoords(DT)   #reverse Ycoords
  DT[,biomass_cohort:=biomass*stems]     #create cohort_biomass
  DT <- DT[,.(biomass_cohort=sum(biomass_cohort)),by=list(species, cell, xcoord, ycoord)]  #
  DT[species=="abiealba", biomass_aa:=biomass_cohort]
  DT[species=="piceabie", biomass_pa:=biomass_cohort]
  DT <- DT[,.(biomass_cell=sum(biomass_cohort), biomass_aa=sum(biomass_aa, na.rm=T), biomass_pa=sum(biomass_pa, na.rm=T) ), by=.(cell, xcoord, ycoord)]
  DT[,biomass_perc:=((biomass_aa/biomass_cell)*100)]
  as.data.frame(DT)
}

## Apply mask to raster!
mask_apply <- function (x, y, mask=mask1) {
  cell_mask<-rowColFromCell(mask, which(mask[]==1))-1
  out1 <- x[-which(paste(x$row, x$col) %in% paste(cell_mask[,"row"], cell_mask[,"col"])),]
  out1 <- rbind(out1, y[which(paste(y$row, y$col) %in% paste(cell_mask[,"row"], cell_mask[,"col"])),])
  out1 <- out1[order(out1$row, out1$col, -out1$age),] #Reorder data.frame
  #row.names(out1) <- seq_along(out1[,1])
  print("mask apply sequence finished")
  out1
  #View(out1)
}
