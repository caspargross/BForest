##############################
## FRONT DISPERSAL (NEWEST) ##
##############################

library(xtable)  # Export to Latex
require(parallel)


############## 
#Functions:
##############

# create landscape
dis_landscape_frnew <- function (alt) {
  ## Flat Gradient Model Landscape
  gk_projection<-CRS("+init=epsg:31467") # GK Zone 3 (Black Forest)
  nr <-80  #corresponds to 25000m
  nc <- 20   # corresponds to 25000m
  res <- 25   #25m grid cells
  ex <- extent(0, nc*res, 0, nr*res)
  dem_dis <- raster(nrows=nr, ncols=nc, ex)
  projection(dem_dis) <- gk_projection
  dem_dis[] <- alt
  maps25_frnew <- create_LandClim_Maps(dem_dis, no_aspect=T)
}

# create .xml for different timesteps / initfiles
changexml <- function (x) {
  ctl <- readLines("Data/Landclim/ctl_frnew.xml")
  string <- paste("\t\t<decades>",x,"</decades></full>")
  ctl[64] <- string
  writeLines(ctl, "Data/Landclim/ctl_frnew.xml")
  # change line no 64
}


### Create map list

list_alt <- as.list (seq (400, 1600, 200))
maps_list_frnew <- lapply (list_alt, dis_landscape_frnew)


## Create Simulations
for (i in seq_along (maps_list_frnew)) { create_inputdir (paste("dis_frnew_fs-pa", i,  sep=""),
                                                    climpath="Data/DWD/climate_feldberg.dat",
                                                    species=c("fagusilv", "piceabie"), ex=F,
                                                    LandClimRasterStack=maps_list_frnew[[i]],
                                                    inputfile=F,
                                                    ctlfile="Data/Landclim/ctl_bforest_50.xml",
                                                    landtypefile="Data/Landclim/landtype.xml") }

for (i in seq_along (maps_list_frnew)) { create_inputdir (paste("dis_frnew_aa", i,  sep=""),
                                                    climpath="Data/DWD/climate_feldberg.dat",
                                                    species=c("abiealba"), ex=F,
                                                    LandClimRasterStack=maps_list_frnew[[i]],
                                                    inputfile=F,
                                                    ctlfile="Data/Landclim/ctl_bforest_50.xml",
                                                    landtypefile="Data/Landclim/landtype.xml") }



## Run LandClim Model
lapply(paste("dis_frnew_aa", seq_along(maps_list_frnew),  sep=""), function(x) run_landclim_model(x, ctl_file="ctl_bforest_50.xml"))
lapply(paste("dis_frnew_fs-pa", seq_along(maps_list_frnew),  sep=""), function(x) run_landclim_model(x, ctl_file="ctl_bforest_50.xml"))


## Create Inputfiles
initlist_fspaa <- as.list( paste("Simulations/dis_frnew_fs-pa", seq_along(maps_list_frnew), "/Output/fullOut_50.csv", sep=""))
initlist_aa <- as.list( paste("Simulations/dis_frnew_aa",  seq_along(maps_list_frnew), "/Output/fullOut_50.csv", sep=""))


### Front Dispersal
#############
#Parameter Initialisation
n <- length(maps_list_frnew)
th_y <- list(4) ## Starting Front: 100m
th_y <- rep(th_y, n)
t <- 200  ## Number of Decades 
initlist2 <- initlist_aa

# Run Model
for (t in seq(t)){
outputlist<- as.list( paste("Data/Init_State/frnew/dis_frnew_alt",  1:n , "step_", t, ".csv", sep=""))

mask_frnew <- list(raster(extent(maps_list_frnew[[1]]), resolution=res(maps_list_frnew[[1]])))
mask_frnew <- rep(mask_frnew, n)
for (i in 1:n){ mask_frnew[[i]][1:th_y[[i]],] <- 1 }

## Combine Input +mask
init1 <- lapply(initlist_fspaa, function (x) fread(x, integer64="numeric"))
init2 <- lapply(initlist2, function (x) fread(x, integer64="numeric")) 
init1 <- lapply(init1, function(x) rev_ycoordsDT(x, extent(mask_frnew[[1]])))
init2 <- lapply(init2, function(x) rev_ycoordsDT(x, extent(mask_frnew[[1]])))
out <- list()
out <- mapply(mask_apply, x=init1, y=init2, mask=mask_frnew,  SIMPLIFY=F, USE.NAMES = F )
for (i in 1:length(out)) {write.table(out[[i]], file=outputlist[[i]], row.names=F, col.names=T, sep=",", quote=F)}

for (i in seq_along(maps_list_frnew))
  create_inputdir (paste("dis_frnew/dis_frnew_alt", i, "_step", t, sep=""),
                 climpath=paste("Data/DWD/climate_feldberg_",k,".dat", sep=""),
                 species=c("abiealba", "fagusilv", "piceabie"), ex=F,
                 LandClimRasterStack=maps_list_frnew[[i]],
                 inputfile=outputlist[[i]],
                 ctlfile="Data/Landclim/ctl_frnew.xml",
                 landtypefile="Data/Landclim/landtype.xml",
                 clim_rename=T)

mclapply(as.list(paste("dis_frnew/dis_frnew_alt", seq_along(maps_list_frnew), "_step", t, sep="")), function(x) run_landclim_model(x, ctl_file="ctl_frnew.xml"), mc.cores=3)
frfiles <- as.list(paste("Simulations/dis_frnew/dis_frnew_alt", seq_along(maps_list_frnew), "_step",t,"/Output/fullOut_1.csv",  sep=""))
#print(frfiles)
initlist2 <- frfiles
frfiles <- lapply(frfiles, fread)
frfiles <- lapply(frfiles, function(x) rev_ycoordsDT(x, aui=extent(maps_list_frnew[[2]])))
if(plot==T) print(levelplot(out2rasterDT(frfiles[[4]])))

foo <- function(x) {  # Function to find the front of established seed trees. 
  DT <- x
  DT[,bio_cohort:=biomass*stems,]
  th_y <- DT[species=="abiealba" & age>70 , max(ycoord)] +160    ## IMPORTANT STEP: Cutoff Distance = Max range of esdtablished trees (>70) + maximal seedling distance
  th_y <- ceiling(th_y/25)+1  # round for 25m patch cell
  th_y
}

th_y <- lapply(frfiles, foo)
print(th_y[1:n])
}



# LOAD ALL RESULT FILES
stats_frnew <- data.table(year=NA, elevation=NA, yfront=NA)
y_origin <- 100   ## Starting Front: 100m
for (t in 1:t) {
print(paste("Timestep: ", t))
frnew_res <- as.list(paste("Simulations/dis_frnew/dis_frnew_alt", 1:n, "_step", t, "/Output/fullOut_1.csv" ,sep=""))
frnew_res <- lapply(frnew_res, fread)
frnew_res <- lapply(frnew_res, rev_ycoordsDT, aui=extent(maps_list_frnew[[2]]))


for (i in 1:length(frnew_res)) {
  DT <- frnew_res[[i]]
  DT[,bio_cohort:=biomass*stems,]
  DT <- DT[,.( bio_tot=sum(bio_cohort)), by=.(species, cell, xcoord, ycoord, elevation, age)]
  setkey(DT,species)
  DT <- DT["abiealba"]
  setkey(DT, age)
  DT <- DT[age>=70]
  DT[,dist:= abs(quantile(ycoord, probs=0.95) - y_origin),]
  DT <- DT[, .(year=unique(t*10), elevation=unique(elevation), yfront=max(dist) )]
  stats_frnew <- rbind(stats_frnew, DT)
  }  
}
stats_frnew <- stats_frnew[-1]
stats_frnew$elevation <- as.factor(stats_frnew$elevation) 
stats_frnew[,rm_yfront:=runmean(yfront, 5), by=elevation]

### NICE PRINTS
ele <- c(1,2,3,4,5,6,7)
levels(ele) <- c("400", "600", "800", "1000", "1200", "1400", "1600")
dec <- c(10,50,100,150)
levels(dec) <- c("100 yr", "500 yr", "1000 yr", "1500 yr")
gr <- expand.grid(ele,dec)
frnew_plot <- as.list(paste("Simulations/dis_frnew/dis_frnew_alt", gr[,1], "_step", gr[,2], "/Output/fullOut_1.csv" ,sep=""))
frnew_plot <- lapply(frnew_plot, fread)
frnew_plot <- lapply(frnew_plot, rev_ycoordsDT, aui=extent(maps_list_frnew[[2]]))
pl <- data.table(cell=NA, xcoord=NA, ycoord=NA, species=NA, bio_tot=NA, ratio_aa=NA, sumbio_tot=NA, ele=NA, dec=NA)
for (i in 1:length(frnew_plot)) {
  DT <- frnew_plot[[i]]
  DT[,bio_cohort:=biomass*stems,]
  DT <- DT[,.( bio_tot=sum(bio_cohort)), by=.(species, cell, xcoord, ycoord, elevation)]
  DT[,bio_aa:=as.double(NA),]  # Create new column with abies alba biomass
  DT[species == "abiealba", bio_aa:=bio_tot, by=.(cell)]
  DT[, sumbio_tot := sum(bio_tot), by=cell]
  DT[, ratio_aa := na.omit(bio_aa) / sumbio_tot, by=cell ] #Calculate aa/total biomass ratio
  DT <- DT[, .SD[which.max(bio_tot), .(xcoord, ycoord, species, bio_tot, ratio_aa, sumbio_tot)], by=.(cell)]
  DT[,ele:=gr[i,1],]
  DT[,dec:=gr[i,2],]
  pl <- rbind(pl, DT)
}  

pl <- pl[-1]
pl[, ele_factor:=factor(pl$ele, levels=c(7,6,5,4,3,2,1), labels=c("1600", "1400", "1200", "1000", "800", "600", "400")), ]


### Plot 1: Images
frnew <- ggplot(pl, aes(x=ycoord, y=xcoord)) +
  theme_cas() +
  theme(legend.position = "bottom") +
  coord_fixed() +
  geom_tile(aes(fill=species),  data= pl[species == "fagusilv"], fill="grey80")+
  geom_tile(aes(fill=species),  data= pl[species == "piceabie"], fill="wheat1")+
  geom_tile(aes(fill=ratio_aa)) +
  scale_fill_gradient(low="lightgreen", high="darkgreen", limits=c(0, 1), na.value="transparent")+
  facet_grid(ele_factor ~ dec) +
  labs(x= "Length (m)", y= "Width (m)", fill="Biomass ratio \n A. alba / total") +
  geom_vline(aes(xintercept=yfront), data=stats_frnew[year==dec  & elevation==dec])+
  scale_y_continuous(breaks=seq(0, 400, 200))  

frnew

### Linear model for the propagation speeds:
lm_models <- stats_frnew[yfront<1700]

models <- data.frame(elevation=character(0), coef=character(0), rsq=character(0))
for (i in 1:7) { 
  mod <- stats_frnew[as.integer(elevation) == i & yfront < 1700, lm(yfront ~ 0 + year), ]
  coef <- round(mod$coef,3)
  rsq <- round(summary(mod)$r.squared, 3)
  ele <- seq(400,1600,200)[i]
  models <- rbind(models, data.frame(elevation=ele, coef=coef[[1]], rsq=rsq) )
}

models$elevation <- as.factor(models$elevation)

### Plot 2: Distances
frnew_dist <- ggplot(stats_frnew[yfront < 1700], aes(x=year, y=yfront, col=elevation, group=elevation)) +
  #geom_line()+
  geom_point(pch=3) +
  theme_cas_big() +
  geom_abline(aes(slope=coef, col=elevation), data=models)+
  labs(x ="Time (yr)", y="Distance covered by A. alba front (m)", col="Elevation \n (m a.s.l.)") 
  #geom_abline(intercept=0, slope=1.6, lty=2)+
  #geom_abline(intercept=0, slope=0.5, lty=3) 
    #scale_colour_brewer(palette="YlGnBu")
 #(values= terrain.colors(10, alpha = 1))
frnew_dist

### Export figures as .pdf

pdf(file="Figures/dis_front_lm.pdf", width=8, height=6)
print(frnew_dist)
dev.off()

pdf(file="Figures/dis_front_map.pdf", width=10, height=10)
print(frnew)
dev.off()

### Export table to latex

write(print(xtable(models), booktabs=T), file="lm_table.txt")

