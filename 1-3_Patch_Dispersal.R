###########################
## Dispersal Dynamics    ##
##   Patch Dispersal     ##
##     Area: Model BF    ##
###########################

setwd("/home/caspar/Bachelor_Thesis/BForest/")
library(ModellingTools)
library(sp)
library(rgdal)
library(raster)
library(classInt)  #nicer plot colour gradients

dis_landscape_p <- function (alt) {
  ## Flat Gradient Model Landscape
  gk_projection<-CRS("+init=epsg:31467") # GK Zone 3 (Black Forest)
  nr <-200  #corresponds to 5000m
  nc <- 200   # corresponds to 5000m
  res <- 25   #25m grid cells
  ex <- extent(0, nc*res, 0, nr*res)
  dem_dis <- raster(nrows=nr, ncols=nc, ex)
  projection(dem_dis) <- gk_projection
  dem_dis[] <- alt
  maps25_dis <- create_LandClim_Maps(dem_dis, no_aspect=T)
}

### Create map list

list_alt_p <- as.list (seq (400, 1600, 200))
maps_list_p <- lapply (list_alt_p, dis_landscape_p)


## Create Simulations
for (i in seq_along (maps_list_p)) { create_inputdir (paste("dis_patch_init_fs-pa_", i, sep=""),
                                                    climpath="Data/DWD/climate_feldberg.dat",
                                                    species=c("fagusilv", "piceabie"), ex=F,
                                                    LandClimRasterStack=maps_list_p[[i]],
                                                    inputfile=F,
                                                    ctlfile="Data/Landclim/ctl_bforest_50.xml",
                                                    landtypefile="Data/Landclim/landtype.xml") }

for (i in seq_along (maps_list_p)) { create_inputdir (paste("dis_patch_init_aa_", i, sep=""),
                                                    climpath="Data/DWD/climate_feldberg.dat",
                                                    species=c("abiealba"), ex=F,
                                                    LandClimRasterStack=maps_list_p[[i]],
                                                    inputfile=F,
                                                    ctlfile="Data/Landclim/ctl_bforest_50.xml",
                                                    landtypefile="Data/Landclim/landtype.xml") }

## Run LandClim Model for Background World
mclapply(as.list(paste("dis_patch_init_fs-pa_", seq_along(maps_list_p), sep="")), function(x) run_landclim_model(x, ctl_file="ctl_bforest_50.xml"), mc.cores=3)
mclapply(as.list(paste("dis_patch_init_aa_", seq_along(maps_list_p), sep="")), function(x) run_landclim_model(x, ctl_file="ctl_bforest_50.xml"), mc.cores=3)


## Check Background World
input_files <- as.list(paste("Simulations/dis_patch_init_fs-pa_", seq_along(maps_list_p),"/Output/fullOut_50.csv", sep=""))
input_files <- lapply(input_files, fread)
input_files <- lapply(input_files, function(x) rev_ycoordsDT(x, extent(maps_list_p[[1]])))

##
#### Create Inputfiles
mask <- raster(extent(maps_list_p[[1]]), resolution=res(maps_list_p[[1]]))

m1 <- (create_mask_p(1,1,20,20, mask, p=TRUE)) # 1x400
m2 <- (create_mask_p(1,2,20,10, mask, p=TRUE)) # 2x200
m3 <- (create_mask_p(2,2,10,10, mask, p=TRUE)) # 4x100
m4.1 <- (create_mask_p(3,3,6,7, mask, p=TRUE)) # 8x50  # 6 6x7 patches --> 252 Cells
m4.2 <- (create_mask_p(1,3,7,7, mask, p=TRUE)) # 8x50  # 3 7x7 patches --> 147 Cells | Total: 399 Cells
m4 <- m4.1 | m4.2 # Combine rasters
m5 <- (create_mask_p(4,4,5,5, mask, p=TRUE))   # 16x25
m6 <- (create_mask_p(5,4,4,5, mask, p=TRUE))   # 20x20
m7 <- (create_mask_p(5,5,4,4, mask, p=TRUE))   # 25x16
m8.1 <- (create_mask_p(7,7,3,3, mask, p=TRUE)) # 50x8  # 49 3x3 patches --> 441 Cells
m8.2 <- (create_mask_p(7,1,2,3, mask, p=TRUE, vshift=0)) # remove 7 2*3 = 42 Patches | Total: 399 Cells
m8 <- m8.1 - !is.na(m8.2)
m8[m8[]==0] <- NA
m9 <- (create_mask_p(10,10,2,2, mask, p=TRUE)) # 100x4
m10 <-(create_mask_p(10,20,2,1, mask, p=TRUE)) # 200x2
m11 <-(create_mask_p(20,20,1,1, mask, p=TRUE)) # 400x1


sink("Logs/log_dis_patch.txt")
for (k in 3) {
print(paste("Patch pattern no", k, "startedt at:", Sys.time()))
  
initlist1_p <- as.list( paste("Simulations/dis_patch_init_fs-pa_", seq_along(maps_list_p), "/Output/fullOut_50.csv", sep=""))
initlist2_p <- as.list( paste("Simulations/dis_patch_init_aa_", seq_along(maps_list_p), "/Output/fullOut_50.csv", sep=""))
outputlist_p<- as.list( paste("Data/Init_State/dis_patch_init_alt_", seq_along(maps_list_p),"_m", k, ".csv", sep=""))

print("Loading Files [OK]")

combine_input (initlist1_p,  
               initlist2_p, 
               outputlist_p,
               1, 1, 2, 100,  
               extent(maps_list_p[[1]]),
               ma=eval(as.name(paste("m",k, sep=""))))

print("Creating Input File  [OK]")
## Check Input Files

#input_files <- as.list(paste("Data/Init_State/dis_patch_init_alt_", 1:8,"_m11.csv", sep=""))
#input_files <- lapply(input_files, fread)
#input_files <- lapply(input_files, function(x) rev_ycoordsDT(x, extent(maps_list_p[[1]])))
#levelplot(out2rasterDT(input_files[[4]], var="species"))

## Run Model with Dispersal

for (i in seq_along (maps_list_p)) { create_inputdir (paste("dis_patch_alt_", i,"_m",k, sep=""),
                                                    climpath="Data/DWD/climate_feldberg.dat",
                                                    species=c("abiealba", "fagusilv", "piceabie"), ex=F,
                                                    LandClimRasterStack=maps_list_p[[i]],
                                                    paste("Data/Init_State/dis_patch_init_alt_", i,"_m", k, ".csv", sep=""),
                                                    ctlfile="Data/Landclim/ctl_bforest_dis_5000.xml",
                                                    landtypefile="Data/Landclim/landtype.xml") }
print(paste("Created Input Directory for Pattern", k, "Elevation", i, "[OK]"))
mclapply(as.list(paste("dis_patch_alt_", seq_along(maps_list), "_m", k, sep="")), function(x) run_landclim_model(x, ctl_file="ctl_bforest_dis_5000.xml"), mc.cores=3)
print(paste("Finished Model Calculation Pattern", k, "Elevation", i, "[OK]"))
}
print("SUCCESS!")
sink()



## Check Dispersal Output files
output_files <- as.list(paste("Simulations/dis_patch_alt_1_m1/Output/fullOut_", seq(5,500,5), ".csv", sep=""))
output_files <- lapply(output_files, fread)
output_files <- mclapply(output_files, function(x) rev_ycoordsDT(x,  aui_rev=extent(maps_list_p[[1]])))
for (i in 1:length(output_files)) {
  png(filename=paste("Animate/dis_patch_alt1_m1_dec", i,".png", sep=""))
  print(levelplot(out2rasterDT(output_files[[25]], var="species")))
  dev.off()
}


######### Calculation of Output-Analysis 
stats_p <- data.table(mask=NA, altitude=NA, year=NA, ratio_bio_aa=NA, mass_bio_aa=NA, n_cohorts=NA)  #Data.table with stats
plot <- F
myTheme <- rasterTheme(region=brewer.pal(9, "BuGn"))
for (m in 1:11) {
  
  for (j in seq_along (maps_list_p)) {  ## Different Alitudes!
  
    list_results_p <- as.list (paste("Simulations/dis_patch_alt_", j,"_m",m,"/Output/fullOut_", seq(5, 500, 5), ".csv", sep=""))
    dis_results_p <- lapply(list_results_p, fread)
    dis_results_p <- mclapply(dis_results_p, function(x) rev_ycoordsDT(x, extent(maps_list_p[[1]])),  mc.cores=3)
    print(paste("Loaded map",j, "Patch Pattern no:", m))
     
  
    for (i in 1:length(dis_results_p)) { ## Different Decades !
    
      PDT<-dis_results_p[[i]]
      PDT[,bio_cohort:=biomass*stems] #Calculate cohort biomas
      PDT <- PDT[,.(bio_cohort=sum(bio_cohort), age=max(age)),by=list(species, cell, xcoord, ycoord)]
      setkey(PDT,species)
      PDT[,bio_aa:=as.double(NA),]  ## Create empty columnt (NA) in double type
      PDT["abiealba", bio_aa:=bio_cohort, by=.(cell)]
      PDT2 <- PDT[,.(bio_cell=sum(bio_cohort), bio_aa), by=.(xcoord, ycoord, cell)]
      PDT2[is.na(PDT2)] = 0  # replace NA with 0
      PDT2 <- PDT2[,.(bio_cell=mean(bio_cell), bio_aa=max(bio_aa), bio_quo=max((bio_aa/bio_cell)*100)), by=.(xcoord, ycoord, cell)]
    
      stats_p <- rbind(stats_p, list(m, list_alt_p[[j]], i*50, PDT2[,round((sum(bio_aa)/sum(bio_cell))*100,2)], PDT2[,sum(bio_aa),], nrow(PDT["abiealba", .(age=max(age)) , by="cell"][age>=70,,])))
      print(paste("finished decade", i))
    
      if (plot == T) {
        png(filename=paste("Animate/Patch/dis_patch_alt",j,"_",m,"_dec", i,".png", sep=""))
        print(levelplot(out2rasterDT(PDT2, var="bio_aa"), main=paste("Year", i*50, "Ele", j, "Mask", m ), par.settings=myTheme))
        dev.off()
      
      }
    }
  }
}    
    
stats_p <- stats_p[-1,]
save(stats_p, file="Data/Results/results_patch_dispersal.RData")
load("Data/Results/results_patch_dispersal.RData")
#write.table(stats_p, "Data/Results/results_patch_dispersal.txt", sep="\t", row.names=F)
#stats_p <- fread("Data/Results/results_patch_dispersal_old.txt")
#save(stats_p, file="Data/Results/patch_disp.RData")
setnames(stats_p, "altitude", "elevation")
stats_p$mask <- as.integer(stats_p$mask)
stats_p$elevation <- as.integer(stats_p$elevation)
stats_p$elevation <- seq(400, 1600, 200)[stats_p$elevation]
#stats_p <- stats_p[altitude!=1800,,]
setkey(stats_p, mask, elevation)

## Make Running mean for Elevation 400,
stats_p[, ele_factor:=factor(stats_p$elevation, levels=c("1600", "1400", "1200", "1000", "800", "600", "400")), ] ## Reverse Levels for plot order
stats_p[elevation==400, ratio_bio_aa_org:=ratio_bio_aa, ]
stats_p[elevation==400, ratio_bio_aa:=runmean(ratio_bio_aa, 4), by=.(mask)]
stats_p[, th:=0.75*max(ratio_bio_aa), by=.(elevation)]


## Calculate the year when biomass reaches 0.75 of total biomass
th_time <- stats_p[, .(max_ratio=(0.75*max(ratio_bio_aa)), year, ratio_bio_aa, mask ), by=.(elevation)]
setkey(th_time, mask)
th_time[,low_ratio:=shift(ratio_bio_aa, 1, type="lag")]
th_time <- th_time[ratio_bio_aa>=max_ratio, .(th_year=min(year), year, upper_limit=(min(year)),  up_ratio=ratio_bio_aa, low_ratio) , by=.(mask, elevation, max_ratio)]
th_time[, .( up_year=min(year), low_ratio=low_ratio[1], up_ratio=up_ratio[1], th_ratio=max_ratio[1]), by=.(mask, elevation)]
th_time <- th_time[, .(low_year=min(year)-50, up_year=min(year), low_ratio=low_ratio[1], up_ratio=up_ratio[1], th_ratio=max_ratio[1]), by=.(mask, elevation)]

th_time[,m:=((up_ratio-low_ratio)/50),]
th_time[,x:=((th_ratio-low_ratio)/m),]
th_time[,th_year:=round(low_year+x),]
th_time[,mean_year:=mean(th_year), by=.(elevation)]
th_time[,th_year01:=range01(th_year),]
th_time[, ele_factor:=factor(th_time$elevation, levels=c("1600", "1400", "1200", "1000", "800", "600", "400")), ]
th_time

### Plot the biomass ratios
library(ggplot2)
distplot <- ggplot(stats_p, aes(x=year, group=mask, col = as.factor(mask))) +
  theme_cas_big() +
  theme(legend.position = "bottom") +
  geom_hline(aes(yintercept=th), col="grey20", lty=2, show_guide = TRUE)+
  geom_line(aes(y=ratio_bio_aa)) +
  geom_rug(mapping=aes(x=th_year, group=mask), data=th_time, sides="t", size=1)+
  labs( y="Biomass percentage Abies alba / Total", x="Year", col="Mask") + 
  facet_grid(ele_factor ~ .) 
distplot


## PLot Heatmap
# create colour breaks based on bclust 
brks <- classIntervals(th_time$th_year, 10)
cas_palette <- colorRampPalette(c("darkgreen", "palegreen", "white"))
hm <- ggplot (th_time, aes(x=as.factor(mask), y=as.factor(elevation))) +
  geom_tile(aes(fill=th_year), colour="grey80") +
  coord_fixed()+
  theme_cas()+
  #scale_fill_gradientn(colours=cas_palette(10), breaks=round(brks$brks)) +
  scale_fill_gradientn(colours=cas_palette(10), guide = "colourbar") +
  labs( x="Initial configuration mask", y="Elevation (m a.s.l.)", fill="Threshold \n year")
hm

## PLot Dispersal Maps + GIF
dec_p <- c(5, 50, 100, 200)
levels(dec_p) <- c("50 yr", "200 yr", "500 yr", "1500 yr")
ele_p <- c(1, 2, 4, 6, 7)
levels(ele_p) <- c("600", "1000",  "1400" )
gr_p <- expand.grid(dec_p, ele_p)

dis_p_files <- as.list(paste("Simulations/dis_patch_alt_", gr_p[,2], "_m3/Output/fullOut_", gr_p[,1], ".csv" ,sep=""))
dis_p_files <- lapply(dis_p_files, fread)
dis_p_files <- lapply(dis_p_files, function(x) rev_ycoordsDT(x, extent(maps_list_p[[1]])))
pl_p <- data.table(cell=NA, xcoord=NA, ycoord=NA, species=NA, bio_tot=NA, ratio_aa=NA, sumbio_tot=NA, ele=NA, dec=NA)
for (i in 1:length(dis_p_files)) {
  DT <- dis_p_files[[i]]
  print(paste("Load Step", i))
  DT[,bio_cohort:=biomass*stems,]
  DT <- DT[,.( bio_tot=sum(bio_cohort)), by=.(species, cell, xcoord, ycoord, elevation)]
  DT[,bio_aa:=as.double(NA),]  # Create new column with abies alba biomass
  DT[species == "abiealba", bio_aa:=bio_tot, by=.(cell)]
  DT[, sumbio_tot := sum(bio_tot), by=cell]
  DT[, ratio_aa := na.omit(bio_aa) / sumbio_tot, by=cell ] #Calculate aa/total biomass ratio
  DT <- DT[, .SD[which.max(bio_tot), .(xcoord, ycoord, species, bio_tot, ratio_aa, sumbio_tot)], by=.(cell)]
  DT[,ele:=gr_p[i,1],]
  DT[,dec:=gr_p[i,2],]
  pl_p <- rbind(pl_p, DT)
  print(paste("Finish Step", i))
}  

pl_p <- pl_p[-1]

pnew <- ggplot(pl_p, aes(x=ycoord, y=xcoord)) +
  theme_cas_big() +
  theme(legend.position = "bottom") +
  coord_fixed() +
  geom_tile(aes(fill=species),  data= pl_p[species == "fagusilv"], fill="grey80")+
  geom_tile(aes(fill=species),  data= pl_p[species == "piceabie"], fill="wheat1")+
  geom_tile(aes(fill=ratio_aa)) +
  scale_fill_gradient(low="lightgreen", high="darkgreen", limits=c(0, 1), na.value="transparent")+
  facet_grid(dec ~ ele) +
  labs(x= "Length (m)", y= "Width (m)", fill="Biomass ratio \n A. alba / total") +
  scale_x_continuous(breaks=seq(0, 4000, 2000)) +
  scale_y_continuous(breaks=seq(0, 4000, 2000)) 

#### PLOT RESULTS
pdf(file="Figures/dis_patch_hm.pdf", width=8, height=6)
print(hm)
dev.off()

pdf(file="Figures/dis_patch_map.pdf", width=10, height=14)
print(pnew)
dev.off()
