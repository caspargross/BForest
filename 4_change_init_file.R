##########################
# Manipulate Input Files #
#                        #
###########################


type1 <- "mono_piceabie_dec50"
init1 <- read.csv(paste("Data/Init_State/", type, ".csv", sep=""))

type2 <-"mono_abiealba_dec50"
init2 <- read.csv(paste("Data/Init_State/", type2, ".csv", sep=""))

### Reverse and fix YCoordinates and Row Numbers
init1 <- rev_ycoords(init1)
init2 <- rev_ycoords(init2)

plot(out2raster(init1))
plot(out2raster(init1, var="biomass_cohort"))
plot(out2raster(init2))
plot(out2raster(init2, var="biomass_cohort"))


## Create Mask (see 0_functions file)
mask1 <- create_mask( npatch_col=3, npatch_row=4, patch_width=2, patch_length=2, outputfile=init1)
plot(mask1)

## Apply Mask, change Piceabie to Abiealba in Original List
cell_mask<-rowColFromCell(mask1, which(mask1[]==1))

out1 <- init1[-which(paste(init1$row, init1$col) %in% paste(cell_mask[,"row"], cell_mask[,"col"])),]
out1 <- rbind(out1, init2[which(paste(init2$row, init2$col) %in% paste(cell_mask[,"row"], cell_mask[,"col"])),])
out1 <- out1[order(out1$row, out1$col, -out1$age),] #Reorder data.frame
row.names(out1) <- seq_along(out1[,1])
View(out1)

### Check if successfull!
r_out1 <- out2raster2(out1, var="species")
levelplot(r_out1)   ## SUCCESS!

## Write table to /Data Folder
View(out1)
write.table(out1, file="Data/Init_State/3x4_2x2_mono50.csv", row.names=F, col.names=T, sep=",", quote=F)













####### Problemsuche
in1 <- raster("Simulations/3br/Input/dem.asc")
out1 <- read.csv("Simulations/3simple/Output/fullOut_4.csv")
ex_out1 <- extent(c(min(out1$xcoord), max(out1$xcoord), min(out1$ycoord), max(out1$ycoord)))
plot(dem_gk)
plot(extent(in1), add=T, col="red")
plot(r, add=T, col="blue")

