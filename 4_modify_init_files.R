##########################
# Manipulate Input Files #
#                        #
###########################

combine_input <- function (filelist1, filelist2, outputfile, npatch_col=3, npatch_row=4, patch_width=2, patch_length=2, aui=aui ) {
require(data.table)
require(bit64)
require(raster)

print("step1")
init1 <- lapply(filelist1, fread)

print("step2")
init2 <- lapply(filelist2, fread) 
print("step3")
#type1 <- "mono_piceabie_dec50"
#init1 <- read.csv(paste("Data/Init_State/", type, ".csv", sep=""))

#type2 <-"mono_abiealba_dec50"
#init2 <- read.csv(paste("Data/Init_State/", type2, ".csv", sep=""))

### Reverse and fix YCoordinates and Row Numbers
init1 <- lapply(init1, rev_ycoords)
print("step4")
init2 <- lapply(init2, rev_ycoords)
print("step5")
plot(out2raster(out[[1]], var="species"))
#plot(out2raster(init2 [[3]], var="species"))


## Create Mask (see 0_functions file)
mask1 <- create_mask(npatch_col, npatch_row, patch_width, patch_length, outputfile=init2[[1]])
plot(mask1)
print("step6")
## Apply Mask, change Piceabie to Abiealba in Original List
out <- mapply(mask_apply, x=init1, y=init2, SIMPLIFY=F)


### Check if successfull!
#r_out1 <- out2raster(out[[1]], var="species")
#levelplot(r_out1)   ## SUCCESS!
print("step7")
## Write table to /Data Folder
lapply(out, function(x) write.table(x, file=outputfile , row.names=F, col.names=T, sep=",", quote=F))

}
