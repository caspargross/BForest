##########################
# Manipulate Input Files #
#                        #
###########################

combine_input <- function (filelist1, filelist2, outputfile, npatch_col, npatch_row, patch_width, patch_length, aui ) {
require(data.table)
require(bit64)
require(raster)

print(aui)
print("step1")
init1 <- lapply(filelist1, function (x) fread(x, integer64="numeric"))

print("step2")
init2 <- lapply(filelist2, function (x) fread(x, integer64="numeric")) 
print("step3")
#type1 <- "mono_piceabie_dec50"
#init1 <- read.csv(paste("Data/Init_State/", type, ".csv", sep=""))

#type2 <-"mono_abiealba_dec50"
#init2 <- read.csv(paste("Data/Init_State/", type2, ".csv", sep=""))

### Reverse and fix YCoordinates and Row Numbers
init1 <- lapply(init1, function(x) rev_ycoordsDT(x, aui))
print("step4")

init2 <- lapply(init2, function(x) rev_ycoordsDT(x, aui))
print("step5")

plot(out2rasterDT(init1 [[3]], var="species"))
#init1 <- lapply(init1, as.data.frame)
#init2 <- lapply(init2, as.data.frame)
#print(class(init1[[1]]))
#print(class(init2[[1]]))
## Create Mask (see 0_functions file)
mask1 <- create_mask(npatch_col, npatch_row, patch_width, patch_length, outputfile=init2[[1]])
print(plot(mask1))
print("step6")

## Apply Mask, change Piceabie to Abiealba in Original List
out <- list
out <- mapply(mask_apply, x=init1, y=init2, MoreArgs=list(mask=mask1), SIMPLIFY=F, USE.NAMES = F )

print(class(out))
print(length(out))
print(class(out[[1]]))
### Check if successfull!
#print(levelplot(out2rasterDT(out[[4]], var="age")))   ## SUCCESS!
print("step7")
## Write table to /Data Folder
for (i in 1:length(out)) {write.table(out[[i]], file=outputfile[[i]], row.names=F, col.names=T, sep=",", quote=F)}
}
