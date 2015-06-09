###########################
## Make Beautiful Images ##
##                       ##
###########################

library(rasterVis)
library()


## Make Timeline

decades <- 1:50
simfolder <- "Simulations/disp_3x4_2x2/Output/"

file_list <- paste(simfolder,"fullOut_",decades,".csv", sep="")

out <- lapply(file_list, function(x) read.csv(x))
 lapply(out, rev_ycoords)
out_pic <- lapply(out, function(x) out2raster(x, var="species"))
 
for (i in seq_along(out_pic)) {
  print(i)
  png(file=paste("Animate/pic",i,".png", sep=""), width=600, height=800)
  print(levelplot(out_pic[[i]]))   
  dev.off()
}
 
cr_ani<- function (pic=out_pic) {
  for (i in seq_along(out_pic)) {
   print(i)
   #png(file=paste("Animate/pic",i,".png", sep=""), width=600, height=800)
   print(levelplot(out_pic[[i]]))   
   #dev.off()
 }
}

saveHTML(cr_ani())