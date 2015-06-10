###########################
## Make Beautiful Images ##
##                       ##
###########################

library(rasterVis)


## Make Timeline Animation

decades <- 1:50
simfolder <- "Simulations/disp_3x4_2x2/Output/"

file_list <- paste(simfolder,"fullOut_",decades,".csv", sep="")
out <- lapply(file_list, biomass_dt)
out_pic <- lapply(out, function(x) out2raster(x, var="biomass_perc"))

casparTheme<-rasterTheme(pch=19, cex=0.7, region=rev(brewer.pal(9, 'BuGn')))
for (i in seq_along(out_pic)) {
  print(i)
  png(file=paste("Animate/pic",i,".png", sep=""), width=600, height=800)
  print(levelplot(out_pic[[i]], margin=F, main=paste(i*10, "Years"), colorkey=list(mypalette), xlab="Longitude", ylab="Latitude", par.settings=casparTheme))  
  dev.off()
}


out <- lapply(file_list, function(x) read.csv(x))

cr_ani<- function (pic=out_pic) {
  for (i in seq_along(out_pic)) {
   print(i)
   #png(file=paste("Animate/pic",i,".png", sep=""), width=600, height=800)
   print(levelplot(out_pic[[i]]))   
   #dev.off()
 }
}

saveHTML(cr_ani())
saveVideo(cr_ani(), ffmpeg="avconv")
saveLatex(cr_ani())