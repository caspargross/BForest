################################
## Perform analysis on output ##
################################

## Use Data.Table format for data manipulation
library(data.table)
library(reshape)

#### Aggregate Output with additional values!

DT <- fread("Simulations/disp_3x4_2x2/Output/fullOut_50.csv")
DT <- rev_ycoords(DT)   #reverse Ycoords
      DT[,biomass_cohort:=biomass*stems]     #create cohort_biomass
DT <- DT[,.(biomass_cohort=sum(biomass_cohort)),by=list(species, cell, xcoord, ycoord)]  #
      DT[species=="abiealba", biomass_aa:=biomass_cohort]
      DT[species=="piceabie", biomass_pa:=biomass_cohort]
DT <- DT[,.(biomass_cell=sum(biomass_cohort), biomass_aa=sum(biomass_aa, na.rm=T), biomass_pa=sum(biomass_pa, na.rm=T) ), by=.(cell, xcoord, ycoord)]
      DT[,biomass_perc:=((biomass_aa/biomass_cell)*100)]



##### As function:



plot(out2raster(data.frame(DT), var="biomass_perc"))
out2raster(data.frame(DT), var="biomass_perc")

class(DT$biomass_quo)
View(DT)

head(melt(agg, id=c("cell", "biomass_cohort"), masure=species))
DT[species=="abiealba",.(biomass_cohort=sum(biomass_cohort)),by=list(species, )]
agg<-DT[,.(biomass_cohort=sum(biomass_cohort)),by=list(species)]


data[, list(metric=sum(metric)), by=list(a,b)]

f_out[c(""),1]
f_out$biomass_cohort <- f_out$biomass * f_out$stems


dcast(DT, c("species", "col"), "biomass_cohort")


f_agg <- aggregate(f_out$biomass_cohort, list(species=f_out$species, cell=f_out$cell, xcoord=f_out$xcoord, ycoord=f_out$ycoord  ), sum)

f_agg <- sapply(f_agg, foo)
foo <- function (x) {
  if (x["species"]=="abiealba") {
  x["species"]<-"abiesalbabum"
  }              
}

f_agg[,"x"] <- f_agg[,"x"] / (ifelse(duplicated(f_agg$cell),f_agg$x,1))
f_quo <- subset(f_agg, species=="piceabie", select=c("xcoord", "ycoord", "x"))
f_quo[,"x"] <- f_quo[,"x"] / (ifelse(duplicated(f_agg$cell),f_agg$x,1))

plot(out2raster(f_quo, var="x"))



head(f_agg[f_agg$species=="piceabie", "x"]/(ifelse(duplicated(f_agg$cell),f_agg$x,1)))

head(f_agg[duplicated(f_agg$cell),"x"])
head(subset(f_agg, species=="abiealba")[,"x"])
head(subset(f_agg, species=="piceabie")[,"x"])

foo <- function (x) {}

duplicated(f_agg$cell)



f_agg <- aggregate(f_agg$x, list(cell=f_out$cell, xcoord=f_out$xcoord, ycoord=f_out$ycoord  ), mean)

View(f_agg)

