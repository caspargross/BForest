#############################################################################
### Climate data from DWD Weather Station (Source: www.dwd.de/stationsliste)
#############################################################################
library(ModellingTools)

datapath <- paste(getwd(), "/Data/DWD/", sep="")

### DWD Climate Station Feldberg.
dat <- read.table("~/Bachelor_Thesis/BForest/Data/DWD/produkt_monat_Monatswerte_19450101_20131231_01346.txt", sep=";", dec=".", header=T)

### DWD Climate Station Wendelstein
dat <- read.csv("~/Bachelor_Thesis/BForest/Data/DWD/produkt_monat_Monatswerte_19510101_20120920_05467.txt", sep=";")
names(dat)[2]<-"Mess_Datum_Beginn"
names(dat)[3]<-"Mess_Datum_Ende"

head(dat)

dat <- dat[,c("Mess_Datum_Beginn", "Mess_Datum_Ende","LUFTTEMPERATUR", "NIEDERSCHLAGSHOEHE")]
dat$Mess_Datum_Beginn <- as.Date(as.character(dat$Mess_Datum_Beginn), format = "%Y%m%d")

ym <- format(dat$Mess_Datum_Beginn, "%Y%m")
table(table(ym)==1) #Check for duplicate entries?

## Cleaning for Feldberg
dat <- dat[as.numeric(format(dat$Mess_Datum_Beginn, "%Y"))>1945, ] #Data newer than 1945
ym <- format(dat$Mess_Datum_Beginn, "%Y%m")
table(table(ym)==1) #Check for duplicate entries?

## Cleaning for Wendelstein
dat <- dat[as.numeric(format(dat$Mess_Datum_Beginn, "%Y"))<2012, ] #Data newer than 1945
ym <- format(dat$Mess_Datum_Beginn, "%Y%m")
table(table(ym)==1) #Check for duplicate entries

### Make climate matrix
head(dat)

dat$Y <- as.numeric(format(dat$Mess_Datum_Beginn, "%Y"))
dat$m <- as.numeric(format(dat$Mess_Datum_Beginn, "%m"))

head(dat)
dat <- na.omit(dat)

temp <- reshape(dat[,c("Y", "m", "LUFTTEMPERATUR")],  timevar= "m", idvar="Y", direction = "wide") 
rownames(temp) <- temp[,"Y"]
temp <- temp[,-1]
colnames(temp) <- 1:12


prec <- reshape(dat[,c("Y", "m", "NIEDERSCHLAGSHOEHE")],  timevar= "m", idvar="Y",direction = "wide") 
rownames(prec) <- prec[,"Y"]
prec <- prec[,-1]
colnames(prec) <- 1:12
temp
prec

clim <- cbind(temp, prec)

### Shuffle climate
csamp <- sample(1:nrow(clim), 5000, replace=T)
clim <-clim[csamp,]


##### FOR FELDBERG
row.names(clim) <- seq(nrow(clim))
write.table(clim, paste(datapath,"clim_feldberg.txt", sep=""), sep=" ", dec=".", quote = FALSE, col.names=F,)
header <- readLines(paste(datapath, "climate_template_feldberg.dat", sep=""))
header <- header[1:14]
writeLines(header, paste(datapath,"climate_feldberg.dat", sep=""))
write.table(clim, paste(datapath,"climate_feldberg.dat", sep=""), append=T, row.names=T, col.names=F, quote=F)



#### FOR WENDELSTEIN
row.names(clim) <- seq(nrow(clim))
write.table(clim, paste(datapath,"clim_wendelstein.txt", sep=""), sep=" ", dec=".", quote = FALSE, col.names=F,)
header <- readLines(paste(datapath, "climate_template_wendelstein.dat", sep=""))
header <- header[1:14]
writeLines(header, paste(datapath,"climate_wendelstein.dat", sep=""))
write.table(clim, paste(datapath,"climate_wendelstein.dat", sep=""), append=T, row.names=T, col.names=F, quote=F)


## Function to create random samples of weather (pseudo-random)
random_weather_bf <- function (file, yr=5000, c=clim) {
 
    csamp <- sample(1:nrow(c), yr, replace=T)
    c <- c[csamp,]
    row.names(c) <- seq(nrow(c))
    write.table(c, "Data/DWD/clim_feldberg.txt", sep=" ", dec=".", quote = FALSE, col.names=F,)
    header <- readLines("Data/DWD/climate_template_feldberg.dat")
    header <- header[1:14]
    writeLines(header, file)
    write.table(c,  file, append=T, row.names=T, col.names=F, quote=F)
  } 

random_weather_bf 


