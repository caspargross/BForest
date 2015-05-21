#########################################################################
### Climate data from DWD Weather Station (Source: www.dwd.de/stationsliste)
#########################################################################
library(ModellingTools)

datapath <- paste(getwd(), "/Data/DWD/", sep="")

### DWD Climate Station Feldberg.
dat <- read.table(paste(datapath, "produkt_monat_Monatswerte_19450101_20131231_01346.txt", sep=""), sep=";", dec=".", header=T)
head(dat)

dat <- dat[,c("Mess_Datum_Beginn", "Mess_Datum_Ende","LUFTTEMPERATUR", "NIEDERSCHLAGSHOEHE")]
dat$Mess_Datum_Beginn <- as.Date(as.character(dat$Mess_Datum_Beginn), format = "%Y%m%d")

ym <- format(dat$Mess_Datum_Beginn, "%Y%m")
table(table(ym)==1) #Check for duplicate entries?

dat <- dat[as.numeric(format(dat$Mess_Datum_Beginn, "%Y"))>1945, ] #Data newer than 1945
ym <- format(dat$Mess_Datum_Beginn, "%Y%m")
table(table(ym)==1) #Check for duplicate entries?

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
row.names(clim) <- seq(nrow(clim))
write.table(clim, paste(datapath,"clim_feldberg.txt", sep=""), sep=" ", dec=".", quote = FALSE)


#### Save Climate as .dat with header
header <- readLines(paste(datapath, "climate_template.dat", sep=""))
header <- header[1:14]
header[2] <- "48 #latitude; used for cacluating drought index in model bugmann#"  ## LATITUDE OF STATION
header[3] <-   "1490 #meter a.s.l. #"            ## Altitude of Station
writeLines(header, paste(datapath,"climate_feldberg.dat", sep=""))
write.table(clim, paste(datapath,"climate_feldberg.dat", sep=""), append=T, row.names=T, col.names=F)



