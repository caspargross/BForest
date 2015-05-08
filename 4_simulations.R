### LandClim Simulations.
### Klara Dolos 2015

### Simulation with all species ####
sim_name<-"model_simple"

run.lanclim.model<-function(sim_name, ctl_file="ctl_berchtesgaden.xml", lcpath="/home/caspar/Bachelor_Thesis/RCode/Data/Landclim/LandClim"){
oldwd <- getwd()
setwd(paste("Simulations/",sim_var,"/", sep=""))
dir.create("Output")   ### Erstellt das Verzeichnis, in das die Ergebnisse geschrieben werden. Es muss so heißen, wie im LandClim-Controll-File angegeben.
inputpath <-  paste(oldwd, "/Input/", sep="")

setwd(inputpath) ### Damit der "system"-Befehl funktioniert, muss in R das working directory der Ordner "Input" der entsprechenden Simulation sein.
system(paste(lcpath,ctl_file, sep=" "))   ### Die Simulation wird ausgeführt und schreibt die Ergebnisse ins Verzeichnis ../Output.
setwd(oldwd) ### Zurücksetzen des Working Directories.
moveOutput(newname = "Output")
}


