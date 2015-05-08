#########################
### R Script Bachelor ###
###     Caspar        ###
#########################

setwd("/home/caspar/Bachelor_Thesis/RCode/")
source("0_functions.R")




# Monocultures
##############

### Abies Alba
create.inputdir("mono_abiealba", species="abiealba")
run.landclim.model("model_simple")






function(file) {
  doc <- xmlTreeParse(file)
  daten <- t(xmlSApply(xmlRoot(doc), function(x) xmlSApply(x, xmlValue)))
  rownames(daten) <- NULL
  daten <- data.frame(daten)
  tmpfile <- file()
  write.table(daten, tmpfile)
  daten <- read.table(tmpfile)
  daten
  #  temp
}

doc <- xmlTreeParse("Input/species.xml")
doc
daten <- t(xmlSApply(xmlRoot(doc), function(x) xmlSApply(x, xmlValue)))
rownames(daten) <- NULL
daten <- data.frame(daten)
daten2<-daten[2:3,]
xmlSApply(xmlRoot(doc), function(x) xmlSApply(x, xmlValue))<-daten

doc2 <- xmlParse("Input/species.xml")
write.xml(daten2, file="species_test.xml")


write.species.xml(daten2, file="species.text.xml")