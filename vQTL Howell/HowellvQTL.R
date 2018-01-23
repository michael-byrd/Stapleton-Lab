  #####NAM RIL vQTL#####
library("tidyverse")
library("readr")
setwd("C:/Users/Thomas/Documents/GitHub/Stapleton-Lab/vQTL Howell")
dat = read_csv(file = "Heat Stress.csv")
zvals = dat$taxa_NAM_RIL_ID
start = Sys.time()
gbs = read.table(file = "Howell_scrubbed_Z_to_SNPs.txt", header = TRUE, stringsAsFactors = FALSE)
timetaken = Sys.time() - start;timetaken
gbsnames = substr(as.character(gbs$Taxa),1,9)
fullnames = vector()
org  = lapply(1:length(zvals),function(x){
  gbs[which(gbsnames == zvals[x]),1]
})
scrubbed = sapply(org,function(ex){
  cln = gregexpr(":",ex)[[1]]
  bad = c(cln[1],cln[length(cln)])
  paste(substr(ex,1,bad[1]),substr(ex,bad[2]+1,nchar(as.character(ex))), sep = "")
})
write.csv(scrubbed, file = "ScrubbedNames.csv", row.names = FALSE, col.names = FALSE)

#Test1, just assume to use the first of the availabe Z vals with multiple Taxa
used = (unlist(lapply(org, function(x){
  x[1]
})))
usedLoc = sapply(1:length(used), function(x){
  which(gbs$Taxa == used[x])
})
usedDat = gbs[usedLoc,]
usedDat$Taxa = substr(usedDat$Taxa, 1, 9)
precrossObj = cbind(dat[,2:3],usedDat[,-1])
precrossObj = rbind(colnames(precrossObj), precrossObj)
colnames(precrossObj) = NA
crossobj = read.cross(precrossObj)
1+3