  #####NAM RIL vQTL#####
setwd("C:/Users/Thomas/Documents/GitHub/Stapleton-Lab/vQTL Howell")
dat = read.csv(file = "Heat Stress.csv")
zvals = as.character(dat[,1])
gbs = read.table(file = "Howell_scrubbed_Z_to_SNPs.txt", header = TRUE)
gbsnames = as.character(gbs[,3])
fullnames = vector()
org  = sapply(1:length(zvals),function(x){
  as.character(gbs[which(gbsnames == zvals[x]),1])
})
org = unlist(org)
ex = org[1]
cln = gregexpr(":",ex)[[1]]
bad = c(cln[1],cln[length(cln)])
rem = substr(ex)
paste(substr(ex,1,bad[1]),substr(ex,bad[2]+1,nchar(ex)), sep = "")

scrubbed = sapply(org,function(ex){
  cln = gregexpr(":",ex)[[1]]
  bad = c(cln[1],cln[length(cln)])
  paste(substr(ex,1,bad[1]),substr(ex,bad[2]+1,nchar(ex)), sep = "")
})
write.csv(scrubbed, file = "ScrubbedNames.csv", row.names = FALSE, col.names = FALSE)