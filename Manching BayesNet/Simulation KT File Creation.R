kt <- read.csv(file = "C:/Users/Thomas/Documents/GitHub/Stapleton-Lab/Manching\ BayesNet/SNPData.csv")
ndat <- read.csv(file = "C:/Users/Thomas/Documents/GitHub/Stapleton-Lab/Manching\ BayesNet/newData.csv")
dat <- read.csv(file = "C:/Users/Thomas/Documents/GitHub/Stapleton-Lab/Manching\ BayesNet/ManchingScrubbed.csv")
#dataset with new height data
gdat <- dat
gdat[-1:-2,1] <- ndat$Height

write.csv(gdat, file = "C:/Users/Thomas/Documents/GitHub/Stapleton-Lab/Manching\ BayesNet/SimulatedResponse.csv")

#making a key with sorted NPV names and column numbers
key <- match(as.character(colnames(dat)[-c(1:4)]), as.character(kt$SNP.Name))
keep <- as.character(kt$SNP.Name)[key]
gcol <- which(!is.na(keep)) + 4
gnames <- colnames(gdat[,gcol])
gdf <- data.frame(gnames,gcol)
colnames(gdf) <- c("SNP Name","SNP Column")
write.csv(gdf, file = "C:/Users/Thomas/Documents/GitHub/Stapleton-Lab/Manching\ BayesNet/SimulatedKTSNPs.csv",
          row.names = F)
