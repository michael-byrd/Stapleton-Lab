library("tidyverse")
fin <- read.csv(file = "C://Users/Thomas/Documents/Github/Stapleton-Lab/vQTL\ Howell/HowellvQTL_Ratio_LOD,Pvals,EffectSizes.csv")
fin <- fin[,-1]
fin[1:6,1:6]
org <- read_csv(file = "C://Users/Thomas/Documents/Github/Stapleton-Lab/vQTL\ Howell/Howell-Cross-ObjectC3.csv")
org[1:6,1:6]
newdat<- as.numeric(org[2,-c(1,2)])
plot(newdat)
plot(fin$Position..cM.)
table(fin$Chromosome)
table(as.numeric(org[1,-c(1,2)]))
unique(as.numeric(org[1,-c(1,2)]))
#18446-19956 on should be Chr 10 not chr 2
co1 <- read_csv(file = "C://Users/Thomas/Documents/Github/Stapleton-Lab/vQTL\ Howell/Howell-Cross-Object.csv")
plot(as.numierc(co1))

unique(as.numeric(org[1,18446:19956]))
org[1,18446:19956] <- 10
org[1:2,1:2] <- ""
write_csv(org,path = "C://Users/Thomas/Documents/Github/Stapleton-Lab/vQTL\ Howell/Howell-Cross-ObjectC3.csv")

ratio <- unlist(map(org[-1:-2,1], as.numeric))/unlist(map(org[-1:-2,2], as.numeric))
names(ratio) <- c()
org[-1:-2,2] <- ratio
org <- org[,-1]
colnames(org)[1] <- "Ratio"
org[1:6,1:6]
write_csv(org,path = "C://Users/Thomas/Documents/Github/Stapleton-Lab/vQTL\ Howell/Howell-Cross-Object-Ratio.csv")
