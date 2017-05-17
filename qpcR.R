#qpcR logistic PCR curve fitting and PCR ratios
#download UCN-Q data as UCN-Q.csv for the FAM data
#download UCN-Q data as UCN-QH.csv for the HEX data
library(qpcR)
reps <- read.csv(file="reps.csv",header = TRUE, sep = ",")
## Simple l4 fit of F1.1 of the reps dataset.
m1 <- pcrfit(reps, 1, 2, l4)
plot(m1)
#in pcr fit data syntax is pcrfit(data, column with cycle data, colum with florescense data, model)

#trying with our data
#ucnf stands for ucn- fam sheet data
ucnf<- read.csv(file="UCN-Q.csv",header=TRUE, sep = ",")
#mufl4 stands for Model of No-Hormone Fam Logistic 4-paramter model
munfl4<- pcrfit(ucnf, 2, 6, l4)
plot(munfl4)

#trying with unc-q Hex data sheet
ucnh<- read.csv(file="UCN-QH.csv", header = TRUE, sep = ",") # as above but hex data
munhl4<-pcrfit(ucnh, 2, 6, l4)
plot(munhl4)

#now working with replicate models
#example with reps data
m2 <- pcrfit(reps, 1, 2:5, l5, weights = "1/error^2")
plot(m2)
#2:5 is chosen as there are 4 replicates in columns 2 through 5 respectively

#with our data
mufl5<- pcrfit(ucnf, 2, c(6,12,18,24,30), l5, weights = "1/error^2")
plot(mufl5)
#this data has like 3 out of 5 sets which don't "take off"

muhl5<- pcrfit(ucnh, 2, c(6,12,18,24,30), l5, weights = "1/error^2")
plot(muhl5)
#it seems like D8 did not grow, and the max diff in C8 is low as well
plot(x=ucnh[,2], y=ucnh[,6])
plot(x=ucnh[,2], y=ucnh[,30])
