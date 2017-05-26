#Cleaner qpcR workflow
library(qpcR)
#data without repetitions (1 curve fit)

#input example data
reps <- read.csv(file="reps.csv",header = TRUE, sep = ",")
## Simple l4 fit of F1.1 of the reps dataset.
m1 <- pcrfit(reps, 1, 2, l4)
plot(m1)

#trying with our data
#ucnh stands for ucn- fam sheet data
ucnh<- read.csv(file="UCN-Q.csv",header=TRUE, sep = ",")
#mufl4 stands for Model of No-Hormone Fam Logistic 4-paramter model
ucnhl4<- pcrfit(ucnh, 2, 6, l4)
pcrGOF(ucnhl4)
plot(ucnhl4)

#same process with hex data
ucnh<- read.csv(file="UCN-QH.csv", header = TRUE, sep = ",") # as above but hex data
ucnhl4<-pcrfit(ucnh, 2, 6, l4)
pcrGOF(ucnhl4)
plot(ucnhl4)

#now with replicate trials

#exmaple data
mlist<-modlist(reps, 1, 2:5, l5, weights= "1/error^2")
KODlist<- KOD(mlist, method = "uni1")
is.outlier(KODlist)
#We got all false, so we continue normally
m2 <- pcrfit(reps, 1, 2:5, l5, weights = "1/error^2")
pcrGOF(m2)
plot(m2)
#R^2 value is .9997 which is very good
#2:5 is chosen as there are 4 replicates in columns 2 through 5 respectively


#with our data

#UCN FAM data
ucnfrmod<- modlist(ucnf, 2, c(6*1:5), l5, weights = "1/error^2")
ucnfrKOD<- KOD(ucnfrmod, method = "uni1")
is.outlier(ucnfrKOD)
#We got all false, so we continue normally
ucnfrfit<- pcrfit(ucnf, 2, c(6*1:5), l5, weights = "1/error^2")
plot(ucnfrfit)
pcrGOF(ucnfrfit)
#this R^2 value is .162, which is awful

#UCN HEX data
ucnhrmod<- modlist(ucnh, 2, c(6*1:5), l5, weights = "1/error^2")
ucnhrKOD<- KOD(ucnhrmod, method = "uni1")
is.outlier(ucnhrKOD)
#We got all false, so we continue normally
ucnhrfit<- pcrfit(ucnh, 2, c(6*1:5), l5, weights = "1/error^2")
plot(ucnhrfit)
pcrGOF(ucnhrfit)
#this R^2 value is .387 which is better, but bad

#Now to test the arithmetically adjusted data
aaucnh<- read.csv(file="UCNnhHexArithAdjust.csv",header=TRUE, sep = ",")
aaucnhmod<- modlist(aaucnh, cyc = 1, fluo = c(7:11), model = l5)
aaucnhKOD<- KOD(aaucnhmod, method = "uni1", plot = TRUE)
is.outlier(aaucnhKOD)
aaucnhfit<-pcrfit(aaucnh, cyc = 1, fluo = c(7,9,10,11), model = l5)
plot(aaucnhfit)
pcrGOF(aaucnhfit)
plot(aaucnhmod)
#R^2 is .548 which is indeterminant

#With Geometrically Adjusted Hex data
gaucnhmod<- modlist(aaucnh, cyc = 1, fluo = c(12:16), model = l5)
gaucnhKOD<- KOD(gaucnhmod, method = "uni1", plot = TRUE)
is.outlier(gaucnhKOD)
gaucnhfit<-pcrfit(aaucnh, cyc = 1, fluo = c(12,13,14,16), model = l5)
plot(gaucnhfit)
pcrGOF(gaucnhfit)
#R^2 is .805, nice