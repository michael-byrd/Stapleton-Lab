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
goodness<-pcrGOF(ucnhl4);print(goodness$Rsq.ad)
plot(ucnhl4)

#now with replicate trials

#exmaple data
mlist<-modlist(reps, 1, 2:5, l5, weights= "1/error^2")
KODlist<- KOD(mlist, method = "uni1")
is.outlier(KODlist)
#We got all false, so we continue normally
m2 <- pcrfit(reps, 1, 2:5, l5, weights = "1/error^2")
mchoice(m2)
pcrGOF(m2)
plot(m2)
#R^2 value is .9997 which is very good
#2:5 is chosen as there are 4 replicates in columns 2 through 5 respectively


#with our data

#UCN FAM data
ucnfrmod<- modlist(ucnf, 2, c(5,11,17,23,29), l5, weights = "1/error^2")
ucnfrKOD<- KOD(ucnfrmod, method = "uni1")
is.outlier(ucnfrKOD)
#We got all false, so we continue normally
ucnfrfit<- pcrfit(ucnf, 2, c(5,11,17,23,29), l5, weights = "1/error^2")
plot(ucnfrfit)
pcrGOF(ucnfrfit)
#this R^2 value is .162, which is awful

#UCN HEX data
ucnhrmod<- modlist(ucnh, 2, c(5,11,17,23,29), l5, weights = "1/error^2")
ucnhrKOD<- KOD(ucnhrmod, method = "uni2")
is.outlier(ucnhrKOD)
#We got all false, so we continue normally
ucnhrfit<- pcrfit(ucnh, 2, c(5,11,17,23,29), l5, weights = "1/error^2")
plot(ucnhrfit)
goodness<-pcrGOF(ucnhrfit);print(goodness$Rsq.ad)
#this R^2 value is .387 which is better, but bad

#Now to test the arithmetically adjusted data
aaucnh<- read.csv(file="UCNnhHexArithAdjust.csv",header=TRUE, sep = ",")
aaucnhmod<- modlist(aaucnh, cyc = 1, fluo = c(7:11), model = l5)
aaucnhKOD<- KOD(aaucnhmod, method = "uni1", plot = TRUE)
is.outlier(aaucnhKOD)
aaucnhfit<-pcrfit(aaucnh, cyc = 1, fluo = c(7,9,10,11), model = l5)
plot(aaucnhfit)
goodness<-pcrGOF(aaucnhfit);print(goodness$Rsq.ad)
plot(aaucnhmod)
#R^2 is .536 which is indeterminant

#With Geometrically Adjusted Hex data
gaucnhmod<- modlist(aaucnh, cyc = 1, fluo = c(12:16), model = l5)
gaucnhKOD<- KOD(gaucnhmod, method = "uni1", plot = TRUE)
is.outlier(gaucnhKOD)
gaucnhfit<-pcrfit(aaucnh, cyc = 1, fluo = c(12,13,14,16), model = l5)
m<-mselect(gaucnhfit,do.all = TRUE)
summary(m)
plot(m)
goodness<-pcrGOF(m);print(goodness$Rsq.ad)
efficiency(m)
#R^2 is .806,nice

#with Logarithmically Adjusted Hex data
laucnhmod<- modlist(aaucnh, cyc = 1, fluo = c(17:21), model = l5)
laucnhKOD<- KOD(laucnhmod, method = "uni1", plot = TRUE)
is.outlier(laucnhKOD)
laucnhfit<-pcrfit(aaucnh, cyc = 1, fluo = c(17,18,19,21), model = l5)
plot(laucnhfit)
goodness<-pcrGOF(laucnhfit);print(goodness$Rsq.ad)
#R2 of .686

#Now trying to use Geom Data with pcrbatch
m<-pcrbatch(aaucnh, cyc = 1, fluo = 12:16, model = l4, check = "uni2", remove= "none", norm = TRUE, plot = TRUE)
mselect(m, do.all = TRUE)

#First attempt at ratiocalc (first merge ucnf, ucnh)
ucn<-read.csv(file="UCN.csv",header=TRUE, sep = ",")
m<-pcrbatch(ucn, cyc = 1, fluo = 2:11, model = b4, check = "uni2",group = NULL, remove= "none", norm = TRUE, plot = TRUE)
rat<-ratiocalc(data = m,group = c(rep("gc",5),rep("gs",5)), which.eff = "sli",type.eff = "mean.single", which.cp = "cpD2")
rat$summary
#testing indiviual curve fits
fit<-pcrfit(data = ucn, cyc = 1, fluo = 7:11, model = l4)
best<-mselect(fit, do.all = TRUE)
#b4 is best for ucnf .08 to l4, l4 is best for ucnh .3 to b4
plot(best)
