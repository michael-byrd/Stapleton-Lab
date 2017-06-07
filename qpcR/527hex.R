#5/27 Hex analysis
#Neill 5/25/17 Data Analysis
library(qpcR)
h517<- read.csv(file = "neillh517.csv", header = TRUE, sep = ",")

#Analysis of 3000 level Green Data

fluoh<- c((1+2*1:41),(12*1:8),(1+12*1:8),86:94)
fluoh<-unique(fluoh)
fluoh<-sort(fluoh)
fluog<-c(rep(c(rep("gs",5),rep("rs",2)),6),rep("gs",5),rep("rs",6),rep("rc",5),"rs","rs")
#new analysis
h517mod<- modlist(h517, cyc = 1, fluo = fluoh, model = l4)
h517pcr<- pcrfit(h517, cyc = 1, fluo = fluoh, model = l4)
best<-mselect(h517pcr, do.all = TRUE)
# best model is l4
h517KOD<- KOD(h517mod, method = "uni1", plot = TRUE)
is.outlier(h517KOD)
#Outliers are A4, B4, B6, D11
m<-pcrbatch(h517, cyc = 1, fluo = fluoh, model = l4, remove= "none", plot = TRUE)
rat<-ratiocalc(data = m, group = fluog, which.eff = "sli",type.eff = "individual", which.cp = "cpD2")
rat$summary
#testing indiviual curve fits
fit<-pcrfit(data = h517, cyc = 1, fluo = fluoh, model = l4)
best<-mselect(fit, do.all = TRUE)
plot(best)