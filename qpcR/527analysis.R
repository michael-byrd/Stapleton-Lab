#Neill 5/25/17 Data Analysis
library(qpcR)
f517<- read.csv(file = "neillf517.csv", header = TRUE, sep = ",")

#Analysis of 3000 level Green Data

fluof<- c(1+2*1:41)
fluof<- fluof
cat(fluof)
fluog<- c(rep("rs",2), rep("gs",5), rep("rs", 3), rep("rc",5))

#new analysis
f517mod<- modlist(f517, cyc = 1, fluo = fluof, model = l5)
f517KOD<- KOD(f517mod, plot = TRUE)
is.outlier(f517KOD)
#so we remove the A11 data
fluof<- fluof[c(-1,-9)]
fluog<- fluog[c(-1,-9)]
m<-pcrbatch(f517, cyc = 1, fluo = fluof, model = l4, remove= "none", plot = TRUE)
rat<-ratiocalc(data = m, group = fluog, which.eff = "sli",type.eff = "mean.single", which.cp = "cpD2")
rat$summary
#testing indiviual curve fits
fit<-pcrfit(data = f517, cyc = 1, fluo = fluof, model = l4)
best<-mselect(fit, do.all = TRUE)
plot(best)