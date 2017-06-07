#Trying the 527 HEx data with ratiobatch instead:
library(qpcR)
h517<- read.csv(file = "neillh517.csv", header = TRUE, sep = ",")

#Analysis of 3000 level Green Data

fluoh<- c(2:94,96,97)
fluog<-c(rep(c("g1s1","g2s1"),5),rep("r1s3",2),
         rep(c("g1s2","g2s2"),5),rep("r1s4",2),
         rep(c("g1s3","g2s3"),5),rep("r1s5",2),
         rep(c("g1s4","g2s4"),5),rep("r1s6",2),
         rep(c("g1s5","g2s5"),5),rep("r1s3",2),
         rep(c("g1s6","g2s6"),5),rep("r1s4",2),
         rep(c("g1c1","g2c1"),5),rep("r1s5",2),
         "r1s3","r1s4","r1s5","r1s6",rep("r1c1",5),rep("r1s6",2))
h517mod<- modlist(h517, cyc = 1, fluo = fluoh, model = l4)
h517pcr<- pcrfit(h517, cyc = 1, fluo = fluoh, model = l4)
best<-mselect(h517pcr, do.all = TRUE)
#still best to use l4
h517KOD<- KOD(h517mod, plot = TRUE)
is.outlier(h517KOD)
#A4, B4, B6, and D11 are outliers, and will be removed in ratiobatch()
m<-pcrbatch(h517, cyc = 1, fluo = fluoh, model = l4, remove= "none", plot = TRUE)
rat<-ratiobatch(data = m, group = fluog)

#maybe the issue is with unequal number of treatment and control samples
#to remedy this, try C1 to C9 odd compared to c(A11,A12,E11,E12,H1)
fluot<- c(26,28,30,32,34,12,13,60,61,86)
h517mod<- modlist(h517, cyc = 1, fluo = fluot, model = l4)
h517pcr<- pcrfit(h517, cyc = 1, fluo = fluot, model = l4)
best<-mselect(h517pcr, do.all = TRUE)
#still best to use l4
h517KOD<- KOD(h517mod, plot = TRUE)
is.outlier(h517KOD)
#A4, B4, B6, and D11 are outliers, and will be removed in ratiobatch()
m<-pcrbatch(h517, cyc = 1, fluo = fluot, model = l4, remove= "KOD", plot = TRUE)
#only ones left are C1:9 odd, A11,12 E11,H1
fluog<- c(rep("gs",5),rep("rs",4))
rat<-ratiocalc(data = m, group = fluog)
#this doesn't work either