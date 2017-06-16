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
h517mod<- modlist(h517, cyc = 1, fluo = fluoh, model = l4, check ="uni1")
h517pcr<- pcrfit(h517, cyc = 1, fluo = fluoh, model = l4)
best<-mselect(h517pcr, do.all = TRUE)
#still best to use l4
h517KOD<- KOD(h517mod, plot = TRUE, method = "uni1")
is.outlier(h517KOD)
#A4, B4, B6, and D11 are outliers, and will be removed in ratiobatch()
m<-pcrbatch(h517, cyc = 1, fluo = fluoh, model = l4, remove= "none",check = "uni1", plot = TRUE)
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


#trying to run this with only 16 trials in any group (based on paper involving REST software)
fluoa<- c(12,13,
          24,25,
          26:33,36,37,
          38:45,48,49,
          50:57,60,61,
          62:69,72,73,
          74:81,84,85,
          90:93,96,97)
fluog<- c(rep("r1s1",2),
          rep("r1s2",2),
          rep(c("g1s1","g2s1"),4),rep("r1s3",2),
          rep(c("g1s2","g2s2"),4),rep("r1s4",2),
          rep(c("g1s3","g2s3"),4),rep("r1s1",2),
          rep(c("g1s4","g2s4"),4),rep("r1s2",2),
          rep(c("g1c1","g2c1"),4),rep("r1s3",2),
          rep("r1c1",4),rep("r1s4",2))
h517mod<- modlist(h517, cyc = 1, fluo = fluoa, model = l4)
h517pcr<- pcrfit(h517, cyc = 1, fluo = fluoa, model = l4)
best<-mselect(h517pcr, do.all = TRUE)
#still best to use l4
h517KOD<- KOD(h517mod, plot = TRUE)
is.outlier(h517KOD)
#A4, B4, B6, and D11 are outliers, and will be removed in ratiobatch()
m<-pcrbatch(h517, cyc = 1, fluo = fluoa, model = l4, remove= "none", plot = TRUE)
rat<-ratiocalc(data = m, group = fluog)


#maybe trying just one gene at a time
#first with the IRE
fluob<- c(12,13,
          24,25,
          c(24+2*1:4),36,37,
          c(36+2*1:4),48,49,
          c(48+2*1:4),60,61,
          c(60+2*1:4),72,73,
          84,55,
          96,97)
fluog<- c(rep("r1s1",2),
          rep("r1s2",2),
          rep("g1s1",4),rep("r1s3",2),
          rep("g1s2",4),rep("r1s4",2),
          rep("g1s3",4),rep("r1s1",2),
          rep("g1s4",4),rep("r1s2",2),
          rep("r1s3",2),
          rep("r1s4",2))
h517mod<- modlist(h517, cyc = 1, fluo = fluob, model = l4)
h517pcr<- pcrfit(h517, cyc = 1, fluo = fluob, model = l4)
best<-mselect(h517pcr, do.all = TRUE)
#still best to use l4
h517KOD<- KOD(h517mod, plot = TRUE)
is.outlier(h517KOD)
#A4, B4, B6, and D11 are outliers, and will be removed in ratiobatch()
m<-pcrbatch(h517, cyc = 1, fluo = fluob, model = l4, remove= "none", norm = TRUE, plot = TRUE)
rat<-ratiocalc(data = m, group = fluog)
