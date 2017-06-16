#Correctly Running FAM/HEX Data
h517<- read.csv(file = "neillh517.csv", header = TRUE, sep = ",")
fluoh<- c( 1+2*1:5, #A
          13+2*1:5, #B
          25+2*1:5, #C
          37+2*1:5, #D
          49+2*1:5, #E
          61+2*1:5, #F
          73+2*1:5) #G
fluog<- c(rep("g1s1",5),
          rep("g1s2",5),
          rep("g1s3",5),
          rep("g1s4",5),
          rep("g1s5",5),
          rep("g1s6",5),
          rep("g1c1",5))
h517mod<- modlist(h517, cyc = 1, fluo = fluoh, model = l4, check ="uni1")
replist(h517mod,gl(7,5), names = "first")
h517pcr<- pcrfit(h517, cyc = 1, fluo = fluoh, model = l4)
best<-mselect(h517pcr, do.all = TRUE)
#still best to use l4
h517KOD<- KOD(h517mod, plot = TRUE, method = "multi3")
is.outlier(h517KOD)
m<-pcrbatch(h517, cyc = 1, fluo = fluoh, model = l4, remove= "none",check = "uni1", plot = TRUE)
rat<- ratiobatch(m,fluog)
rat

#individual curve plotting
fl <- 65
mdl<-pcrfit(h517,1,fl,l4)
efficiency(mdl)


fluof<- sapply(fluoh, function(x) x-1)

f517mod<- modlist(f517, cyc = 1, fluo = fluof, model = l4, check ="uni1")
replist(f517mod,gl(7,5), names = "first")
f517pcr<- pcrfit(f517, cyc = 1, fluo = fluof, model = l4)
efficiency(f517pcr, plot = TRUE, type = "cpD1")
best<-mselect(f517pcr, do.all = TRUE)
#still best to use l4
f517KOD<- KOD(f517mod, plot = TRUE, method = "uni1")
is.outlier(f517KOD)
m<-pcrbatch(f517, cyc = 1, fluo = fluof, model = l4, remove= "none",check = "uni1", plot = TRUE)
rat<- ratiobatch(m,fluog)
rat

#individual curve plotting
fl <- 65
mdl<-pcrfit(h517,1,fl,l4)
efficiency(mdl)
