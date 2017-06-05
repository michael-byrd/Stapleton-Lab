#Just checking if the curves look good
library(qpcR)
f517<- read.csv(file = "neillf517.csv", header = TRUE, sep = ",")

#I'd like to look at cells A to G and 1 to 10 as experimental cells
#         A       B        C       D       E       F       G
fluof<- c(2:11 , 14:23 , 26:35 , 38:47 , 50:59 , 62:71 , 74:83)
f517mod<- modlist(f517, cyc = 1, fluo = fluof, model = l5)
plot(f517mod)
f517KOD<- KOD(f517mod, method = "uni1", plot = TRUE)
is.outlier(f517KOD)
#some of these have asterisks, but still aren't labeled as outliers
m1<-pcrfit(data = f517, cyc = 1, fluo =9+36, model = l4)
plot(m1)
m<-pcrbatch(f517, cyc = 1, fluo = fluof, model = l4, check = "uni1", remove= "none", plot = TRUE)

#too large to plot
#how about to go one column at a time
fluof<- c(2,14,26,38,50,62,74) # 1st row
f517mod<- modlist(f517, cyc = 1, fluo = fluof, model = l5)
plot(f517mod)
f517KOD<- KOD(f517mod, method = "uni1", plot = TRUE)
is.outlier(f517KOD)
m<-pcrbatch(f517, cyc = 1, fluo = fluof, model = l4, check = "uni1", remove= "none", plot = TRUE)
