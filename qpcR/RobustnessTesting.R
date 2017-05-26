#Testing the power of different parameter combinations
#beginning with pcrfit runthrough
upper<-10
mid<-6
steep<-2
randc<-.05
#run this through LogCurveRepPower.R
curveseta<-curveset
curvesetamod<- modlist(curveseta, cyc = 1, fluo = c(2:13), model = l5)
curvesetaKOD<- KOD(curvesetamod, method = "uni1", plot = TRUE)
is.outlier(curvesetaKOD)
#none are false
curvesetafit<-pcrfit(curveseta, cyc = 1, fluo = c(2:13), model = l5)
plot(curvesetafit)
pcrGOF(curvesetafit)
#R^2 is .99, nice


#making adding one additional curve and retesting
added<-curve(upper,mid,steep,randc)
curveseta <- cbind(curveseta, added[,2])
#rerunning analysis
curvesetamod<- modlist(curveseta, cyc = 1, fluo = c(2:14), model = l5)
curvesetaKOD<- KOD(curvesetamod, method = "uni1", plot = TRUE)
is.outlier(curvesetaKOD)
#none are false
curvesetafit<-pcrfit(curveseta, cyc = 1, fluo = c(2:14), model = l5)
plot(curvesetafit)
pcrGOF(curvesetafit)

#now throw in another curve with different parameters
upper<-13
mid<-6
steep<-2
randc<-.05
#This drops the R^2 value to .976 which is still very strong

#how about we try it again with one that barely took off
upper<-1
mid<-6
steep<-2
randc<-.05
#This drops the R^2 to .85, which is more significant

#how about one that never takes off
upper<-0.1
mid<-6
steep<-2
randc<-.05