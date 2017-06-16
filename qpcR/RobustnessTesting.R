#Testing the power of different parameter combinations
#beginning with pcrfit runthrough
upper<-20
mid<-20
steep<-.35
randc<-.05
numrep<-35
#run this through LogCurveRepPower.R
curveseta<-curvesetgen(upper,mid,steep,randc,numrep)
curvesetamod<- modlist(curveseta, cyc = 1, fluo = 2:(numrep+1), model = l5, weights = "1/error^2")
curvesetaKOD<- KOD(curvesetamod, method = "uni2", plot = TRUE)
is.outlier(curvesetaKOD)
#none are false
curvesetafit<-pcrfit(curveseta, cyc = 1, fluo = 2:(numrep+1), model = l5, weights = "1/error^2")
plot(curvesetafit)
goodness<-pcrGOF(curvesetafit)
print(goodness$Rsq)
#R^2 is .99, nice


#making adding one additional curve and retesting
upper<-.5
mid<-20
steep<-.35
randc<-.05
added<-curve(upper,mid,steep,randc)
curveseta <- cbind(curveseta, added[,2])
#rerunning analysis
curvesetamod<- modlist(curveseta, cyc = 1, fluo = 2:(numrep+2), model = l5, weights = "1/error^2")
curvesetaKOD<- KOD(curvesetamod, method = "uni2", plot = TRUE)
is.outlier(curvesetaKOD)
#none are false
curvesetafit<-pcrfit(curveseta, cyc = 1, fluo = 2:(numrep+2), model = l5, weights = "1/error^2")
plot(curvesetafit)
goodness<-pcrGOF(curvesetafit)
print(goodness$Rsq)
