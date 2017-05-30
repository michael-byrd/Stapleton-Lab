#Generating Known-Truth Logistic Curves
#Altering LogCurveGeneration to be selfcontained

#parameters
upper<-12
mid<-6
steep<-1
randc<-.08
numrep<- 12
#defined curve function
curve<- function(upper,mid,steep, randc){
#Generating signal
#Making x values from -6*steepness+mid to 6*steepness+mid
x<-c(mid+(1/steep)*(-6:6))
#signal from the logistic model
sig<- function(x,steep = 1, mid = 3, upper = 3) {
  upper/(1+exp(-steep*(x-mid)))
}
#generating y values
#signal
signal<-c(sig(x,steep, mid, upper))
#Generating noise
noiseval<-sapply(signal, function(signal) rnorm(n = 1,mean = 0,sd = signal*randc))
#Adding noise and signal to produce stochastic curve
y<-c(signal+noiseval)
#preparing column names
curvesdf<- data.frame(x,y)
return(curvesdf)
}
#generate whole set of curve
curvesetgen<-function(upper,mid,steep,randc,numrep){
#base case set up
curvesdf<-curve(upper,mid,steep,randc)
curveset<-curvesdf
i=2

#recursion for upper number of trials
for(i in 2:numrep){
  #increase count of i
  i<- i+1
  #curve dataframe for next individual curve
  curvesdf<- curve(upper,mid,steep,randc)
  #add fluo column to existing dataframe
  curveset <- cbind(curveset, curvesdf[,2])
  #rename columns
  j=1
  colnames(curveset)<- "Cycle"
  for (j in 2:i) {
    colnames(curveset)[j] <-paste("Fluorescense", j-1)
  }
}
return(curveset)
}
#final data
curveset<-curvesetgen(upper,mid,steep, randc,numrep)