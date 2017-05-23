#Generating Known-Truth Logistic Curves
#importing known truth values
ktf<-read.csv(file = "KTF.csv", header = TRUE, sep = ",")
#translating Excel values into R values
upper<-ktf[1,1]
mid<-ktf[1,2]
steep<-ktf[1,3]
randc<-ktf[1,4]
#Generating signal

#Making x values from -6*steepness+mid to 6*steepness+mid
x<-c(mid+steep*(-6:6))

#signal from the logistic model
sig<- function(x,steep = 1, mid = 3, upper = 3) {
  upper/(1+exp(-steep*(x-mid)))
}
#generating y values
signal<-c(sig(x,steep, mid, upper))

#Generating noise
noise <- function(input = 1, randc = .05){
  rnorm(1,0, input*randc)
}
noiseval<-sapply(signal, noise)
y<-c(signal+noiseval)
plot(y)