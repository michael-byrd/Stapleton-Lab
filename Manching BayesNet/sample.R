######
library("qtl")
library("vqtl")
dat <- read.csv(file = "C:/Users/Thomas/Documents/GitHub/Stapleton-Lab/Manching\ BayesNet/SimulatedResponse.csv")
dat[1:6,1:6]
dat[1:2,1:4] <- ""
dat[1:6,1:6]
m <- dim(dat)[1]
n <- dim(dat)[2]
