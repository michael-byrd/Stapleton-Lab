
setwd("~/StapletonLab/Thomas_Stapleton_Lab_Forked/Stapleton-Lab/Manching BayesNet")
dat <- read.csv(file = "SimulatedResponse.csv")

small_dat <- dat[sample(4:6675, 5),]




small_dat <- rbind(dat[1:3,],small_dat)

small_dat[1,]

write.csv(small_dat, "small_dat.csv")
