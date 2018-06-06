
setwd("~/StapletonLab/Thomas_Stapleton_Lab_Forked/Stapleton-Lab/Manching BayesNet")
dat <- read_csv(file = "SimulatedResponse.csv")

small_dat <- dat[sample(4:6675, 100),]




small_dat <- rbind(dat[1:3,],small_dat)

small_dat[1,]

write_csv(small_dat, "small_dat.csv")
