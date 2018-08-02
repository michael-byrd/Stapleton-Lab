setwd("/Users/mbyrd/StapletonLab/Thomas/Stapleton-Lab/Manching BayesNet/JulyWork")
dat <- read.csv(file = "ManchingScrubbed.csv")

small_dat <- dat[sample(4:6675, 500),]

small_dat <- rbind(dat[1:3,],small_dat)

small_dat[1,]

# col_samp <- small_dat[, sample(5:3239, 100)]
# 
# 
# small_dat <- cbind(small_dat[1:4],col_samp)

write.csv(small_dat, "small_manching.csv")

