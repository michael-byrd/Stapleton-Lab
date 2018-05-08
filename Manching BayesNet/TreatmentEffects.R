library(tidyverse)
dat <- read.csv(file = "C://Users/Thomas/Documents/Github/Stapleton-Lab/Manching\ BayesNet/ManchingScrubbed.csv")
dat[1:6,1:6]
as.factor(dat[,2])
mod <- lm(Height ~ Low.Water + Low.Nitrogen + Pathogen, data = dat)
summary(mod)
max(dat$Height[-1:-2])
hist(dat$Height[-1:-2], breaks = 50)
plot(aov(mod))

