library(qtl)
library(vqtl)
library(purrr)
library(readr)
library(dplyr)
setwd("/Users/mbyrd/StapletonLab/Thomas/Stapleton-Lab/Manching BayesNet/JulyWork")

dat <- read_rds("manching_cross.rds")

outv <- read_rds("outv.rds")

p1 <- mean_var_plot_model_based(cross = dat, phenotype.name = "Height", genotype.names = c("AA","BB"), focal.groups = "mmp177c")

plot(p1)

min