#####Multiple Stress vQTL#####
library(qtl)
library(vqtl)
cross = read.cross(format = "csv", file = "C:/Users/Thomas/Documents/GitHub/Stapleton-Lab/Manching BayesNet/ManchingScrubbed.csv")
gc()
cgp = calc.genoprob(cross = cross)
scanv = scanonevar(cross = cgp, mean.formula = Height ~ mean.QTL.add + Low.Water + Low.Nitrogen + Pathogen,
                   var.formula =  ~ var.QTL.add + Low.Water + Low.Nitrogen + Pathogen)
