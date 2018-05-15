#####Howell Visualization
library('tidyverse')
library('qtl')
library('vqtl')
setwd('C://Users/Thomas/Documents/Github/Stapleton-Lab/vQTL\ Howell/')
dat <- read.csv(file = "VarLabels.csv")
sig <- which(dat$QTL_IDvarP_fdrTPM !="")
dat$Variance_P_Value_F9[sig]

repr <- unlist(map(1:10, function(x){
  keep <- which(dat$QTL_IDvarP_fdrTPM == x)
  repr <- which(dat$Variance_P_Value_F9[keep] == min(dat$Variance_P_Value_F9[keep]))
  return(keep[repr][1])
}))
as.character(dat$SNP_Name[repr])
cross <- read.cross(file = "Howell-Cross-Object-Ratio.csv")
names <- read.csv(file = "HowellvQTL_Ratio_LOD,Pvals,EffectSizes.csv")
names <- names$SNP.Name
pdf(file = "RatioSNPs.pdf")
map(1:10, function(x){
  mean_var_plot_model_based(cross = cross, phenotype.name = 'Ratio',
                            focal.groups = as.character(names[repr][x]),
                            genotype.names = c("AA","BB"))
})
dev.off()
