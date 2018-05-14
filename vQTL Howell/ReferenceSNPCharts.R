#####Howell Visualization
library('tidyverse')
library('qtl')
library('vqtl')
dat <- read.csv(file = "C://Users/Thomas/Documents/Github/Stapleton-Lab/vQTL\ Howell/HowellvQTL_Ratio_LOD,Pvals,EffectSizes.csv")
dat <- dat[,-1]
plot(dat$Mean.P.Value)
ggplot(data = dat) +
  geom_point(aes(x = Position..cM., y = Mean.P.Value)) + 
  geom_hline(yintercept = .05, color = 'red')

ggplot(data = dat) +
  geom_point(aes(x = Position..cM., y = Variance.P.Value)) +
  geom_hline(yintercept = .05, color = 'red')

ggplot(data = dat) +
  geom_point(aes(x = Position..cM., y = Joint.P.Value)) +
  geom_hline(yintercept = .05, color = 'red')

?smooth
nas <- which(is.na(dat$Joint.P.Value))

yval <- smooth(dat$Joint.P.Value[-nas], kind = '3')

ggplot() +
  geom_point(aes(x = dat$Position..cM.[-nas], y = as.numeric(yval))) +
  geom_hline(yintercept = .05, color = 'red')

sig <- which(yval < .05)
repr <- c()
which(yval[sig][1:4] == min(yval[sig][1:4]))
repr <- c(repr,1)
which(yval[sig][5:6] == min(yval[sig][5:6]))
repr <- c(repr,5)
repr <- c(repr,7)
which(yval[sig][8:9] == min(yval[sig][8:9]))
repr <- c(repr,8)
which(yval[sig][10:11] == min(yval[sig][10:11]))
repr <- c(repr,10)
which(yval[sig][12:15] == min(yval[sig][12:15]))
repr <- c(repr,12)
which(yval[sig][16:17] == min(yval[sig][16:17]))
repr <- c(repr,16)

sig[repr]

cross <- read.cross(file = "C://Users/Thomas/Documents/Github/Stapleton-Lab/vQTL\ Howell/Howell-Cross-Object-Ratio.csv")
pdf(file = "C://Users/Thomas/Documents/Github/Stapleton-Lab/vQTL\ Howell/RatioSNPs.pdf")
map(sig[repr], function(x){
  mean_var_plot_model_based(cross = cross, phenotype.name = 'Ratio',
                            focal.groups = as.character(dat$SNP.Name[-nas][x]),
                            genotype.names = c("AA","BB"))
})
dev.off()
