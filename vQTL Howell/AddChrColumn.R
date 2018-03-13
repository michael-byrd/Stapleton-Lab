library("tidyverse")
#####Issue with only 9 chromosomes showing up#####
outp <- read_csv(file = "HowellvQTL_Ratio_LOD,Pvals,EffectSizes.csv")
pos <- outp$`Position (cM)`
dfr <- pos[-1]- pos[-length(pos)]
fnew <- c(1,which(dfr < 0), length(pos)+1) 
chr <-unlist(sapply(1:(length(fnew)-1), function(x){
  rep(x,fnew[x+1]-fnew[x])
}))
nms <- colnames(outp)
outp <-cbind(outp[,1:3],chr, outp[,4:ncol(outp)])
colnames(outp)[4] = "Chromosome"
write_csv(outp,path = "HowellvQTL_Ratio_LOD,Pvals,EffectSizes.csv")

#C3 Object still has 10 chromosomes#
#Maybe it has 
par(mfrow = c(2,1))
orgr <- read_csv(file = "Howell-Cross-Object-Ratio.csv")
posr <- as.numeric(orgr[2,-1])
plot(posr)
which(posr == max(posr))

vqtlr <- read_csv(file = "HowellvQTL_Ratio_LOD,Pvals,EffectSizes.csv")
posvr <- as.numeric(vqtlr$`Position (cM)`)
plot(posvr)
which(posvr == max(posvr[which(!is.na(posvr))]))