library("tidyverse")
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
