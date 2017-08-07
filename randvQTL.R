install.packages("vqtl")
install.packages("qtl")
library("qtl")
library("vqtl")
# random <-read.cross(dir = "~/PCR\ data/", file = "Random2.csv")
# random <- drop.nullmarkers(random)
# #scan with variance
# random <- calc.genoprob(random)
# routv <-scanonevar(random)
# random= sim.geno(random)
# #trying to script fitqtls for all the markers
# varlist<- paste("c",1:10, sep = "")
# values <- rep(1,10)
# df<-data.frame(varlist,values)
# pval <- rep(2,3235)
# names <- rep("empty",3235)
# names[1:531] <- colnames(random$geno$`1`$data)
# names[532:879] <- colnames(random$geno$`2`$data)
# names[880:1263] <- colnames(random$geno$`3`$data)
# names[1264:1583] <- colnames(random$geno$`4`$data)
# names[1584:1908] <- colnames(random$geno$`5`$data)
# names[1909:2205] <- colnames(random$geno$`6`$data)
# names[2206:2475] <- colnames(random$geno$`7`$data)
# names[2476:2734] <- colnames(random$geno$`8`$data)
# names[2735:3003] <- colnames(random$geno$`9`$data)
# names[3004:3235] <- colnames(random$geno$`10`$data)
# for (x in 1:531) {
#   for (y in 1:10) {
#     if (length(random$geno[[y]]$data[1,]) >= x) {
#       df[y,2] <- x
#     }
#     else {
#       df[y,2] <- 1
#     }
#   }
#   poslist<- c(routv[colnames(random$geno$`1`$data)[df[1,2]],2],
#               routv[colnames(random$geno$`2`$data)[df[2,2]],2],
#               routv[colnames(random$geno$`3`$data)[df[3,2]],2],
#               routv[colnames(random$geno$`4`$data)[df[4,2]],2],
#               routv[colnames(random$geno$`5`$data)[df[5,2]],2],
#               routv[colnames(random$geno$`6`$data)[df[6,2]],2],
#               routv[colnames(random$geno$`7`$data)[df[7,2]],2],
#               routv[colnames(random$geno$`8`$data)[df[8,2]],2],
#               routv[colnames(random$geno$`9`$data)[df[9,2]],2],
#               routv[colnames(random$geno$`10`$data)[df[10,2]],2])
#   qtl <- makeqtl(random, chr = 1:10, pos = poslist, what = "prob")
#   fit <- fitqtl(random, qtl = qtl, get.ests = TRUE)
#   sum <- summary(fit)
#   pval[x] <- sum$result.drop[1,7]
#   pval[df[2,2]+531] <- sum$result.drop[2,7]
#   pval[df[3,2]+879] <- sum$result.drop[3,7]
#   pval[df[4,2]+1263] <- sum$result.drop[4,7]
#   pval[df[5,2]+1583] <- sum$result.drop[5,7]
#   pval[df[6,2]+1908] <- sum$result.drop[6,7]
#   pval[df[7,2]+2205] <- sum$result.drop[7,7]
#   pval[df[8,2]+2475] <- sum$result.drop[8,7]
#   pval[df[9,2]+2734] <- sum$result.drop[9,7]
#   pval[df[10,2]+3003] <- sum$result.drop[10,7]
# }
# allpos<-routv[,2]
# pval <- data.frame(names,allpos,pval)
# pval
# write.csv(pval, file = "RandomPvalues.csv")
# max(pval[,3])
# 
# #Creating a Manhattan plot
# #   Specifically using the qqman package
# CHR = rep(1,length(pval[,1]))
# plotframe <- data.frame(SNP = pval[,1], BP = pval[,2], P = pval[,3], CHR = CHR)
# manhattan(plotframe)

#trying to actually run this correctly
random <-read.cross(dir = "~/PCR\ data/", file = "Random2.csv")
random <- drop.nullmarkers(random)
scan <- scanonevar(cross = random,
                   mean.formula = height.in. ~ mean.QTL.add + mean.QTL.dom,
                   var.formula = ~ var.QTL.add + var.QTL.dom,
                   chrs = 1)
perm = scanonevar.perm(scan,100, n.cores = 2)

#now doing it will all 10

scanall <- scanonevar(cross = random,
                   mean.formula = height.in. ~ mean.QTL.add + mean.QTL.dom,
                   var.formula = ~ var.QTL.add + var.QTL.dom)
perm = scanonevar.perm(scanall,100, n.cores = 2)

#finding important SNPs for mean, variance or both
#current functionality:
#Inputs, perm, threshold, type
#type legend
# mean = mean.empir.p
# var = var.empir.p
# both = joint.empir.p
SNP <- function(perm, type, thresh) {
  values <- perm$result[[type]]
  thresh = thresh
  newv = NULL
  for (x in values) {
    if (x <= thresh) {
      newv = c(newv,x)
    }
  } 
  return(newv)
}

#finding significant joint SNPS for different thresholds
sigj<-SNP(perm = perm, type = "joint.empir.p", thresh = .001)

#matching this for SNP index
SNP.match <- function(perm, sig, type) {
  markers = NULL
  for (x in sig) {
    markers = c(markers,match(x, perm$result[[type]]))
  }
  return(markers)
} 
#continuing example
matchj <- SNP.match(perm = perm, sig = sigj, type = "joint.empir.p")

#matching indeces to names
SNP.name <- function(perm, match, type) {
  snp.names = NULL
  for (x in match) {
    snp.names = c(snp.names, perm$result$loc.name[x])
  }
  return(snp.names)
}
#example
namesj <- SNP.name(perm = perm, match = matchj, type = "joint.empir.p")

#there seems to be issues with it adding markers we don't have
#let's try getting rid of these extra ones
SNP.scrub <- function(names) {
  scrubbed.names = namesj
  bad = grep("_",names)
  return(scrubbed.names[-bad])
}
scrub <- SNP.scrub(names = namesj)
#trying to plot these SNPs
mean_var_plot_model_based(cross = random,
                                phenotype.name = "height.in.",
                                genotype.names = c("AA","BB"),
                                focal.groups = scrub[2])
#exporting plots
x = 1 
for (x in 1:length(scrub)) {
  pdf(file = paste("VQTLplots/",scrub[x],".pdf",sep = ""))
  print(mean_var_plot_model_based(cross = random,
                                   phenotype.name = "height.in.",
                                   genotype.names = c("AA","BB"),
                                   focal.groups = scrub[x])
  )
  dev.off()
}
dev.cur()
