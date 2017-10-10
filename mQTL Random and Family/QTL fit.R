dat = read.csv(file = "Github/Stapleton-Lab/mQTL\ Random\ and\ Family/Input\ Files/Family.csv", header = FALSE, stringsAsFactors = FALSE)
family <-read.cross(dir = "~/PCR\ data/", file = "family.csv")
sigmark <- read.csv(file ="Github/Stapleton-Lab/mQTL\ Random\ and\ Family/Input\ Files/adj_jointP_fam.csv", header = TRUE)
family <- drop.nullmarkers(family)
#scan with variance
family <- calc.genoprob(family)
foutv <- scanonevar(cross = family,
                    mean.formula = PlantHeight ~ mean.QTL.add,
                    var.formula = ~ var.QTL.add,
                    chrs = 1)
fam1 = NULL
family = sim.geno(family)
sigmrk1 = sigmark[-which(match(sigmark$jQTL_ID, "")== 1),]
sig1 = sigmrk1[which(sigmrk1$jQTL_ID == "q1"),3]
fam1 = makeqtl(cross = family, chr = rep(1,length(sig1)), pos = sig1)

sigmrk1 = sigmark[-which(match(sigmark$jQTL_ID, "")== 1),]
sig1 = sigmrk1[which(sigmrk1$jQTL_ID == "q1"),3]
fam1 = makeqtl(cross = family, chr = rep(1,length(sig1)), pos = sig1)

qtlmaker = function(cross,sigmark, chr, index){
  qx = paste("q",index, sep = "")
  famx = paste("fam",index, sep = "")
  sigmrk = sigmark[-which(match(sigmark$jQTL_ID,"") == 1),]
  sig = sigmrk[which(sigmrk1$jQTL_ID == qx),]
  sigpos = sig[,3]
  names = as.character(sig$SNP_Names)
  
  famx = makeqtl(cross = cross, chr = chr, pos = sigpos)
  return(famx)
}
vect = (1:40)[-c(4,18)]
#####Matching chromosome values to QTL indices#####
lapply(vect, function(x){
  assign(paste("famj",x, sep = ""),qtlmaker(cross = family, sigmark = sigmark, chrvect[x], index= x))
})