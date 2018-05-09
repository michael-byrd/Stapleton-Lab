library("qtl")
library("vqtl")
setwd("C://Users/Thomas/Documents")
dat = read.csv(file = "Github/Stapleton-Lab/mQTL\ Random\ and\ Family/Input\ Files/Family.csv", header = FALSE, stringsAsFactors = FALSE)
family <-read.cross(dir = "~/PCR\ data/", file = "family.csv")
sigmark <- read.csv(file ="Github/Stapleton-Lab/mQTL\ Random\ and\ Family/Input\ Files/adj_jointP_fam.csv", header = TRUE)
family <- drop.nullmarkers(family)
#scan with variance
family <- calc.genoprob(family)
foutv <- scanonevar(cross = family,
                    mean.formula = PlantHeight ~ mean.QTL.add,
                    var.formula = ~ var.QTL.add)
fam1 = NULL
family = sim.geno(family)
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
  
  famx = makeqtl(cross = cross, chr = rep(chr,length(sigpos)), pos = sigpos)
  return(famx)
}
vect = (1:40)[-c(4,18)]
chrvect = unlist(lapply(vect,function(x){
  row = match(paste("q",x,sep = ""),sigmark$jQTL_ID)
  snp = as.character(sigmark[row,2])
  col = match(snp,dat[1,])
  return(as.numeric(dat[2,col]))
}))
#####Matching chromosome values to QTL indices#####
qtls = list()
for (x in 1:length(vect)){
  assign(paste("famj",vect[x], sep = ""),qtlmaker(cross = family, sigmark = sigmark, chr = chrvect[x], index= vect[x]))
  qtls = c(qtls,qtlmaker(cross = family, sigmark = sigmark, chr = chrvect[x], index= vect[x]))
}
class(qtls)
qtls = lapply(1:length(vect), function(x){
  fitqtl(cross = family, qtl = qtlmaker(cross = family, sigmark = sigmark, chr = chrvect[x], index = vect[x]))
})
extrpos <- function(x){
  as.numeric(substr(x,regexpr("@",x)[1]+1, nchar(x)))
}
x5 = extrpos(rownames(qtls[[1]]$result.drop)[which(qtls[[1]]$result.drop[,7] == min(qtls[[1]]$result.drop[,7]))])
x6 = extrpos(x = rownames(qtls[[1]]$result.drop)[1])
x7 = extrpos(x = rownames(qtls[[1]]$result.drop)[dim(qtls[[1]]$result.drop)[1]])
lod = unlist(lapply(1:length(vect), function(x){
  qtls[[x]]$result.full[1,4]
}))
r2 = unlist(lapply(1:length(vect), function(x){
  qtls[[x]]$result.full[1,5]
}))
x5 = unlist(lapply(1:length(vect), function(x){
  unique(extrpos(rownames(qtls[[x]]$result.drop)[which(qtls[[x]]$result.drop[,7] == min(qtls[[x]]$result.drop[,7]))]))
}))
x6 = unlist(lapply(1:length(vect), function(x){
  extrpos(x = rownames(qtls[[x]]$result.drop)[1])
  }))
x7 = unlist(lapply(1:length(vect), function(x){
  extrpos(x = rownames(qtls[[x]]$result.drop)[dim(qtls[[x]]$result.drop)[1]])
}))

l= length(vect)
chrvect = lapply(chrvect,function(x){
  paste("chr",x,sep ="")
})

rdf = cbind(paste("QTL",vect, sep = ""), rep("Trait 1", l), rep("ID:0",l),
            rep("P1",l),rep("Y1", l), chrvect,
            rep("linkage_group\ 1", l), round(lod,2), round(r2),
            x5,x6,x7)
colnames(rdf) = c("QTL", "Trait", "ID",
                  "Place", "Year", "Chromosome",
                  "Linkage Group", "LOD Score", "R Squared",
                  "Most Likely QTL Location", "QTL Start Location", "QTL Stop Location")

write.table(rdf, file = "Github/Stapleton-Lab/mQTL\ Random\ and\ Family/Famj_QTLs.txt",
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
#Now for the random data
dat = read.csv(file = "Github/Stapleton-Lab/mQTL\ Random\ and\ Family/Input\ Files/Random2.csv", header = FALSE, stringsAsFactors = FALSE)
random <-read.cross(dir = "Github/Stapleton-Lab/mQTL\ Random\ and\ Family/Input\ Files", file = "Random2.csv")
sigmark <- read.csv(file ="Github/Stapleton-Lab/mQTL\ Random\ and\ Family/Input\ Files/adj_jointP_random.csv", header = TRUE)
random <- drop.nullmarkers(random)

#scan with variance
random <- calc.genoprob(random)
routv <- scanonevar(cross = random,
                    mean.formula = height.in. ~ mean.QTL.add,
                    var.formula = ~ var.QTL.add)
ran1 = NULL
random = sim.geno(random)
sigmrk1 = sigmark[-which(match(sigmark$qtl_id, "")== 1),]
sig1 = sigmrk1$Position__cM_[which(sigmrk1$qtl_id == "q1")]
rand1 = makeqtl(cross = random, chr = rep(1,length(sig1)), pos = sig1)

qtlmaker = function(cross,sigmark, chr, index){
  qx = paste("q",index, sep = "")
  famx = paste("fam",index, sep = "")
  sigmrk = sigmark[-which(match(sigmark$qtl_id,"") == 1),]
  sig = sigmrk[which(sigmrk1$qtl_id == qx),]
  sigpos = sig[,2]
  names = as.character(sig$SNP_Names)
  
  famx = makeqtl(cross = cross, chr = rep(chr,length(sigpos)), pos = sigpos)
  return(famx)
}
vect <- 1:6
chrvect = unlist(lapply(vect,function(x){
  row = which(sigmark$qtl_id == paste("q",x,sep = ""))
  snp = as.character(sigmark[row,1])
  col = match(snp,dat[1,])
  return(as.numeric(dat[2,col])[1])
}))
#####Matching chromosome values to QTL indices#####
qtls = list()
for (x in 1:length(vect)){
  assign(paste("randj",vect[x], sep = ""),qtlmaker(cross = random, sigmark = sigmark, chr = chrvect[x], index= vect[x]))
  qtls = c(qtls,qtlmaker(cross = random, sigmark = sigmark, chr = chrvect[x], index= vect[x]))
}
class(qtls)
qtls = lapply(1:length(vect), function(x){
  fitqtl(cross = random, qtl = qtlmaker(cross = random, sigmark = sigmark, chr = chrvect[x], index = vect[x]))
})
extrpos <- function(x){
  as.numeric(substr(x,regexpr("@",x)[1]+1, nchar(x)))
}
x5 = extrpos(rownames(qtls[[1]]$result.drop)[which(qtls[[1]]$result.drop[,7] == min(qtls[[1]]$result.drop[,7]))])
x6 = extrpos(x = rownames(qtls[[1]]$result.drop)[1])
x7 = extrpos(x = rownames(qtls[[1]]$result.drop)[dim(qtls[[1]]$result.drop)[1]])
lod = unlist(lapply(1:length(vect), function(x){
  qtls[[x]]$result.full[1,4]
}))
r2 = unlist(lapply(1:length(vect), function(x){
  qtls[[x]]$result.full[1,5]
}))
x5 = unlist(lapply(1:length(vect), function(x){
  unique(extrpos(rownames(qtls[[x]]$result.drop)[which(qtls[[x]]$result.drop[,7] == min(qtls[[x]]$result.drop[,7]))]))
}))
x6 = unlist(lapply(1:length(vect), function(x){
  extrpos(x = rownames(qtls[[x]]$result.drop)[1])
}))
x7 = unlist(lapply(1:length(vect), function(x){
  extrpos(x = rownames(qtls[[x]]$result.drop)[dim(qtls[[x]]$result.drop)[1]])
}))

x5 <- c(x5[1:4],4158.68, x5[5])
x6 <- c(x6[1:4],4158.68, x6[5])
x7 <- c(x7[1:4],4158.68, x7[5])

l= length(vect)


rdf = cbind(paste("QTL",vect, sep = ""), rep("Trait 1", l), rep("ID:0",l),
            rep("P1",l),rep("Y1", l), chrvect,
            rep("linkage_group_one", l), lod, r2,
            x5,x6,x7)
colnames(rdf) = c("QTL", "Trait", "ID",
                  "Place", "Year", "Chromosome",
                  "Linkage Group", "LOD Score", "R Squared",
                  "Most Likely QTL Location", "QTL Start Location", "QTL Stop Location")

write.table(rdf, file = "Github/Stapleton-Lab/mQTL\ Random\ and\ Family/Randj1_QTLs.txt",
            sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)

