#Manching Data Analysis
#edited all the Mo/M0/MO to be MO for consistancy across the sheets

#read in the data
library("qtl")
library("vqtl")
geno<-read.csv("C:/Users/Thomas/Documents/GitHub/Stapleton-Lab/vQTL IBM and Manching/Data/IBM94markerset08seq.csv")
randp<- read.csv("C:/Users/Thomas/Documents/GitHub/Stapleton-Lab/vQTL IBM and Manching/Data/RandIBM642006PyearData.csv")
famp<- read.csv("C:/Users/Thomas/Documents/GitHub/Stapleton-Lab/vQTL IBM and Manching/Data/Pyear2006fielddata.csv")
#tidying data with bad genotypes (MO062, MO075, and MO066)
bad = NULL
for (w in 1:length(randp$Line)){
  if (as.character(randp$Line[w]) == "MO075" | as.character(randp$Line[w]) == "MO062"| as.character(randp$Line[w]) == "MO066"){
    bad = c(bad,w)
  }
}
randp = randp[-bad,]
#start forming a cross object
dat = data.frame(randp$height.in.,randp$NumTasselBranches,randp$TasselBranchAngle.deg.,randp$Line)
for(x in 1:3235){
  dat = cbind(dat,"0")
}
colnames(dat) <- NULL
dat <- data.frame(lapply(dat, as.character), stringsAsFactors=FALSE)
#stiching genotypes to phenotypes
for (x in 1:dim(dat)[1]){
a <- as.character(dat[x,4])
if (0 == length(which(colnames(geno) == a))){
  next
}
vect = as.character((geno[,which(colnames(geno) == a)]))
for (q in 1:3235){
  dat[x,q+4] = vect[q]
}
}
beep()
#looking for columns that did not populate
empty = which(dat[,10] == "0")
nonpop = NULL
for (x in 1:length(empty)){
nonpop = c(nonpop,dat[empty[x],2])
}
nonpop

#all these columns are either B73 or have not Line as such, we discard them
randp = randp[-empty,]
#now we have a solid cross object
inc = as.numeric(geno$incre_new)
colnames(dat) = c("Height","Number of Branches","Branch Angle","Line",names)
#we want to add the increment values
dat = dat[,-4]
vect = c(rep("",3),inc)
dat = rbind(vect,dat)
write.csv(dat, file = "Github/Stapleton-Lab/vQTL IBM and Manching/cross.csv")
#then delete the first column
#####With the cross#####
cross = read.cross(dir = "Github/Stapleton-Lab/vQTL IBM and Manching/", file= "cross.csv")cross <- drop.nullmarkers(cross)
#scan with variance
cross <- calc.genoprob(cross)
outv <- scanonevar(cross = cross,
                    mean.formula = Height ~ mean.QTL.add,
                    var.formula = ~ var.QTL.add)
cross= sim.geno(cross)
library("dplyr")
effect.sizes = function (cross, phenotype.name, focal.groups = NULL, nuisance.groups = NULL, 
                         genotype.names = c("AA", "AB", "BB"), xlim = NULL, ylim = NULL, 
                         title = paste(phenotype.name, "by", paste(focal.groups, 
                                                                   collapse = ", ")), draw_ribbons = TRUE, se_line_size = 1, 
                         point_size = 1) 
{
  indiv.mean.estim <- indiv.mean.lb <- indiv.mean.ub <- "fake_global_for_CRAN"
  indiv.sd.estim <- indiv.sd.lb <- indiv.sd.ub <- "fake_global_for_CRAN"
  group.mean.estim <- group.mean.ub <- group.mean.lb <- "fake_global_for_CRAN"
  group.sd.estim <- group.sd.ub <- group.sd.lb <- "fake_global_for_CRAN"
  modeling.df <- dplyr::data_frame(placeholder = rep(NA, qtl::nind(cross)))
  modeling.df[[phenotype.name]] <- cross[["pheno"]][[phenotype.name]]
  marker.names <- c(focal.groups[focal.groups %in% colnames(qtl::pull.geno(cross = cross))], 
                    nuisance.groups[nuisance.groups %in% colnames(qtl::pull.geno(cross = cross))])
  phen.names <- c(focal.groups[focal.groups %in% colnames(qtl::pull.pheno(cross = cross))], 
                  nuisance.groups[nuisance.groups %in% colnames(qtl::pull.pheno(cross = cross))])
  for (marker.name in marker.names) {
    modeling.df[[marker.name]] <- factor(x = qtl::pull.geno(cross = cross)[, 
                                                                           marker.name], labels = genotype.names)
  }
  for (phen.name in phen.names) {
    modeling.df[[phen.name]] <- factor(qtl::pull.pheno(cross = cross)[[phen.name]])
  }
  modeling.df[["placeholder"]] <- NULL
  covar.form <- paste(focal.groups, collapse = "+")
  if (!is.null(nuisance.groups)) {
    covar.form <- paste(covar.form, "+", paste(nuisance.groups, 
                                               collapse = "+"))
  }
  mean.form <- paste(phenotype.name, "~", covar.form)
  var.form <- paste("~", covar.form)
  dglm.fit <- dglm::dglm(formula = stats::formula(mean.form), 
                         dformula = stats::formula(var.form), data = modeling.df)
  mean.pred <- stats::predict(dglm.fit, se.fit = TRUE)
  mean.estim <- mean.pred$fit
  mean.se <- mean.pred$se.fit
  sd.pred <- stats::predict(dglm.fit$dispersion.fit, se.fit = TRUE)
  sd.estim <- sd.pred$fit/sd.pred$residual.scale
  sd.se <- sd.pred$se.fit
  indiv.prediction.tbl <- dplyr::bind_cols(stats::na.omit(modeling.df), 
                                           dplyr::data_frame(indiv.mean.estim = mean.estim, indiv.mean.lb = mean.estim - 
                                                               mean.se, indiv.mean.ub = mean.estim + mean.se, indiv.sd.estim = exp(sd.estim), 
                                                             indiv.sd.lb = exp(sd.estim - sd.se), indiv.sd.ub = exp(sd.estim + 
                                                                                                                      sd.se)))
  group.prediction.tbl <- indiv.prediction.tbl %>% dplyr::group_by_(.dots = c(focal.groups)) %>% 
    dplyr::summarise(group.mean.estim = mean(indiv.mean.estim), 
                     group.mean.lb = mean(indiv.mean.lb), group.mean.ub = mean(indiv.mean.ub), 
                     group.sd.estim = mean(indiv.sd.estim), group.sd.lb = mean(indiv.sd.lb), 
                     group.sd.ub = mean(indiv.sd.ub))
  return(group.prediction.tbl)
}
#this works
#this is an example of the first one
debug(effect.sizes)
sizes = effect.sizes(cross = cross,
                     phenotype.name = "Height",
                     genotype.names = c("AA","BB"),
                     focal.groups = outv$result$loc.name[1])
sizes

sizedf <- data.frame(NULL)

y = 1:length(outv$result$loc.name)
y = y[-c(458,2482,2483)]
#effect sizes can not be computed for these 3 SNPs
for (x in y){
  tempm =  effect.sizes(cross = cross,
                        phenotype.name = "Height",
                        genotype.names = c("AA","BB"),
                        focal.groups = outv$result$loc.name[x])
  tempv = c(tempm[1,2:7],tempm[2,2:7])
  sizedf = rbind(sizedf,tempv)
}
undebug(effect.sizes)
effect.sizes(cross = cross,
             phenotype.name = "height.in.",
             genotype.names = c("AA","BB"),
             focal.groups = routv$result$loc.name[470])


routvdf<- data.frame(routv$result$loc.name,
                     routv$result$pos,
                     routv$result$mean.lod,
                     routv$result$mean.asymp.p,
                     routv$result$var.lod,
                     routv$result$var.asymp.p,
                     routv$result$joint.lod,
                     routv$result$joint.asymp.p)
#dropping the SNPs whose effect sizes could not be computed
routvdf = routvdf[-c(458,2482,2483),]
routvdf = cbind(routvdf,rsizedf)
colnames(routvdf) = c("SNP Names",
                      "Position (cM)",
                      "Mean LOD",
                      "Mean P Value",
                      "Variance LOD",
                      "Variance P Value",
                      "Joint LOD",
                      "Joint P Value",
                      "A Mean Est",
                      "A Mean Lower Bound",
                      "A Mean Upper Bound",
                      "A Standard Deviation Est",
                      "A Standard Deviation Lower Bound",
                      "A Standard Deviation Upper Bound",
                      "B Mean Est",
                      "B Mean Lower Bound",
                      "B Mean Upper Bound",
                      "B Standard Deviation Est",
                      "B Standard Deviation Lower Bound",
                      "B Standard Deviation Upper Bound")

write.csv(routvdf, file = "crossvQTL_LOD,Pvals,EffectSizes.csv")
