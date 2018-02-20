#####Multiple Stress vQTL#####
#to be run in stampede2
library(tidyverse)
library(qtl);
library(vqtl);
setwd("C:/Users/Thomas/Documents/Github/Stapleton-Lab/Manching BayesNet")
#####Constructing a  Small scale sample set#####
fullobj <- read_csv(file = "ManchingScrubbed.csv")
set.seed(569612)
srows <- c(1,2,sample(3:dim(fullobj)[1], 500))
scols <- c(1:4,sample(5:dim(fullobj)[2], 500))
sampobj <- fullobj[srows,scols]
sampobj[1:2,1:4] = ""
colnames(sampobj) <- make.names(colnames(sampobj))
write_csv(sampobj, path = "MSSample.csv")
#####Small Scale Analysis#####
cross = read.cross(format = "csv", file = "MSSample.csv");
cross = drop.nullmarkers(cross);
cgp = calc.genoprob(cross = cross);
scanv = scanonevar(cross = cgp, mean.formula = Height ~ Low.Water + Low.Nitrogen + Pathogen +  mean.QTL.add ,
                   var.formula =  ~ Low.Water + Low.Nitrogen + Pathogen + var.QTL.add);
library("dplyr");
library("dglm")
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
};
#set up a vector to run the function on
y = 1:length(scanv$result$loc.name);
#effect sizes can not be computed for these 3 SNPs so we remove them from the vector
#populating a dataframe with effect size estimates
sizedf = sapply(y, function(x){
  tempm =  effect.sizes(cross = cgp,
                        phenotype.name = "Height",
                        genotype.names = c("AA","BB"),
                        focal.groups = scanv$result$loc.name[x])
  tempv = c(tempm[1,2:7],tempm[2,2:7])
  return(unlist(tempv))
});
#gathering data from the initial scan
scanvdf<- data.frame(scanv$result$loc.name,
                     scanv$result$pos,
                     scanv$result$mean.lod,
                     scanv$result$mean.asymp.p,
                     scanv$result$var.lod,
                     scanv$result$var.asymp.p,
                     scanv$result$joint.lod,
                     scanv$result$joint.asymp.p)
#combining both 
scanvdf = cbind(scanvdf,t(sizedf))
colnames(scanvdf) = c("SNP Name",
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
write.csv(scanvdf, file = "MultipleStress_Sample_vQTL_Height_LOD,Pvals,EffectSizes.csv")
#####Full Sized Analysis####
cross = read.cross(format = "csv", file = "ManchingScrubbed_GenoOnly.csv");
cross = drop.nullmarkers(cross);
evf = read_csv(file = "ManchingScrubbed_Environment.csv");
cross[['pheno']][['LowWater']] = factor(evf$`Low Water`[-1:-2]);
cross[['pheno']][['LowNitrogen']] = factor(evf$`Low Nitrogen`[-1:-2]);
cross[['pheno']][['Pathogen']] = factor(evf$Pathogen[-1:-2]);
gc();
cgp = calc.genoprob(cross = cross);
gc();
scanv = scanonevar(cross = cgp, mean.formula = Height ~ LowWater + LowNitrogen + Pathogen +  mean.QTL.add ,
                   var.formula =  ~ LowWater + LowNitrogen + Pathogen + var.QTL.add);
library("dplyr");
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
};
#set up a vector to run the function on
 y = 1:length(scanv$result$loc.name);
nnames <- sapply(scanv$result$loc.name, function(x){
   if (substr(x,1,3) == "chr"){
     return(NA)
   } else return(x)
})
nnames <- nnames[which(!is.na(nnames))]
#effect sizes can not be computed for these 3 SNPs so we remove them from the vector
#populating a dataframe with effect size estimates
sizedf = sapply(1:length(nnames), function(x){
  print(x)
  tryCatch({
    tempm =  effect.sizes(cross = cgp,
                        phenotype.name = "Height",
                        genotype.names = c("AA","BB"),
                        focal.groups = nnames[x])
    }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")},
    finally = {
    tempv = c(tempm[1,2:7],tempm[2,2:7])
  return(unlist(tempv))
  })
});
  #gathering data from the initial scan
scanvdf<- data.frame(scanv$result$loc.name,
                     scanv$result$pos,
                     scanv$result$mean.lod,
                     scanv$result$mean.asymp.p,
                     scanv$result$var.lod,
                     scanv$result$var.asymp.p,
                     scanv$result$joint.lod,
                     scanv$result$joint.asymp.p)
#combining both 
scanvdf = cbind(scanvdf[which(scanv$result$loc.name %in% nnames),],t(sizedf))
colnames(scanvdf) = c("SNP Name",
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
write_csv(scanvdf, path = "MultipleStressvQTL_Height_LOD,Pvals,EffectSizes1.csv")
