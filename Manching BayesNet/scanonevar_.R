function (modeling.df, genoprob.df, loc.info.df, scan.types, 
          scan.formulae, return.covar.effects) 
{
  loc.name <- "fake_global_for_CRAN"
  result <- initialize.scanonevar.result_(loc.info.df = loc.info.df, 
                                          scan.types = scan.types, scan.formulae = scan.formulae, 
                                          return.covar.effects = return.covar.effects)
  mean.df <- sum(grepl(pattern = "mean.QTL", x = labels(stats::terms(scan.formulae[["mean.alt.formula"]]))))
  var.df <- sum(grepl(pattern = "var.QTL", x = labels(stats::terms(scan.formulae[["var.alt.formula"]]))))
  if ("joint" %in% scan.types) {
    joint.null.fit <- fit.model.00_(formulae = scan.formulae, 
                                    df = modeling.df)
  }
  for (loc.idx in 1:nrow(result)) {
    this.loc.name <- result[["loc.name"]][loc.idx]
    loc.genoprobs <- dplyr::filter(.data = genoprob.df, 
                                   loc.name == this.loc.name)
    this.loc.modeling.df <- make.loc.specific.modeling.df(general.modeling.df = modeling.df, 
                                                          loc.genoprobs = loc.genoprobs, model.formulae = scan.formulae)
    if (all(this.loc.modeling.df[["mean.QTL.dom"]] == 0)) {
      this.loc.mean.df <- mean.df - 1
    }
    else {
      this.loc.mean.df <- mean.df
    }
    if (all(this.loc.modeling.df[["var.QTL.dom"]] == 0)) {
      this.loc.var.df <- var.df - 1
    }
    else {
      this.loc.var.df <- var.df
    }
    alternative.fit <- fit.model.mv_(formulae = scan.formulae, 
                                     df = this.loc.modeling.df)
    if (all(return.covar.effects)) {
      if (identical(alternative.fit, NA)) {
        if (loc.idx == 1) {
          stop("Cant fit model on locus 1.  Due to programming weirdness, cant return effect estimates.")
        }
        coef_mtx <- rbind(coef_mtx, rep(NA, ncol(coef_mtx)))
      }
      else {
        mean_coef_mtx <- stats::coef(summary(alternative.fit))
        mean_ests <- mean_ses <- rep(NA, length(stats::coef(alternative.fit)))
        mean_ests[!is.na(stats::coef(alternative.fit))] <- mean_coef_mtx[, 
                                                                         "Estimate"]
        names(mean_ests) <- paste0("mef_", names(stats::coef(alternative.fit)))
        mean_ses[!is.na(stats::coef(alternative.fit))] <- mean_coef_mtx[, 
                                                                        "Std. Error"]
        names(mean_ses) <- paste0("mse_", names(stats::coef(alternative.fit)))
        disp_fit <- alternative.fit$dispersion.fit
        disp_coef_mtx <- stats::coef(summary(disp_fit))
        disp_ests <- disp_ses <- rep(NA, length(stats::coef(disp_fit)))
        disp_ests[!is.na(stats::coef(disp_fit))] <- disp_coef_mtx[, 
                                                                  "Estimate"]
        names(disp_ests) <- paste0("vef_", names(stats::coef(disp_fit)))
        disp_ses[!is.na(stats::coef(disp_fit))] <- disp_coef_mtx[, 
                                                                 "Std. Error"]
        names(disp_ses) <- paste0("vse_", names(stats::coef(disp_fit)))
        if (loc.idx == 1) {
          coef_mtx <- matrix(data = c(mean_ests, mean_ses, 
                                      disp_ests, disp_ses), nrow = 1)
        }
        else {
          coef_mtx <- rbind(coef_mtx, c(mean_ests, mean_ses, 
                                        disp_ests, disp_ses))
        }
      }
    }
    if ("mean" %in% scan.types) {
      mean.null.fit <- fit.model.0v_(formulae = scan.formulae, 
                                     df = this.loc.modeling.df)
      result[["mean.lod"]][loc.idx] <- LOD(alt = alternative.fit, 
                                           null = mean.null.fit)
      result[["mean.asymp.p"]][loc.idx] <- stats::pchisq(q = result[["mean.lod"]][loc.idx], 
                                                         df = this.loc.mean.df, lower.tail = FALSE)
    }
    if ("var" %in% scan.types) {
      var.null.fit <- fit.model.m0_(formulae = scan.formulae, 
                                    df = this.loc.modeling.df)
      result[["var.lod"]][loc.idx] <- LOD(alt = alternative.fit, 
                                          null = var.null.fit)
      result[["var.asymp.p"]][loc.idx] <- stats::pchisq(q = result[["var.lod"]][loc.idx], 
                                                        df = this.loc.var.df, lower.tail = FALSE)
    }
    if ("joint" %in% scan.types) {
      result[["joint.lod"]][loc.idx] <- LOD(alt = alternative.fit, 
                                            null = joint.null.fit)
      result[["joint.asymp.p"]][loc.idx] <- stats::pchisq(q = result[["joint.lod"]][loc.idx], 
                                                          df = this.loc.mean.df + this.loc.var.df, lower.tail = FALSE)
    }
  }
  if (return.covar.effects) {
    result <- cbind(result, coef_mtx)
  }
  class(result) <- c("scanonevar_result", class(result))
  return(result)
}