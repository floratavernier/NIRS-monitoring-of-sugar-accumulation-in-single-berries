# Wold criterion for determining the number of LVs
optLV <- function(mod, crit) {
  if(crit == 0) {
    which.min( mod$validation$PRESS)
  } else {
    R <- vector(mode = "numeric", length = mod$ncomp)
    for(i in 1:mod$ncomp) {
      R[i] <- (mod$validation$PRESS[i + 1] / mod$validation$PRESS[i])
    }
    if(max(R, na.rm = TRUE) < crit) {mod$ncomp} else {which(R >= crit)[1]}
  }
}

# Definition of outlier according to the distribution of residuals from a linear regression between obs and pred
outlier <- function(mod, threshold = 0, crit){
  residReg <- residuals(lm(mod$validation$pred[,, optLV(mod, crit)] ~ mod$model$trait))
  if(threshold == 0){
    return(names(boxplot(residReg, plot = FALSE)$out))
  } else {
    sigma <- vector(mode = 'numeric', length = length(residReg))
    for (i in 1:length(residReg)) {sigma[i] <- sd(residReg[-i])}
    zscore <- abs(residReg) / sigma
    p_value <- 1 - pnorm(zscore)
    return(which(p_value <= threshold))
  }
}

# Identification of outliers in a loop according to a threshold
dropOutliersLOO <- function(dat, trait, maxcomp = 20, threshold = 0, crit = 1, maxsteps = 20){
  outliers <- list()
  #step 1
  dataSet <- data.frame("echantillon" = dat$echantillon, "trait" = dat[, which(colnames(dat) == trait)], "NIRS" = I(dat$NIRS))
  modele <- mvr(trait ~ NIRS, data = dataSet, ncomp = maxcomp, validation = "LOO")
  outliers[[1]] <- as.vector(dataSet[outlier(mod = modele, threshold = threshold, crit = crit), ]$echantillon)
  cat(1)
  #step 2
  dataSet <- dataSet[!dataSet$echantillon %in% outliers[[1]], ]
  modele <-  mvr(trait ~ NIRS, data = dataSet, ncomp = maxcomp, validation = "LOO")
  outliers[[2]] <- c(outliers[[1]], as.vector(dataSet[outlier(mod = modele, threshold = threshold, crit = crit), ]$echantillon))
  cat(2)
  #steps 3 to i
  for (i in 3:maxsteps){
    if (identical(outliers[[i - 2]], outliers[[i - 1]])) break else {
      dataSet <- dataSet[!dataSet$echantillon %in% outliers[[i-1]], ]
      modele <-  mvr(trait ~ NIRS, data = dataSet, ncomp = maxcomp, validation = "LOO")
      outliers[[i]] <- c(outliers[[i-1]], as.vector(dataSet[outlier(mod = modele, threshold = threshold, crit = crit), ]$echantillon))
    }
    cat(i)
  }
  cat("\n")
  return(outliers[[length(outliers)]])
}

# Repeated cross-validation
MCCV <- function(dat, trait, maxcomp = 20, fold = 5, iter = 50, crit = 1){
  ncomp <- min(ncol(dat$NIRS), maxcomp)
  dataSet <- data.frame("echantillon" = dat$echantillon, "trait" = dat[, which(colnames(dat) == trait)], "NIRS" = I(dat$NIRS))
  modele <- mvr(trait ~ NIRS, data = dataSet, ncomp = ncomp)
  R2.train <- 1 - (colSums(modele$residuals[, 1, ] ^ 2) / (var(dataSet$trait) * (nrow(dataSet) - 1)))
  MC.PRESS.mat <- matrix(NA, nrow = ncomp, ncol = iter)
  predMat <- array(NA, dim = c(nrow(dataSet), ncomp, iter))
  dimnames(predMat)[[1]] <- rownames(dataSet)
  dimnames(predMat)[[2]] <- paste0(1:ncomp, "_comp")
  dimnames(predMat)[[3]] <- paste0("iter_", 1:iter)
  for (i in 1:iter){
    CV <- crossval(modele, segments = fold, data = dataSet, segment.type = "random")
    MC.PRESS.mat[, i] <- CV$validation$PRESS
    predMat[, , i] <- CV$validation$pred[, 1, ]
    rm(CV)
    cat(i)
  }
  cat('\n')
  MC.R2.CV.mat <- 1 - (MC.PRESS.mat / (var(dataSet$trait) * (nrow(dataSet) - 1)))
  MC.RMSEP.mat <- sqrt(MC.PRESS.mat / nrow(dataSet))
  MC.RPD.mat <- sd(dataSet$trait) / MC.RMSEP.mat
  if(crit == 0){
    optComp <- which.min(rowMeans(MC.PRESS.mat))
  } else {
    R <- vector(mode = "numeric", length = ncomp)
    for(i in 1:ncomp) {R[i] <- (rowMeans(MC.PRESS.mat)[i + 1] / rowMeans(MC.PRESS.mat)[i])}
    if(max(R, na.rm = TRUE) >= crit){optComp = which(R >= crit)[1]} else {optComp = length(R)}
  }
  R2.MCCV.mean <- rowMeans(MC.R2.CV.mat)
  R2.MCCV.sd <- apply(MC.R2.CV.mat, 1, sd)
  RMSE.MCCV.mean <- rowMeans(MC.RMSEP.mat)
  RMSE.MCCV.sd <- apply(MC.RMSEP.mat, 1, sd)
  RPD.MCCV.mean <- rowMeans(MC.RPD.mat)
  RPD.MCCV.sd <- apply(MC.RPD.mat, 1, sd)
  if(iter == 1){
    predMatOptComp <- predMat[, optComp, ]
  } else {
    predMatOptComp <- rowMeans(predMat[, optComp, ])
  }
  output <- c("nbcomp" = optComp[1], "R2.train" = R2.train[optComp],
              "R2.MCCV.mean" = R2.MCCV.mean[optComp], "R2.MCCV.sd" = R2.MCCV.sd[optComp],
              "RMSE.MCCV.mean" = RMSE.MCCV.mean[optComp], "RMSE.MCCV.sd" = RMSE.MCCV.sd[optComp],
              "RPD.MCCV.mean" = RPD.MCCV.mean[optComp], "RPD.MCCV.sd" = RPD.MCCV.sd[optComp])
  return(list("model" = modele, "output" = output, "predicted" = predMatOptComp))
}

# Plot MCCV model
plotMCCV<-function(MCCV, ...){
  plot(MCCV$model$model$trait, MCCV$predicted, xlab = "observed", ylab = "predicted", ...)
  abline(0,1, col = "black")
  legend("topleft", bty = "n", c(
    paste("nb_comp = ", MCCV$output[1], sep = ""),
    paste("R2_train = ", round(MCCV$output[2], 2) , sep = ""),
    paste("R2_cv = ", round(MCCV$output[3], 2), sep = ""),
    paste("RMSE_cv = ", round(MCCV$output[5], 2), sep = ""),
    # paste("RPD_cv = ", round(MCCV$output[7], 2), sep = ""),
    paste("nobs = ", length(MCCV$predicted), sep = "")), text.col = "black")
}

# Calibration wrapper
calibrationSPIR <- function(calibSetList, chem, maxcomps = 20, LOO = FALSE, folds = 4, iterMCCV = 50,
                            wold = 1, thresh = NULL, maxstepOutliers = 20) {
  if(is.null(thresh)){
    outliers <- NULL
    if(isTRUE(LOO)) {
      iterMCCV = 1
      folds = unique(sapply(calibSetList, nrow))
    }
    calibMCCV <- list()
    for (i in 1:length(calibSetList)){
      calibMCCV[[i]] <- MCCV(dat = calibSetList[[i]], trait = chem, maxcomp = maxcomps, fold = folds, iter = iterMCCV, crit = wold)
    }
    names(calibMCCV) <- names(calibSetList)
    # calibMCCV <- lapply(calibSetList, function(x) {
    #   MCCV(dat = x, trait = chem, maxcomp = maxcomps, fold = folds, iter = iterMCCV, crit = wold)
    # })
    outcalib <- data.frame(t(sapply(calibMCCV,function(x){x$output})),
                           "n" = sapply(calibSetList, nrow),
                           "nb.outliers" = rep(0, length(calibSetList)))
  }
  else {
    outliers <- lapply(calibSetList, function(x){
      dropOutliersLOO(dat = x, trait = chem, maxcomp = maxcomps, threshold = thresh, crit = wold, maxsteps = maxstepOutliers)
    })
    calibFiltList <- calibSetList
    for (i in 1:length(calibFiltList)){
      calibFiltList[[i]] <- subset(calibSetList[[i]], ! echantillon %in% outliers[[i]])
    }
    if(isTRUE(LOO)){
      iterMCCV = 1
      folds = unique(sapply(calibFiltList, nrow))
    }
    calibMCCV <- list()
    for (i in 1:length(calibFiltList)){
      calibMCCV[[i]] <- MCCV(dat = calibFiltList[[i]], trait = chem, maxcomp = maxcomps, fold = folds, iter = iterMCCV, crit = wold)
    }
    names(calibMCCV) <- names(calibFiltList)
    # calibMCCV <- lapply(calibFiltList, function(x) {
    #   MCCV(dat = x, trait = chem, maxcomp = maxcomps, fold = folds, iter = iterMCCV, crit = wold)
    # })
    outcalib <- data.frame(t(sapply(calibMCCV,function(x){x$output})),
                           "n" = sapply(calibFiltList, nrow),
                           "nb.outliers" = sapply(outliers, length))
  }
  output <- list("chem" = chem,
                 "thresh" = thresh,
                 "outliers" = outliers,
                 "outcalib" = outcalib)
  return(output)
}

# Plot of a model
stat_pred <- function(dataSet, obs, pred) {
  dataSet_ok <- na.omit(dataSet)
  observed_data <- dataSet_ok[, colnames(dataSet_ok) == obs]
  predicted_data <- dataSet_ok[, colnames(dataSet_ok) == pred]
  PRESS <- sum((observed_data - predicted_data) ^ 2)
  TSS <- sum((observed_data - mean(observed_data)) ^ 2)
  R2 <- 1 - (PRESS / TSS)
  RMSE <- sqrt(PRESS / nrow(dataSet_ok))
  RPD <- sd(observed_data) / RMSE
  modlin <- lm(predicted_data ~ observed_data)
  return(list(
    "output" = c(
      "R2_reg" = summary(modlin)$r.sq,
      "R2_val" = R2,
      "RMSE" = RMSE,
      "RPD" = RPD,
      "n" = nrow(dataSet_ok)
    ),
    "modlin" = modlin
  ))
}